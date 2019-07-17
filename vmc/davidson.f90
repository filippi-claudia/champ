!> \namespace Davidson eigensolver
!> \author Felipe Zapata
!> The current implementation uses a general  davidson algorithm, meaning
!> that it compute all the eigenvalues simultaneusly using a variable size block approach.
!> The family of Davidson algorithm only differ in the way that the correction
!> vector is computed.
!> Computed pairs of eigenvalues/eigenvectors are deflated using algorithm
!> described at: https://doi.org/10.1023/A:101919970

module davidson_free

  use numeric_kinds, only: dp
  use lapack_wrapper, only: lapack_generalized_eigensolver, lapack_matmul, lapack_matrix_vector, &
       lapack_qr, lapack_solver
  use array_utils, only: concatenate, generate_preconditioner, norm
  implicit none

  type davidson_parameters
     INTEGER :: nparm
     INTEGER :: nparm_max
     INTEGER :: lowest
     INTEGER :: max_dim_sub
     INTEGER :: mvec
  end type davidson_parameters
  
  !> \private
  private
  !> \public
  public :: davidson_parameters, generalized_eigensolver

  
contains

  subroutine generalized_eigensolver(fun_mtx_gemv, eigenvalues, ritz_vectors, nparm, nparm_max, &
       lowest, mvec, method, max_iters, tolerance, iters, max_dim_sub, fun_stx_gemv, idtask)
    !> \brief use a pair of functions fun_mtx and fun_stx to compute on the fly the matrices to solve
    !>  the general eigenvalue problem
    !> The current implementation uses a general  davidson algorithm, meaning
    !> that it compute all the eigenvalues simultaneusly using a block approach.
    !> The family of Davidson algorithm only differ in the way that the correction
    !> vector is computed.
    
    !> \param[in] fun_mtx_gemv: Function to apply the matrix to a buncof vectors
    !> \param[in, opt] fun_stx_gemv: (optional) function to apply the pencil to a bunch of vectors.
    !> \param[out] eigenvalues Computed eigenvalues
    !> \param[inout] ritz_vectors approximation to the eigenvectors
    !> \param nparm[in] Leading dimension of the matrix to diagonalize
    !> \param nparm_max[in] Maximum dimension of the matrix to diagonalize
    !> \param[in] lowest Number of lowest eigenvalues/ritz_vectors to compute
    !> \param[in] method Method to compute the correction vector. Available
    !> methods are,
    !>    DPR: Diagonal-Preconditioned-Residue
    !>    GJD: Generalized Jacobi Davidson
    !> \param[in] max_iters: Maximum number of iterations
    !> \param[in] tolerance norm-2 error of the eigenvalues
    !> \param[in] method: Method to compute the correction vectors
    !> \param[in, opt] max_dim_sub: maximum dimension of the subspace search   
    !> \param[out] iters: Number of iterations until convergence
    !> \return eigenvalues and ritz_vectors of the matrix `mtx`

    implicit none

    include 'mpif.h'

    ! input/output variable
    integer, intent(in) :: nparm, nparm_max, mvec, lowest, idtask
    real(dp), dimension(lowest), intent(out) :: eigenvalues
    real(dp), dimension(:, :), intent(inout) :: ritz_vectors
    integer, intent(in) :: max_iters
    integer, intent(in), optional :: max_dim_sub
    real(dp), intent(in) :: tolerance
    character(len=*), intent(in) :: method
    integer, intent(out) :: iters
    
    ! Function to compute the target matrix on the fly
    interface

     function fun_mtx_gemv(parameters, input_vect) result(output_vect)
       !> \brief Function to compute the action of the hamiltonian on the fly
       !> \param[in] dimension of the arrays to compute the action of the hamiltonian
       !> \param[in] input_vec Array to project
       !> \return Projected matrix

       use numeric_kinds, only: dp
       import :: davidson_parameters
       type(davidson_parameters) :: parameters
       real (dp), dimension(:,:), intent(in) :: input_vect
       real (dp), dimension(size(input_vect,1),size(input_vect,2)) :: output_vect
       
     end function fun_mtx_gemv
     
     function fun_stx_gemv(parameters, input_vect) result(output_vect)
       !> \brief Fucntion to compute the optional stx matrix on the fly
       !> \param[in] dimension of the arrays to compute the action of the hamiltonian
       !> \param[in] input_vec Array to project
       !> \return Projected matrix
       
       use numeric_kinds, only: dp
       import :: davidson_parameters
       type(davidson_parameters) :: parameters
       real (dp), dimension(:,:), intent(in) :: input_vect
       real (dp), dimension(size(input_vect,1), size(input_vect,2)) :: output_vect
       
       end function fun_stx_gemv

    end interface
    
    !local variables
    integer :: dim_mtx, dim_sub, max_dim, i, j, ier
    
    ! ! Basis of subspace of approximants
    real(dp), dimension(size(ritz_vectors, 1),1) :: guess, rs
    real(dp), dimension(:), allocatable :: diag_mtx, diag_stx
    real(dp), dimension(lowest):: errors
    
    ! ! Working arrays
    real(dp), dimension(:), allocatable :: eigenvalues_sub
    real(dp), dimension(:, :), allocatable :: correction, eigenvectors_sub, mtx_proj, stx_proj, V
    ! Arrays dimension
    type(davidson_parameters) :: parameters
    
    ! Iteration subpsace dimension
    dim_sub = lowest * 2
    
    ! maximum dimension of the basis for the subspace
    if (present(max_dim_sub)) then
       max_dim  = max_dim_sub
    else
       max_dim = lowest * 10
    endif
    
    ! dimension of the matrix
    dim_mtx = size(ritz_vectors, 1)
    parameters = davidson_parameters(nparm, nparm_max, lowest, max_dim_sub, mvec)

    ! extract the diagonals of the matrices

    write(6,'(''DAV: Compute diagonals of S and H'')')

    !! Diagonal of the arrays
    allocate(diag_mtx(dim_sub))
    allocate(diag_stx(dim_sub))
    diag_mtx = extract_diagonal_free(fun_mtx_gemv, parameters, dim_sub)
    diag_stx = extract_diagonal_free(fun_stx_gemv, parameters, dim_sub)
    
    ! 1. Variables initialization
    ! Select the initial ortogonal subspace based on lowest elements
    ! of the diagonal of the matrix
    V = ritz_vectors(1:nparm, 1:lowest) ! Initial orthonormal basis
    call lapack_qr(V)
    
    ! 2. Generate subspace matrix problem by projecting into V
    write(6,'(''DAV: Setup subspace problem'')')

    mtx_proj = lapack_matmul('T', 'N', V, fun_mtx_gemv(parameters, V))
    stx_proj = lapack_matmul('T', 'N', V, fun_stx_gemv(parameters, V))

    write(6,'(''DAV: Begin iterations'')')

    ! Outer loop block Davidson schema
    outer_loop: do i=1, max_iters
      write(6,'(''DAV: Davidson iteration: '', I10)') i

       ! IF IDTASK=0
       if (idtask.eq.0) then

       ! 3. compute the eigenvalues and their corresponding ritz_vectors
       ! for the projected matrix using lapack
       call check_deallocate_matrix(eigenvectors_sub)
       
       if (allocated(eigenvalues_sub)) then
          deallocate(eigenvalues_sub)
       end if
       
       allocate(eigenvalues_sub(size(mtx_proj, 1)))
       allocate(eigenvectors_sub(size(mtx_proj, 1), size(mtx_proj, 2)))
       
       write(6,'(''DAV: enter lapack_generalized_eigensolver'')') 
       call lapack_generalized_eigensolver(mtx_proj, eigenvalues_sub, eigenvectors_sub, stx_proj)
       write(6,'(''DAV: exit lapack_generalized_eigensolver'')') 
       
       ! 4. Check for convergence
       ritz_vectors = lapack_matmul('N', 'N', V, eigenvectors_sub(:, :lowest))
       do j=1,lowest
          guess = eigenvalues_sub(j) * fun_stx_gemv(parameters, reshape(ritz_vectors(:, j),(/dim_mtx,1/) ) )
          rs = fun_mtx_gemv(parameters, reshape(ritz_vectors(:, j), (/dim_mtx,1/))) - guess
          errors(j) = norm(reshape(rs,(/dim_mtx/)))
       end do
       write(6,'(''DAV: eigv'',1000f12.5)') (eigenvalues_sub(j),j=1,lowest)

      ! ENDIF IDTASK
       endif
       
       call MPI_BCAST(errors,lowest,MPI_REAL8,0,MPI_COMM_WORLD,ier)

       if (all(errors < tolerance)) then
          iters = i
          exit
       end if
       
       ! 5. Add the correction vectors to the current basis
       if (size(V, 2) <= max_dim) then
          
          ! append correction to the current basis
          call check_deallocate_matrix(correction)
          allocate(correction(size(ritz_vectors, 1), size(V, 2)))
          
          correction = compute_DPR_free(fun_mtx_gemv, fun_stx_gemv, parameters, lowest, V, &
               eigenvalues_sub, eigenvectors_sub, diag_mtx, diag_stx)
          
          ! IF IDTASK=0
          if(idtask.eq.0) then

          ! 6. Increase Basis size
          call concatenate(V, correction)
          
          ! 7. Orthogonalize basis
          call lapack_qr(V)

          ! ENDIF IDTASK
          endif

          ! 8. Update the the projection
          call update_projection_free(fun_mtx_gemv, parameters, lowest, V, mtx_proj)
          call update_projection_free(fun_stx_gemv, parameters, lowest, V, stx_proj)
          
       else
          
          ! 6. Otherwise reduce the basis of the subspace to the current correction
          V = lapack_matmul('N', 'N', V, eigenvectors_sub(:, :dim_sub))
          
          ! we refresh the projected matrices
          mtx_proj = lapack_matmul('T', 'N', V, fun_mtx_gemv(parameters, V))
          stx_proj = lapack_matmul('T', 'N', V, fun_stx_gemv(parameters, V))

       end if
       
    end do outer_loop
    
    !  8. Check convergence
    if (i > max_iters / dim_sub) then
       print *, "Warning: Algorithm did not converge!!"
    end if
    
    ! Select the lowest eigenvalues and their corresponding ritz_vectors
    ! They are sort in increasing order
    eigenvalues = eigenvalues_sub(:lowest)
    
    ! Free memory
    call check_deallocate_matrix(correction)
    deallocate(eigenvalues_sub, eigenvectors_sub, V, mtx_proj, diag_mtx, diag_stx)
    
    ! free optional matrix
    call check_deallocate_matrix(stx_proj)
    
  end subroutine generalized_eigensolver
  

  function compute_DPR_free(fun_mtx_gemv, fun_stx_gemv, parameters, lowest, V, eigenvalues, eigenvectors, &
    diag_mtx, diag_stx) result(correction)

    !> compute the correction vector using the DPR method for a matrix free diagonalization
    !> See correction_methods submodule for the implementations
    !> \param[in] fun_mtx: function to compute matrix
    !> \param[in] fun_stx: function to compute the matrix for the generalized case
    !> \param[in] V: Basis of the iteration subspace
    !> \param[in] eigenvalues: of the reduce problem
    !> \param[in] eigenvectors: of the reduce problem
    !> \return correction matrix
    integer, intent(in) :: lowest
    
    real(dp), dimension(:), intent(in) :: eigenvalues
    real(dp), dimension(:, :), intent(in) :: V, eigenvectors
    real(dp), dimension(:), intent(in) :: diag_mtx, diag_stx
    
    
    ! Function to compute the target matrix on the fly
    interface
     function fun_mtx_gemv(parameters, input_vect) result(output_vect)
       !> \brief Function to compute the action of the hamiltonian on the fly
       !> \param[in] dimension of the arrays to compute the action of the hamiltonian
       !> \param[in] input_vec Array to project
       !> \return Projected matrix

       use numeric_kinds, only: dp
       import :: davidson_parameters
       type(davidson_parameters) :: parameters
       real (dp), dimension(:,:), intent(in) :: input_vect
       real (dp), dimension(size(input_vect,1),size(input_vect,2)) :: output_vect
       
     end function fun_mtx_gemv
     
     function fun_stx_gemv(parameters, input_vect) result(output_vect)
       !> \brief Function to compute the action of the overlap on the fly
       !> \param[in] dimension of the arrays to compute the action of the hamiltonian
       !> \param[in] input_vec Array to project
       !> \return Projected matrix
       
       use numeric_kinds, only: dp
       import :: davidson_parameters
       type(davidson_parameters) :: parameters
       real (dp), dimension(:,:), intent(in) :: input_vect
       real (dp), dimension(size(input_vect,1),size(input_vect,2)) :: output_vect
       
       end function fun_stx_gemv      
       
    end interface
    
    ! local variables
    !real(dp), dimension(size(V, 1),1) :: vector
    real(dp), dimension(size(V, 1), size(V, 2)) :: correction, vectors
    real(dp), dimension(size(V, 2),size(V, 2)) :: diag_eigenvalues
    type(davidson_parameters) :: parameters
    integer :: ii, j
    integer :: dim
    
    dim = size(V,1)
    
    ! create a diagonal matrix of eigenvalues
    diag_eigenvalues = 0E0
    do ii =1, size(V,2)
       diag_eigenvalues(ii,ii) = eigenvalues(ii)
    end do
    
    ! create all the ritz vectors
    vectors = lapack_matmul('N','N', V, eigenvectors)
    
    ! initialize the correction vectors
    correction = fun_mtx_gemv(parameters, vectors) - lapack_matmul( &
         'N','N', fun_stx_gemv(parameters, vectors), diag_eigenvalues)
!   
! Pablo 
!
! A threshold in the denomitor of the correction is imposed.  
    do j=1, size(V, 2)
       do ii=1,size(correction,1)
          if ((eigenvalues(j) * diag_stx(ii)  - diag_mtx(ii)).gt.0.001_dp) then 
            correction(ii, j) = correction(ii, j) / (eigenvalues(j) * diag_stx(ii)  - diag_mtx(ii))
          else
            correction(ii,j) =0.0_dp
          endif
       end do
    end do
!
! End Pablo
!
!*! Original
!*!    do j=1, size(V, 2)
!*!       do ii=1,size(correction,1)
!*!          correction(ii, j) = correction(ii, j) / (eigenvalues(j) * diag_stx(ii)  - diag_mtx(ii))
!*!       end do
!*!    end do

  end function compute_DPR_free

 function extract_diagonal_free(fun_mtx_gemv, parameters, dim) result(out)
   !> \brief extract the diagonal of the matrix
   !> \param parameters: dimensions of the arrays
   implicit none
   
   type(davidson_parameters) :: parameters
   integer, intent(in) :: dim
   real(dp), dimension(dim) :: out
   
   interface
     function fun_mtx_gemv(parameters, input_vect) result(output_vect)
       !> \brief Function to compute the action of the hamiltonian on the fly
       !> \param[in] dimension of the arrays to compute the action of the hamiltonian
       !> \param[in] input_vec Array to project
       !> \return Projected matrix

       use numeric_kinds, only: dp
       import :: davidson_parameters
       type(davidson_parameters) :: parameters
       real (dp), dimension(:,:), intent(in) :: input_vect
       real (dp), dimension(size(input_vect,1),size(input_vect,2)) :: output_vect
       
     end function fun_mtx_gemv      
   end interface
   
   ! local variable
   integer :: ii
   real(dp), dimension(parameters%nparm, 1) :: tmp_array

   write(6,'(''DAV: extracting diagonal'')')
   
   do ii = 1, dim
      tmp_array = 0E0
      tmp_array(ii,1) = 1.0
      tmp_array = fun_mtx_gemv(parameters, tmp_array)
      out(ii) = tmp_array(ii,1)
   end do
   
 end function extract_diagonal_free

 
 subroutine update_projection_free(fun_mtx_gemv, parameters, lowest, V, mtx_proj)
   !> \brief update the projected matrices
   !> \param A: full matrix
   !> \param V: projector
   !> \param mtx_proj: projected matrix
   
   implicit none
   integer, intent(in) :: lowest
   real(dp), dimension(:, :), intent(in) :: V
   real(dp), dimension(:, :), intent(inout), allocatable :: mtx_proj
   real(dp), dimension(:, :), allocatable :: tmp_array
   type(davidson_parameters) :: parameters
   
   interface
     function fun_mtx_gemv(parameters, input_vect) result(output_vect)
       !> \brief Function to compute the action of the hamiltonian on the fly
       !> \param[in] dimension of the arrays to compute the action of the hamiltonian
       !> \param[in] input_vec Array to project
       !> \return Projected matrix

       use numeric_kinds, only: dp
       import :: davidson_parameters
       type(davidson_parameters) :: parameters
       real (dp), dimension(:,:), intent(in) :: input_vect
       real (dp), dimension(size(input_vect,1),size(input_vect,2)) :: output_vect

     end function fun_mtx_gemv

   end interface
   
   ! local variables
   integer :: nvec, old_dim
   
   ! dimension of the matrices
   nvec = size(V,2)
   old_dim = size(mtx_proj,1)    
!
! Pablo makes the full projetion:
!
    deallocate(mtx_proj)
    allocate(mtx_proj(nvec,nvec))
    allocate(tmp_array(size(V,1), size(V,2)))
    tmp_array=fun_mtx_gemv(parameters, V)
    mtx_proj=lapack_matmul('T', 'N', V, fun_mtx_gemv(parameters, V)) 
!
! End Pablo
!  
!*! Original
!*!   ! move to temporal array
!*!   allocate(tmp_array(nvec, nvec))
!*!   tmp_array(:old_dim, :old_dim) = mtx_proj
!*!   tmp_array(:,old_dim+1:) = lapack_matmul('T', 'N', V, fun_mtx_gemv(parameters, V(:, old_dim+1:)))
!*!   tmp_array( old_dim+1:,:old_dim ) = transpose(tmp_array(:old_dim, old_dim+1:))
!*!   
!*!   ! Move to new expanded matrix
!*!   deallocate(mtx_proj)
!*!   call move_alloc(tmp_array, mtx_proj)
   
  end subroutine update_projection_free  
    
  subroutine check_deallocate_matrix(mtx)
    !> deallocate a matrix if allocated
    real(dp), dimension(:, :), allocatable, intent(inout) ::  mtx
    
    if (allocated(mtx)) then
       deallocate(mtx)
    end if
    
  end subroutine check_deallocate_matrix  

end module davidson_free
