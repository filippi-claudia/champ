!> \namespace Davidson eigensolver
!> \author Felipe Zapata
!> The current implementation uses a general  davidson algorithm, meaning
!> that it compute all the eigenvalues simultaneusly using a variable size block approach.
!> The family of Davidson algorithm only differ in the way that the correction
!> vector is computed.
!> Computed pairs of eigenvalues/eigenvectors are deflated using algorithm
!> described at: https://doi.org/10.1023/A:101919970


module davidson_dense
  !> Submodule containing the implementation of the Davidson diagonalization method
  !> for dense matrices
  use numeric_kinds, only: dp
  use lapack_wrapper, only: lapack_generalized_eigensolver, lapack_matmul, lapack_matrix_vector, &
       lapack_qr, lapack_solver, lapack_sort
  use array_utils, only: concatenate, diagonal, eye, generate_preconditioner, norm

  implicit none
  
  !> \private
  private
  !> \public
  public :: generalized_eigensolver_dense

  interface
     module function compute_correction_generalized_dense(mtx, V, eigenvalues, eigenvectors, method, stx) &
          result(correction)
       !> compute the correction vector using a given `method` for the Davidson algorithm
       !> See correction_methods submodule for the implementations
       !> \param[in] mtx: Original matrix
       !> \param[in] stx: Matrix to compute the general eigenvalue problem
       !> \param[in] V: Basis of the iteration subspace
       !> \param[in] eigenvalues: of the reduce problem
       !> \param[in] eigenvectors: of the reduce problem
       !> \param[in] method: name of the method to compute the correction
       
       real(dp), dimension(:), intent(in) :: eigenvalues
       real(dp), dimension(:, :), intent(in) :: mtx, V, eigenvectors
       real(dp), dimension(:, :), intent(in), optional :: stx
       character(len=*), optional, intent(in) :: method
       real(dp), dimension(size(mtx, 1), size(V, 2)) :: correction
       
     end function compute_correction_generalized_dense

  end interface
  
contains

    subroutine generalized_eigensolver_dense(mtx, eigenvalues, ritz_vectors, lowest, method, max_iters, &
        tolerance, iters, max_dim_sub, stx)
     !> Implementation storing in memory the initial densed matrix mtx.
      
    !> \param[in] mtx: Matrix to diagonalize
    !> \param[in, opt] Optional matrix to solve the general eigenvalue problem:
    !> \f$ mtx \lambda = V stx \lambda \f$
    !> \param[out] eigenvalues Computed eigenvalues
    !> \param[out] ritz_vectors approximation to the eigenvectors
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
    ! input/output variable
    integer, intent(in) :: lowest
    real(dp), dimension(:, :), intent(in) :: mtx
    real(dp), dimension(:, :), intent(in), optional :: stx
    real(dp), dimension(lowest), intent(out) :: eigenvalues
    real(dp), dimension(:, :), intent(out) :: ritz_vectors
    integer, intent(in) :: max_iters
    integer, intent(in), optional :: max_dim_sub
    real(dp), intent(in) :: tolerance
    character(len=*), intent(in) :: method
    integer, intent(out) :: iters
    
    !local variables
    integer :: i, j, dim_sub, max_dim
    integer :: n_converged ! Number of converged eigenvalue/eigenvector pairs
    
    ! Basis of subspace of approximants
    real(dp), dimension(size(mtx, 1)) :: guess, rs
    real(dp), dimension(lowest):: errors

    ! Working arrays
    real(dp), dimension(:), allocatable :: eigenvalues_sub
    real(dp), dimension(:, :), allocatable :: correction, eigenvectors_sub, mtx_proj, stx_proj, V

    ! Diagonal matrix
    real(dp), dimension(size(mtx, 1)) :: d
    
    ! generalize problem
    logical :: gev 

    ! indices of the eigenvalues/eigenvectors pair that have not converged
    logical, dimension(lowest) :: has_converged

    ! Iteration subpsace dimension
    dim_sub = lowest * 2

    ! Initial number of converged eigenvalue/eigenvector pairs
    n_converged = 0
    has_converged = .False.
    
    ! maximum dimension of the basis for the subspace
    if (present(max_dim_sub)) then
       max_dim  = max_dim_sub
    else
       max_dim = lowest * 10
    endif

    ! generalied problem
    gev = present(stx)

    ! 1. Variables initialization
    ! Select the initial ortogonal subspace based on lowest elements
    ! of the diagonal of the matrix
    d = diagonal(mtx)
    V = generate_preconditioner(d, dim_sub)

   ! 2. Generate subpace matrix problem by projecting into V
   mtx_proj = lapack_matmul('T', 'N', V, lapack_matmul('N', 'N', mtx, V))

   if(gev) then
    stx_proj = lapack_matmul('T', 'N', V, lapack_matmul('N', 'N', stx, V))
   end if

    ! ! Outer loop block Davidson schema
    outer_loop: do i=1, max_iters

       ! 3. compute the eigenvalues and their corresponding ritz_vectors
       ! for the projected matrix using lapack
       call check_deallocate_matrix(eigenvectors_sub)

       if (allocated(eigenvalues_sub)) then
          deallocate(eigenvalues_sub)
       end if

       allocate(eigenvalues_sub(size(mtx_proj, 1)))
       allocate(eigenvectors_sub(size(mtx_proj, 1), size(mtx_proj, 2)))


       if (gev) then
        call lapack_generalized_eigensolver(mtx_proj, eigenvalues_sub, eigenvectors_sub, stx_proj)
       else
        call lapack_generalized_eigensolver(mtx_proj, eigenvalues_sub, eigenvectors_sub)
       end if

       ! 4. Check for convergence
       ritz_vectors = lapack_matmul('N', 'N', V, eigenvectors_sub(:, :lowest))
       do j=1,lowest
          if(gev) then
            guess = eigenvalues_sub(j) * lapack_matrix_vector('N',stx,ritz_vectors(:, j))
          else
            guess = eigenvalues_sub(j) * ritz_vectors(:, j)
          end if
          rs = lapack_matrix_vector('N', mtx, ritz_vectors(:, j)) - guess
          errors(j) = norm(rs)
          ! Check which eigenvalues has converged
          if (errors(j) < tolerance) then
             has_converged(j) = .true.
          end if
       end do
       
       ! Count converged pairs of eigenvalues/eigenvectors
       n_converged = n_converged + count(errors < tolerance)
       
       if (all(has_converged)) then
          iters = i
          exit
       end if
       
       ! 5. Add the correction vectors to the current basis
       if (size(V, 2) <= max_dim) then

          ! append correction to the current basis
          call check_deallocate_matrix(correction)
          allocate(correction(size(mtx, 1), size(V, 2)))

          if(gev) then
            correction = compute_correction_generalized_dense(mtx, V, eigenvalues_sub, eigenvectors_sub, method, stx)
          else
            correction = compute_correction_generalized_dense(mtx, V, eigenvalues_sub, eigenvectors_sub, method)
          end if


          ! 6. Increase Basis size
          call concatenate(V, correction)
       
          ! 7. Orthogonalize basis
          call lapack_qr(V)
          
       else

          ! 6. Otherwise reduce the basis of the subspace to the current correction
          V = lapack_matmul('N', 'N', V, eigenvectors_sub(:, :dim_sub))
       end if

       ! we refresh the projected matrices
       mtx_proj = lapack_matmul('T', 'N', V, lapack_matmul('N', 'N', mtx, V))
       
       if(gev) then
          stx_proj = lapack_matmul('T', 'N', V, lapack_matmul('N', 'N', stx, V))
       end if
      
    end do outer_loop

    !  8. Check convergence
    if (i > max_iters) then
       iters = i
       print *, "Warning: Algorithm did not converge!!"
    end if
    
    ! Select the lowest eigenvalues and their corresponding ritz_vectors
    ! They are sort in increasing order
    eigenvalues = eigenvalues_sub(:lowest)
    
    ! Free memory
    call check_deallocate_matrix(correction)
    deallocate(eigenvalues_sub, eigenvectors_sub, V, mtx_proj)

    ! free optional matrix
    if (gev) then
       call check_deallocate_matrix(stx_proj)
    endif
    
  end subroutine generalized_eigensolver_dense

  subroutine update_projection_dense(A, V, A_proj)
    !> \brief update the projected matrices
    !> \param A: full matrix
    !> \param V: projector
    !> \param A_proj: projected matrix

    implicit none
    real(dp), dimension(:, :), intent(in) :: A
    real(dp), dimension(:, :), intent(in) :: V
    real(dp), dimension(:, :), intent(inout), allocatable :: A_proj
    real(dp), dimension(:, :), allocatable :: tmp_array

    ! local variables
    integer :: nvec, old_dim

    ! dimension of the matrices
    nvec = size(V,2)
    old_dim = size(A_proj,1)    

    ! move to temporal array
    allocate(tmp_array(nvec, nvec))
    tmp_array(:old_dim, :old_dim) = A_proj
    tmp_array(:,old_dim+1:) = lapack_matmul('T', 'N', V, lapack_matmul('N', 'N', A, V(:, old_dim+1:)))
    tmp_array( old_dim+1:,:old_dim ) = transpose(tmp_array(:old_dim, old_dim+1:))
    
    ! Move to new expanded matrix
    deallocate(A_proj)
    call move_alloc(tmp_array, A_proj)
 
  end subroutine update_projection_dense

  subroutine check_deallocate_matrix(mtx)
    !> deallocate a matrix if allocated
    real(dp), dimension(:, :), allocatable, intent(inout) ::  mtx
    
    if (allocated(mtx)) then
       deallocate(mtx)
    end if
    
  end subroutine check_deallocate_matrix  
  
end module davidson_dense

submodule (davidson_dense) correction_methods_generalized_dense
  !> submodule containing the implementations of different kind
  !> algorithms to compute the correction vectors for the Davidson's diagonalization

  implicit none
  
contains

  module function compute_correction_generalized_dense(mtx, V, eigenvalues, eigenvectors, method, stx) &
       result(correction)
    !> see interface in davidson module
    real(dp), dimension(:), intent(in) :: eigenvalues
    real(dp), dimension(:, :), intent(in) :: mtx, V, eigenvectors
    real(dp), dimension(:, :), intent(in), optional :: stx
    character(len=*), optional,intent(in) :: method
    logical :: gev 

    ! local variables
    character(len=10) :: opt 
    real(dp), dimension(size(mtx, 1), size(V, 2)) :: correction

    !check optional arguments
    gev = present(stx)
    opt="DPR"
    if (present(method)) opt=trim(method)
    
    select case (method)
    case ("DPR")
      if(gev) then
       correction = compute_DPR_generalized_dense(mtx, V, eigenvalues, eigenvectors, stx)
      else
        correction = compute_DPR_generalized_dense(mtx, V, eigenvalues, eigenvectors)
      end if
    case ("GJD")
      if(gev) then
       correction = compute_GJD_generalized_dense(mtx, V, eigenvalues, eigenvectors, stx)
      else
        correction = compute_GJD_generalized_dense(mtx, V, eigenvalues, eigenvectors)
      end if
    end select
    
  end function compute_correction_generalized_dense

  function compute_DPR_generalized_dense(mtx, V, eigenvalues, eigenvectors, stx) result(correction)
    !> compute Diagonal-Preconditioned-Residue (DPR) correction
    real(dp), dimension(:), intent(in) :: eigenvalues
    real(dp), dimension(:, :), intent(in) :: mtx, V, eigenvectors
    real(dp), dimension(:, :), intent(in), optional ::  stx 
    real(dp), dimension(size(mtx, 1), size(V, 2)) :: correction
    
    ! local variables
    integer :: ii,j, m
    real(dp), dimension(size(mtx, 1), size(mtx, 2)) :: diag, arr
    real(dp), dimension(size(mtx, 1)) :: vec
    logical :: gev
    
    ! shape of matrix
    m = size(mtx, 1)
    gev = (present(stx))

    do j=1, size(V, 2)
       if(gev) then
          diag = eigenvalues(j) * stx
       else
          diag = eye(m , m, eigenvalues(j))
       end if
       arr = mtx - diag
       vec = lapack_matrix_vector('N', V, eigenvectors(:, j))
      
       correction(:, j) = lapack_matrix_vector('N', arr, vec) 

       do ii=1,size(correction,1)
          if (gev) then
             correction(ii, j) = correction(ii, j) / (eigenvalues(j) * stx(ii,ii) - mtx(ii, ii))
           else
              correction(ii, j) = correction(ii, j) / (eigenvalues(j)  - mtx(ii, ii))
           endif
        end do
    end do

  end function compute_DPR_generalized_dense

  function compute_GJD_generalized_dense(mtx, V, eigenvalues, eigenvectors, stx) result(correction)
    !> Compute the Generalized Jacobi Davidson (GJD) correction
    
    real(dp), dimension(:), intent(in) :: eigenvalues
    real(dp), dimension(:, :), intent(in) :: mtx, V, eigenvectors
    real(dp), dimension(:, :), intent(in), optional :: stx
    real(dp), dimension(size(mtx, 1), size(V, 2)) :: correction

    ! local variables
    integer :: k, m
    logical :: gev
    real(dp), dimension(size(mtx, 1), 1) :: rs
    real(dp), dimension(size(mtx, 1), size(V, 2)) :: ritz_vectors
    real(dp), dimension(size(mtx, 1), size(mtx, 2)) :: arr, xs, ys
    real(dp), dimension(size(mtx, 1), 1) :: brr

    ! Diagonal matrix
    m = size(mtx, 1)
    ritz_vectors = lapack_matmul('N', 'N', V, eigenvectors)

    gev = present(stx)

    do k=1, size(V, 2)
       rs(:, 1) = ritz_vectors(:, k)
       xs = eye(m, m) - lapack_matmul('N', 'T', rs, rs)
       if(gev) then
        ys = mtx - eigenvalues(k)*stx
       else
         ys = substract_from_diagonal(mtx, eigenvalues(k))
       end if
       arr = lapack_matmul('N', 'N', xs, lapack_matmul('N', 'N', ys, xs))
       brr = -rs
       
       call lapack_solver(arr, brr)
       correction(:, k) = brr(:, 1)
    end do
    
  end function compute_GJD_generalized_dense

  function substract_from_diagonal(mtx, alpha) result(arr)
    !> susbstract an scalar from the diagonal of a matrix
    !> \param mtx: square matrix
    !> \param alpha: scalar to substract
    real(dp), dimension(:, :), intent(in) :: mtx
    real(dp), dimension(size(mtx, 1), size(mtx, 2)) :: arr
    real(dp), intent(in) :: alpha
    integer :: i

    arr = mtx
    do i=1,size(mtx, 1)
       arr(i, i) = arr(i, i) - alpha
    end do
    
  end function substract_from_diagonal
  
end submodule correction_methods_generalized_dense
