!> \namespace Davidson eigensolver
!> \author Felipe Zapata
!> The current implementation uses a general  davidson algorithm, meaning
!> that it compute all the eigenvalues simultaneusly using a variable size block approach.
!> The family of Davidson algorithm only differ in the way that the correction
!> vector is computed.
!> Computed pairs of eigenvalues/eigenvectors are deflated using algorithm
!> described at: https://doi.org/10.1023/A:101919970
!>
!> Revised by P. Lopez-Tarifa@NLeSC(2019)
module davidson_free

  use numeric_kinds, only: dp
  use lapack_wrapper, only: lapack_generalized_eigensolver, lapack_matmul, lapack_matrix_vector, &
       lapack_qr, lapack_solver
  use array_utils, only: concatenate, generate_preconditioner, norm, write_matrix, write_vector, & 
                        eye, check_deallocate_matrices, check_deallocate_matrix
  implicit none

  type davidson_parameters
     INTEGER :: nparm
     INTEGER :: nparm_max
     INTEGER :: lowest
     INTEGER :: nvecx 
     INTEGER :: basis_size 
  end type davidson_parameters

  !> \private
  private
  !> \public
  public :: generalized_eigensolver_free, davidson_parameters, die
  
contains

  subroutine generalized_eigensolver_free(fun_mtx_gemv, eigenvalues, ritz_vectors, nparm, nparm_max, &
       lowest, nvecx, method, max_iters, tolerance, iters, fun_stx_gemv, nproc, idtask)
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
    integer, intent(in) :: nparm, nparm_max, nvecx, lowest, nproc, idtask
    real(dp), dimension(lowest), intent(out) :: eigenvalues
    real(dp), dimension(:, :), allocatable, intent(out) :: ritz_vectors
    integer, intent(in) :: max_iters
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
!    
    !local variables
    integer :: dim_sub, max_size_basis, i, j, ier
!    
    ! ! Basis of subspace of approximants
    real(dp), dimension(:), allocatable :: diag_mtx, diag_stx
    real(dp), dimension(:,:), allocatable :: guess, rs 
    real(dp), dimension(lowest):: errors
    ! generalize problem
    
    ! ! Working arrays
    real(dp), dimension(:), allocatable :: eigenvalues_sub
    real(dp), dimension(:,:), allocatable :: lambda       ! eigenvalues_sub in a diagonal matrix 
    real(dp), dimension(:, :), allocatable :: correction, eigenvectors_sub, mtx_proj, stx_proj, V
    real(dp), dimension(:, :), allocatable :: mtxV, stxV 
    real(dp), dimension(nparm, 1) :: xs, gs
    real(dp), dimension(:), allocatable :: d 
    ! Arrays dimension
    type(davidson_parameters) :: parameters
    
    ! Iteration subpsace dimension
    dim_sub = lowest  * 2

    ! Lapack qr safety check 
    if (nvecx > nparm) then 
      if( idtask == 1) call die('DAV: nvecx > nparm, increase nparm or decrese lin_nvecx')
    endif

    ! dimension of the matrix
    parameters = davidson_parameters(nparm, nparm_max, lowest, nvecx, dim_sub) 

    ! 1. Variables initialization
    ! extract the diagonals of the matrices

    write(6,'(''DAV: Compute diagonals of S and H'')')

    !! Diagonal of the arrays
    allocate(diag_mtx(parameters%nparm))
    allocate(diag_stx(parameters%nparm))
    allocate(d(dim_sub))
    allocate(guess(parameters%nparm,1), rs(parameters%nparm,1))

    diag_mtx= extract_diagonal_free(fun_mtx_gemv, parameters, parameters%nparm)
    diag_stx= extract_diagonal_free(fun_stx_gemv, parameters, parameters%nparm)

    d= 0.0_dp
    do i= 1, dim_sub 
      xs= 0.0_dp 
      xs( i, 1)= 1.0_dp 
      gs= fun_mtx_gemv( parameters, xs)
      d( i)= gs( i,1)
    enddo

    if (nproc > 1) call MPI_BCAST( d, dim_sub, MPI_REAL8, 0, MPI_COMM_WORLD, ier)
 
    ! 2.  Select the initial ortogonal subspace based on lowest elements
    !     of the diagonal of the matrix

    V= generate_preconditioner( d( 1: dim_sub), dim_sub, nparm) ! Initial orthonormal basis
    
    if( idtask== 0) write(6,'(''DAV: Setup subspace problem'')')

    ! 3. Outer loop block Davidson schema

    outer_loop: do i= 1, max_iters

      if( idtask== 0) write(6,'(''DAV: Davidson iteration: '', I10)') i

      !! Array deallocation/allocation
      call check_deallocate_matrices(mtx_proj, stx_proj, lambda, eigenvectors_sub, ritz_vectors, mtxV, stxV)
      if( allocated( eigenvalues_sub)) then
          deallocate( eigenvalues_sub)
      end if
      allocate( mtxV( parameters%nparm, parameters%basis_size), stxV(parameters%nparm, parameters%basis_size)) 
      allocate( eigenvalues_sub( parameters%basis_size)) 
      allocate( lambda( parameters%basis_size, parameters%basis_size)) 
      allocate( eigenvectors_sub( parameters%basis_size, parameters%basis_size)) 

      ! 4. Projection of H and S matrices

      mtxV= fun_mtx_gemv( parameters, V)
      stxV= fun_stx_gemv( parameters, V)

      if (nproc > 1) then
        call MPI_BCAST( mtxV, parameters%nparm* parameters%basis_size, MPI_REAL8, 0, MPI_COMM_WORLD, ier)
        call MPI_BCAST( stxV, parameters%nparm* parameters%basis_size, MPI_REAL8, 0, MPI_COMM_WORLD, ier)
      endif

      mtx_proj= lapack_matmul( 'T', 'N', V, mtxV)
      stx_proj= lapack_matmul( 'T', 'N', V, stxV)

      !! IF IDTASK=0
      if( idtask.eq.0) then

        ! 5.Compute the eigenvalues and their corresponding ritz_vectors
        ! for the projected matrix using lapack
        
        write( 6,'(''DAV: enter lapack_generalized_eigensolver'')') 
        call lapack_generalized_eigensolver( mtx_proj, eigenvalues_sub, eigenvectors_sub, stx_proj)
        write( 6,'(''DAV: exit lapack_generalized_eigensolver'')') 
        write( 6,'(''DAV: eigv'',1000f12.5)') (eigenvalues_sub(j),j=1,lowest)
    
        ! 6. Construction of lambda matrix (a squared one with eigenvalues_sub in the diagonal)

        lambda= eye( parameters%basis_size, parameters%basis_size) 
        do j= 1, parameters%basis_size 
          lambda( j, j)= eigenvalues_sub( j)
        enddo
 
        ! 7. Residue calculation  

        rs = lapack_matmul('N', 'N', stxV, eigenvectors_sub)
        guess =  lapack_matmul('N', 'N', rs, lambda)  
        deallocate(rs)
        rs =  lapack_matmul('N', 'N', mtxV, eigenvectors_sub) - guess 
        do j=1,lowest
          errors(j) = norm(reshape(rs(:,j),(/parameters%nparm/)))
        end do

      !! ENDIF IDTASK
      endif

      if (nproc > 1) then
        call MPI_BCAST(errors,lowest,MPI_REAL8,0,MPI_COMM_WORLD,ier)
        call MPI_BCAST(eigenvalues_sub, parameters%basis_size, MPI_REAL8, 0, MPI_COMM_WORLD,ier)
        call MPI_BCAST(eigenvectors_sub, parameters%basis_size*parameters%basis_size, & 
                      MPI_REAL8, 0, MPI_COMM_WORLD,ier)
      endif       

      ritz_vectors = lapack_matmul('N', 'N', V, eigenvectors_sub)

      ! 8. Check for convergence

      if (all(errors < tolerance)) then
        iters = i
        write( 6, '(''DAV: roots are converged'')') 
        eigenvalues = eigenvalues_sub(:lowest)
        exit outer_loop
      end if

      ! 7. Add the correction vectors to the current basis

      if(( parameters%basis_size <= nvecx).and.( 2*parameters%basis_size< nparm)) then

        ! append correction to the current basis
        call check_deallocate_matrix(correction)
        allocate(correction(size(ritz_vectors, 1), size(V, 2)))

        select case( method)
        case( "DPR")
          if( idtask == 1)  write(6,'(''DAV: Diagonal-Preconditioned-Residue (DPR)'')')
          correction= compute_DPR_free( mtxV, stxV, parameters, eigenvalues_sub,             &
                                        eigenvectors_sub, diag_mtx, diag_stx)
        case( "GJD")
          if( idtask == 1)  write(6,'(''DAV: Diagonal-Preconditioned-Residue (GJD)'')')
          correction= compute_GJD_free( parameters, ritz_vectors, rs, eigenvectors_sub,      &
                                        eigenvalues_sub)
        end select

        ! 8. Increase Basis size

        call concatenate( V, correction)
           
        ! IF IDTASK=0
        if( idtask .eq. 0) then

          ! 9. Orthogonalize basis

          call lapack_qr( V)

        ! ENDIF IDTASK
        endif

        if (nproc > 1) then
          call MPI_BCAST( V, size( V, 1)* size( V, 2), MPI_REAL8, 0, MPI_COMM_WORLD, ier)
        endif 

      else

        ! 10. Otherwise reduce the basis of the subspace to the current correction
        V = lapack_matmul('N', 'N', V, eigenvectors_sub(:, :dim_sub))

      end if
    
    ! Update basis size
      parameters%basis_size = size( V, 2) 

    end do outer_loop

    !  11. Check convergence
    if (i > max_iters ) then
       iters = i
       print *, "Warning: Algorithm did not converge!!"
    end if

    ! Select the lowest eigenvalues and their corresponding ritz_vectors
    ! They are sort in increasing order
    eigenvalues = eigenvalues_sub(:lowest)
    
    ! Free memory
    call check_deallocate_matrix(correction)
    deallocate(eigenvalues_sub, eigenvectors_sub, mtx_proj, diag_mtx, diag_stx, d)
    deallocate(V, mtxV, stxV, guess, rs, lambda)
    
    ! free optional matrix
    call check_deallocate_matrix(stx_proj)
    
  end subroutine generalized_eigensolver_free
!  
  subroutine die(msg)
  !> Subroutine that dies the calculation raising an errror message
  !
  character msg*(*)
  integer ierr
  include 'mpif.h'

  write(6,'(''Fatal error: '',a)') msg
  call mpi_abort(MPI_COMM_WORLD,0,ierr)

  end subroutine

  function compute_DPR_free(mtxV, stxV, parameters, eigenvalues, eigenvectors,    &
                            diag_mtx, diag_stx) result(correction)

    !> compute the correction vector using the DPR method for a matrix free diagonalization
    !> See correction_methods submodule for the implementations
    !> \param[in] fun_mtx: function to compute matrix
    !> \param[in] fun_stx: function to compute the matrix for the generalized case
    !> \param[in] V: Basis of the iteration subspace
    !> \param[in] eigenvalues: of the reduce problem
    !> \param[in] eigenvectors: of the reduce problem
    !> \return correction matrix
    !
    use array_utils, only: eye
    !
    
    real(dp), dimension(:), intent(in) :: eigenvalues
    real(dp), dimension(:, :), intent(in) ::  eigenvectors, mtxV, stxV
    real(dp), dimension(:), intent(in) :: diag_mtx, diag_stx
    
    ! local variables
    !real(dp), dimension(size(V, 1),1) :: vector
    real(dp), dimension(size(mtxV, 1), size(mtxV, 2)) :: correction
    real(dp), dimension(size(mtxV, 1), size(mtxV, 2)) :: proj_mtx, proj_stx
    real(dp), dimension(size(mtxV, 1),size(mtxV, 1)) :: diag
    type(davidson_parameters) :: parameters
    integer :: ii, j
    integer :: m
    
    ! calculate the correction vectors
    m= size( mtxV, 1) 

    ! computed the projected matrices
    proj_mtx = lapack_matmul('N', 'N', mtxV, eigenvectors)
    proj_stx = lapack_matmul('N', 'N', stxV, eigenvectors)
    diag = 0.0_dp

    do j = 1, size( mtxV, 2)
     diag= eye( m , m, eigenvalues( j))
     correction( :, j)= proj_mtx( :, j)- lapack_matrix_vector( 'N', diag, proj_stx( :, j))

     do ii= 1, size( correction, 1)
       correction( ii, j)= correction( ii, j)/( eigenvalues( j)* diag_stx( ii)- diag_mtx( ii))
     end do

    end do

  end function compute_DPR_free

  function compute_GJD_free( parameters, ritz_vectors, residues, eigenvectors, & 
             eigenvalues) result( correction)

    !> Compute the correction vector using the GJD method for a matrix free
    !> diagonalization. We follow the notation of:
    !> I. Sabzevari, A. Mahajan and S. Sharma,  arXiv:1908.04423 (2019)
    !>
    !> \param[in] ritz_vectors: ritz_vectors.
    !> \param[in] residues: residue vectors.
    !> \param[in] parameter: davidson_parameters type.
    !> \param[in] eigenvectors. 
    !> \param[in] eigenvalues. 
    !> \return correction matrix

    use array_utils, only: eye

    type( davidson_parameters)               :: parameters
    real( dp), dimension( :, :), intent( in) :: ritz_vectors
    real( dp), dimension( :, :), intent( in) :: residues
    real( dp), dimension( :, :), intent( in) :: eigenvectors
    real( dp), dimension( :),    intent( in) :: eigenvalues 
    !
    ! local variables
    !
    real( dp), dimension( parameters%nparm, parameters%basis_size) :: correction
    integer :: k, m
    logical :: gev
    real( dp), dimension( :, :), allocatable   ::  F
    real( dp), dimension( parameters%nparm, 1) :: brr

    do k= 1, parameters%basis_size 

      F= fun_F_matrix( ritz_vectors, parameters, k, eigenvalues( k))  
      brr( :, 1) = -residues(:,k)
      call lapack_solver( F, brr)
      correction( :, k)= brr( :, 1)

    end do

     ! Deallocate
     deallocate( F)

  end function compute_GJD_free

  function fun_F_matrix( ritz_vectors, parameters, eigen_index, eigenvalue) &
           result( F_matrix)
    !> \brief Function that computes the F matrix: 
    !> F= ubut*( A- theta* B)* uubt
    !> in a pseudo-free way for a given engenvalue. 
    !> 
    !> ritz_vectors( in) :: ritz_vectors.  
    !> parameters( in)   :: array_sizes  
    !> eigen_index( in)  :: index of the passing eingenvalue.
    !> eigenvalue( in)   :: eigen_index eigenvalue.  

    use array_utils, only: eye

    real( dp), dimension( :, :), intent( in) :: ritz_vectors 
    type( davidson_parameters) :: parameters
    integer   :: eigen_index 
    real( dp) :: eigenvalue 

    interface

     function fun_mtx_gemv( parameters, input_vect) result( output_vect)
       !> \brief Function to compute the action of the hamiltonian on the fly
       !> \param[in] dimension of the arrays to compute the action of the
       !             hamiltonian
       !> \param[in] input_vec Array to project
       !> \return Projected matrix
       use numeric_kinds, only: dp
       import                                   :: davidson_parameters
       type( davidson_parameters)               :: parameters
       real( dp), dimension( :, :), intent( in) :: input_vect
       real( dp), dimension( size( input_vect, 1), size( input_vect, 2)) :: output_vect
     end function fun_mtx_gemv

     function fun_stx_gemv( parameters, input_vect) result( output_vect)
       !> \brief Fucntion to compute the optional stx matrix on the fly
       !> \param[in] dimension of the arrays to compute the action of the
       !             hamiltonian
       !> \param[in] input_vec Array to project
       !> \return Projected matrix
       use numeric_kinds, only: dp
       import                                   :: davidson_parameters
       type( davidson_parameters)               :: parameters
       real( dp), dimension( :, :), intent( in) :: input_vect
       real( dp), dimension( size( input_vect, 1), size( input_vect, 2)) :: output_vect
     end function fun_stx_gemv

    end interface
 
    real( dp), dimension( parameters%nparm, parameters%nparm) :: F_matrix, lambda
    real( dp), dimension( parameters%nparm, 1) :: ritz_tmp
    real( dp), dimension( :, :), allocatable :: ys 
    real( dp), dimension( parameters%nparm, parameters%nparm) :: ubut, uubt

    ritz_tmp( :, 1)= ritz_vectors( :, eigen_index)

    lambda= eye( parameters%nparm, parameters%nparm, eigenvalue)

    ubut= eye( parameters%nparm, parameters%nparm)- &
          lapack_matmul( 'N', 'T', fun_stx_gemv( parameters, ritz_tmp), ritz_tmp)

    uubt= eye( parameters%nparm, parameters%nparm)- &
          lapack_matmul( 'N', 'T', ritz_tmp, fun_stx_gemv( parameters, ritz_tmp)) 

    ys = lapack_matmul( 'N', 'N', lambda, fun_stx_gemv( parameters, uubt)) 

    F_matrix= lapack_matmul( 'N', 'N', ubut, fun_mtx_gemv( parameters, uubt)) - &
              lapack_matmul( 'N', 'N', ubut, ys) 

  end function fun_F_matrix

  function extract_diagonal_free( fun_mtx_gemv, parameters, dim) result( out)
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

!    write(6,'(''DAV: extracting diagonal'')')
    
    do ii = 1, dim
       tmp_array = 0E0
       tmp_array(ii,1) = 1.0
       tmp_array = fun_mtx_gemv(parameters, tmp_array)
       out(ii) = tmp_array(ii,1)
    end do
    
  end function extract_diagonal_free
 
end module davidson_free

module davidson_dense
  !> Submodule containing the implementation of the Davidson diagonalization method
  !> for dense matrices
  use numeric_kinds, only: dp
  use lapack_wrapper, only: lapack_generalized_eigensolver, lapack_matmul, lapack_matrix_vector, &
       lapack_qr, lapack_solver, lapack_sort
  use array_utils, only: concatenate, diagonal, eye, generate_preconditioner, norm, check_deallocate_matrices, &
                         check_deallocate_matrix
  use davidson_free, only: die

  implicit none

  type davidson_parameters
     INTEGER :: nparm
     INTEGER :: nparm_max
     INTEGER :: lowest
     INTEGER :: nvecx
  end type davidson_parameters
  
  !> \private
  private
  !> \public
  public :: generalized_eigensolver_dense

  interface
     module function compute_correction_generalized_dense(mtx, V, eigenvalues, eigenvectors, method, rs, ritz_vectors, stx) &
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
       real(dp), dimension(:, :), intent(in), optional :: rs 
       real(dp), dimension(:, :), intent(in), optional :: ritz_vectors 
       real(dp), dimension(:, :), intent(in), optional :: stx
       character(len=*), optional, intent(in) :: method
       real(dp), dimension(size(mtx, 1), size(V, 2)) :: correction
       
     end function compute_correction_generalized_dense

  end interface
  
contains

  subroutine generalized_eigensolver_dense(mtx, eigenvalues, ritz_vectors, nparm, nparm_max, &
             lowest, nvecx, method, max_iters, tolerance, iters, stx, nproc, idtask)
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

    include 'mpif.h'

    ! input/output variable
    integer, intent(in) :: nparm, nparm_max, nvecx, lowest, nproc, idtask 
    real(dp), dimension(:, :), intent(in) :: mtx
    real(dp), dimension(:, :), intent(in), optional :: stx
    real(dp), dimension(lowest), intent(out) :: eigenvalues
    real(dp), dimension(:, :), allocatable, intent(out) :: ritz_vectors
    integer, intent(in) :: max_iters
    real(dp), intent(in) :: tolerance
    character(len=*), intent(in) :: method
    integer, intent(out) :: iters
    
    !local variables
    integer :: i, j, dim_sub, max_dim, ier
    
    ! Basis of subspace of approximants
    real(dp), dimension( :, :), allocatable:: guess, rs  
    real(dp), dimension(lowest):: errors

    ! Working arrays
    real(dp), dimension(:), allocatable :: eigenvalues_sub
    real(dp), dimension(:,:), allocatable :: lambda       ! eigenvalues_sub in a diagonal matrix
    real(dp), dimension(:, :), allocatable :: correction, eigenvectors_sub, mtx_proj, stx_proj, V
    ! Arrays dimension
    type(davidson_parameters) :: parameters

    ! Diagonal matrix
    real(dp), dimension(size(mtx, 1)) :: d
    
    ! generalize problem
    logical :: gev 

    ! Iteration subpsace dimension
    dim_sub = lowest  * 2

    ! Lapack qr safety check
    if (nvecx > nparm) then
      if( idtask == 1) call die('DAV: nvecx > nparm, increase nparm or decrese lin_nvecx')
    endif

    ! dimension of the matrix
    parameters = davidson_parameters(nparm, nparm_max, lowest, nvecx)

    
    ! Generalied problem
    gev = present(stx)
!
    ! 1. Variables initialization
    ! extract the diagonals of the matrices

    write(6,'(''DAV: Compute diagonals of H'')')  
    d = diagonal(mtx)

    if (nproc > 1) call MPI_BCAST( d, dim_sub, MPI_REAL8, 0, MPI_COMM_WORLD, ier)

    ! 2. Select the initial ortogonal subspace based on lowest elements
    !    of the diagonal of the matrix

    V = generate_preconditioner(d(1:dim_sub), dim_sub, size(mtx,1))

    write(6,'(''DAV: Setup subspace problem'')')

!   
   ! 3. Outer loop block Davidson schema
    outer_loop: do i= 1, max_iters
      if( idtask == 0) write(6,'(''DAV: Davidson iteration: '', I10)') i

      !! Array deallocation/allocation
      call check_deallocate_matrices( mtx_proj, stx_proj, lambda, eigenvectors_sub, ritz_vectors)
      if( allocated( eigenvalues_sub)) then
         deallocate( eigenvalues_sub)
      end if
      allocate( eigenvalues_sub( size( V, 2)))
      allocate( lambda( size( V, 2), size( V, 2))) 
      allocate( eigenvectors_sub( size( V, 2), size( V, 2)))

      ! 4. Projection of H and S matrices

      mtx_proj= lapack_matmul( 'T', 'N', V, lapack_matmul( 'N', 'N', mtx, V))

      if( gev) then
        stx_proj= lapack_matmul( 'T', 'N', V, lapack_matmul( 'N', 'N', stx, V))
      end if

      if (nproc > 1) then
        call MPI_BCAST( mtx_proj, size( V, 1)* size( V, 2), MPI_REAL8, 0, MPI_COMM_WORLD, ier)
        if( gev) & 
        call MPI_BCAST( stx_proj, size( V, 1)* size( V, 2), MPI_REAL8, 0, MPI_COMM_WORLD, ier)
      endif

      !! IF IDTASK=0
      if( idtask.eq. 0) then

        ! 5. compute the eigenvalues and their corresponding ritz_vectors
        ! for the projected matrix using lapack

        write(6,'(''DAV: enter lapack_generalized_eigensolver'')')
        if( gev) then
          call lapack_generalized_eigensolver( mtx_proj, eigenvalues_sub, eigenvectors_sub, stx_proj)
        else
          call lapack_generalized_eigensolver( mtx_proj, eigenvalues_sub, eigenvectors_sub)
        end if
        write(6,'(''DAV: exit lapack_generalized_eigensolver'')')
        write(6,'(''DAV: eigv'',1000f12.5)') ( eigenvalues_sub( j), j= 1, lowest)

        ! 6. Construction of lambda matrix (a squared one with eigenvalues_sub
        ! in the diagonal)

        lambda= eye( size( V, 2), size( V, 2))
        do j= 1, size(V,2)
          lambda( j, j)= eigenvalues_sub( j)
        enddo
      !! ENDIF IDTASK
      endif

      if (nproc > 1) then
        call MPI_BCAST(eigenvalues_sub, size( mtx_proj, 1), MPI_REAL8, 0, MPI_COMM_WORLD,ier)
        call MPI_BCAST(eigenvectors_sub, size( mtx_proj, 1)* size( mtx_proj, 2), MPI_REAL8, 0, MPI_COMM_WORLD,ier)
      endif       

      ritz_vectors = lapack_matmul('N', 'N', V, eigenvectors_sub(:, :))

      ! 7. Residue calculation

      !! IF IDTASK=0

      if( idtask .eq. 0) then
        if( gev) then
          rs= lapack_matmul('N', 'N', stx, ritz_vectors) 
         else 
          rs = ritz_vectors 
        endif
        guess =  lapack_matmul('N', 'N', rs, lambda)
        deallocate(rs)
        rs =  lapack_matmul('N', 'N', mtx, ritz_vectors) - guess
        do j=1,lowest
          errors( j) = norm(reshape( rs( :, j), (/ parameters%nparm/)))
        end do

      !! ENDIF IDTASK
      endif

      if (nproc > 1) then
        call MPI_BCAST(errors, lowest, MPI_REAL8, 0, MPI_COMM_WORLD,ier)
      endif       

      ! 8. Check for convergence

      if (all(errors < tolerance)) then
        iters = i
        write( 6, '(''DAV: roots are converged'')')
        eigenvalues = eigenvalues_sub(:lowest)
        exit outer_loop
      end if

      ! 7. Add the correction vectors to the current basis

      if(( size( V, 2) <= nvecx).and.( 2* size( V, 2)< nparm)) then

          ! append correction to the current basis
          call check_deallocate_matrix( correction)
          allocate( correction( size( mtx, 1), size( V, 2)))

          if( gev) then
            correction= compute_correction_generalized_dense( mtx, V, eigenvalues_sub, eigenvectors_sub, method, rs, ritz_vectors, stx)
          else
            correction= compute_correction_generalized_dense( mtx, V, eigenvalues_sub, eigenvectors_sub, method)
          end if

          ! 8. Increase Basis size
          call concatenate(V, correction)
       
          ! IF IDTASK=0
          if(idtask.eq.0) then

            ! 9. Orthogonalize basis

            call lapack_qr(V)

          ! ENDIF IDTASK
          endif

        if (nproc > 1) then
          call MPI_BCAST(V, size(V,1)*size(V,2), MPI_REAL8, 0, MPI_COMM_WORLD,ier)
        endif
          
      else

        ! 10. Otherwise reduce the basis of the subspace to the current correction

        V = lapack_matmul('N', 'N', V, eigenvectors_sub(:, :dim_sub))

      end if

    end do outer_loop

    !  11. Check convergence
    if (i > max_iters) then
       iters = i
       print *, "Warning: Algorithm did not converge!!"
    end if

    ! Select the lowest eigenvalues and their corresponding ritz_vectors
    ! They are sort in increasing order
    eigenvalues= eigenvalues_sub( :lowest)

    ! Free memory
    call check_deallocate_matrix( correction)
    deallocate( eigenvalues_sub, eigenvectors_sub, mtx_proj, V, guess, rs, lambda)

    ! free optional matrix
    if (gev) then
       call check_deallocate_matrix(stx_proj)
    endif

  end subroutine generalized_eigensolver_dense
  
end module davidson_dense

submodule (davidson_dense) correction_methods_generalized_dense
  !> submodule containing the implementations of different kind
  !> algorithms to compute the correction vectors for the Davidson's diagonalization

  implicit none
  
contains

  module function compute_correction_generalized_dense(mtx, V, eigenvalues, eigenvectors, method, rs, ritz_vectors, stx) &
       result(correction)
    !> see interface in davidson module
    real(dp), dimension(:), intent(in) :: eigenvalues
    real(dp), dimension(:, :), intent(in) :: mtx, V, eigenvectors
    real(dp), dimension(:, :), intent(in), optional :: ritz_vectors 
    real(dp), dimension(:, :), intent(in), optional :: rs 
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
       correction = compute_GJD_generalized_dense(mtx, size(V,2), eigenvalues,  rs, ritz_vectors, stx)
      else
        correction = compute_GJD_generalized_dense(mtx, size(V,2), eigenvalues, residues=rs, ritz_vectors=ritz_vectors)
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
    integer :: ii, j, m
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
               correction(ii, j) = correction(ii, j) / (eigenvalues(j)  * stx(ii,ii)  - mtx(ii, ii))
           else
              correction(ii, j) = correction(ii, j) / (eigenvalues(j)  - mtx(ii, ii))
          endif
        end do
    end do

  end function compute_DPR_generalized_dense

  function compute_GJD_generalized_dense(mtx, basis_size, eigenvalues, residues, ritz_vectors, stx) &
           result(correction)
    !> Compute the Generalized Jacobi Davidson (GJD) correction
    
    real(dp), dimension( :), intent( in) :: eigenvalues
    integer, intent( in):: basis_size
    real( dp), dimension( :, :), intent( in) :: mtx
    real( dp), dimension( :, :), intent( in), optional :: residues 
    real( dp), dimension( :, :), intent( in), optional :: ritz_vectors 
    real( dp), dimension( :, :), intent( in), optional :: stx
    real( dp), dimension( size(mtx, 1), basis_size) :: correction

    ! local variables
    integer :: k, m
    logical :: gev
    real( dp), dimension( size( mtx, 1), 1) :: ritz_tmp
    real( dp), dimension( size( mtx, 1), size( mtx, 2)) :: F, ubut, uubt, ys
    real( dp), dimension( size( mtx, 1), 1) :: brr

    ! Diagonal matrix
    m= size(mtx, 1)

    gev= present(stx)

    do k= 1, basis_size 
       ritz_tmp( :, 1)= ritz_vectors( :, k)

       if( gev) then
         ys= mtx- eigenvalues( k)* stx
       else
         ys= substract_from_diagonal( mtx, eigenvalues( k))
       end if

! Following I. Sabzevari, A. Mahajan and S. Sharma, arXiv:1908.04423 (2019)
       ubut= eye( m, m)- lapack_matmul( 'N', 'T', lapack_matmul( 'N', 'N', stx, ritz_tmp), ritz_tmp)
       uubt= eye( m, m)- lapack_matmul( 'N', 'T', ritz_tmp, lapack_matmul( 'N', 'N', stx, ritz_tmp))

       F= lapack_matmul( 'N', 'N', ubut, lapack_matmul( 'N', 'N', ys, uubt))
       brr( :, 1) = -residues(:,k)

       call lapack_solver( F, brr)
       correction( :, k)= brr( :, 1)
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

module davidson

  use numeric_kinds, only: dp
  use lapack_wrapper, only: lapack_generalized_eigensolver, lapack_matmul, lapack_matrix_vector, &
       lapack_qr, lapack_solver
  use array_utils, only: concatenate, eye, norm
  use davidson_dense, only: generalized_eigensolver_dense
  use davidson_free, only: generalized_eigensolver_free
  implicit none

!
  !> \private
  private
  !> \public
  public :: generalized_eigensolver

  interface generalized_eigensolver
 !> \brief Solve a (general) eigenvalue problem using different types of Davidson algorithms.

 !> \param[in] mtx: Matrix to diagonalize
 !> \param[in, opt] stx: Optional matrix for the general eigenvalue problem:
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

     procedure generalized_eigensolver_dense
     procedure generalized_eigensolver_free
  end interface generalized_eigensolver
!
end module davidson
!> \namespace Davidson eigensolver
