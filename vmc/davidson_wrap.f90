SUBROUTINE davidson_wrap( nparm, nparmx, nvec, nvecx, mvec, eigenvectors, ethr, &
                    eigenvalues, btype, notcnv, dav_iter, ipr, idtask, free)
  !----------------------------------------------------------------------------
  !
  ! ... iterative solution of the eigenvalue problem:
  !
  ! ... ( H - e S ) * v = 0
  !
  ! ... where H is an hermitean operator, e is a real scalar,
  ! ... S is an symmetric matrix, v is a real vector
  ! ... (real wavefunctions with only half plane waves stored)
  ! 
  use numeric_kinds, only: dp
  use davidson, only: generalized_eigensolver
  use davidson_free, only: davidson_parameters
  use array_utils, only: eye

  IMPLICIT NONE

  !> \param npram dimension of the matrix to be diagonalized
  !> \param nparmx leading dimension of matrix eigenvectors
  !> \param nvec integer number of searched low-lying roots
  !> \param nvecx maximum dimension of the reduced basis set
  !> \param eigenvectors   contains the refined estimates of the eigenvectors
  !> \param ethr  energy threshold for convergence: root improvement is stopped,
  !> \param btype band type ( 1 = occupied, 0 = empty )
  !> \param eigenvalues Eigenvalues
  !> \param dav_iter integer  number of iterations performed
  !> \param notcnv number of unconverged roots
  !> \param notcnv number of unconverged roots
  !> \param free free version (.true.) or dens version (.false.)  


  REAL(dp), dimension(nparmx, nvec),  INTENT(INOUT) :: eigenvectors
  REAL(dp), dimension(nvec), INTENT(OUT) :: eigenvalues
  REAL(dp), INTENT(IN) :: ethr

  INTEGER, INTENT(IN) :: nparm, nparmx, nvec, nvecx, mvec, ipr, idtask
  INTEGER, dimension(nvec), INTENT(IN) :: btype
  INTEGER, INTENT(OUT) :: dav_iter, notcnv
  LOGICAL, INTENT(IN)  :: free 
  ! local variables
  integer :: i
  real(dp), dimension(nparm, nparm) :: mtx, stx
  real(dp), dimension(nparmx, nparmx) :: psi
  real(dp), dimension(nparm, nvec) :: ritz_vectors
  real(dp), dimension(:, :), allocatable :: hpsi, spsi
  
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

  notcnv=0 !Not used in davidson_wrap
    
  ! Allocate variables
  IF ( nvec > nvecx / 2 ) CALL fatal_error( 'regter: nvecx is too small')
  is_free_or_dens: IF (free) then   
!
      call generalized_eigensolver(fun_mtx_gemv, eigenvalues, eigenvectors, nparm, &
             nparmx,  nvec, mvec, "DPR", 100, ethr, dav_iter, nvecx, fun_stx_gemv, idtask)
!
      do i=1,size(eigenvalues)
        print *, "eigenvalue ", i, " : ", eigenvalues(i)
      end do
!
    ELSEIF (.not.free) then 
!
      ! Allocate Arrays to compute H ans S
      psi = eye(nparmx, nparmx, 1.0_dp)
      allocate(hpsi(nparmx, nparmx))
      allocate(spsi(nparmx, nparmx))
      hpsi = 0.0_dp
      spsi = 0.0_dp
!
      call h_psi_lin_d(nparm, nparm, psi, hpsi)
      call s_psi_lin_d(nparm, nparm, psi, spsi)
      
      mtx(1:nparm, 1:nparm) = hpsi(1:nparm, 1:nparm)
      stx(1:nparm, 1:nparm) = spsi(1:nparm, 1:nparm)
!     
      call generalized_eigensolver(mtx, eigenvalues, ritz_vectors, nvec, &
           "DPR", 100, ethr, dav_iter, nvecx, stx)
!
      eigenvectors(1:nparm,1:nvec) = ritz_vectors
!
      do i=1,size(eigenvalues)
         print *, "eigenvalue ", i, " : ", eigenvalues(i)
      end do
!
      notcnv = 0
      dav_iter = 0
!
      deallocate(hpsi, spsi)
 
  ENDIF is_free_or_dens
  
END SUBROUTINE davidson_wrap

function fun_mtx_gemv(parameters,  input_vect) result(output_vect)
  !> \brief Function to compute the action of the hamiltonian on the fly
  !> \param[in] dimension of the arrays to compute the action of the hamiltonian
  !> \param[in] input_vec Array to project
  !> \return Projected matrix

  use numeric_kinds, only: dp
  use davidson_free, only: davidson_parameters

  type(davidson_parameters) :: parameters
  real (dp), dimension(:,:), intent(in) :: input_vect
  real (dp), dimension(size(input_vect,1),size(input_vect,2)) :: output_vect
  real(dp), dimension(:, :), allocatable :: psi, hpsi
  
  allocate(psi(parameters%nparm_max, parameters%mvec))
  allocate(hpsi(parameters%nparm_max, parameters%mvec))
  psi = 0.0_dp
  psi(1:size(input_vect,1),1:size(input_vect,2)) = input_vect
  
  call h_psi_lin_d(parameters%nparm, size(input_vect,2), psi, hpsi)

  output_vect = hpsi(1:size(input_vect,1),1:size(input_vect,2))
  deallocate(psi, hpsi)
  
end function fun_mtx_gemv

function fun_stx_gemv(parameters, input_vect) result(output_vect)
  !> \brief Fucntion to compute the optional stx matrix on the fly
  !> \param[in] dimension of the arrays to compute the action of the hamiltonian
  !> \param[in] input_vec Array to project
  !> \return Projected matrix
  
  use numeric_kinds, only: dp
  use davidson_free, only: davidson_parameters
  type(davidson_parameters) :: parameters
  real (dp), dimension(:,:), intent(in) :: input_vect
  real (dp), dimension(size(input_vect,1),size(input_vect,2)) :: output_vect
  real(dp), dimension(:, :), allocatable :: psi, spsi

  allocate(psi(parameters%nparm_max, parameters%mvec))
  allocate(spsi(parameters%nparm_max, parameters%mvec))
  psi = 0.0_dp
  psi(1:size(input_vect,1),1:size(input_vect,2)) = input_vect
  
  call s_psi_lin_d(parameters%nparm, size(input_vect,2), psi, spsi)

  output_vect = spsi(1:size(input_vect,1),1:size(input_vect,2))
  deallocate(psi, spsi)
  
end function fun_stx_gemv
