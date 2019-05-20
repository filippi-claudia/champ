SUBROUTINE davidson_wrap( nparm, nparmx, nvec, nvecx, eigenvectors, ethr, &
                    eigenvalues, btype, notcnv, dav_iter, ipr, idtask)
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


  REAL(dp), dimension(nparmx, nvec),  INTENT(INOUT) :: eigenvectors
  REAL(dp), dimension(nvec), INTENT(OUT) :: eigenvalues
  REAL(dp), INTENT(IN) :: ethr

  INTEGER, INTENT(IN) :: nparm, nparmx, nvec, nvecx, ipr, idtask
  INTEGER, dimension(nvec), INTENT(IN) :: btype
  INTEGER, INTENT(OUT) :: dav_iter, notcnv

  
  ! ... LOCAL variables
  REAL(dp), ALLOCATABLE :: psi(:,:), hpsi(:,:), spsi(:,:)
  INTEGER :: ierr ! Error info

    ! Function to compute the target matrix on the fly
    interface

       function apply_mtx_to_vect(nparm, lowest, input_vect) result(output_vect)
         !> \brief Function to compute the optional mtx on the fly
         !> \param lowest Number of lowest eigenvalues
         !> \param[in] i column/row to compute from mtx
         !> \return vec column/row from mtx
         use numeric_kinds, only: dp
         integer, intent(in) :: nparm, lowest 
         real (dp), dimension(:,:), intent(in) :: input_vect
         real (dp), dimension(size(input_vect,1),size(input_vect,2)) :: output_vect

       end function apply_mtx_to_vect
       
       function apply_stx_to_vect(nparm, lowest, input_vect) result(output_vect)
         !> \brief Fucntion to compute the optional stx matrix on the fly
         !> \param[in] i column/row to compute from stx
         !> \param vec column/row from stx
         use numeric_kinds, only: dp
         integer, intent(in) :: nparm, lowest 
         real(dp), dimension(:,:), intent(in) :: input_vect
         real (dp), dimension(size(input_vect,1),size(input_vect,2)) :: output_vect
         
       end function apply_stx_to_vect

    end interface

  
  ! Allocate variables
  IF ( nvec > nvecx / 2 ) CALL fatal_error( 'regter: nvecx is too small')

  ALLOCATE( psi(  nparmx, nvecx ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL fatal_error( 'davidson_wrap: cannot allocate psi ')
  ALLOCATE( hpsi( nparmx, nvecx ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL fatal_error( 'davidson_wrap: cannot allocate hpsi ')
  !
  ALLOCATE( spsi( nparmx, nvecx ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL fatal_error( ' davidson_wrap: cannot allocate spsi ')
  
  ! Initialize
  spsi = 0.0_dp
  hpsi = 0.0_dp
  psi  = 0.0_dp
  psi(:,1:nvec) = eigenvectors(:,1:nvec)
 
  notcnv=0 !Not used in davidson_wrap
  
  call generalized_eigensolver(apply_mtx_to_vect, eigenvalues, eigenvectors, nparm, nvec, &
       "DPR", 100, ethr, dav_iter, nvecx, apply_stx_to_vect, idtask)
  
END SUBROUTINE davidson_wrap

function apply_mtx_to_vect(nparm, nvec, psi) result(hpsi)
  !> \brief interface to compute the product between the Hamiltonian and the eigenvectors
  !> \param nvec Number of pair eigenvalues/eigenvectors to compute
  !> \return the result of applying the Hamiltonian to the trial eigenvectors
  
  use numeric_kinds, only: dp

  implicit none
  
  ! IO variables
  integer, intent(in) :: nparm, nvec
  real(dp), dimension(:, :), intent(in) :: psi

  real(dp), dimension(size(psi, 1), size(psi, 2)) :: hpsi

  ! nparmx = size(psi, 1)
  call h_psi_lin_d(nparm, nvec, psi, hpsi) 

end function apply_mtx_to_vect


function apply_stx_to_vect(nparm, nvec, psi) result(spsi)
  !> \brief interface to compute the product between the Overlap and the eigenvectors
  !> \param nvec Number of pair eigenvalues/eigenvectors to compute
  !> \return the result of applying the Hamiltonian to the trial eigenvectors
  
  use numeric_kinds, only: dp
  
  implicit none
  
  ! IO variables
  integer, intent(in) :: nparm, nvec
  real(dp), dimension(:, :), intent(in) :: psi

  real(dp), dimension(size(psi, 1), size(psi, 2)) :: spsi
  ! nparm = size(psi, 1)
  CALL s_psi_lin_d(nparm, nvec, psi, spsi)
  
end function apply_stx_to_vect

