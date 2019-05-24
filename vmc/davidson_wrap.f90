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

  use array_utils, only: eye
  use numeric_kinds, only: dp
  use davidson_dense, only: generalized_eigensolver_dense

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

  ! local variables
  integer :: i
  real(dp), dimension(nparm, nparm) :: mtx, stx
  real(dp), dimension(nparmx, nparmx) :: psi

  ! Allocate Arrays to compute H ans S
  psi = eye(nparm, nparm, 1.0_dp)
  call h_psi_lin_d(nparm, nparm, psi, mtx)
  call s_psi_lin_d(nparm, nparm, psi, stx)

  call write_matrix("H.txt", mtx, nparm)
  call write_matrix("S.txt", stx, nparm)
     
  call generalized_eigensolver_dense(mtx, eigenvalues, eigenvectors, nvec, &
       "DPR", 100, ethr, dav_iter, nvecx, stx)

  do i=1,size(eigenvalues)
     print *, "eigenvalue ", i, " : ", eigenvalues(i)
  end do
  
  notcnv = 0

END SUBROUTINE davidson_wrap


subroutine write_matrix(path_file, mtx, dim)
  !> Write matrix to path_file
  use numeric_kinds, only: dp
  integer, INTENT(in) :: dim
  character(len=*), intent(in) :: path_file
  real(dp), dimension(dim, dim), intent(in) :: mtx
  integer :: i, j

  open(unit=2, file=path_file, status="REPLACE")
  print *, "size: ", size(mtx, 1), size(mtx, 2)
  do i=1, dim
     do j=1, dim
        write(2, *) mtx(i, j)
     end do
  end do
  close(2)

end subroutine write_matrix
