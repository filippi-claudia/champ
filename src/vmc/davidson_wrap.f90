!>----------------------------------------------------------------------------
!>
!> ... iterative solution of the eigenvalue problem:
!>
!> ... ( H - e S ) * v = 0
!>
!> ... where H is an hermitean operator, e is a real scalar,
!> ... S is an symmetric matrix, v is a real vector
!> ... (real wavefunctions with only half plane waves stored)
!>
module davidson_wrap_mod
contains
SUBROUTINE davidson_wrap(nparm, nparmx, nvec, nvecx, mvec, eigenvectors, ethr, &
                         eigenvalues, btype, notcnv, dav_iter, ipr, method)

      use array_utils, only: eye,write_matrix,write_vector
      use contrl_file, only: errunit,ounit
      use davidson, only: davidson_parameters,fun_mtx_gemv,fun_stx_gemv
      use davidson, only: generalized_eigensolver
      use error,   only: fatal_error
      use mpi
      use optwf_lin_dav_extra, only: s_psi_lin_d
      use precision_kinds, only: dp

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

    INTEGER, INTENT(IN) :: nparm, nparmx, nvec, nvecx, mvec, ipr
    REAL(dp), dimension(nparmx, nvec), INTENT(INOUT) :: eigenvectors
    REAL(dp), dimension(nvec), INTENT(OUT) :: eigenvalues
    REAL(dp), INTENT(IN) :: ethr

    INTEGER, dimension(nvec), INTENT(IN) :: btype
    INTEGER, INTENT(OUT) :: dav_iter, notcnv
    character(len=*), intent(in) :: method

    ! local variables
    integer :: i, ierr
    integer :: nproc, idtask
    real(dp), dimension(nparm, nparm) :: mtx, stx
    real(dp), dimension(nparmx, nparmx) :: psi
    real(dp), dimension(:, :), allocatable :: hpsi, spsi, ritz_vectors


    allocate (ritz_vectors(nparm, nvec))

    call mpi_comm_rank(MPI_COMM_WORLD, idtask, ierr)
    call mpi_comm_size(MPI_COMM_WORLD, nproc, ierr)
    write (ounit, '(''DAV: idtask      : '', I10)') idtask
    notcnv = 0 !Not used in davidson_wrap

    ! Allocate variables
    IF (nvec > nvecx/2) CALL fatal_error('regter: nvecx is too small')
    call generalized_eigensolver(fun_mtx_gemv, eigenvalues, ritz_vectors, nparm, &
                                 nparmx, nvec, nvecx, method, 200, ethr, dav_iter, fun_stx_gemv, nproc, idtask)

    if (idtask == 0) then
        do i = 1, size(eigenvalues)
            write(ounit,*) "(DAV) Eigenvalue ", i, " : ", eigenvalues(i)
        end do
    endif

    eigenvectors(1:nparm, 1:nvec) = ritz_vectors

    deallocate (ritz_vectors)
END SUBROUTINE davidson_wrap

end module
