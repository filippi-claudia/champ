SUBROUTINE davidson_wrap(nparm, nparmx, nvec, nvecx, mvec, eigenvectors, ethr, &
                         eigenvalues, btype, notcnv, dav_iter, ipr, method)
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
    use precision_kinds, only: dp
    use davidson, only: generalized_eigensolver
    use davidson, only: davidson_parameters
    use array_utils, only: eye, write_matrix, write_vector
    use mpi
    use contrl_file,    only: ounit, errunit

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

    REAL(dp), dimension(nparmx, nvec), INTENT(INOUT) :: eigenvectors
    REAL(dp), dimension(nvec), INTENT(OUT) :: eigenvalues
    REAL(dp), INTENT(IN) :: ethr

    INTEGER, INTENT(IN) :: nparm, nparmx, nvec, nvecx, mvec, ipr
    INTEGER, dimension(nvec), INTENT(IN) :: btype
    INTEGER, INTENT(OUT) :: dav_iter, notcnv
    character(len=*), intent(in) :: method

    ! local variables
    integer :: i, ierr
    integer :: nproc, idtask
    real(dp), dimension(nparm, nparm) :: mtx, stx
    real(dp), dimension(nparmx, nparmx) :: psi
    real(dp), dimension(:, :), allocatable :: hpsi, spsi, ritz_vectors

    ! Function to compute the target matrix on the fly

    interface
        function fun_mtx_gemv(parameters, input_vect) result(output_vect)
            !> Brief Function to compute the action of the hamiltonian on the fly
            !> \param[in] dimension of the arrays to compute the action of the hamiltonian
            !> \param[in] input_vec Array to project
            !> return Projected matrix

            use precision_kinds, only: dp
            import :: davidson_parameters
            type(davidson_parameters) :: parameters
            real(dp), dimension(:, :), intent(in) :: input_vect
            real(dp), dimension(size(input_vect, 1), size(input_vect, 2)) :: output_vect

        end function fun_mtx_gemv

        function fun_stx_gemv(parameters, input_vect) result(output_vect)
            !> Brief Function to compute the action of the overlap on the fly
            !> \param[in] dimension of the arrays to compute the action of the hamiltonian
            !> \param[in] input_vec Array to project
            !> return Projected matrix

            use precision_kinds, only: dp
            import :: davidson_parameters
            type(davidson_parameters) :: parameters
            real(dp), dimension(:, :), intent(in) :: input_vect
            real(dp), dimension(size(input_vect, 1), size(input_vect, 2)) :: output_vect

        end function fun_stx_gemv

    end interface

    allocate (ritz_vectors(nparm, nvec), source=0.0_dp)

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

function fun_mtx_gemv(parameters, input_vect) result(output_vect)
    !> Brief Function to compute the action of the hamiltonian on the fly
    !> \param[in] dimension of the arrays to compute the action of the hamiltonian
    !> \param[in] input_vec Array to project
    !> return Projected matrix

    use precision_kinds, only: dp
    use davidson, only: davidson_parameters

    type(davidson_parameters) :: parameters
    real(dp), dimension(:, :), intent(in) :: input_vect
    real(dp), dimension(size(input_vect, 1), size(input_vect, 2)) :: output_vect
    real(dp), dimension(:, :), allocatable :: psi, hpsi

    allocate (psi(parameters%nparm_max, 2*parameters%nvecx), source=0.0_dp)
    allocate (hpsi(parameters%nparm_max, 2*parameters%nvecx), source=0.0_dp)

    psi = 0.0_dp
    psi(1:size(input_vect, 1), 1:size(input_vect, 2)) = input_vect

    call h_psi_lin_d(parameters%nparm, size(input_vect, 2), psi, hpsi)

    output_vect = hpsi(1:size(input_vect, 1), 1:size(input_vect, 2))
    deallocate (psi, hpsi)

end function fun_mtx_gemv

function fun_stx_gemv(parameters, input_vect) result(output_vect)
    !> Brief Fucntion to compute the optional stx matrix on the fly
    !> \param[in] dimension of the arrays to compute the action of the hamiltonian
    !> \param[in] input_vec Array to project
    !> return Projected matrix

    use precision_kinds, only: dp
    use davidson, only: davidson_parameters

    type(davidson_parameters) :: parameters
    real(dp), dimension(:, :), intent(in) :: input_vect
    real(dp), dimension(size(input_vect, 1), size(input_vect, 2)) :: output_vect
    real(dp), dimension(:, :), allocatable :: psi, spsi

    allocate (psi(parameters%nparm_max, 2*parameters%nvecx), source=0.0_dp)
    allocate (spsi(parameters%nparm_max, 2*parameters%nvecx), source=0.0_dp)

    psi = 0.0_dp
    psi(1:size(input_vect, 1), 1:size(input_vect, 2)) = input_vect

    call s_psi_lin_d(parameters%nparm, size(input_vect, 2), psi, spsi)

    output_vect = spsi(1:size(input_vect, 1), 1:size(input_vect, 2))
    deallocate (psi, spsi)

end function fun_stx_gemv
