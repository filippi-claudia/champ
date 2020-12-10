!! The current implementation uses a general  davidson algorithm, meaning
!! that it compute all the eigenvalues simultaneusly using a variable size block approach.
!! The family of Davidson algorithm only differ in the way that the correction
!! vector is computed.
!! Computed pairs of eigenvalues/eigenvectors are deflated using algorithm
!! described at: https://doi.org/10.1023/A:101919970
!!
!! Authors: Felipe Zapata and revised by P. Lopez-Tarifa NLeSC(2019)
!>
!> rief Solves the Davidson diagonalisation without storing matrices in memory.
!>
!> Matrices mtx and stx are calculated on the fly using the fun_mtx_gemv and fun_stx_gemv
!> functions.
!>
!> uthor Felipe Zapata and P. Lopez-Tarifa NLeSC(2019)
!>
!> \param[in] mtx: Matrix to diagonalize
!> \param[inout] stx: Optional matrix for the general eigenvalue problem:
!> $ mtx \lambda = V stx \lambda $
!> \param[out] eigenvalues Computed eigenvalues
!> \param[out] ritz_vectors Approximation to the eigenvectors
!> \param[in] lowest Number of lowest eigenvalues/ritz_vectors to compute
!> \param[in] method Method to compute the correction vector. Available
!> methods are:
!>    - DPR: Diagonal-Preconditioned-Residue.
!>    - GJD: Generalized Jacobi Davidsoni.
!> \param[in] max_iters: Maximum number of iterations.
!> \param[in] tolerance Norm**2 error of the eigenvalues.
!> \param[in] method: Method to compute the correction vectors.
!> \param[out] iters: Number of iterations until convergence.
!>
eturn eigenvalues and ritz_vectors of the matrix.
module davidson
    use numeric_kinds, only: dp
    use lapack_wrapper, only: lapack_generalized_eigensolver, lapack_matmul, lapack_matrix_vector, &
                              lapack_qr, lapack_solver
    use array_utils, only: concatenate, initialize_subspace, norm, write_matrix, write_vector, &
                           eye, check_deallocate_matrix, check_deallocate_vector, modified_gram_schmidt, diag_mat
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
    public :: generalized_eigensolver, davidson_parameters, die

contains

    subroutine generalized_eigensolver(fun_mtx_gemv, eigenvalues, eigenvectors, nparm, nparm_max, &
                                       lowest, nvecx, method, max_iters, tolerance, iters, fun_stx_gemv, nproc, idtask)
        !> rief use a pair of functions fun_mtx and fun_stx to compute on the fly the matrices to solve
        !>  the general eigenvalue problem
        !> The current implementation uses a general  davidson algorithm, meaning
        !> that it compute all the eigenvalues simultaneusly using a block approach.
        !> The family of Davidson algorithm only differ in the way that the correction
        !> vector is computed.

        !> \param[in] fun_mtx_gemv: Function to apply the matrix to a buncof vectors
        !> \param[in, opt] fun_stx_gemv: (optional) function to apply the pencil to a bunch of vectors.
        !> \param[out] eigenvalues Computed eigenvalues
        !> \param[inout] eigenvectors approximation to the eigenvectors
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
        !> \param[out] iters: Number of iterations until convergence
        !>
        eturn eigenvalues and ritz_vectors of the matrix

        implicit none

        include 'mpif.h'

        ! input/output variable
        integer, intent(in) :: nparm, nparm_max, nvecx, lowest, nproc, idtask
        real(dp), dimension(lowest), intent(out) :: eigenvalues
        real(dp), dimension(:, :), allocatable, intent(out) :: eigenvectors
        integer, intent(in) :: max_iters
        real(dp), intent(in) :: tolerance
        character(len=*), intent(in) :: method
        integer, intent(out) :: iters

        ! Function to compute the target matrix on the fly
        interface

            function fun_mtx_gemv(parameters, input_vect) result(output_vect)
                !> rief Function to compute the action of the hamiltonian on the fly
                !> \param[in] dimension of the arrays to compute the action of the hamiltonian
                !> \param[in] input_vec Array to project
                !>
                eturn Projected matrix

                use precision_kinds, only: dp
                import :: davidson_parameters
                type(davidson_parameters) :: parameters
                real(dp), dimension(:, :), intent(in) :: input_vect
                real(dp), dimension(size(input_vect, 1), size(input_vect, 2)) :: output_vect

            end function fun_mtx_gemv

            function fun_stx_gemv(parameters, input_vect) result(output_vect)
                !> rief Fucntion to compute the optional stx matrix on the fly
                !> \param[in] dimension of the arrays to compute the action of the hamiltonian
                !> \param[in] input_vec Array to project
                !>
                eturn Projected matrix

                use precision_kinds, only: dp
                import :: davidson_parameters
                type(davidson_parameters) :: parameters
                real(dp), dimension(:, :), intent(in) :: input_vect
                real(dp), dimension(size(input_vect, 1), size(input_vect, 2)) :: output_vect

            end function fun_stx_gemv

        end interface

        ! Local variables
        integer :: i, j, ier
        integer :: init_subspace_size, max_size_basis, size_update

        ! Basis of subspace of approximants
        real(dp), dimension(:), allocatable :: diag_mtx, diag_stx, diag_mtx_cpy
        real(dp), dimension(:, :), allocatable :: residues
        real(dp), dimension(lowest):: errors

        ! Working arrays
        real(dp), dimension(:), allocatable :: eigenvalues_sub
        real(dp), dimension(:, :), allocatable :: ritz_vectors
        real(dp), dimension(:, :), allocatable :: correction, eigenvectors_sub, mtx_proj, stx_proj, V
        real(dp), dimension(:, :), allocatable :: mtxV, stxV

        ! tmp arrays for the vectorsi residue calculation
        real(dp), dimension(:, :), allocatable :: lambda              ! eigenvalues_sub in a diagonal matrix
        real(dp), dimension(:, :), allocatable :: tmp_res_array       ! tmp array for vectorized res calculation

        ! Arrays dimension
        type(davidson_parameters) :: parameters

        ! Indices of the eigenvalues/eigenvectors pair that have not converged
        logical, dimension(lowest) :: has_converged
        logical :: update_proj
        integer :: n_converged ! Number of converged eigenvalue/eigenvector pairs
        integer :: sizeV ! size of V for broadcasting
        logical, parameter :: use_gs_ortho = .true.! which orthogonalization method to use gs/qr
        integer :: not_cnv

        ! Iteration subpsace dimension
        init_subspace_size = lowest*2

        ! number of correction vectors appended to V at each iteration
        size_update = lowest*2

        ! Make sure nvecx is lower than nparm
        if (nvecx > nparm) then
            if (idtask == 1) call die('DAV: nvecx > nparm, increase nparm or decrease lin_nvecx')
        endif

        ! Dimension of the matrix
        parameters = davidson_parameters(nparm, nparm_max, lowest, nvecx, init_subspace_size)

        ! Initial number of converged eigenvalue/eigenvector pairs
        n_converged = 0
        has_converged = .false.
        update_proj = .false.

        ! Diagonal of the arrays
        allocate (diag_mtx(parameters%nparm))
        allocate (diag_mtx_cpy(parameters%nparm))
        allocate (diag_stx(parameters%nparm))

        if (idtask == 0) call store_diag_hs(parameters%nparm, diag_mtx, diag_stx)

        ! why ?
        ! wouldn't it be faster to have all the procs computing that
        ! instead of master computes and then broadcast ?
        ! Needed for initalizing V so we need it
        if (nproc > 1) then
            call MPI_BCAST(diag_mtx, parameters%nparm, MPI_REAL8, 0, MPI_COMM_WORLD, ier)
            call MPI_BCAST(diag_stx, parameters%nparm, MPI_REAL8, 0, MPI_COMM_WORLD, ier)
        endif

        ! Select the initial ortogonal subspace based on lowest elements
        ! of the diagonal of the matrix.
        diag_mtx_cpy = diag_mtx
        V = initialize_subspace(diag_mtx_cpy, init_subspace_size, nparm) ! Initial orthonormal basis
        deallocate (diag_mtx_cpy)

        if (idtask == 0) write (6, '(''DAV: Setup subspace problem'')')

        ! Warning we reset the diag
        ! if (idtask==0) then
        !   do i=1, parameters%nparm
        !           diag_mtx(i) = 1.0
        !           diag_stx(i) = 0.0
        !   end do
        ! end if

        ! allocate mtxV and stxV
        allocate (mtxV(parameters%nparm, parameters%basis_size))
        allocate (stxV(parameters%nparm, parameters%basis_size))

        ! Calculation of HV and SV
        ! Only the master has the correct matrix
        ! nut only the master needs it
        mtxV = fun_mtx_gemv(parameters, V)
        stxV = fun_stx_gemv(parameters, V)

        ! allocate eigenvalues/vectors
        ! allocate(eigenvalues(parameters%lowest))
        allocate (eigenvectors(parameters%nparm, parameters%lowest))

        if (idtask == 0) then

            select case (method)
            case ("DPR")
                write (6, '(''DAV: Diagonal-Preconditioned-Residue (DPR)'')')
            case ("GJD")
                write (6, '(''DAV: Generalized Jacobi-Davidson (GJD)'')')
            end select

            write (6, '(''DAV: tolerance         : '', E10.3)') tolerance
            write (6, '(''DAV: Number eigenvalues : '', I10)') parameters%lowest
            write (6, '(''DAV: Update size        : '', I10)') size_update
            write (6, '(''DAV: Max basis size     : '', I10)') nvecx
            if (use_gs_ortho) then
                write (6, '(''Modified Gram-Schmidt orthogonalization with projection update'')')
            else
                write (6, '(''QR orthogonalization with full projection'')')
            end if
        endif

        ! 3. Outer loop block Davidson
        outer_loop: do i = 1, max_iters

            ! do most of the calculation on the master
            ! we could try to use openMP here via the lapack routines
            if (idtask == 0) then

                write (6, '(''DAV: -----------------------------'')')
                write (6, '(''DAV: Iteration: '', I10)') i
                write (6, '(''DAV: Basis size: '', I10)') parameters%basis_size

                ! needed if we want to vectorix the residue calculation
                call check_deallocate_matrix(lambda)
                call check_deallocate_matrix(tmp_res_array)
                ! allocate( lambda(size_update, size_update ))
                ! allocate( tmp_res_array(parameters%nparm, size_update ))

                ! reallocate eigenpairs of the small system
                call check_deallocate_vector(eigenvalues_sub)
                call check_deallocate_matrix(eigenvectors_sub)
                allocate (eigenvalues_sub(parameters%basis_size))
                allocate (eigenvectors_sub(parameters%basis_size, parameters%basis_size))

                ! deallocate the corection/residues
                call check_deallocate_matrix(residues)
                call check_deallocate_matrix(correction)
                allocate (residues(parameters%nparm, size_update))
                allocate (correction(parameters%nparm, size_update))

                ! deallocate ritz vectors
                call check_deallocate_matrix(ritz_vectors)

                ! update the projected matrices in the small subspace
                if (update_proj) then

                    ! update the projected matrices
                    call update_projection(V, mtxV, mtx_proj)
                    call update_projection(V, stxV, stx_proj)

                    ! recompute it from scratch when restarting
                else

                    ! Array deallocation/allocation.
                    call check_deallocate_matrix(mtx_proj)
                    call check_deallocate_matrix(stx_proj)

                    ! recompute the projected matrix
                    mtx_proj = lapack_matmul('T', 'N', V, mtxV)
                    stx_proj = lapack_matmul('T', 'N', V, stxV)

                end if

                ! Solve the small eigenvalue problem
                call lapack_generalized_eigensolver(mtx_proj, eigenvalues_sub, eigenvectors_sub, stx_proj)
                write (6, '(''DAV: eigv'',1000f12.5)') (eigenvalues_sub(j), j=1, parameters%lowest)

                ! Compute the necessary ritz vectors
                ritz_vectors = lapack_matmul('N', 'N', V, eigenvectors_sub(:, :size_update))

                ! Residue calculation (vectorized)
                lambda = diag_mat(eigenvalues_sub(:size_update))
                tmp_res_array = lapack_matmul('N', 'N', lapack_matmul('N', 'N', stxV, eigenvectors_sub(:, :size_update)), lambda)
                residues = lapack_matmul('N', 'N', mtxV, eigenvectors_sub(:, :size_update)) - tmp_res_array

                ! Check which eigenvalues has converged
                errors = norm2(residues(:, :parameters%lowest), 1)
                do j = 1, parameters%lowest
                    if (errors(j) < tolerance) has_converged(j) = .true.
                end do
                write (6, '(''DAV: resd'',1000f12.5)') (errors(j), j=1, parameters%lowest)

                not_cnv = count(.not. has_converged(:))
                write (6, '(''DAV: Root not yet converged     : '', I10)') not_cnv

                ! Append correction vectors
                if (parameters%basis_size + size_update <= nvecx) then

                    ! compute the correction vectors
                    select case (method)
                    case ("DPR")
                        correction = compute_DPR(residues, parameters, eigenvalues_sub, diag_mtx, diag_stx, has_converged)
                    case ("GJD")
                        correction = compute_GJD_free(parameters, ritz_vectors, residues, eigenvectors_sub, eigenvalues_sub)
                    end select

                    ! Add the correction vectors to the current basis.
                    call concatenate(V, correction(:, :not_cnv))

                    ! Orthogonalize basis using modified GS
                    if (use_gs_ortho) then
                        call modified_gram_schmidt(V, parameters%basis_size + 1)
                        update_proj = .true.
                    else
                        call lapack_qr(V)
                        update_proj = .false.
                    end if

                    ! Restart.
                else

                    write (6, '(''DAV: --- Restart ---'')')
                    V = ritz_vectors(:, :init_subspace_size)
                    update_proj = .false.

                end if

      !! ENDIF IDTASK
            end if

            ! Check for convergence
            ! all the procs need to know when to exit
            call MPI_BCAST(has_converged, parameters%lowest, MPI_LOGICAL, 0, MPI_COMM_WORLD, ier)
            if (all(has_converged)) then
                iters = i
                if (idtask == 0) then
                    write (6, '(''DAV: roots are converged'')')
                endif
                exit outer_loop
            end if

            ! broadcast the basis vector
            if (nproc > 1) then
                sizeV = size(V, 2)
                call MPI_BCAST(sizeV, 1, MPI_INT, 0, MPI_COMM_WORLD, ier)
                if (idtask > 0) then
                    call check_deallocate_matrix(V)
                    allocate (V(parameters%nparm, sizeV))
                end if
                call MPI_BCAST(V, parameters%nparm*sizeV, MPI_REAL8, 0, MPI_COMM_WORLD, ier)
            endif

            ! broadcast update proj
            call MPI_BCAST(update_proj, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ier)

            ! update mtxV and stxV
            if (update_proj) then

                call check_deallocate_matrix(tmp_res_array)
                tmp_res_array = fun_mtx_gemv(parameters, V(:, parameters%basis_size + 1:))
                call concatenate(mtxV, tmp_res_array)

                call check_deallocate_matrix(tmp_res_array)
                tmp_res_array = fun_stx_gemv(parameters, V(:, parameters%basis_size + 1:))
                call concatenate(stxV, tmp_res_array)

                ! recompute mtxV and stxV
            else

                call check_deallocate_matrix(mtxV)
                mtxV = fun_mtx_gemv(parameters, V)

                call check_deallocate_matrix(stxV)
                stxV = fun_stx_gemv(parameters, V)

            end if

            ! Update basis size
            parameters%basis_size = size(V, 2)

        end do outer_loop

        ! Master store eigenpairs
        if (idtask == 0) then

            ! if we didnt converge
            if (i > max_iters) then
                write (6, '(''DAV: Warning Davidson did not converge'')')
            else
                write (6, '(''DAV: Davidson converged'')')
            end if

            ! store eigenpairs
            eigenvalues = eigenvalues_sub(:parameters%lowest)
            eigenvectors = ritz_vectors(:, :parameters%lowest)
            iters = i

        end if

        ! Broadcast solutions to all slaves
        if (nproc > 1) then

            if (idtask == 0) then
                write (6, '(''DAV: Broadcasting solutions'')')
            end if
            call MPI_BCAST(iters, 1, MPI_INT, 0, MPI_COMM_WORLD, ier)
            call MPI_BCAST(eigenvalues, parameters%lowest, MPI_REAL8, 0, MPI_COMM_WORLD, ier)
            call MPI_BCAST(eigenvectors, parameters%nparm*parameters%lowest, MPI_REAL8, 0, MPI_COMM_WORLD, ier)

        endif

        ! Free memory
        deallocate (V, mtxV, stxV)

        if (idtask == 0) then
            call check_deallocate_matrix(correction)
            deallocate (eigenvalues_sub, ritz_vectors)
            deallocate (eigenvectors_sub)
            deallocate (diag_mtx, diag_stx)
            deallocate (residues)
            deallocate (lambda, tmp_res_array)
            deallocate (mtx_proj)
            call check_deallocate_matrix(stx_proj)
        endif

        if (idtask == 0) then
            write (6, '(''DAV: Exiting Davidson'')')
        endif

    end subroutine generalized_eigensolver
!
    subroutine update_projection(V, mtxV, mtx_proj)
        !> update the projected matrices
        !> \param mtxV: full matrix x projector
        !> \param V: projector
        !> \param mtx_proj: projected matrix

        implicit none
        real(dp), dimension(:, :), intent(in) :: mtxV
        real(dp), dimension(:, :), intent(in) :: V
        real(dp), dimension(:, :), intent(inout), allocatable :: mtx_proj
        real(dp), dimension(:, :), allocatable :: tmp_array

        ! local variables
        integer :: nvec, old_dim

        ! dimension of the matrices
        nvec = size(mtxV, 2)
        old_dim = size(mtx_proj, 1)

        ! move to temporal array
        allocate (tmp_array(nvec, nvec))
        tmp_array(:old_dim, :old_dim) = mtx_proj
        tmp_array(:, old_dim + 1:) = lapack_matmul('T', 'N', V, mtxV(:, old_dim + 1:))
        tmp_array(old_dim + 1:, :old_dim) = transpose(tmp_array(:old_dim, old_dim + 1:))

        ! Move to new expanded matrix
        deallocate (mtx_proj)
        call move_alloc(tmp_array, mtx_proj)

    end subroutine update_projection

    subroutine die(msg)
        !> Subroutine that dies the calculation raising an errror message
        !
        character msg*(*)
        integer ierr
        include 'mpif.h'

        write (6, '(''Fatal error: '',a)') msg
        call mpi_abort(MPI_COMM_WORLD, 0, ierr)

    end subroutine

    function compute_DPR(residues, parameters, eigenvalues, diag_mtx, diag_stx, has_converged) &
        result(correction)

        !> compute the correction vector using the DPR method for a matrix free diagonalization
        !> See correction_methods submodule for the implementations
        !> \param[in] fun_mtx: function to compute matrix
        !> \param[in] fun_stx: function to compute the matrix for the generalized case
        !> \param[in] V: Basis of the iteration subspace
        !> \param[in] eigenvalues: of the reduce problem
        !> \param[in] eigenvectors: of the reduce problem
        !>
        eturn correction matrix
        !
        use array_utils, only: eye
        !
        real(dp), dimension(:, :), intent(in) :: residues
        real(dp), dimension(:), intent(in) :: eigenvalues
        real(dp), dimension(:), intent(in) :: diag_mtx, diag_stx
        logical, dimension(:), intent(in) :: has_converged

        ! local variables
        type(davidson_parameters) :: parameters

        real(dp), dimension(parameters%nparm, size(residues, 2)) :: correction
        integer :: ii, j, k

        j = 1
        do k = 1, size(residues, 2)
            if (.not. has_converged(k)) then
                correction(:, j) = residues(:, k)

                do ii = 1, size(correction, 1)
                    correction(ii, j) = correction(ii, j)/(eigenvalues(k)*diag_stx(ii) - diag_mtx(ii))
                end do
                j = j + 1
            endif
        end do

    end function compute_DPR

    function compute_GJD_free(parameters, ritz_vectors, residues, eigenvectors, &
                              eigenvalues) result(correction)

        !> Compute the correction vector using the GJD method for a matrix free
        !> diagonalization. We follow the notation of:
        !> I. Sabzevari, A. Mahajan and S. Sharma,  arXiv:1908.04423 (2019)
        !>
        !> \param[in] ritz_vectors: ritz_vectors.
        !> \param[in] residues: residue vectors.
        !> \param[in] parameter: davidson_parameters type.
        !> \param[in] eigenvectors.
        !> \param[in] eigenvalues.
        !>
        eturn correction matrix

        use array_utils, only: eye

        type(davidson_parameters)               :: parameters
        real(dp), dimension(:, :), intent(in) :: ritz_vectors
        real(dp), dimension(:, :), intent(in) :: residues
        real(dp), dimension(:, :), intent(in) :: eigenvectors
        real(dp), dimension(:), intent(in) :: eigenvalues
        !
        ! local variables
        !
        real(dp), dimension(parameters%nparm, parameters%basis_size) :: correction
        integer :: k, m
        logical :: gev
        real(dp), dimension(:, :), allocatable   ::  F
        real(dp), dimension(parameters%nparm, 1) :: brr

        do k = 1, parameters%basis_size

            F = fun_F_matrix(ritz_vectors, parameters, k, eigenvalues(k))
            call write_matrix("F.txt", F)
            brr(:, 1) = -residues(:, k)
            call lapack_solver(F, brr)
            call write_matrix("brr.txt", brr)
            correction(:, k) = brr(:, 1)

        end do

        ! Deallocate
        deallocate (F)

    end function compute_GJD_free

    function fun_F_matrix(ritz_vectors, parameters, eigen_index, eigenvalue) &
        result(F_matrix)
        !> rief Function that computes the F matrix:
        !> F= ubut*( A- theta* B)* uubt
        !> in a pseudo-free way for a given engenvalue.
        !>
        !> ritz_vectors( in) :: ritz_vectors.
        !> parameters( in)   :: array_sizes
        !> eigen_index( in)  :: index of the passing eingenvalue.
        !> eigenvalue( in)   :: eigen_index eigenvalue.

        use array_utils, only: eye

        real(dp), dimension(:, :), intent(in) :: ritz_vectors
        type(davidson_parameters) :: parameters
        integer   :: eigen_index
        real(dp) :: eigenvalue

        interface

            function fun_mtx_gemv(parameters, input_vect) result(output_vect)
                !> rief Function to compute the action of the hamiltonian on the fly
                !> \param[in] dimension of the arrays to compute the action of the
                !             hamiltonian
                !> \param[in] input_vec Array to project
                !>
                eturn Projected matrix
                use precision_kinds, only: dp
                import                                   :: davidson_parameters
                type(davidson_parameters)               :: parameters
                real(dp), dimension(:, :), intent(in) :: input_vect
                real(dp), dimension(size(input_vect, 1), size(input_vect, 2)) :: output_vect
            end function fun_mtx_gemv

            function fun_stx_gemv(parameters, input_vect) result(output_vect)
                !> rief Fucntion to compute the optional stx matrix on the fly
                !> \param[in] dimension of the arrays to compute the action of the
                !             hamiltonian
                !> \param[in] input_vec Array to project
                !>
                eturn Projected matrix
                use precision_kinds, only: dp
                import                                   :: davidson_parameters
                type(davidson_parameters)               :: parameters
                real(dp), dimension(:, :), intent(in) :: input_vect
                real(dp), dimension(size(input_vect, 1), size(input_vect, 2)) :: output_vect
            end function fun_stx_gemv

        end interface

        real(dp), dimension(parameters%nparm, parameters%nparm) :: F_matrix, lambda
        real(dp), dimension(parameters%nparm, 1) :: ritz_tmp
        real(dp), dimension(:, :), allocatable :: ys
        real(dp), dimension(parameters%nparm, parameters%nparm) :: ubut, uubt

        ritz_tmp(:, 1) = ritz_vectors(:, eigen_index)

        lambda = eye(parameters%nparm, parameters%nparm, eigenvalue)

        ubut = eye(parameters%nparm, parameters%nparm) - &
               lapack_matmul('N', 'T', fun_stx_gemv(parameters, ritz_tmp), ritz_tmp)

        uubt = eye(parameters%nparm, parameters%nparm) - &
               lapack_matmul('N', 'T', ritz_tmp, fun_stx_gemv(parameters, ritz_tmp))

        ys = lapack_matmul('N', 'N', lambda, fun_stx_gemv(parameters, uubt))

        F_matrix = lapack_matmul('N', 'N', ubut, fun_mtx_gemv(parameters, uubt)) - &
                   lapack_matmul('N', 'N', ubut, ys)

    end function fun_F_matrix

end module davidson
