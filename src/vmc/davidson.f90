!! The current implementation uses a general  davidson algorithm, meaning
!! that it compute all the eigenvalues simultaneusly using a variable size block approach.
!! The family of Davidson algorithm only differ in the way that the correction
!! vector is computed.
!! Computed pairs of eigenvalues/eigenvectors are deflated using algorithm
!! described at: https://doi.org/10.1023/A:101919970
!!
!! @Authors: Felipe Zapata and revised by P. Lopez-Tarifa NLeSC(2019)
!>
!> brief Solves the Davidson diagonalisation without storing matrices in memory.
!>
!> Matrices mtx and stx are calculated on the fly using the fun_mtx_gemv and fun_stx_gemv
!> functions.
!>
!> @author Felipe Zapata and P. Lopez-Tarifa NLeSC(2019)
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
!> return eigenvalues and ritz_vectors of the matrix.
module davidson
      use array_utils, only: check_deallocate_matrix
      use array_utils, only: check_deallocate_vector,concatenate
      use array_utils, only: diag_mat,eye,initialize_subspace
      use array_utils, only: modified_gram_schmidt,norm,write_matrix
      use array_utils, only: write_vector
      use contrl_file, only: errunit,ounit
      use lapack_wrapper, only: lapack_generalized_eigensolver
      use lapack_wrapper, only: lapack_matmul,lapack_matrix_vector
      use lapack_wrapper, only: lapack_qr,lapack_solver
      use precision_kinds, only: dp
    implicit none

    type davidson_parameters
        INTEGER :: nparm
        INTEGER :: nparm_max
        INTEGER :: lowest
        INTEGER :: nvecx
        INTEGER :: basis_size
    end type davidson_parameters

contains

    subroutine generalized_eigensolver(eigenvalues, eigenvectors, nparm, nparm_max, &
                                       lowest, nvecx, method, max_iters, tolerance, iters, nproc, idtask)
        !> brief use a pair of functions fun_mtx and fun_stx to compute on the fly the matrices to solve
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
        !> return eigenvalues and ritz_vectors of the matrix

        use mpi
        use store_diag_hs_mod, only: store_diag_hs

        implicit none


        ! input/output variable
        integer, intent(in) :: nparm, nparm_max, nvecx, lowest, nproc, idtask
        real(dp), dimension(lowest), intent(out) :: eigenvalues
        real(dp), dimension(nparm, lowest), intent(out) :: eigenvectors
        integer, intent(in) :: max_iters
        real(dp), intent(in) :: tolerance
        character(len=*), intent(in) :: method
        integer, intent(out) :: iters

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
        real(dp), dimension(:, :), allocatable :: mtx, identity

        ! tmp arrays for the vectorsi residue calculation
        real(dp), dimension(:, :), allocatable :: lambda              ! eigenvalues_sub in a diagonal matrix
        real(dp), dimension(:, :), allocatable :: tmp_res_array       ! tmp array for vectorized res calculation
        real(dp), dimension(:, :), allocatable :: tmp_array           ! tmp array

        ! Arrays dimension
        type(davidson_parameters) :: parameters

        ! Indices of the eigenvalues/eigenvectors pair that have not converged
        logical, dimension(lowest*2) :: has_converged
        logical :: update_proj
        integer :: n_converged ! Number of converged eigenvalue/eigenvector pairs
        integer :: sizeV ! size of V for broadcasting
        logical, parameter :: use_gs_ortho = .true. ! which orthogonalization method to use gs/qr
        logical, parameter :: bool_export_matrices = .false. ! turn true if you want to export matrices
        integer :: not_cnv
        integer :: ii, jj, n 

        ! Iteration subpsace dimension
        init_subspace_size = lowest*2

        ! Allocate our stuff
        allocate(V(nparm, init_subspace_size))

        ! number of correction vectors appended to V at each iteration
        size_update = lowest*2

        ! Make sure nvecx is lower than nparm
        if (nvecx > nparm) then
            call die('DAV: nvecx > nparm, increase nparm or decrease lin_nvecx')
        endif

        ! Dimension of the matrix
        parameters = davidson_parameters(nparm, nparm_max, lowest, nvecx, init_subspace_size)

        ! export the matrices
        if (bool_export_matrices) then
            call export_matrices(parameters)
        endif

        ! Initial number of converged eigenvalue/eigenvector pairs
        n_converged = 0
        has_converged = .false.
        update_proj = .false.

        ! Diagonal of the arrays
        allocate (diag_mtx(parameters%nparm))
        allocate (diag_mtx_cpy(parameters%nparm))
        allocate (diag_stx(parameters%nparm))

        if (idtask == 0) call store_diag_hs(parameters%nparm, diag_mtx, diag_stx)

        ! mtx is incorrect so use a identity multiplication to get the diagonal
        allocate (mtx(parameters%nparm, parameters%nparm))
        allocate (identity(parameters%nparm, parameters%nparm))
        call eye(identity)
        call fun_mtx_gemv(parameters, identity, mtx)
        do i=1,nparm
          diag_mtx( i)=mtx(i,i)
        enddo
        deallocate (mtx, identity)

        if (nproc > 1) then
            call MPI_BCAST(diag_mtx, parameters%nparm, MPI_REAL8, 0, MPI_COMM_WORLD, ier)
            call MPI_BCAST(diag_stx, parameters%nparm, MPI_REAL8, 0, MPI_COMM_WORLD, ier)
        endif

        ! Select the initial ortogonal subspace based on lowest elements
        ! of the diagonal of the matrix.
        diag_mtx_cpy = diag_mtx
        call initialize_subspace(diag_mtx_cpy(1:init_subspace_size), init_subspace_size, nparm, V) ! Initial orthonormal basis
        deallocate (diag_mtx_cpy)

        if (idtask == 0) write(ounit, '(''DAV: Setup subspace problem'')')

        ! allocate mtxV and stxV
        allocate (mtxV(parameters%nparm, parameters%basis_size))
        allocate (stxV(parameters%nparm, parameters%basis_size))

        ! Calculation of HV and SV
        ! after the call only
        ! the master has the correct matrix
        ! but only the master needs it
        call fun_mtx_gemv(parameters, V, mtxV)
        call fun_stx_gemv(parameters, V, stxV)

        if (idtask == 0) then

            select case (method)
            case ("DPR")
                write(ounit, '(''DAV: Diagonal-Preconditioned-Residue (DPR)'')')
            case ("GJD")
                write(ounit, '(''DAV: Generalized Jacobi-Davidson (GJD)'')')
            end select

            write(ounit, '(''DAV: tolerance         : '', E10.3)') tolerance
            write(ounit, '(''DAV: Number eigenvalues : '', I10)') parameters%lowest
            write(ounit, '(''DAV: Update size        : '', I10)') size_update
            write(ounit, '(''DAV: Max basis size     : '', I10)') nvecx
            if (use_gs_ortho) then
                write(ounit, '(''DAV: Modified Gram-Schmidt orthogonalization with projection update'')')
            else
                write(ounit, '(''DAV: QR orthogonalization with full projection'')')
            end if
        endif

        ! 3. Outer loop block Davidson
        outer: do i = 1, max_iters

            ! do most of the calculation on the master
            ! we could try to use openMP here via the lapack routines
            if (idtask == 0) then

                write(ounit, '(''DAV: -----------------------------'')')
                write(ounit, '(''DAV: Iteration: '', I10)') i
                write(ounit, '(''DAV: Basis size: '', I10)') parameters%basis_size

                ! needed if we want to vectorix the residue calculation
                if (allocated(lambda)) deallocate(lambda)
                allocate(lambda(size_update,size_update))

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


                ! update the projected matrices in the small subspace
                if (update_proj) then

                    ! update the projected matrices
                    call update_projection(V, mtxV, mtx_proj)
                    call update_projection(V, stxV, stx_proj)

                    ! recompute it from scratch when restarting
                else

                    ! Array deallocation/allocation.
                    if (allocated(mtx_proj)) deallocate(mtx_proj)
                    allocate(mtx_proj(size(V,2),size(mtxV,2)))
                    if (allocated(stx_proj)) deallocate(stx_proj)
                    allocate(stx_proj(size(V,2),size(stxV,2)))

                    ! recompute the projected matrix
                    call lapack_matmul('T', 'N', V, mtxV, mtx_proj)
                    call lapack_matmul('T', 'N', V, stxV, stx_proj)

                end if

                ! Solve theV small eigenvalue problem
                call lapack_generalized_eigensolver(mtx_proj, eigenvalues_sub, eigenvectors_sub, stx_proj)
                write(ounit, '(''DAV: eigv'',1000d12.5)') (eigenvalues_sub(j), j=1, parameters%lowest)

                ! Compute the necessary ritz vectors
                if (allocated(ritz_vectors)) deallocate(ritz_vectors)
                allocate(ritz_vectors(size(V,1),size_update))
                call lapack_matmul('N', 'N', V, eigenvectors_sub(:, :size_update), ritz_vectors)

                ! Residue calculation (vectorized)
                call diag_mat(eigenvalues_sub(:size_update), lambda)
                if (allocated(tmp_array)) deallocate(tmp_array)
                allocate(tmp_array(size(stxV,1),size_update))
                if (allocated(tmp_res_array)) deallocate(tmp_res_array)
                allocate(tmp_res_array(size(stxV,1),size_update))

                call lapack_matmul('N', 'N', stxV, eigenvectors_sub(:, :size_update),tmp_array)
                call lapack_matmul('N', 'N', tmp_array, lambda, tmp_res_array)
                call lapack_matmul('N', 'N', mtxV, eigenvectors_sub(:, :size_update),residues)
                residues = residues - tmp_res_array

                ! Check which eigenvalues has converged
                errors = norm2(residues(:, :parameters%lowest), 1)
                not_cnv = 0
                do j = 1, parameters%lowest
                    if (errors(j) < tolerance) then
                        has_converged(j) = .true.
                    else
                        not_cnv = not_cnv + 1 
                    endif
                end do
                write(ounit, '(''DAV: resd'',1000d12.5)') (errors(j), j=1, parameters%lowest)
                write(ounit, '(''DAV: Root not yet converged     : '', I10)') not_cnv

                ! exit the loop if all eigenvalues have been found
                if(not_cnv == 0)    exit

                ! Append correction vectors
                if (parameters%basis_size + size_update <= nvecx) then

                    ! compute the correction vectors
                    select case (method)
                    case ("DPR")
                        call compute_DPR(residues, parameters, eigenvalues_sub, diag_mtx, diag_stx, has_converged, correction)
                    case ("GJD")
                        call compute_GJD_free(parameters, ritz_vectors, residues, eigenvectors_sub, eigenvalues_sub, correction)
                    end select

                    ! Add the correction vectors to the current basis.

                    call concatenate(V, correction(:, :not_cnv))

                    ! Orthogonalize basis using modified GS
                    if (use_gs_ortho) then
                        call modified_gram_schmidt(V, parameters%basis_size+1)
                        update_proj = .true.
                    else
                        call lapack_qr(V)
                        update_proj = .false.
                    end if

                    ! Restart.
                else

                    write(ounit, '(''DAV: --- Restart ---'')')
                    V = ritz_vectors(:, :init_subspace_size)
                    do n=1,size(V,2) ! normalize V
                      V(:,n) = V(:,n)/sqrt(sum(V(:,n)**2))
                    end do
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
                    write(ounit, '(''DAV: roots are converged'')')
                endif
                exit outer
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

                if (allocated(tmp_res_array)) deallocate(tmp_res_array)
                allocate(tmp_res_array(size(V,1),size(V,2)-parameters%basis_size))
                call fun_mtx_gemv(parameters, V(:, parameters%basis_size + 1:), tmp_res_array)
                call concatenate(mtxV, tmp_res_array)

                if (allocated(tmp_res_array)) deallocate(tmp_res_array)
                allocate(tmp_res_array(size(V,1),size(V,2)-parameters%basis_size))
                call fun_stx_gemv(parameters, V(:, parameters%basis_size + 1:), tmp_res_array)
                call concatenate(stxV, tmp_res_array)

                ! recompute mtxV and stxV
            else

                if (allocated(mtxV)) deallocate(mtxV)
                allocate(mtxV(size(V,1),size(V,2)))
                call fun_mtx_gemv(parameters, V, mtxV)

                if (allocated(stxV)) deallocate(stxV)
                allocate(stxV(size(V,1),size(V,2)))
                call fun_stx_gemv(parameters, V, stxV)

            end if

            ! Update basis size
            parameters%basis_size = size(V, 2)

        end do outer

        ! Master store eigenpairs
        if (idtask == 0) then

            ! if we didnt converge
            if (i > max_iters) then
                write(ounit, '(''DAV: Warning Davidson did not converge'')')
            else
                write(ounit, '(''DAV: Davidson converged'')')
            end if

            ! store eigenpairs
            eigenvalues = eigenvalues_sub(:parameters%lowest)
            eigenvectors = ritz_vectors(:, :parameters%lowest)
            iters = i

        end if

        ! Broadcast solutions to all slaves
        if (nproc > 1) then

            if (idtask == 0) then
                write(ounit, '(''DAV: Broadcasting solutions'')')
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
            deallocate (residues)
            deallocate (lambda, tmp_res_array)
            deallocate (mtx_proj)
            call check_deallocate_matrix(stx_proj)
        endif

        deallocate (diag_mtx, diag_stx)

        if (idtask == 0) then
            write(ounit, '(''DAV: Exiting Davidson'')')
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
        call lapack_matmul('T', 'N', V, mtxV(:, old_dim + 1:), tmp_array(:, old_dim + 1:))
        tmp_array(old_dim + 1:, :old_dim) = transpose(tmp_array(:old_dim, old_dim + 1:))

        ! Move to new expanded matrix
        call move_alloc(tmp_array, mtx_proj)

    end subroutine update_projection

    subroutine die(msg)
        !> Subroutine that dies the calculation raising an errror message
        !
        use mpi
        character msg*(*)
        integer ierr

        write(ounit, '(''Fatal error: '',a)') msg
        call mpi_abort(MPI_COMM_WORLD, 0, ierr)

    end subroutine

    subroutine compute_DPR(residues, parameters, eigenvalues, diag_mtx, diag_stx, has_converged, correction)
        !> compute the correction vector using the DPR method for a matrix free diagonalization
        !> See correction_methods submodule for the implementations
        !> \param[in] fun_mtx: function to compute matrix
        !> \param[in] fun_stx: function to compute the matrix for the generalized case
        !> \param[in] V: Basis of the iteration subspace
        !> \param[in] eigenvalues: of the reduce problem
        !> \param[in] eigenvectors: of the reduce problem
        !> return correction matrix
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

    end subroutine compute_DPR

    subroutine compute_GJD_free(parameters, ritz_vectors, residues, eigenvectors, &
                              eigenvalues, correction)

        !> Compute the correction vector using the GJD method for a matrix free
        !> diagonalization. We follow the notation of:
        !> I. Sabzevari, A. Mahajan and S. Sharma,  arXiv:1908.04423 (2019)
        !>
        !> \param[in] ritz_vectors: ritz_vectors.
        !> \param[in] residues: residue vectors.
        !> \param[in] parameter: davidson_parameters type.
        !> \param[in] eigenvectors.
        !> \param[in] eigenvalues.
        !> return correction matrix

        use array_utils, only: eye

        type(davidson_parameters)               :: parameters
        real(dp), dimension(:, :), intent(in) :: ritz_vectors
        real(dp), dimension(:, :), intent(in) :: residues
        real(dp), dimension(:, :), intent(in) :: eigenvectors
        real(dp), dimension(:), intent(in) :: eigenvalues
        !
        ! local variables
        !
        real(dp), dimension(:, :) :: correction
        integer :: k, m
        logical :: gev
        real(dp), dimension(parameters%nparm, parameters%nparm) ::  F
        real(dp), dimension(parameters%nparm, 1) :: brr

        do k = 1, parameters%basis_size

            call fun_F_matrix(ritz_vectors, parameters, k, eigenvalues(k), F)
            ! call write_matrix("F.txt", F)
            brr(:, 1) = -residues(:, k)
            call lapack_solver(F, brr)
            ! call write_matrix("brr.txt", brr)
            correction(:, k) = brr(:, 1)

        end do

    end subroutine compute_GJD_free

    subroutine fun_F_matrix(ritz_vectors, parameters, eigen_index, eigenvalue, F_matrix)
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
        real(dp), dimension(:, :) :: F_matrix
        type(davidson_parameters) :: parameters
        integer   :: eigen_index
        real(dp) :: eigenvalue

        real(dp), dimension(parameters%nparm, parameters%nparm) :: lambda
        real(dp), dimension(parameters%nparm, 1) :: ritz_tmp, stx
        real(dp), dimension(parameters%nparm, parameters%nparm) :: I, ubut, uubt, tmp, ys, mtx

        ritz_tmp(:, 1) = ritz_vectors(:, eigen_index)

        call eye(lambda, eigenvalue)

        call eye(I)
        call fun_stx_gemv(parameters, ritz_tmp, stx)
        call lapack_matmul('N', 'T', stx, ritz_tmp, ubut)
        ubut = I - ubut

        call lapack_matmul('N', 'T', ritz_tmp, stx, uubt)
        uubt = I - uubt
        
        call fun_stx_gemv(parameters, uubt, tmp)
        call lapack_matmul('N', 'N', lambda, tmp, ys)

        call fun_mtx_gemv(parameters, uubt, mtx)
        call lapack_matmul('N', 'N', ubut, mtx, F_matrix)
        call lapack_matmul('N', 'N', ubut, ys, tmp)
        F_matrix = F_matrix - tmp

    end subroutine fun_F_matrix

    subroutine export_matrices(parameters)
        !> export the H and S matrix to file for comparison with numpy
        !> \param[in] parameters: structure containing the davidson parameters

        use array_utils, only: write_matrix, eye
        type(davidson_parameters)               :: parameters
        real(dp), dimension(:, :), allocatable :: I
        real(dp), dimension(:, :), allocatable :: mtx
        real(dp), dimension(:, :), allocatable :: stx
        integer n, m

        allocate (mtx(parameters%nparm, parameters%nparm))
        allocate (stx(parameters%nparm, parameters%nparm))
        allocate (I(parameters%nparm, parameters%nparm))

        call eye(I)

        call fun_mtx_gemv(parameters, I, mtx)
        call fun_stx_gemv(parameters, I, stx)

        write(*,*) "Full STX"
        do n = 1,parameters%nparm
          do m = 1,parameters%nparm
            write(*,'("  ",E10.4)',advance='no') stx(m,n)
          end do
          write(*,*)
        end do
        write(*,*) "Full MTX"
        do n = 1,parameters%nparm
          do m = 1,parameters%nparm
            write(*,'("  ",E10.4)',advance='no') mtx(m,n)
          end do
          write(*,*)
        end do
        ! call write_matrix("H.dat", mtx)
        ! call write_matrix("S.dat", stx)

        deallocate (mtx)
        deallocate (stx)
        deallocate (I)

    end subroutine export_matrices

    subroutine fun_mtx_gemv(parameters, input_vect, output_vect)
        !> Brief Function to compute the action of the hamiltonian on the fly
        !> \param[in] dimension of the arrays to compute the action of the hamiltonian
        !> \param[in] input_vec Array to project
        !> return Projected matrix

        use precision_kinds, only: dp
        use optwf_lin_dav_extra,  only: h_psi_lin_d

        type(davidson_parameters) :: parameters
        real(dp), dimension(:, :), intent(in) :: input_vect
        real(dp), dimension(:, :) :: output_vect
        real(dp), dimension(:, :), allocatable :: psi, hpsi

        allocate (psi(parameters%nparm_max, 2*parameters%nvecx), source=0.0_dp)
        allocate (hpsi(parameters%nparm_max, 2*parameters%nvecx), source=0.0_dp)

        psi = 0.0_dp
        psi(1:size(input_vect, 1), 1:size(input_vect, 2)) = input_vect
        write(*,*) size(input_vect, 2)
        call h_psi_lin_d(parameters%nparm, size(input_vect, 2), psi, hpsi)

        output_vect = hpsi(1:size(input_vect, 1), 1:size(input_vect, 2))
        deallocate (psi, hpsi)

    end subroutine fun_mtx_gemv

    subroutine fun_stx_gemv(parameters, input_vect, output_vect)
        !> Brief Fucntion to compute the optional stx matrix on the fly
        !> \param[in] dimension of the arrays to compute the action of the hamiltonian
        !> \param[in] input_vec Array to project
        !> return Projected matrix

        use precision_kinds, only: dp
        use optwf_lin_dav_extra,  only: s_psi_lin_d

        type(davidson_parameters) :: parameters
        real(dp), dimension(:, :), intent(in) :: input_vect
        real(dp), dimension(:, :) :: output_vect
        real(dp), dimension(:, :), allocatable :: psi, spsi

        allocate (psi(parameters%nparm_max, 2*parameters%nvecx), source=0.0_dp)
        allocate (spsi(parameters%nparm_max, 2*parameters%nvecx), source=0.0_dp)

        psi = 0.0_dp
        psi(1:size(input_vect, 1), 1:size(input_vect, 2)) = input_vect

        call s_psi_lin_d(parameters%nparm, size(input_vect, 2), psi, spsi)

        output_vect = spsi(1:size(input_vect, 1), 1:size(input_vect, 2))
        ! output_vect = input_vect
        deallocate (psi, spsi)

    end subroutine fun_stx_gemv

end module davidson
