module lapack_wrapper

      use precision_kinds, only: dp
  implicit none
  private
  interface!LAPACK interface
!*  -- LAPACK driver routine (version 3.1) --
     SUBROUTINE DSYSV( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, WORK, LWORK, INFO )
        CHARACTER          UPLO
        INTEGER            INFO, LDA, LDB, LWORK, N, NRHS
        INTEGER            IPIV( * )
        DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), WORK( * ) 
      END SUBROUTINE
      SUBROUTINE DGEMV(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
        DOUBLE PRECISION ALPHA,BETA
        INTEGER INCX,INCY,LDA,M,N
        CHARACTER TRANS
        DOUBLE PRECISION A(LDA,*),X(*),Y(*)
      END SUBROUTINE
      SUBROUTINE DGEQRF( M, N, A, LDA, TAU, WORK, LWORK, INFO )
        INTEGER            INFO, LDA, LWORK, M, N
        DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
      END SUBROUTINE
      SUBROUTINE DSYEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO )
        CHARACTER          JOBZ, UPLO
        INTEGER            INFO, LDA, LWORK, N
        DOUBLE PRECISION   A( LDA, * ), W( * ), WORK( * )
      END SUBROUTINE
      SUBROUTINE DSYGV( ITYPE, JOBZ, UPLO, N, A, LDA, B, LDB, W, WORK, LWORK, INFO )
        CHARACTER          JOBZ, UPLO
        INTEGER            INFO, ITYPE, LDA, LDB, LWORK, N
        DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), W( * ), WORK( * )
      END SUBROUTINE
      SUBROUTINE DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
        DOUBLE PRECISION ALPHA,BETA
        INTEGER K,LDA,LDB,LDC,M,N
        CHARACTER TRANSA,TRANSB
        DOUBLE PRECISION A(LDA,*),B(LDB,*),C(LDC,*)
      END SUBROUTINE
      SUBROUTINE DLASRT( ID, N, D, INFO )
        CHARACTER          ID
        INTEGER            INFO, N
        DOUBLE PRECISION   D( * )
      END SUBROUTINE
      SUBROUTINE DSYGVX( ITYPE, JOBZ, RANGE, UPLO, N, A, LDA, B, LDB, &
                         VL, VU, IL, IU, ABSTOL, M, W, Z, LDZ, WORK,  &
                         LWORK, IWORK, IFAIL, INFO )
        CHARACTER          JOBZ, RANGE, UPLO
        INTEGER            IL, INFO, ITYPE, IU, LDA, LDB, LDZ, LWORK, M, N
        DOUBLE PRECISION   ABSTOL, VL, VU
        INTEGER            IFAIL( * ), IWORK( * )
        DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), W( * ), WORK( * ), Z( LDZ, * )
      END SUBROUTINE
      SUBROUTINE DORGQR( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
        INTEGER            INFO, K, LDA, LWORK, M, N
        DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
      END SUBROUTINE
  end interface

  public :: lapack_generalized_eigensolver, lapack_generalized_eigensolver_lowest, &
       lapack_matmul, lapack_matrix_vector, lapack_qr, lapack_solver, lapack_sort

contains

    subroutine lapack_generalized_eigensolver(mtx, eigenvalues, eigenvectors, stx)
    !> Call the DSYGV subroutine lapack to compute ALL the eigenvalues
    !> and corresponding eigenvectors of mtx
    !> \param mtx: Matrix to diaogonalize
    !> \param stx: Overlap Matrix to diaogonalize
    !> \param eigenvalues: lowest eigenvalues
    !> \param eigenvectors: corresponding eigenvectors
    !> \return eigenvalues/eigenvectors

    ! input/output
    implicit none
    real(dp), dimension(:, :), intent(in) :: mtx
    real(dp), dimension(:, :), intent(in), optional :: stx
    real(dp), dimension(:), intent(inout) :: eigenvalues
    real(dp), dimension(:, :), intent(inout) :: eigenvectors

    real(dp), dimension(:, :), allocatable :: mtx_copy
    real(dp), dimension(:, :), allocatable :: stx_copy


    ! Local variables
    integer :: dim, info, lwork, itype = 1, n
    logical :: gev

     ! ALL the eigenvalues of the subpace (re, im)
    real(dp), dimension(size(mtx, 1)) :: eigenvalues_work
    real(dp), dimension(:), allocatable :: work ! workspace, see lapack documentation

    ! ! dimension of the guess space
    dim = size(mtx, 1)

    gev = present(stx)

    ! local copy of the matrices
    allocate(mtx_copy(dim,dim))
    mtx_copy = mtx

    if (gev) then
      allocate(stx_copy(dim,dim))
      stx_copy = stx
    end if

    ! Query size of the optimal workspace
    allocate(work(1))

    if (gev) then
      call DSYGV(itype,"V", "U", dim, mtx_copy, dim, stx_copy, dim, eigenvalues_work, work, -1, info)
      call check_lapack_call(info, "DSYGV")
    else
      call DSYEV("V", "U", dim, mtx_copy, dim, eigenvalues_work, work, -1, info)
      call check_lapack_call(info, "DSYEV")
    end if

    ! Allocate memory for the workspace
    lwork = max(1, int(work(1)))
    deallocate(work)
    allocate(work(lwork))

    ! Compute Eigenvalues
    if (gev) then
      call DSYGV(itype,"V", "U", dim, mtx_copy, dim, stx_copy, dim, eigenvalues_work, work, lwork, info)
      call check_lapack_call(info, "DSYGV")
    else
      call DSYEV("V", "U", dim, mtx_copy, dim, eigenvalues_work, work, lwork, info)
      call check_lapack_call(info, "DSYEV")
    end if

    ! Sort the eigenvalues and eigenvectors of the basis
    eigenvalues = eigenvalues_work
    eigenvectors = mtx_copy

    ! release memory
    deallocate(work)
    deallocate(mtx_copy)
    if (gev) then
      deallocate(stx_copy)
    end if

  end subroutine lapack_generalized_eigensolver

    subroutine lapack_generalized_eigensolver_lowest(mtx, stx, eigenvalues, eigenvectors, lowest)
    !> Call the DSYGVX subroutine lapack to compute the lowest eigenvalues
    !> and corresponding eigenvectors of mtx
    !> \param mtx: Matrix to diaogonalize
    !> \param stx: Overlap Matrix to diaogonalize
    !> \param eigenvalues: lowest eigenvalues
    !> \param eigenvectors: corresponding eigenvectors
    !> \param lowest: number of lowest eigenvalues/eigenvectors pairs to compute
    !> \return eigenvalues/eigenvectors

    ! input/output
    implicit none
    real(dp), dimension(:, :), intent(in) :: mtx
    real(dp), dimension(:, :), intent(in) :: stx
    real(dp), dimension(:), intent(out) :: eigenvalues
    real(dp), dimension(:, :), intent(out) :: eigenvectors
    integer :: lowest

    ! Local variables
    real(dp), dimension(:, :), allocatable :: mtx_copy
    real(dp), dimension(:, :), allocatable :: stx_copy

    integer :: dim, info,  lwork, m, itype = 1
    real(dp) :: vl, vu, abstol

    ! Failed eigenvalues
    integer, dimension(size(mtx, 1)) :: ifail


     ! ALL the eigenvalues of the subpace (re, im)
    real(dp), dimension(size(mtx, 1)) :: eigenvalues_work
    real(dp), dimension(size(mtx, 1), lowest) :: eigenvectors_work
    real(dp), dimension(:), allocatable :: work ! workspace, see lapack documentation
    integer, dimension(:), allocatable :: iwork ! workspace, see lapack documentation

    ! ! dimension of the guess space
    dim = size(mtx, 1)

    ! local copy of the matrices
    allocate(mtx_copy(dim,dim))
    mtx_copy = mtx

    allocate(stx_copy(dim,dim))
    stx_copy = stx

    ! LAPACK SAYS: If range = 'A' or 'I', vl and vu are not referenced
    vl = 0.0_dp
    vu = 0.0_dp

    ! Absolute tolerance

    ! Query size of the optimal workspace
    allocate(work(1))
    allocate(iwork(1))

    call DSYGVX(itype,"V", "I", "U", dim, mtx_copy, dim, stx_copy, dim, vl, vu, &
         1, lowest, abstol, m, eigenvalues_work, eigenvectors_work, &
         dim, work, -1, iwork, ifail, info)
    call check_lapack_call(info, "DSYGVX")

    ! Allocate memory for the workspace
    lwork = max(1, int(work(1)))
    deallocate(work, iwork)
    allocate(work(lwork))
    allocate(iwork(lwork))

    ! Compute Eigenvalues

    call DSYGVX(itype,"V", "I", "U", dim, mtx_copy, dim, stx_copy, dim, vl, vu, &
         1, lowest, abstol, m, eigenvalues_work, eigenvectors_work, &
         dim, work, lwork, iwork, ifail, info)

    call check_lapack_call(info, "DSYGVX")

    ! Copy the eigenvalues and eigenvectors
    eigenvalues = eigenvalues_work(1:lowest)
    eigenvectors = eigenvectors_work(1:dim, 1:lowest)

    ! release memory
    deallocate(work, iwork, mtx_copy, stx_copy)

  end subroutine lapack_generalized_eigensolver_lowest

  subroutine lapack_qr(basis)
    !> Orthoghonalize the basis using the QR factorization.
    !> QR factorization of the M-by-N (M>N) matrx A=Q*R in the form where
    !> Q is square M-by-M matrix and R is an upper triangular M-by-N matrix.
    !> The equality A=Q*R can be re-written also as a product Q1*R1 where Q1
    !> is a rectangular M-by-N submatrix of the matrix Q and R1 is M-by-M
    !> submatrix of the R. Let us note that columns of Q1 are orthonormal
    !> (they are orthogonal to each other and have norms equal to 1).
    !> The equality A=Q1*R1 can be treated as every column of A is a linear
    !> combination of Q1 columns, i.e. they span the same linear space.
    !> In other words, columns of Q1 is the result of ortogonalization of columns A.
    !> DGEQRF does not not compute Q directly, DORGQR must be call subsequently.

    !> \param basis
    !> \return orthogonal basis

    implicit none
    real(dp), dimension(:, :), contiguous, intent(inout) :: basis
    real(dp), dimension(:), allocatable :: work ! workspace, see lapack documentation
    real(dp), dimension(size(basis, 2)) :: tau ! see DGEQRF documentation
    integer :: info, lwork, m, n

    ! Matrix shape
    m = size(basis, 1)
    n = size(basis, 2)

    ! 1. Call the QR decomposition
    ! 1.1 Query size of the workspace (Check lapack documentation)
    allocate(work(1))
    call DGEQRF(m, n, basis, m, tau, work, -1, info)
    call check_lapack_call(info, "DGEQRF")

    ! 1.2 Allocate memory for the workspace
    lwork = max(1, int(work(1)))
    deallocate(work)
    allocate(work(lwork))

    ! 1.3 Call QR factorization
    call DGEQRF(m, n, basis, m, tau, work, lwork, info)
    call check_lapack_call(info, "DGEQRF")
    deallocate(work)

    ! 2. Generates an orthonormal matrix
    ! 2.1 Query size of the workspace (Check lapack documentation)
    allocate(work(1))
    call DORGQR(m, n, min(m, n), basis, m, tau, work, -1, info)
    call check_lapack_call(info, "DORGQR")

    ! 2.2 Allocate memory fo the workspace
    lwork = max(1, int(work(1)))
    deallocate(work)
    allocate(work(lwork))

    ! 2.3 compute the matrix Q
    call DORGQR(m, n, min(m, n), basis, m, tau, work, lwork, info)
    call check_lapack_call(info, "DORGQR")

    ! release memory
    deallocate(work)

  end subroutine lapack_qr

  subroutine lapack_solver(arr, brr)
    !> Call lapack DSYSV subroutine to solve a AX=B Linear system
    !> \param arr: matrix with the coefficients of the linear system
    !> \param brr: Vector with the constant terms
    !> \returns: Solution vector X (overwriten brr)

    implicit none

    real(dp), dimension(:, :), contiguous, intent(inout) :: arr, brr

    ! local variables
    real(dp), dimension(:), allocatable :: work
    integer :: n, info, lwork
    integer, dimension(size(arr, 1)) :: ipiv

    n = size(arr, 1)

    ! query spacework size
    allocate(work(1))
    call DSYSV("U", n, 1, arr, n, ipiv, brr, n, work, -1, info)
    call check_lapack_call(info, "DSYSV")

    ! Allocate memory fo the workspace
    lwork = max(1, int(work(1)))
    deallocate(work)
    allocate(work(lwork))

    ! run linear solver
    call DSYSV("U", n, 1, arr, n, ipiv, brr, n, work, lwork, info)
    ! If the diagonalization fails due to a singular value try to recover
    ! by replacing the 0 value with a tiny one
    if (info > 0) then
       arr(info, info) = tiny(arr(1, 1))
       call DSYSV("U", n, 1, arr, n, ipiv, brr, n, work, lwork, info)
       call check_lapack_call(info, "DSYSV")
    end if

    deallocate(work)

  end subroutine lapack_solver

    subroutine lapack_matmul(transA, transB, arr, brr, mtx, alpha)
    !> perform the matrix multiplication alpha * arr ^ (transA) * brr ^ (transB)
    !> see Lapack DGEMM for further details
    !> \param transA: 'T' transpose A, 'N' do not tranpose
    !> \param transB: 'T' transpose B, 'N' do not tranpose
    !> \param arr: first matrix to multiply
    !> \param brr: second matrix
    !> \param alpha: optional scalar number
    !> \return matrix multiplication

    implicit none

    character(len=1), intent(in) :: transA, transB
    real(dp), dimension(:, :), contiguous, intent(in) :: arr, brr
    real(dp), optional, intent(in) :: alpha
    real(dp), dimension(:, :), contiguous :: mtx

    ! local variables
    real(dp) :: x
    integer :: m, n, k, lda, ldb
    x = 1.d0

    ! check optional variable
    if (present(alpha)) x=alpha

    if (transA == 'T') then
       k = size(arr, 1)
       m = size(arr, 2)
       lda = k
    else
       k = size(arr, 2)
       m = size(arr, 1)
       lda = m
    end if

    if (transB == 'T') then
       n = size(brr, 1)
       ldb = n
    else
       n = size(brr, 2)
       ldb = k
    end if

    ! resulting array
    mtx = 0.d0

    call DGEMM(transA, transB, m, n, k, x, arr, lda, brr, ldb, 0.d0, mtx, m)

  end subroutine lapack_matmul

  subroutine lapack_matrix_vector(transA, mtx, vector, rs, alpha)
    !> perform the Matrix vector multiplication alpha * mtx ^ (transA) * vector
    !> see DGEMV for details
    !> \param transA: 'T' transpose A; 'N' do not transpose
    !> \param mtx: matrix to multiply
    !> \param vector: vector to multiply
    !> \param alpha: optional scalar value
    !> \return resulting vector

    implicit none

    character(len=1), intent(in) :: transA
    real(dp), dimension(:, :), contiguous, intent(in) :: mtx
    real(dp), dimension(:), contiguous, intent(in) :: vector
    real(dp), optional, intent(in) :: alpha
    real(dp), dimension(:), contiguous :: rs

    ! local variable
    integer :: m, n
    real(dp) :: scalar
    scalar = 1.d0

    ! check optional variable
    if (present(alpha)) scalar=alpha

    ! number of row of mtx
    m = size(mtx, 1)
    n = size(mtx, 2)

    rs = 0.d0

    call DGEMV(transA, m, n, scalar, mtx, m, vector, 1, 0.d0, rs, 1)

  end subroutine lapack_matrix_vector


  subroutine lapack_sort(id, vector,keys)
    !> \brief sort a vector of douboles
    !> \param vector Array to sort
    !> \param keys Indices of the sorted array
    !> \param id Sorted either `I` increasing or `D` decreasing
    real(dp), dimension(:), contiguous, intent(inout) :: vector
    character(len=1), intent(in) :: id
    integer, dimension(:) :: keys

    ! local variables
    real(dp), dimension(size(vector)) :: xs
    integer :: i, j, info
    xs = vector

    call DLASRT(id, size(vector), vector, info)
    call check_lapack_call(info, "DLASRT")

    do i=1,size(vector)
       do j=1, size(vector)
          if (abs(vector(j) - xs(i)) < 1e-16) then
             keys(i) = j
          end if
       end do
    end do

  end subroutine lapack_sort


  subroutine check_lapack_call(info, name)
    !> Check if a subroutine finishes sucessfully
    !> \param info: Termination signal
    !> \param name: Name of the subroutine
    use contrl_file,    	only: iunit, ounit, errunit
    use error,   only: fatal_error
    implicit none
    integer :: info
    character(len=*), intent(in) :: name

    if (info /= 0) then
       print *, "call to subroutine: ", name, " has failed!"
       print *, "info: ", info
       call fatal_error("Stopping in check_lapack_call subroutine")
    end if

  end subroutine check_lapack_call

end module lapack_wrapper
