!> Utility module for array manipulation
!> @author Viktor Azizi
!> @author Alice Cuzzocrea
!> @author Ravindra Shinde
!> @author Nico Renaud
!> @date 24/05/2022
module array_utils

      use lapack_wrapper, only: lapack_matrix_vector,lapack_sort
      use precision_kinds, only: dp
    implicit none

    !> \private
    private
    !> \public
    public :: concatenate, diagonal, eye, generate_diagonal_dominant, norm, &
              initialize_subspace, write_matrix, write_vector, check_deallocate_vector, &
              check_deallocate_matrix, modified_gram_schmidt, diag_mat
    public :: unique_elements, unique_string_elements

contains

    !> Create an identity matrix
    pure function eye(m, n, alpha)
        !> \param m: number of rows
        !> \param n: number of colums
        !> \param alpha: optional diagonal value
        !>
        !> return matrix of size n x m
        integer, intent(in) :: n, m
        real(dp), dimension(m, n) :: eye
        real(dp), intent(in), optional :: alpha

        !local variable
        integer :: i, j, min_dim
        real(dp) :: x

        ! check optional values
        x = 1.d0
        if (present(alpha)) x = alpha

        ! that's sooooo costly
        do i = 1, m
            do j = 1, n
                if (i /= j) then
                    eye(i, j) = 0.d0
                else
                    eye(i, i) = x
                end if
            end do
        end do

        ! why not that ?
        !  eye = 0.0_dp
        !  min_dim = min(m,n)
        !  do i=1, min_dim
        !    eye(i,i) = x
        !  end do

    end function eye

    !> Create a diagnoal matrix from a vector
    pure function diag_mat(vec)
        real(dp), dimension(:), intent(in) :: vec
        real(dp), dimension(size(vec, 1), size(vec, 1)) :: diag_mat

        integer :: i
        diag_mat = 0.0_dp
        do i = 1, size(vec, 1)
            diag_mat(i, i) = vec(i)
        end do

    end function diag_mat

    !> compute the norm-2 of a vector
    pure function norm(vector)
        real(dp), dimension(:), intent(in) :: vector
        real(dp) :: norm

        norm = sqrt(sum(vector**2.d0))

    end function norm

    !> Concatenate two matrices
    subroutine concatenate(arr, brr)
        !> \param arr: first array
        !> \param brr: second array
        !>
        !> return arr concatenate brr(overwrites arr)

        real(dp), dimension(:, :), intent(inout), allocatable :: arr
        real(dp), dimension(:, :), intent(in) :: brr
        real(dp), dimension(:, :), allocatable :: tmp_array
        integer :: new_dim, dim_cols, dim_rows

        ! dimension
        dim_rows = size(arr, 1)
        dim_cols = size(arr, 2)

        ! Number of columns of the new matrix
        new_dim = dim_cols + size(brr, 2)

        ! move to temporal array
        allocate (tmp_array(dim_rows, new_dim))
        tmp_array(:, :dim_cols) = arr

        ! Move to new expanded matrix
        deallocate (arr)
        call move_alloc(tmp_array, arr)

        arr(:, dim_cols + 1:) = brr

    end subroutine concatenate

    !> Generate a diagonal dominant square matrix of dimension m
    function generate_diagonal_dominant(m, sparsity, diag_val) result(arr)
        !> \param m dimension of the matrix
        !> \param sparsity magnitude order of the off-diagonal values

        integer, intent(in) :: m ! size of the square matrix
        real(dp), optional :: diag_val
        integer :: i, j
        real(dp) :: sparsity
        real(dp), dimension(m, m) :: arr
        call random_number(arr)

        arr = arr*sparsity
        do j = 1, m
            do i = 1, m
                if (i > j) then
                    arr(i, j) = arr(j, i)
                else if (i == j) then
                    if (present(diag_val)) then
                        arr(i, i) = diag_val
                    else
                        arr(i, i) = i
                    end if
                end if
            end do
        end do

    end function generate_diagonal_dominant

    !> return the diagonal of a matrix
    function diagonal(matrix)
        !> \param matrix: input matrix
        real(dp), dimension(:, :), intent(in) :: matrix
        real(dp), dimension(size(matrix, 1)) :: diagonal

        ! local variables
        integer :: i, j, m

        ! dimension of the matrix
        m = size(matrix, 1)

        do i = 1, m
            do j = 1, m
                if (i == j) then
                    diagonal(i) = matrix(i, j)
                end if
            end do
        end do

    end function diagonal

    !> Brief generates a diagonal preconditioner for .
    function initialize_subspace(diag, dim_sub, dim_base) result(precond)
        !> return diagonal matrix

        ! input variable
        real(dp), dimension(:), intent(inout) :: diag
        integer, intent(in) :: dim_sub, dim_base

        ! local variables
        real(dp), dimension(dim_base, dim_sub) :: precond
        integer, dimension(size(diag)) :: keys
        integer :: i, k

        ! sort diagonal
        keys = lapack_sort('I', diag)
        ! Fill matrix with zeros
        precond = 0.0_dp

        ! Add one depending on the order of the matrix diagonal
        do i = 1, dim_sub
            k = search_key(keys, i)
            precond(k, i) = 1.d0
        end do

    end function initialize_subspace

    !> Brief use modifed gram-schmidt orthogonalization on mat
    !> Brief nstart is the index of the first vector to orthogonalize
    subroutine modified_gram_schmidt(mat, nstart)
        ! input
        real(dp), dimension(:, :), intent(inout) :: mat
        integer, optional, intent(in) :: nstart

        integer :: i
        integer :: nrows, ncols
        integer :: idx_start
        real(dp), dimension(:), allocatable :: tmp_array

        idx_start = 1
        if (present(nstart)) idx_start = nstart

        nrows = size(mat, 1)
        ncols = size(mat, 2)

        allocate (tmp_array(nrows))

        do i = idx_start, ncols

            tmp_array = 0.0_dp
            tmp_array = lapack_matrix_vector('T', mat(:, :i - 1), mat(:, i))
            mat(:, i) = mat(:, i) - lapack_matrix_vector('N', mat(:, :i - 1), tmp_array)
            mat(:, i) = mat(:, i)/norm(mat(:, i))

        end do

        deallocate (tmp_array)

    end subroutine modified_gram_schmidt

    !> Brief Search for a given index  in a vector
    function search_key(keys, i) result(k)
        !> \param keys Vector of index
        !> \param i Index to search for
        !>
        !> return index of i inside keys

        integer, dimension(:), intent(in) :: keys
        integer, intent(in) :: i
        integer :: j, k

        k = -1
        do j = 1, size(keys)
            if (keys(j) == i) then
                k = j
                exit
            end if
        end do

    end function search_key

    !> Write matrix to path_file
    subroutine write_matrix(path_file, mtx)
        character(len=*), intent(in) :: path_file
        real(dp), dimension(:, :), intent(in) :: mtx
        integer :: i, j

        open (unit=314, file=path_file, status="REPLACE")
        do i = 1, size(mtx, 1)
            do j = 1, size(mtx, 2)
                write (314, *) mtx(i, j)
            end do
        end do
        close (314)

    end subroutine write_matrix

    !> Write vector to path_file
    subroutine write_vector(path_file, vector)
        character(len=*), intent(in) :: path_file
        real(dp), dimension(:), intent(in) :: vector
        integer :: i

        open (unit=314, file=path_file, status="REPLACE")
        do i = 1, size(vector)
            write (314, *) vector(i)
        end do
        close (314)

    end subroutine write_vector

    !> deallocate a matrix if allocated
    subroutine check_deallocate_matrix(mtx)

        real(dp), dimension(:, :), allocatable, intent(inout) ::  mtx

        if (allocated(mtx)) then
            deallocate (mtx)
        end if

    end subroutine check_deallocate_matrix

    !> deallocate a vector if allocated
    subroutine check_deallocate_vector(vec)


        real(dp), dimension(:), allocatable, intent(inout) ::  vec

        if (allocated(vec)) then
            deallocate (vec)
        end if

    end subroutine check_deallocate_vector

    !> Returns the unique elements of a vector
    !> Also returns the number of unique elements, their indices in the original vector.
    subroutine unique_elements(n, arr, res, count, frequency, ind)
        !> \param[in] n: size of the vector
        !> \param[in] arr: vector to be processed
        !> \param[out] res: vector of unique elements
        !> \param[out] count: number of unique elements
        !> \param[out] frequency: frequency of unique elements
        !> \param[out] ind: vector of indices of the unique elements in the original vector
        !> @author: Ravindra Shinde
        !> @email: r.l.shinde@utwente.nl
        !> @date:   24/05/2022

        implicit none
        integer, intent(in)                         :: n
        integer, dimension(n), intent(in)           :: arr    ! The input
        integer, dimension(n), intent(out)          :: res    ! The output
        integer, intent(out)                        :: count                   ! The number of unique elements
        integer, dimension(n), intent(out)          :: frequency ! The output
        integer, dimension(n), intent(out),optional :: ind
        integer                                     :: i,j,k, counter1, counter2

        k = 1
        ind(1) = 1
        frequency(1) = 0
        res(1) = arr(1);
        counter1 = 1

        outer: do i=1,n
            do j=1,counter1
                if (res(j) == arr(i)) then
                    frequency(j) = frequency(j) + 1
                    ind(j+1) = i + 1
                    k = k + 1
                    cycle outer
                end if
            end do
            res(counter1+1) = arr(i)
            frequency(counter1+1) = 1
            ind(k) = i
            k = k + 1
            counter1 = counter1 + 1
        end do outer

        count = counter1
        do i = count + 1, n
            res(i) = 0
            frequency(i) = 0
            ind(i) = 0
        end do

    end subroutine unique_elements

    !> Returns the unique elements of a vector
    !> Also returns the number of unique elements, their indices in the original vector.
    subroutine unique_string_elements(n, arr, res, count)
        !> \param[in] n: size of the vector
        !> \param[in] arr: vector to be processed
        !> \param[out] res: vector of unique elements
        !> \param[out] count: number of unique elements
        !> \param[out] frequency: frequency of unique elements
        !> \param[out] ind: vector of indices of the unique elements in the original vector
        !> @author: Ravindra Shinde
        !> @email: r.l.shinde@utwente.nl
        !> @date:   24/05/2022

        implicit none
        integer, intent(in)                         :: n
        character(len=3), dimension(n), intent(in)  :: arr    ! The input
        character(len=3), dimension(n), intent(out) :: res    ! The output
        integer, intent(out)                        :: count                   ! The number of unique elements
        integer, dimension(n)          :: frequency ! The output
        integer, dimension(n)          :: ind
        integer                                     :: i,j,k, counter1, counter2

        k = 1
        ind(1) = 1
        frequency(1) = 0
        res(1) = arr(1)
        counter1 = 1


        outer: do i=1,n
            do j=1,counter1
                if (res(j) == arr(i)) then
                    frequency(j) = frequency(j) + 1
                    ind(j+1) = i + 1
                    k = k + 1
                    cycle outer
                end if
            end do
            res(counter1+1) = arr(i)
            frequency(counter1+1) = 1
            ind(k) = i
            k = k + 1
            counter1 = counter1 + 1
        end do outer

        count = counter1
        do i = count + 1, n
            res(i) = ""
            frequency(i) = 0
            ind(i) = 0
        end do

    end subroutine unique_string_elements


    !> find "indices", the list of unique numbers in "list"
    !> and return the sorted unique numbers in "sorted"
    subroutine sortedunique(list, n, indices, sorted)
        !> \param[in] list: input list
        !> \param[in] n: size of the list
        !> \param[out] indices: indices of the unique elements
        !> \param[out] sorted: sorted unique elements
        !> @author: Ravindra Shinde
        !> @date: 24/05/2022
        implicit none
        !   find "indices", the list of unique numbers in "list"

        integer :: n, kx, i, nitems
        integer, dimension(n) :: list
        logical, dimension(n) :: mask
        integer, dimension(:), allocatable :: indices, sorted

        mask(1)=.true.

        do kx=n,2,-1
            mask(kx)= .not.(any(list(:kx-1)==list(kx)))
        end do

        if (allocated(indices)) then
            indices=pack([(kx,kx=1,n)],mask)
        else
            allocate(indices, source=pack( [(kx,kx=1,n)],mask) )
        endif

        nitems = size(indices)
        if (allocated(sorted)) then
            do i=1,nitems
                sorted(i) = list(indices(i))
            enddo
        else
            allocate(sorted(nitems))
            do i=1,nitems
                sorted(i) = list(indices(i))
            enddo
        endif
    end subroutine sortedunique


end module array_utils


