module array_utils

      use lapack_wrapper, only: lapack_matrix_vector,lapack_sort
      use precision_kinds, only: dp
    implicit none

contains

   pure subroutine eye(diag, alpha)
        !> Create a matrix with ones in the diagonal and zero everywhere else
        !> \param m: number of rows
        !> \param n: number of colums
        !> \param alpha: optional diagonal value
        !>
        !> return matrix of size n x m
        real(dp), dimension(:, :), intent(inout) :: diag
        real(dp), intent(in), optional :: alpha

        !local variable
        integer :: i
        real(dp) :: x

        ! check optional values
        if (present(alpha)) then 
          x = alpha
        else
          x = 1.d0
        endif

        diag = 0.0_dp

        do i = 1, min(size(diag,1),size(diag,2))
          diag(i, i) = x
        end do

    end subroutine eye

    pure subroutine diag_mat(vec,mat)
        real(dp), dimension(:), intent(in) :: vec
        real(dp), dimension(:,:), intent(out) :: mat

        integer :: i
        mat = 0.0_dp
        do i = 1, size(vec, 1)
            mat(i, i) = vec(i)
        end do

    end subroutine diag_mat

    function norm(vector)
        !> compute the norm-2 of a vector
        real(dp), dimension(:), intent(in) :: vector
        real(dp) :: norm

        norm = sqrt(sum(vector**2.d0))

    end function norm

    subroutine concatenate(arr, brr)

        !> Concatenate two matrices
        !> \param arr: first array
        !> \param brr: second array
        !>
        !> return arr concatenate brr(overwrites arr)

        real(dp), dimension(:, :), allocatable, intent(inout) :: arr
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
        call move_alloc(tmp_array, arr)

        arr(:, dim_cols + 1:) = brr

    end subroutine concatenate

    subroutine generate_diagonal_dominant(m, sparsity, diag_val, arr)
        !> Generate a diagonal dominant square matrix of dimension m
        !> \param m dimension of the matrix
        !> \param sparsity magnitude order of the off-diagonal values
        use random_mod, only: random_dp

        integer, intent(in) :: m ! size of the square matrix
        real(dp), optional :: diag_val
        integer :: i, j
        real(dp) :: sparsity
        real(dp), dimension(:, :) :: arr
        do i = 1,m
          do j = 1,m
            arr(j,i) = random_dp()
          end do
        end do

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

    end subroutine generate_diagonal_dominant

    subroutine diagonal(matrix, vector)
        !> return the diagonal of a matrix
        real(dp), dimension(:, :), intent(in) :: matrix
        real(dp), dimension(:) :: vector

        ! local variables
        integer :: i, m

        ! dimension of the matrix
        m = min(size(matrix, 1),size(matrix,2))

        do i = 1, m
            vector(i) = matrix(i, i)
        end do

    end subroutine diagonal

    subroutine initialize_subspace(diag, dim_sub, dim_base, precond)
        !> Brief generates a diagonal preconditioner for .
        !>
        !> return diagonal matrix

        ! input variable
        real(dp), dimension(:), intent(inout) :: diag
        integer, intent(in) :: dim_sub, dim_base

        ! local variables
        real(dp), dimension(:, :) :: precond
        integer, dimension(:), allocatable :: keys
        integer :: i, k

        allocate(keys(size(diag)))
        keys = 0

        ! sort diagonal
        call lapack_sort('I', diag, keys)
        ! Fill matrix with zeros
        precond = 0.0_dp

        ! Add one depending on the order of the matrix diagonal
        do i = 1, dim_sub
            k = search_key(keys, i)
            precond(k, i) = 1.d0
        end do

        deallocate(keys)
    end subroutine initialize_subspace

    subroutine modified_gram_schmidt(mat, nstart)
        !> Brief use modifed gram-schmidt orthogonalization on mat
        !> Brief nstart is the index of the first vector to orthogonalize

        ! input
        real(dp), dimension(:, :), intent(inout) :: mat
        integer, optional, intent(in) :: nstart

        integer :: i
        integer :: nrows, ncols
        integer :: idx_start
        real(dp), dimension(:), allocatable :: tmp_array, tmp_array2

        idx_start = 1
        if (present(nstart)) idx_start = nstart

        nrows = size(mat, 1)
        ncols = size(mat, 2)

        allocate (tmp_array(nrows))
        allocate (tmp_array2(nrows))

        do i = idx_start, ncols
            call lapack_matrix_vector('T', mat(:, :i - 1), mat(:, i), tmp_array)
            call lapack_matrix_vector('N', mat(:, :i - 1), tmp_array, tmp_array2)
            mat(:, i) = mat(:, i) - tmp_array2 
            mat(:, i) = mat(:, i)/norm(mat(:, i))
        end do

        deallocate (tmp_array)
        deallocate (tmp_array2)

    end subroutine modified_gram_schmidt

    function search_key(keys, i) result(k)
        !> Brief Search for a given index  in a vector
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

    subroutine write_matrix(path_file, mtx)
        !> Write matrix to path_file
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

    subroutine write_vector(path_file, vector)
        !> Write vector to path_file
        character(len=*), intent(in) :: path_file
        real(dp), dimension(:), intent(in) :: vector
        integer :: i

        open (unit=314, file=path_file, status="REPLACE")
        do i = 1, size(vector)
            write (314, *) vector(i)
        end do
        close (314)

    end subroutine write_vector

    subroutine check_deallocate_matrix(mtx)

        !> deallocate a matrix if allocated
        real(dp), dimension(:, :), allocatable, intent(inout) ::  mtx

        if (allocated(mtx)) then
            deallocate (mtx)
        end if

    end subroutine check_deallocate_matrix

    subroutine check_deallocate_vector(vec)

        !> deallocate a matrix if allocated
        real(dp), dimension(:), allocatable, intent(inout) ::  vec

        if (allocated(vec)) then
            deallocate (vec)
        end if

    end subroutine check_deallocate_vector

    subroutine unique_elements(n, arr, res, count, frequency, ind)
        !> Returns the unique elements of a vector
        !> Also returns the number of unique elements, their indices in the original vector.
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
        integer, dimension(:), intent(in)           :: arr    ! The input
        integer, dimension(:), intent(out)          :: res    ! The output
        integer, intent(out)                        :: count                   ! The number of unique elements
        integer, dimension(:), intent(out)          :: frequency ! The output
        integer, dimension(:), intent(out),optional :: ind
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


    subroutine unique_string_elements(n, arr, res, count)
        !> Returns the unique elements of a vector
        !> Also returns the number of unique elements, their indices in the original vector.
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
        character(len=3), dimension(:), intent(in)  :: arr    ! The input
        character(len=3), dimension(:), intent(out) :: res    ! The output
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



    subroutine sortedunique(list, n, indices, sorted)
        implicit none
        !   find "indices", the list of unique numbers in "list"

        integer :: n, kx, i, nitems
        integer, dimension(:) :: list
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


