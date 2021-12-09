module array_resize_utils

    use precision_kinds, only: dp
    implicit none

    !> \private
    private
    !> \public
    public :: resize_array, resize_matrix, resize_tensor

contains

    subroutine resize_array(arr, new_size)
        !> Resize a 1D array
        !> \param arr : matrix to resize
        !> \param new_size : new dimension

        real(dp), dimension(:), allocatable, INTENT(INOUT) :: arr
        integer, INTENT(IN) :: new_size
        integer :: old_dim

        real(dp), dimension(:), allocatable :: tmp_array

        old_dim = size(arr, 1)

        if (old_dim .lt. new_size) then
            allocate (tmp_array(new_size), source=0.0_dp)
            tmp_array(:old_dim) = arr

            deallocate (arr)
            call move_alloc(tmp_array, arr)

        endif

    end subroutine resize_array

    subroutine resize_matrix(mat, new_size, index_dim)
        !> Resize a matrix
        !> \param mat : matrix to resize
        !> \param new_size : new dimension
        !> \index_dim : index of the dim to resize

        real(dp), dimension(:, :), allocatable, INTENT(INOUT) :: mat
        integer, INTENT(IN) :: new_size
        integer, INTENT(IN) :: index_dim
        integer :: old_dim, dim1, dim2

        real(dp), dimension(:, :), allocatable :: tmp_array

        old_dim = size(mat, index_dim)

        if (old_dim .lt. new_size) then

            dim1 = size(mat, 1)
            dim2 = size(mat, 2)

            if (index_dim .eq. 1) then
                dim1 = new_size
                allocate (tmp_array(dim1, dim2), source=0.0_dp)
                tmp_array(:old_dim, :) = mat
            else
                dim2 = new_size
                allocate (tmp_array(dim1, dim2), source=0.0_dp)
                tmp_array(:, :old_dim) = mat
            endif

            deallocate (mat)
            call move_alloc(tmp_array, mat)

        endif

    endsubroutine resize_matrix

    subroutine resize_tensor(mat, new_size, index_dim)
        !> Resize a 3D tensor matrix
        !> \param mat : matrix to resize
        !> \param new_size : new dimension
        !> \index_dim : index of the dim to resize

        real(dp), dimension(:, :, :), allocatable, INTENT(INOUT) :: mat
        integer, INTENT(IN) :: new_size
        integer, INTENT(IN) :: index_dim
        integer :: old_dim, dim1, dim2, dim3

        real(dp), dimension(:, :, :), allocatable :: tmp_array

        old_dim = size(mat, index_dim)

        if (old_dim .lt. new_size) then

            dim1 = size(mat, 1)
            dim2 = size(mat, 2)
            dim3 = size(mat, 3)

            if (index_dim .eq. 1) then
                dim1 = new_size
                allocate (tmp_array(dim1, dim2, dim3), source=0.0_dp)
                tmp_array(:old_dim, :, :) = mat
            elseif (index_dim .eq. 2) then
                dim2 = new_size
                allocate (tmp_array(dim1, dim2, dim3), source=0.0_dp)
                tmp_array(:, :old_dim, :) = mat
            else
                dim3 = new_size
                allocate (tmp_array(dim1, dim2, dim3), source=0.0_dp)
                tmp_array(:, :, :old_dim) = mat
            endif

            deallocate (mat)
            call move_alloc(tmp_array, mat)

        endif

    endsubroutine resize_tensor

end module
