module coefs
    !> need a better name is that the MO ?
    !> if yes can we put it in basis ?
    !> Arguments: coef, nbasis, norb
    use precision_kinds, only: dp

    implicit none

    real(dp), dimension(:, :, :), allocatable :: coef !(MBASIS,norb_tot,MWF)
    integer :: nbasis
    integer :: norb
    integer :: next_max

    private
    public :: coef, nbasis, norb, next_max
!    public :: allocate_coefs
    public :: deallocate_coefs
    save
contains
    ! subroutine allocate_coefs()
    !     use force_mod, only: MWF
    !     use precision_kinds, only: dp
    !     use vmc_mod, only: norb_tot, MBASIS
    !     if (.not. allocated(coef)) allocate (coef(MBASIS, norb_tot, MWF))
    ! end subroutine allocate_coefs

    subroutine deallocate_coefs()
        if (allocated(coef)) deallocate (coef)
    end subroutine deallocate_coefs

    ! subroutine resize_norb_coefs(new_size)
    !     !> Resize the second dimension of coef
    !     !> \param new_size: new size

    !     integer, INTENT(IN) :: new_size
    !     integer :: dim1, dim2, dim3
    !     real(dp), dimension(:, :, :), allocatable :: tmp_array

    !     dim2 = size(coef, 2)

    !     if (dim2 .lt. new_size) then
    !         dim1 = size(coef, 1)
    !         dim3 = size(coef, 3)
    !         allocate (tmp_array(dim1, new_size, dim3))
    !         tmp_array(:, :dim2, :) = coef
    !         deallocate (coef)
    !         call move_alloc(tmp_array, coef)
    !     endif

    ! endsubroutine

end module coefs
