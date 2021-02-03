module coefs
    !> need a better name is that the MO ?
    !> if yes can we put it in basis ?
    !> Arguments: coef, nbasis, norb
    use precision_kinds, only: dp

    real(dp), dimension(:, :, :), allocatable :: coef !(MBASIS,MORB,MWF)
    integer :: nbasis
    integer :: norb

    private
    public :: coef, nbasis, norb
    public :: allocate_coefs, deallocate_coefs
    save
contains
    ! subroutine allocate_coefs()
    !     use force_mod, only: MWF
    !     use precision_kinds, only: dp
    !     use vmc_mod, only: MORB, MBASIS
    !     if (.not. allocated(coef)) allocate (coef(MBASIS, MORB, MWF))
    ! end subroutine allocate_coefs

    subroutine deallocate_coefs()
        if (allocated(coef)) deallocate (coef)
    end subroutine deallocate_coefs

end module coefs
