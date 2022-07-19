module coefs
    !> need a better name is that the MO ?
    !> if yes can we put it in basis ?
    !> Arguments: coef, nbasis, norb

    implicit none

    integer :: nbasis
    integer :: norb
    integer :: next_max

    private
    public :: nbasis, norb, next_max
    save
end module coefs
