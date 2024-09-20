!> Module that defines the system/molecule.
module system
    use precision_kinds, only: dp

    implicit none

    !> Number of electrons in the system.
    integer :: nelec

    !> Number of down electrons.
    integer :: ndn

    !> Number of up electrons.
    integer :: nup

    !> Coordinates of the centers.
    real(dp), dimension(:, :), allocatable :: cent

    !> Nuclear charges.
    real(dp), dimension(:), allocatable :: znuc

    !> Type of the wavefunction.
    integer, dimension(:), allocatable :: iwctype

    !> Number of types of centers.
    integer :: nctype

    !> Number of centers.
    integer :: ncent

    !> Total number of types of centers.
    integer :: nctype_tot

    !> Total number of centers.
    integer :: ncent_tot

    !> Elemental symbol of the centers.
    character(len=2), dimension(:), allocatable :: symbol

    !> Type of the centers.
    character(len=2), dimension(:), allocatable :: atomtyp

    !> New ghost type.
    integer :: newghostype

    !> Number of ghost centers.
    integer :: nghostcent

    private
    public :: nelec, ndn, nup, znuc, cent, iwctype, nctype, ncent, ncent_tot, nctype_tot
    public :: symbol, atomtyp, allocate_atom, deallocate_atom
    public :: newghostype, nghostcent
    save

contains

    !> Allocates memory for the system.
    subroutine allocate_atom()

        !> Allocating memory for centers.
        if (.not. allocated(cent)) allocate (cent(3, ncent_tot))

        !> Allocating memory for nuclear charges.
        if (.not. allocated(znuc)) allocate (znuc(nctype_tot))

        !> Allocating memory for wavefunction type.
        if (.not. allocated(iwctype)) allocate (iwctype(nctype_tot), source=0)

        !> Allocating memory for elemental symbols.
        if (.not. allocated(symbol)) allocate (symbol(ncent_tot))
    end subroutine allocate_atom

    !> Deallocates memory for the system.
    subroutine deallocate_atom()
        if (allocated(iwctype)) deallocate (iwctype)
        if (allocated(znuc)) deallocate (znuc)
        if (allocated(cent)) deallocate (cent)
    end subroutine deallocate_atom

end module system