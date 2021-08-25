module vmc_mod
    !> Arguments:
    use precision_kinds, only: dp

    ! nelec  >= number of electrons
    ! MORB   >= number of orbitals
    ! nbasis >= number of basis functions
    ! ndet   >= number of determinants
    ! ncent_tot  >= number of centers
    ! nctype >= number of center types
    ! nctyp3x = max(3,MCTYPE)

    ! Slater matrices are dimensioned (MELEC/2)**2 assuming
    ! equal numbers of up and down spins. MELEC has to be
    ! correspondingly larger if spin polarized calculations
    ! are attempted.

    ! PLT@eScienceCenter(2020) Moved the parameter here:
    ! "For Jastrow4 NEQSX=2*(MORDJ-1) is sufficient.
    !  For Jastrow3 NEQSX=2*MORDJ should be sufficient.
    !  I am setting NEQSX=6*MORDJ simply because that is how it was for
    !  Jastrow3 for reasons I do not understand."
    !     parameter(NEQSX=2*(MORDJ-1),MTERMS=55)


    real(dp), parameter :: radmax = 10.d0
    integer, parameter :: nrad = 3001
    real(dp), parameter :: delri = (nrad - 1)/radmax

    integer :: MORB
    integer :: nctyp3x
    integer, parameter :: NSPLIN = 1001, MORDJ = 7

    integer :: nmat_dim, nmat_dim2
    integer, parameter :: MORDJ1 = MORDJ + 1

    integer, parameter :: NEQSX = 6*MORDJ, MTERMS = 55
    integer :: ncent3

    integer, parameter :: NCOEF = 5
    integer, parameter :: MEXCIT = 10

    private
    public :: MORB, nctyp3x
    public :: NSPLIN, nrad, MORDJ, MORDJ1, nmat_dim, nmat_dim2
    public :: radmax, delri

    public :: NEQSX, MTERMS

    public :: ncent3, NCOEF, MEXCIT
    public :: set_vmc_size

    save
contains
    subroutine set_vmc_size
        use const, only: nelec
        use atom, only: nctype_tot, ncent_tot
        nmat_dim = nelec*nelec/4
        nmat_dim2 = nelec*(nelec - 1)/2
        nctyp3x = max(3, nctype_tot)
        ncent3 = 3*ncent_tot

    end subroutine set_vmc_size
end module vmc_mod
