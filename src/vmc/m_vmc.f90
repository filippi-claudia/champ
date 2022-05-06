module vmc_mod
    !> Arguments:
    use precision_kinds, only: dp

    ! MELEC  >= number of electrons
    ! MORB   >= number of orbitals
    ! MBASIS >= number of basis functions
    ! MDET   >= number of determinants
    ! MCENT  >= number of centers
    ! MCTYPE >= number of center types
    ! MCTYP3X=max(3,MCTYPE)

    integer :: MMAT_DIM20
    real(dp), parameter :: radmax = 10.d0
    integer, parameter :: nrad = 3001
    real(dp), parameter :: delri = (nrad - 1)/radmax

    integer, parameter :: MELEC = 50, MORB = 550, MBASIS = 550, MDET = 4000, MCENT = 20
    integer, parameter :: MCTYPE = 3
    integer, parameter :: MCTYP3X = 5, NSPLIN = 1001, MORDJ = 7

    integer, parameter :: MMAT_DIM = (MELEC*MELEC)/4, MMAT_DIM2 = (MELEC*(MELEC - 1))/2
    integer, parameter :: MORDJ1 = MORDJ + 1

    integer, parameter :: NEQSX = 6*MORDJ, MTERMS = 55
    integer, parameter :: MCENT3 = 3*MCENT

    integer, parameter :: NCOEF = 5
    integer, parameter :: MEXCIT = 10

    private
    public :: MELEC, MORB, MBASIS, MDET, MCENT, MCTYPE, MCTYP3X
    public :: NSPLIN, nrad, MORDJ, MORDJ1, MMAT_DIM, MMAT_DIM2, MMAT_DIM20
    public :: radmax, delri

    public :: NEQSX, MTERMS
    public :: MCENT3, NCOEF, MEXCIT

    save
end module vmc_mod
