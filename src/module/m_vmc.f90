module vmc_mod
    !> Arguments:
    use precision_kinds, only: dp

    ! nelec      = number of electrons
    ! norb_tot   = number of orbitals read from the orbital file
    ! nbasis     = number of basis functions
    ! ndet       = number of determinants
    ! ncent_tot  = number of centers
    ! nctype     = number of center types
    ! nctyp3x    = max(3,MCTYPE)

    ! Slater matrices are dimensioned (MELEC/2)**2 assuming
    ! equal numbers of up and down spins. MELEC has to be
    ! correspondingly larger if spin polarized calculations
    ! are attempted.

    ! PLT@eScienceCenter(2020) Moved the parameter here:
    ! "For Jastrow4 neqsx=2*(nordj-1) is sufficient.
    !  For Jastrow3 neqsx=2*nordj should be sufficient.
    !  I am setting neqsx=6*nordj simply because that is how it was for
    !  Jastrow3 for reasons I do not understand."
    !     parameter(neqsx=2*(nordj-1),MTERMS=55)


    real(dp), parameter :: radmax = 10.d0
    integer, parameter :: nrad = 3001
    real(dp), parameter :: delri = (nrad - 1)/radmax

    integer :: norb_tot
    integer :: nctyp3x
    integer :: nordj

    integer :: nmat_dim, nmat_dim2
    integer :: nordj1   ! nordj+1

    integer :: neqsx    ! 6*nordj
    integer :: mterms
    integer :: ncent3

    integer, parameter :: NCOEF = 5
    integer, parameter :: MEXCIT = 10

    private
    public :: norb_tot, nctyp3x
    public :: nrad, nordj, nordj1, nmat_dim, nmat_dim2
    public :: radmax, delri

    public :: neqsx, mterms

    public :: ncent3, NCOEF, MEXCIT
    public :: set_vmc_size

    save
contains
    subroutine set_vmc_size
        use atom, only: nctype_tot, ncent_tot

      use system, only: nelec
      use system, only: nup
      use system, only: ndn
        nmat_dim = nup*nup
        nmat_dim2 = nelec*(nelec - 1)/2
        nctyp3x = max(3, nctype_tot)
        ncent3 = 3*ncent_tot

    end subroutine set_vmc_size
end module vmc_mod
