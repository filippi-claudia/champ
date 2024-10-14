!> Module that contains VMC parameters
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


    !> nstoo
    integer, dimension(:), allocatable :: nstoo !(nwftypeorb) allocate later

    !> nstoj
    integer, dimension(:), allocatable :: nstoj !(nwftypejas) allocate later

    !> jtos
    integer, dimension(:, :), allocatable :: jtos  !(nwftypejas,entry) allocate later

    !> otos
    integer, dimension(:, :), allocatable :: otos  !(nwftypeorb,entry) allocate later

    !> stoj
    integer, dimension(:), allocatable :: stoj  !(nstates) allocate later

    !> stoo
    integer, dimension(:), allocatable :: stoo  !(nstates) allocate later

    !> stobjx
    integer, dimension(:), allocatable :: stobjx  !(nstates) allocate later

    !> bjxtoo
    integer, dimension(:), allocatable :: bjxtoo  !(nbjx) allocate later

    !> bjxtoj
    integer, dimension(:), allocatable :: bjxtoj  !(nbjx) allocate later

    !> norb_tot :: total number of orbitals
    integer :: norb_tot

    !> nctyp3x :: max(3, nctype_tot)
    integer :: nctyp3x

    !> nmat_dim :: nup*nup
    integer :: nmat_dim

    !> nmat_dim2 :: nelec*(nelec - 1)/2
    integer :: nmat_dim2

    !> nwftypeorb
    integer :: nwftypeorb

    !> nwftypejas
    integer :: nwftypejas

    !> nstojmax
    integer :: nstojmax

    !> nstoomax
    integer :: nstoomax

    !> nbjx
    integer :: nbjx

    !> nstoo_tot
    integer :: nstoo_tot

    !> nstoj_tot
    integer :: nstoj_tot

    !> extraj
    integer :: extraj

    !> extrao
    integer :: extrao

    !> mterms
    integer :: mterms

    !> ncent3
    integer :: ncent3

    !> @param NCOEF :: 5
    integer, parameter :: NCOEF = 5

    !> @param MEXCIT :: 10
    integer, parameter :: MEXCIT = 10

    private
    public :: norb_tot, nctyp3x
    public :: nmat_dim, nmat_dim2

    public :: mterms, nwftypejas, nwftypeorb, nstojmax, nstoomax, nbjx, nstoj_tot, nstoo_tot
    public :: nstoo, nstoj, jtos, otos, stoj, stoo, stobjx, extrao, extraj, bjxtoo, bjxtoj

    public :: ncent3, NCOEF, MEXCIT
    public :: set_vmc_size

    save
contains
    !> Subroutine that sets the size of the VMC arrays
    subroutine set_vmc_size
      use system,  only: ncent_tot,nctype_tot,ndn,nelec,nup

        nmat_dim = nup*nup
        nmat_dim2 = nelec*(nelec - 1)/2
        nctyp3x = max(3, nctype_tot)
        ncent3 = 3*ncent_tot

    end subroutine set_vmc_size
end module vmc_mod
