!> @brief Module for analytical basis set information.
!> @author CHAMP developers
!> @date 2020
!>
!> @details This module stores the fundamental basis set parameters including
!> exponents (screening constants) and the number of basis functions per angular
!> momentum type for each atomic center. It provides the foundation for constructing
!> Gaussian basis sets used in wavefunction calculations.
!>
!> Key data structures:
!> - zex: Basis function exponents (screening constants)
!> - ns, np, nd, nf, ng: Number of basis functions per angular momentum (s, p, d, f, g)
!> - betaq: Beta parameter for basis function scaling
!>
!> @note This module is primarily used with analytical Gaussian basis sets.
module basis
    use precision_kinds, only: dp

    implicit none

    !> Screening constants (exponents) for each basis function, dimension (nbasis, nctype_tot)
    real(dp), dimension(:, :), allocatable :: zex

    !> Beta coefficient parameter for basis function scaling
    real(dp) :: betaq

    !> Number of s-type basis functions at each center type, dimension (nctype_tot)
    integer, dimension(:), allocatable :: ns

    !> Number of p-type basis functions at each center type, dimension (nctype_tot)
    integer, dimension(:), allocatable :: np

    !> Number of d-type basis functions at each center type, dimension (nctype_tot)
    integer, dimension(:), allocatable :: nd

    !> Number of f-type basis functions at each center type, dimension (nctype_tot)
    integer, dimension(:), allocatable :: nf

    !> Number of g-type basis functions at each center type, dimension (nctype_tot)
    integer, dimension(:), allocatable :: ng

    private
    public :: zex, betaq
    public :: ns, np, nd, nf, ng
    public :: allocate_basis, deallocate_basis
    save
contains
    !> Allocates memory for analytical basis set arrays.
    !>
    !> @details This subroutine is deprecated. Memory allocation for basis set
    !> parameters is now performed during input file parsing when basis set
    !> information is read directly.
    !>
    !> @note Deprecated - allocation occurs during basis set input reading.
    subroutine allocate_basis()
        ! Allocate basis is called while reading the basis set from the input file.
    end subroutine allocate_basis

    !> Deallocates memory for analytical basis set arrays.
    !>
    !> @details This subroutine frees dynamically allocated arrays used for
    !> storing basis set parameters. It releases memory for:
    !> - ng, nf, nd, np, ns: Angular momentum basis function counts
    !> - zex: Basis function exponents (screening constants)
    !>
    !> @note Called at program termination to prevent memory leaks.
    subroutine deallocate_basis()

        if (allocated(ng)) deallocate (ng)
        if (allocated(nf)) deallocate (nf)
        if (allocated(nd)) deallocate (nd)
        if (allocated(np)) deallocate (np)
        if (allocated(ns)) deallocate (ns)
        if (allocated(zex)) deallocate (zex)
    end subroutine deallocate_basis

end module basis

!> @brief Module defining array dimension constants for numerical basis representation.
!> @author CHAMP developers
!> @date 2022
!>
!> @details This module provides compile-time constants that define maximum array
!> dimensions for numerical basis set representations. These parameters are used
!> throughout the code to allocate arrays for storing radial wavefunctions on grids.
!>
!> Constants:
!> - MRWF_PTS: Maximum grid points for radial wavefunction representation
!> - MRWF: Maximum number of radial basis functions per center type
!>
!> @note These are fixed parameters that may need adjustment for very large basis sets.
module numbas_mod

    implicit none

    !> Maximum number of radial grid points for numerical wavefunction representation
    integer, parameter :: MRWF_PTS = 4000

    !> Maximum number of radial basis functions (shells) per center type
    integer, parameter :: MRWF = 200

    private
    public :: MRWF, MRWF_PTS
    save

end module numbas_mod

!> @brief Module for numerical basis function exponent and coefficient storage.
!> @author CHAMP developers
!> @date 2022
!>
!> @details This module stores coefficients for numerical representation of basis
!> functions across multiple geometries. It manages exponent parameters (ae) and
!> expansion coefficients (ce) needed for evaluating numerical basis functions.
!>
!> Key arrays:
!> - ae: Exponent parameters for numerical basis functions
!> - ce: Expansion coefficients for linear combinations
!>
!> @note Used primarily for force calculations with multiple geometries.
module numexp

      use multiple_geo, only: MFORCE
      use numbas_mod, only: MRWF
      use precision_kinds, only: dp
      use vmc_mod, only: NCOEF

    implicit none

    !> Exponent parameters for numerical basis functions, dimension (2, MRWF, nctype_tot, MFORCE)
    real(dp), dimension(:, :, :, :), allocatable :: ae

    !> Expansion coefficients for numerical basis functions, dimension (NCOEF, MRWF, nctype_tot, MFORCE)
    real(dp), dimension(:, :, :, :), allocatable :: ce

    private
    public :: ae, ce
    public :: allocate_numexp, deallocate_numexp
    save
contains
    !> Allocates memory for numerical basis function coefficient arrays.
    !>
    !> @details This subroutine allocates arrays for storing exponent parameters
    !> and expansion coefficients for numerical basis functions across multiple
    !> center types and force calculations. Allocates:
    !> - ae(2, MRWF, nctype_tot, MFORCE): Exponent parameters
    !> - ce(NCOEF, MRWF, nctype_tot, MFORCE): Expansion coefficients
    !>
    !> @note Called during basis set initialization for numerical representations.
    subroutine allocate_numexp()
      use multiple_geo, only: MFORCE
      use numbas_mod, only: MRWF
      use system,  only: nctype_tot
      use vmc_mod, only: NCOEF
        if (.not. allocated(ae)) allocate (ae(2, MRWF, nctype_tot, MFORCE))
        if (.not. allocated(ce)) allocate (ce(NCOEF, MRWF, nctype_tot, MFORCE))
    end subroutine allocate_numexp

    !> Deallocates memory for numerical basis function coefficient arrays.
    !>
    !> @details Frees memory for ce (expansion coefficients) and ae (exponents).
    subroutine deallocate_numexp()
        if (allocated(ce)) deallocate (ce)
        if (allocated(ae)) deallocate (ae)
    end subroutine deallocate_numexp

end module numexp

!> @brief Module for numerical radial basis functions on grids.
!> @author CHAMP developers
!> @date 2022
!>
!> @details This module stores radial basis functions and their derivatives evaluated
!> on numerical grids. It manages the grid parameters, radial function values, and
!> indexing information needed for numerical wavefunction evaluation.
!>
!> Key data structures:
!> - rwf, d2rwf: Radial wavefunctions and their second derivatives on grids
!> - arg, r0, nr: Grid parameters (arguments, starting point, number of points)
!> - igrid, iwrwf: Grid and wavefunction indexing arrays
!> - nrbas, rmaxwf: Basis function counts and maximum radii
!>
!> @note Used for numerical basis set representations and all-electron calculations.
module numbas

      use numbas_mod, only: MRWF,MRWF_PTS
      use precision_kinds, only: dp

    implicit none

    !> Grid argument parameter for each center type, dimension (nctype_tot)
    real(dp), dimension(:), allocatable :: arg

    !> Second derivatives of radial wavefunctions on grid, dimension (MRWF_PTS, MRWF, nctype_tot, nwftype)
    real(dp), dimension(:, :, :, :), allocatable :: d2rwf

    !> Grid type index for each center type, dimension (nctype_tot)
    integer, dimension(:), allocatable :: igrid

    !> Wavefunction index mapping to radial functions, dimension (nbasis, nctype_tot)
    integer, dimension(:, :), allocatable :: iwrwf

    !> Number of radial grid points for each center type, dimension (nctype_tot)
    integer, dimension(:), allocatable :: nr

    !> Number of radial basis functions (shells) per center type, dimension (nctype_tot)
    integer, dimension(:), allocatable :: nrbas

    !> Total number of radial grid points used
    integer :: numr

    !> Starting radius for radial grid per center type, dimension (nctype_tot)
    real(dp), dimension(:), allocatable :: r0

    !> Maximum radius for each radial basis function, dimension (MRWF, nctype_tot)
    real(dp), dimension(:,:), allocatable :: rmaxwf

    !> Radial wavefunctions evaluated on grid, dimension (MRWF_PTS, MRWF, nctype_tot, nwftype)
    real(dp), dimension(:, :, :, :), allocatable :: rwf

    private
    public :: arg, d2rwf, igrid, iwrwf, nr, nrbas, numr, r0, rwf, rmaxwf
    public :: allocate_numbas, deallocate_numbas
    save
contains
    !> Allocates memory for numerical radial basis function arrays.
    !>
    !> @details This subroutine allocates arrays for storing radial wavefunctions
    !> and their properties on numerical grids. It allocates:
    !> - arg(nctype_tot): Grid argument parameters
    !> - d2rwf(MRWF_PTS, MRWF, nctype_tot, nwftype): Second derivatives
    !> - igrid(nctype_tot): Grid type indices (initialized to 0)
    !> - iwrwf(nbasis, nctype_tot): Wavefunction index mapping (initialized to 0)
    !> - nr(nctype_tot): Number of radial points (initialized to 0)
    !> - nrbas(nctype_tot): Number of radial basis functions (initialized to 0)
    !> - r0(nctype_tot): Starting radii
    !> - rmaxwf(MRWF, nctype_tot): Maximum radii (initialized to 0.0)
    !> - rwf(MRWF_PTS, MRWF, nctype_tot, nwftype): Radial wavefunctions
    !>
    !> @note Called during numerical basis set initialization.
    subroutine allocate_numbas()
      use coefs,   only: nbasis
      use multiple_geo, only: nwftype
      use numbas_mod, only: MRWF,MRWF_PTS
      use system,  only: nctype_tot
        if (.not. allocated(arg)) allocate (arg(nctype_tot))
        if (.not. allocated(d2rwf)) allocate (d2rwf(MRWF_PTS, MRWF, nctype_tot, nwftype))
        if (.not. allocated(igrid)) allocate (igrid(nctype_tot), source=0)
        if (.not. allocated(iwrwf)) allocate (iwrwf(nbasis, nctype_tot), source=0)
        if (.not. allocated(nr)) allocate (nr(nctype_tot), source=0)
        if (.not. allocated(nrbas)) allocate (nrbas(nctype_tot), source=0)
        if (.not. allocated(r0)) allocate (r0(nctype_tot))
        if (.not. allocated(rmaxwf)) allocate (rmaxwf(MRWF,nctype_tot), source=0.0d0)
        if (.not. allocated(rwf)) allocate (rwf(MRWF_PTS, MRWF, nctype_tot, nwftype))
    end subroutine allocate_numbas

    !> Deallocates memory for numerical radial basis function arrays.
    !>
    !> @details Frees all arrays used for numerical radial basis function storage.
    subroutine deallocate_numbas()
        if (allocated(rmaxwf)) deallocate (rmaxwf)
        if (allocated(rwf)) deallocate (rwf)
        if (allocated(r0)) deallocate (r0)
        if (allocated(nrbas)) deallocate (nrbas)
        if (allocated(nr)) deallocate (nr)
        if (allocated(iwrwf)) deallocate (iwrwf)
        if (allocated(igrid)) deallocate (igrid)
        if (allocated(d2rwf)) deallocate (d2rwf)
        if (allocated(arg)) deallocate (arg)
    end subroutine deallocate_numbas

end module numbas

!> @brief Module for numerical basis function angular momentum indexing.
!> @author CHAMP developers
!> @date 2022
!>
!> @details This module manages the mapping between basis functions and their
!> angular momentum quantum numbers for numerical basis sets. It stores indices
!> that relate basis functions to their angular momentum types.
!>
!> Key arrays:
!> - iwlbas: Angular momentum (l) indices for each basis function
!> - nbastyp: Number of basis function types per center
!>
!> @note Used for organizing numerical basis functions by angular momentum.
module numbas1

    implicit none

    !> Angular momentum (l) quantum number for each basis function, dimension (nbasis, nctype_tot)
    integer, dimension(:, :), allocatable :: iwlbas

    !> Number of basis function types per center type, dimension (nctype_tot)
    integer, dimension(:), allocatable :: nbastyp

    private
    public :: iwlbas, nbastyp
    public :: allocate_numbas1, deallocate_numbas1
    save
contains
    !> Allocates memory for angular momentum indexing arrays.
    !>
    !> @details This subroutine allocates arrays for mapping basis functions to
    !> their angular momentum quantum numbers. Currently allocates:
    !> - iwlbas(nbasis, nctype_tot): Angular momentum indices (initialized to 0)
    !>
    !> @note nbastyp allocation is commented out and handled elsewhere.
    subroutine allocate_numbas1()
      use coefs,   only: nbasis
      use system,  only: nctype_tot
        if (.not. allocated(iwlbas)) allocate (iwlbas(nbasis, nctype_tot), source=0)
        ! if (.not. allocated(nbastyp)) allocate (nbastyp(nctype_tot))
    end subroutine allocate_numbas1

    !> Deallocates memory for angular momentum indexing arrays.
    !>
    !> @details Frees nbastyp and iwlbas arrays.
    subroutine deallocate_numbas1()
        if (allocated(nbastyp)) deallocate (nbastyp)
        if (allocated(iwlbas)) deallocate (iwlbas)
    end subroutine deallocate_numbas1

end module numbas1

!> @brief Module for basis function range indexing per atomic center.
!> @author CHAMP developers
!> @date 2022
!>
!> @details This module stores the starting and ending indices of basis functions
!> associated with each atomic center. This indexing allows efficient access to
!> basis functions belonging to specific atoms.
!>
!> Key arrays:
!> - ibas0: Starting basis function index for each center
!> - ibas1: Ending basis function index for each center
!>
!> @note Used for partitioning basis functions by atomic center.
module numbas2

    implicit none

    !> Starting basis function index for each atomic center, dimension (ncent_tot)
    integer, dimension(:), allocatable :: ibas0

    !> Ending basis function index for each atomic center, dimension (ncent_tot)
    integer, dimension(:), allocatable :: ibas1

    private
    public :: ibas0, ibas1
    public :: allocate_numbas2, deallocate_numbas2
    save
contains
    !> Allocates memory for basis function range indexing arrays.
    !>
    !> @details This subroutine allocates arrays for storing the starting and
    !> ending indices of basis functions for each atomic center. Allocates:
    !> - ibas0(ncent_tot): Starting indices (initialized to 0)
    !> - ibas1(ncent_tot): Ending indices (initialized to 0)
    !>
    !> @note Used for partitioning basis functions by atomic center.
    subroutine allocate_numbas2()
      use system,  only: ncent_tot
        if (.not. allocated(ibas0)) allocate (ibas0(ncent_tot), source=0)
        if (.not. allocated(ibas1)) allocate (ibas1(ncent_tot), source=0)
    end subroutine allocate_numbas2

    !> Deallocates memory for basis function range indexing arrays.
    !>
    !> @details Frees ibas1 and ibas0 arrays.
    subroutine deallocate_numbas2()
        if (allocated(ibas1)) deallocate (ibas1)
        if (allocated(ibas0)) deallocate (ibas0)
    end subroutine deallocate_numbas2

end module numbas2

!> @brief Master module for coordinated basis set memory management.
!> @author CHAMP developers
!> @date 2022
!>
!> @details This module provides unified allocation and deallocation interfaces
!> for all basis set related modules. It orchestrates memory management across
!> analytical basis (basis), numerical basis (numbas, numbas1, numbas2), and
!> numerical expansion (numexp) modules.
!>
!> Manages allocation for:
!> - numexp: Numerical basis exponents and coefficients
!> - numbas: Radial wavefunctions on grids
!> - numbas1: Angular momentum indexing
!> - numbas2: Basis function range per center
!>
!> @note This is the primary interface for basis set memory management.
module m_basis
contains
!> Allocates memory for all basis set modules.
!>
!> @details This subroutine calls allocation routines from all basis-related
!> modules in the proper order:
!> - allocate_numexp(): Numerical expansion coefficients
!> - allocate_numbas(): Radial wavefunctions and grid data
!> - allocate_numbas1(): Angular momentum indices
!> - allocate_numbas2(): Basis function range indices per center
!>
!> @note Called from allocate_vmc() after system parameters are initialized.
subroutine allocate_m_basis()
      use numbas,  only: allocate_numbas
      use numbas1, only: allocate_numbas1
      use numbas2, only: allocate_numbas2
      use numexp,  only: allocate_numexp

    implicit none

    call allocate_numexp()
    call allocate_numbas()
    call allocate_numbas1()
    call allocate_numbas2()
end subroutine allocate_m_basis

!> Deallocates memory for all basis set modules.
!>
!> @details This subroutine calls deallocation routines from all basis-related
!> modules:
!> - deallocate_basis(): Analytical basis parameters
!> - deallocate_numexp(): Numerical expansion coefficients
!> - deallocate_numbas(): Radial wavefunctions
!> - deallocate_numbas1(): Angular momentum indices
!> - deallocate_numbas2(): Basis function range indices
!>
!> @note Called from deallocate_vmc() at program termination.
subroutine deallocate_m_basis()
      use basis,   only: deallocate_basis
      use numbas,  only: deallocate_numbas
      use numbas1, only: deallocate_numbas1
      use numbas2, only: deallocate_numbas2
      use numexp,  only: deallocate_numexp

    implicit none

    call deallocate_basis()
    call deallocate_numexp()
    call deallocate_numbas()
    call deallocate_numbas1()
    call deallocate_numbas2()
end subroutine deallocate_m_basis
end module

!> @brief Module for TREXIO basis set import and reconstruction.
!> @author CHAMP developers
!> @date 2022
!>
!> @details This module manages basis set information imported from TREXIO files.
!> It stores shell structure, angular momentum information, and indexing arrays
!> needed to reconstruct basis functions from TREXIO format.
!>
!> Key data structures:
!> - Shell information (basis_num_shell, basis_shell_ang_mom)
!> - Spherical harmonic indexing (index_slm, slm_per_l)
!> - Radial function mapping (ao_radial_index)
!> - Per-center basis function counts (num_rad_per_cent, num_ao_per_cent)
!>
!> @note Used when reading basis sets from TREXIO format files.
module m_trexio_basis
      use coefs,   only: nbasis
    implicit none

    !> Number of spherical harmonics per angular momentum: s=1, p=3, d=6, f=10, g=15
    integer, dimension(5)       :: slm_per_l = (/1, 3, 6, 10, 15/)

    !> Total number of shells (contracted basis functions) in the basis set
    integer                     :: basis_num_shell

    !> Spherical harmonic (lm) index for each atomic orbital basis function, dimension (nbasis)
    integer, allocatable        :: index_slm(:)

    !> Radial function index for each atomic orbital basis function, dimension (nbasis)
    integer, allocatable        :: ao_radial_index(:)

    !> Number of radial basis functions per atomic center, dimension (ncent_tot)
    integer, allocatable        :: num_rad_per_cent(:)

    !> Number of atomic orbital basis functions per atomic center, dimension (ncent_tot)
    integer, allocatable        :: num_ao_per_cent(:)

    !> Angular momentum quantum number (l) for each shell, dimension (basis_num_shell)
    integer, allocatable        :: basis_shell_ang_mom(:)

    private
    public :: slm_per_l, index_slm, num_rad_per_cent, num_ao_per_cent
    public :: basis_num_shell, basis_shell_ang_mom
    public :: ao_radial_index
    public :: gnorm

    contains

    !> Calculates normalization constant for Gaussian basis functions.
    !>
    !> @details This function computes the normalization constant for Gaussian
    !> basis functions of the form r^l * exp(-alpha*r^2), ensuring proper
    !> normalization in the radial part. The normalization depends on both
    !> the exponent (alpha) and angular momentum quantum number (l).
    !>
    !> Formulas:
    !> - l=0 (s): N = (2*alpha)^(3/4) * 2 / pi^(1/4)
    !> - l=1 (p): N = (2*alpha)^(5/4) * sqrt(8/3) / pi^(1/4)
    !> - l=2 (d): N = (2*alpha)^(7/4) * sqrt(16/15) / pi^(1/4)
    !> - l=3 (f): N = (2*alpha)^(9/4) * sqrt(32/105) / pi^(1/4)
    !> - l=4 (g): N = (2*alpha)^(11/4) * sqrt(64/945) / pi^(1/4)
    !>
    !> @param[in] exponent Gaussian exponent (alpha)
    !> @param[in] l Angular momentum quantum number (0=s, 1=p, 2=d, 3=f, 4=g)
    !> @return Normalization constant for the Gaussian basis function
    !>
    !> @note Returns 1.0 for unsupported angular momentum values (l > 4).
    double precision function gnorm(exponent, l)
        use precision_kinds,    only: dp
        implicit none
        real(dp), intent (in)       :: exponent
        integer, intent (in)    :: l
        real(dp), parameter     :: pi = 4.0d0*atan(1.0d0)

        gnorm = 1.0d0

        if (l .eq. 0) then
            gnorm = (2.d0*exponent)**(3.d0/4.d0)*2.d0*(1.d0/(pi**(1.d0/4.d0)))
        elseif (l .eq. 1) then
            gnorm = (2.d0*exponent)**(5.d0/4.d0)*dsqrt(8.d0/3.d0)*(1.d0/(pi**(1.d0/4.d0)))
        elseif (l .eq. 2) then
            gnorm = (2.d0*exponent)**(7.d0/4.d0)*dsqrt(16.d0/15.d0)*(1.d0/(pi**(1.d0/4.d0)))
        elseif (l .eq. 3) then
            gnorm = (2.d0*exponent)**(9.d0/4.d0)*dsqrt(32.d0/105.d0)*(1.d0/(pi**(1.d0/4.d0)))
        elseif (l .eq. 4) then
            gnorm = (2.d0*exponent)**(11.d0/4.d0)*dsqrt(64.d0/945.d0)*(1.d0/(pi**(1.d0/4.d0)))
        endif

    end function gnorm

end module
