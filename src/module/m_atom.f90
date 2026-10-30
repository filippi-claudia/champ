!> @brief Module for atomic and molecular system geometry and properties.
!> @author CHAMP developers
!> @date 2022
!>
!> @details This module defines the fundamental properties of the quantum system
!> including atomic positions, nuclear charges, electron counts, and center types.
!> It stores both the actual system parameters and total parameters (including ghost
!> centers if present). The module provides allocation and deallocation routines
!> for dynamic memory management.
!>
!> Key data structures:
!> - Electron configuration (nelec, nup, ndn)
!> - Nuclear geometry (cent, znuc, symbol)
!> - Center type information (nctype, ncent, iwctype, atomtyp)
!> - Ghost center support (nghostcent, newghostype)
!>
!> @note This module must be initialized before any wavefunction or energy calculations.
module system
    use precision_kinds, only: dp

    implicit none

    !> Total number of electrons in the system (up + down)
    integer :: nelec

    !> Number of down-spin (beta) electrons
    integer :: ndn

    !> Number of up-spin (alpha) electrons
    integer :: nup

    !> Cartesian coordinates of atomic centers, dimension (3, ncent_tot) in atomic units (bohr)
    real(dp), dimension(:, :), allocatable :: cent

    !> Nuclear charges for each center type, dimension (nctype_tot)
    real(dp), dimension(:), allocatable :: znuc

    !> Wavefunction type index for each center type, dimension (nctype_tot)
    integer, dimension(:), allocatable :: iwctype

    !> Number of unique center types (excluding ghost centers)
    integer :: nctype

    !> Number of atomic centers (excluding ghost centers)
    integer :: ncent

    !> Total number of unique center types (including ghost centers)
    integer :: nctype_tot

    !> Total number of atomic centers (including ghost centers)
    integer :: ncent_tot

    !> Element symbols for each center, dimension (ncent_tot), e.g., 'H', 'C', 'O'
    character(len=2), dimension(:), allocatable :: symbol

    !> Atom type labels for each center type, dimension (nctype_tot)
    character(len=2), dimension(:), allocatable :: atomtyp

    !> Index for newly added ghost center type
    integer :: newghostype

    !> Number of ghost centers in the system
    integer :: nghostcent

    private
    public :: nelec, ndn, nup, znuc, cent, iwctype, nctype, ncent, ncent_tot, nctype_tot
    public :: symbol, atomtyp, allocate_atom, deallocate_atom
    public :: newghostype, nghostcent
    save

contains

    !> Allocates memory for atomic system geometry and properties.
    !>
    !> @details This subroutine allocates arrays for storing molecular geometry
    !> and nuclear properties. It checks for existing allocations before allocating
    !> to prevent memory leaks. The allocation sizes are based on ncent_tot (total
    !> number of centers including ghosts) and nctype_tot (total number of center types).
    !>
    !> Allocates the following arrays:
    !> - cent(3, ncent_tot): Cartesian coordinates (x, y, z) for all centers
    !> - znuc(nctype_tot): Nuclear charges for each center type
    !> - iwctype(nctype_tot): Wavefunction type indices, initialized to 0
    !> - symbol(ncent_tot): Element symbols for each center
    !>
    !> @note This subroutine must be called after ncent_tot and nctype_tot are set
    !> from input file parsing, typically in the parser module.
    subroutine allocate_atom()

        if (.not. allocated(cent)) allocate (cent(3, ncent_tot))
        if (.not. allocated(znuc)) allocate (znuc(nctype_tot))
        if (.not. allocated(iwctype)) allocate (iwctype(nctype_tot), source=0)
        if (.not. allocated(symbol)) allocate (symbol(ncent_tot))
    end subroutine allocate_atom

    !> Deallocates memory for atomic system geometry and properties.
    !>
    !> @details This subroutine frees all dynamically allocated arrays used for
    !> storing molecular geometry and nuclear properties. It checks for allocation
    !> status before deallocating to avoid runtime errors.
    !>
    !> Deallocates the following arrays:
    !> - iwctype: Wavefunction type indices
    !> - znuc: Nuclear charges
    !> - cent: Atomic center coordinates
    !>
    !> @note Called at program termination or when switching to a different system.
    subroutine deallocate_atom()
        if (allocated(iwctype)) deallocate (iwctype)
        if (allocated(znuc)) deallocate (znuc)
        if (allocated(cent)) deallocate (cent)
    end subroutine deallocate_atom

end module system