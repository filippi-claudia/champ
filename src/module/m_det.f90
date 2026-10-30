!> @brief Module for determinant mapping parameters.
!> @author CHAMP developers
!> @date 2020
!>
!> @details This module stores the number of unique determinant mappings used
!> in multi-determinant wavefunctions. The mapping is used to efficiently store
!> and access determinant information in Configuration State Function (CSF) 
!> expansions.
!>
!> @note This parameter is set during wavefunction initialization and is used
!> by the csfs module for array allocation.
module dets
    implicit none

    !> Number of unique determinant-to-CSF mappings in the wavefunction expansion
    integer :: nmap

    save
end module dets

!> @brief Module for Configuration State Functions (CSFs) and determinant expansions.
!> @author CHAMP developers
!> @date 2020
!>
!> @details This module manages the expansion of the wavefunction in terms of
!> Configuration State Functions (CSFs) and their constituent Slater determinants.
!> It stores CSF coefficients, determinant indices for spin-up and spin-down
!> electrons, and mapping information between determinants and CSFs.
!>
!> Key data structures:
!> - CSF coefficients (ccsf, cxdet) for multiple states and geometries
!> - Determinant occupation indices (iadet, ibdet) for alpha and beta electrons
!> - CSF-determinant mapping (icxdet) and normalization (anormo)
!> - Counts (ncsf, nstates, maxcsf) for dimensioning
!>
!> @note Essential for multi-determinant wavefunctions including MCSCF and CASSCF.
!> Used in both VMC and DMC calculations for multi-reference wavefunctions.
module csfs
      use precision_kinds, only: dp

    implicit none

    !> CSF coefficients for each determinant, state, and wavefunction type, dimension (ndet, MSTATES, nwftype)
    real(dp), dimension(:, :, :), allocatable :: ccsf

    !> Expansion coefficients for determinant-to-CSF mapping, dimension (nmap)
    real(dp), dimension(:), allocatable :: cxdet

    !> Normalization constants for each electronic state, dimension (MSTATES)
    real(dp), dimension(:), allocatable :: anormo

    !> Determinant occupation indices for spin-up (alpha) electrons, dimension (ndet)
    integer, dimension(:), allocatable :: iadet

    !> Determinant occupation indices for spin-down (beta) electrons, dimension (ndet)
    integer, dimension(:), allocatable :: ibdet

    !> Indices for CSF-to-determinant mapping, dimension (nmap)
    integer, dimension(:), allocatable :: icxdet

    !> Maximum CSF index for each state, dimension (MSTATES)
    integer, dimension(:), allocatable :: maxcsf

    !> Total number of Configuration State Functions in the expansion
    integer :: ncsf

    !> Number of electronic states (ground and excited states)
    integer :: nstates

    private
    public   :: ccsf, cxdet, iadet, ibdet, icxdet, maxcsf, ncsf, nstates, anormo
    public :: allocate_csfs, deallocate_csfs
    save
contains
    !> Allocates memory for CSF and determinant expansion arrays.
    !>
    !> @details This subroutine allocates arrays for storing Configuration State
    !> Function coefficients, determinant occupation indices, and mapping information.
    !> The allocation sizes depend on the number of determinants (ndet), states
    !> (MSTATES), wavefunction types (nwftype), and mappings (nmap).
    !>
    !> Allocates the following arrays:
    !> - ccsf(ndet, MSTATES, nwftype): CSF coefficients for all determinants and states
    !> - cxdet(nmap): Determinant-to-CSF expansion coefficients
    !> - anormo(MSTATES): Normalization constants per state
    !> - iadet(ndet): Alpha electron occupation indices
    !> - ibdet(ndet): Beta electron occupation indices
    !> - icxdet(nmap): CSF-determinant mapping indices
    !> - maxcsf(MSTATES): Maximum CSF index per state
    !>
    !> @note Called during wavefunction initialization for multi-determinant calculations.
    subroutine allocate_csfs()
      use dets, only: nmap
      use mstates_mod, only: MSTATES
      use multiple_geo, only: nwftype
      use slater, only: ndet
        if (.not. allocated(ccsf)) allocate (ccsf(ndet, MSTATES, nwftype))
        if (.not. allocated(cxdet)) allocate (cxdet(nmap))
        if (.not. allocated(anormo)) allocate (anormo(MSTATES))
        if (.not. allocated(iadet)) allocate (iadet(ndet))
        if (.not. allocated(ibdet)) allocate (ibdet(ndet))
        if (.not. allocated(icxdet)) allocate (icxdet(nmap))
        if (.not. allocated(maxcsf)) allocate (maxcsf(MSTATES))
    end subroutine allocate_csfs

    !> Deallocates memory for CSF and determinant expansion arrays.
    !>
    !> @details This subroutine frees all arrays used for storing CSF coefficients
    !> and determinant information. Deallocates in reverse allocation order.
    !>
    !> @note Called at program termination to prevent memory leaks.
    subroutine deallocate_csfs()
        if (allocated(icxdet)) deallocate (icxdet)
        if (allocated(maxcsf)) deallocate (maxcsf)
        if (allocated(anormo)) deallocate (anormo)
        if (allocated(ibdet)) deallocate (ibdet)
        if (allocated(iadet)) deallocate (iadet)
        if (allocated(cxdet)) deallocate (cxdet)
        if (allocated(ccsf)) deallocate (ccsf)
    end subroutine deallocate_csfs

end module csfs

