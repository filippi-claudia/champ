!> @brief Module defining maximum array dimensions for external electric field charges.
!> @author CHAMP developers
!> @date 2020
!>
!> @details This module provides compile-time constants for the maximum number
!> of external point charges that can be included in electric field calculations.
!> These charges are used to model external electrostatic potentials acting on
!> the quantum system.
!>
!> @note MCHARGES may need to be increased for systems with many external charges.
module efield_mod

     implicit none

     !> Maximum number of external point charges for electric field calculations
     integer, parameter :: MCHARGES = 100
     private
     public :: MCHARGES
     save
 end module efield_mod

!> @brief Module for external electric field control parameters.
!> @author CHAMP developers
!> @date 2020
!>
!> @details This module stores control parameters for external electric field
!> calculations. It specifies whether an electric field is applied, the screening
!> model to use, and the number of external point charges.
!>
!> Key parameters:
!> - iefield: Switch to enable/disable electric field calculations
!> - iscreen: Screening model selector (0=no screening, >0=screened potential)
!> - ncharges: Actual number of external point charges in the calculation
!>
!> @note These parameters are typically read from the input file during initialization.
 module efield

     implicit none

     !> Flag to enable (1) or disable (0) external electric field calculations
     integer :: iefield

     !> Screening model selector: 0=unscreened, 1=exponential screening, etc.
     integer :: iscreen

     !> Actual number of external point charges (must be ≤ MCHARGES)
     integer :: ncharges

     private
     public :: iscreen, ncharges, iefield
     save
 end module efield

!> @brief Module for external point charge positions, magnitudes, and screening parameters.
!> @author CHAMP developers
!> @date 2022
!>
!> @details This module stores the coordinates, charges, and screening parameters
!> for external point charges used in electric field calculations. These charges
!> create an external electrostatic potential that affects the electronic structure.
!>
!> Key arrays:
!> - xcharge, ycharge, zcharge: Cartesian coordinates of each charge (in bohr)
!> - qcharge: Magnitude of each charge (in atomic units)
!> - ascreen, bscreen: Screening parameters for distance-dependent screening
!>
!> @note All arrays have dimension MCHARGES, but only the first ncharges entries
!> are used. Screening models: exponential screening V(r) = q*exp(-a*r^b)/r.
 module efield_blk
     use efield_mod, only: MCHARGES
     use precision_kinds, only: dp

     implicit none

     !> Screening parameter 'a' for exponential screening function, dimension (MCHARGES)
     real(dp), dimension(:), allocatable :: ascreen

     !> Screening parameter 'b' for exponential screening function, dimension (MCHARGES)
     real(dp), dimension(:), allocatable :: bscreen

     !> Magnitude of each external point charge in atomic units, dimension (MCHARGES)
     real(dp), dimension(:), allocatable :: qcharge

     !> X-coordinate of each external point charge in bohr, dimension (MCHARGES)
     real(dp), dimension(:), allocatable :: xcharge

     !> Y-coordinate of each external point charge in bohr, dimension (MCHARGES)
     real(dp), dimension(:), allocatable :: ycharge

     !> Z-coordinate of each external point charge in bohr, dimension (MCHARGES)
     real(dp), dimension(:), allocatable :: zcharge

     private
     public :: zcharge, bscreen, qcharge, ycharge, xcharge, ascreen
     public :: allocate_efield_blk, deallocate_efield_blk
     save
 contains
     !> Allocates memory for external charge arrays.
     !>
     !> @details This subroutine allocates arrays for storing positions, charges,
     !> and screening parameters for up to MCHARGES external point charges. All
     !> arrays are allocated with fixed dimension MCHARGES.
     !>
     !> Allocates:
     !> - ascreen, bscreen: Screening parameters
     !> - qcharge: Charge magnitudes
     !> - xcharge, ycharge, zcharge: Cartesian coordinates
     !>
     !> @note Called during electric field initialization if iefield is enabled.
     subroutine allocate_efield_blk()
      use efield_mod, only: MCHARGES
         if (.not. allocated(ascreen)) allocate (ascreen(MCHARGES))
         if (.not. allocated(bscreen)) allocate (bscreen(MCHARGES))
         if (.not. allocated(qcharge)) allocate (qcharge(MCHARGES))
         if (.not. allocated(xcharge)) allocate (xcharge(MCHARGES))
         if (.not. allocated(ycharge)) allocate (ycharge(MCHARGES))
         if (.not. allocated(zcharge)) allocate (zcharge(MCHARGES))
     end subroutine allocate_efield_blk

     !> Deallocates memory for external charge arrays.
     !>
     !> @details Frees all arrays used for external point charge storage.
     subroutine deallocate_efield_blk()
         if (allocated(zcharge)) deallocate (zcharge)
         if (allocated(ycharge)) deallocate (ycharge)
         if (allocated(xcharge)) deallocate (xcharge)
         if (allocated(qcharge)) deallocate (qcharge)
         if (allocated(bscreen)) deallocate (bscreen)
         if (allocated(ascreen)) deallocate (ascreen)
     end subroutine deallocate_efield_blk

 end module efield_blk

!> @brief Master module for coordinated electric field memory management.
!> @author CHAMP developers
!> @date 2020
!>
!> @details This module provides unified allocation and deallocation interfaces
!> for all electric field related modules. It currently manages memory for the
!> efield_blk module containing external point charge data.
!>
!> @note This is the primary interface for electric field memory management,
!> called from allocate_vmc() when external electric fields are enabled.
module m_efield
contains
 !> Allocates memory for all electric field modules.
 !>
 !> @details This subroutine calls allocation routines for electric field arrays:
 !> - allocate_efield_blk(): External point charge positions and parameters
 !>
 !> @note Called from allocate_vmc() when iefield is enabled in the input.
 subroutine allocate_m_efield()
      use efield_blk, only: allocate_efield_blk

     implicit none

     call allocate_efield_blk()
 end subroutine allocate_m_efield

 !> Deallocates memory for all electric field modules.
 !>
 !> @details This subroutine calls deallocation routines to free electric field
 !> arrays at program termination.
 !>
 !> @note Called from deallocate_vmc() at program termination.
 subroutine deallocate_m_efield()
      use efield_blk, only: deallocate_efield_blk

     implicit none

     call deallocate_efield_blk()
 end subroutine deallocate_m_efield
end module 
