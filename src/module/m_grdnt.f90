!> @brief Module for storing Hessian matrix in Z-matrix coordinates.
!> @author CHAMP developers
!> @date 2022
!>
!> @details This module stores the Hessian matrix (second derivatives of energy)
!> expressed in Z-matrix (internal) coordinates. The Hessian is used for geometry
!> optimization, transition state searches, and vibrational frequency calculations.
!>
!> Z-matrix coordinates (bonds, angles, dihedrals) are more natural for molecular
!> structures than Cartesian coordinates, particularly for describing molecular
!> vibrations and reaction pathways.
!>
!> The hessian_zmat array stores the second derivatives ∂²E/∂q_i∂q_j where q are
!> internal coordinates.
!>
!> @note Dimension: (3, ncent_tot) for storing Hessian components.
module grdnthes
      use precision_kinds, only: dp

     implicit none

     !> Hessian matrix in Z-matrix (internal) coordinates (3, ncent_tot)
     real(dp), dimension(:, :), allocatable :: hessian_zmat

     private
     public   ::  hessian_zmat
     public :: allocate_grdnthes, deallocate_grdnthes
     save
 contains
     !> @brief Allocate Hessian matrix array in Z-matrix coordinates.
     !> @details Allocates hessian_zmat(3, ncent_tot) for storing second derivatives
     !> of energy with respect to internal coordinates. Array is uninitialized.
     subroutine allocate_grdnthes()
      use system, only: ncent_tot
         if (.not. allocated(hessian_zmat)) allocate (hessian_zmat(3, ncent_tot))
     end subroutine allocate_grdnthes

     !> @brief Deallocate Hessian matrix array.
     !> @details Frees memory for hessian_zmat array.
     subroutine deallocate_grdnthes()
         if (allocated(hessian_zmat)) deallocate (hessian_zmat)
     end subroutine deallocate_grdnthes

 end module grdnthes

!> @brief Module for grdntsmv
!> @author CHAMP developers
!> @date 2022
!>
!> @details This module manages grdntsmv. 
 module grdntsmv
     use multiple_geo, only: MFORCE

     implicit none

     !> igrdaidx (dimensions MFORCE)
     integer, dimension(:), allocatable :: igrdaidx
     
     !> igrdcidx (dimensions MFORCE)
     integer, dimension(:), allocatable :: igrdcidx
     
     !> igrdmv (dimensions (3, MFORCE))
     integer, dimension(:, :), allocatable :: igrdmv

     private
     public :: igrdaidx, igrdcidx, igrdmv
     public :: allocate_grdntsmv, deallocate_grdntsmv
     save
 contains
     !> @brief Allocate various gradient arrays.
     !> @details Allocates igrdaidx(MFORCE) and igrdcidx(MFORCE), and igrdmv(3, ncent_tot) All arrays initialized to zero.
     subroutine allocate_grdntsmv()
      use system, only: ncent_tot
      use multiple_geo, only: MFORCE
         if (.not. allocated(igrdaidx)) allocate (igrdaidx(MFORCE), source=0)
         if (.not. allocated(igrdcidx)) allocate (igrdcidx(MFORCE), source=0)
         if (.not. allocated(igrdmv)) allocate (igrdmv(3, ncent_tot), source=0)
     end subroutine allocate_grdntsmv

     !> @brief Deallocate gradient move arrays.
     !> @details Frees memory for igrdmv, igrdcidx, and igrdaidx arrays.
     subroutine deallocate_grdntsmv()
         if (allocated(igrdmv)) deallocate (igrdmv)
         if (allocated(igrdcidx)) deallocate (igrdcidx)
         if (allocated(igrdaidx)) deallocate (igrdaidx)
     end subroutine deallocate_grdntsmv

 end module grdntsmv

!> @brief Module for gradient calculation parameters and step sizes.
!> @author CHAMP developers
!> @date 2022
!>
!> @details This module stores parameters controlling numerical gradient calculations,
!> including displacement step sizes for different coordinate types. These parameters
!> determine the finite-difference stencil used to compute derivatives numerically.
!>
!> Step size parameters for different coordinate types:
!> - delgrdxyz: Step size for Cartesian coordinate displacements
!> - delgrdbl: Step size for bond length displacements
!> - delgrdba: Step size for bond angle displacements
!> - delgrdda: Step size for dihedral angle displacements
!>
!> Control parameters:
!> - igrdtype: Type of gradient calculation
!> - ngradnts: Total number of gradient components to compute
!>
 module grdntspar
      use precision_kinds, only: dp

     implicit none

     !> Step size for bond angle gradient displacements
     real(dp) :: delgrdba
     
     !> Step size for bond length gradient displacements
     real(dp) :: delgrdbl
     
     !> Step size for dihedral angle gradient displacements
     real(dp) :: delgrdda
     
     !> Step size for Cartesian coordinate gradient displacements
     real(dp) :: delgrdxyz
     
     !> Type of gradient calculation
     integer :: igrdtype
     
     !> Total number of gradient components to compute
     integer :: ngradnts

     private
     public :: delgrdba, delgrdbl, delgrdda, delgrdxyz, igrdtype, ngradnts
     save
 end module grdntspar

!> @brief Master module for coordinating gradient calculation memory management.
!> @author CHAMP developers
!> @date 2022
!>
!> @details This module provides top-level routines for managing memory allocation
!> and deallocation for all gradient-related data structures. It coordinates allocation
!> across multiple sub-modules to ensure consistent initialization and cleanup for
!> gradient, force, and Hessian calculations.
!>
!> Sub-modules managed:
!> - grdnthes: Hessian matrix in Z-matrix coordinates (hessian_zmat)
!> - grdntsmv: Gradient move indices and flags (igrdaidx, igrdcidx, igrdmv)
!> - grdntspar: Gradient step sizes and control parameters (static, no allocation)
!>
module m_grdnt
contains
 !> @brief Allocate all gradient-related data structures.
 !> @details Calls allocation routines for gradient sub-modules:
 !> 1. allocate_grdnthes(): Hessian matrix in Z-matrix coordinates
 !> 2. allocate_grdntsmv(): Gradient move indices and flags
 subroutine allocate_m_grdnt()
      use grdnthes, only: allocate_grdnthes
      use grdntsmv, only: allocate_grdntsmv

     implicit none

     call allocate_grdnthes()
     call allocate_grdntsmv()
 end subroutine allocate_m_grdnt

 !> @brief Deallocate all gradient-related data structures.
 !> @details Calls deallocation routines for gradient sub-modules:
 !> 1. deallocate_grdnthes(): Hessian matrix arrays
 !> 2. deallocate_grdntsmv(): Gradient move arrays
 subroutine deallocate_m_grdnt()
      use grdnthes, only: deallocate_grdnthes
      use grdntsmv, only: deallocate_grdntsmv

     implicit none

     call deallocate_grdnthes()
     call deallocate_grdntsmv()
 end subroutine deallocate_m_grdnt
 end module 
