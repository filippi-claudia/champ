!> @brief Module for state-averaged energy checking and storage.
!> @author CHAMP developers
!> @date 2025
!>
!> @details This module stores energies and their uncertainties for multiple states
!> in state-averaged calculations. Used for comparing and validating energies across
!> different quantum states.
module sa_check
      use mstates_mod, only: MSTATES
      use precision_kinds, only: dp

     implicit none

     !> Energies for all states in state-averaged calculation.
     real(dp), dimension(:), allocatable :: energy_all

     !> Energy uncertainties (errors) for all states.
     real(dp), dimension(:), allocatable :: energy_err_all

     private
     public :: energy_all, energy_err_all
     public :: allocate_sa_check, deallocate_sa_check
     save
 contains

     !> Allocates memory for state-averaged energy arrays.
     subroutine allocate_sa_check()
      use mstates_mod, only: MSTATES
         if (.not. allocated(energy_all)) allocate (energy_all(MSTATES))
         if (.not. allocated(energy_err_all)) allocate (energy_err_all(MSTATES))
     end subroutine allocate_sa_check

     !> Deallocates memory for state-averaged energy arrays.
     subroutine deallocate_sa_check()
         if (allocated(energy_err_all)) deallocate (energy_err_all)
         if (allocated(energy_all)) deallocate (energy_all)
     end subroutine deallocate_sa_check

 end module sa_check

!> @brief Module for state-averaged calculation weights.
!> @author CHAMP developers
!> @date 2025
!>
!> @details This module manages weights for state-averaged calculations, including
!> weight indices and values for each state. The weights determine the contribution
!> of each state to the averaged energy and wavefunction optimization.
 module sa_weights
      use mstates_mod, only: MSTATES
      use precision_kinds, only: dp

     implicit none

     !> Weight indices for each state.
     integer, dimension(:), allocatable :: iweight

     !> Number of weighted states.
     integer :: nweight

     !> Weight values for each state in averaging.
     real(dp), dimension(:), allocatable :: weights

     private
     public :: iweight, nweight, weights
     public :: allocate_sa_weights, deallocate_sa_weights
     save
 contains

     !> Allocates memory for state-averaged weight arrays.
     subroutine allocate_sa_weights()
      use mstates_mod, only: MSTATES
         if (.not. allocated(iweight)) allocate (iweight(MSTATES), source=0)
         if (.not. allocated(weights)) allocate (weights(MSTATES))
     end subroutine allocate_sa_weights

     !> Deallocates memory for state-averaged weight arrays.
     subroutine deallocate_sa_weights()
         if (allocated(weights)) deallocate (weights)
         if (allocated(iweight)) deallocate (iweight)
     end subroutine deallocate_sa_weights

 end module sa_weights

!> @brief Main module for state-averaged calculations.
!> @author CHAMP developers
!> @date 2025
!>
!> @details This module provides a unified interface for allocating and deallocating
!> all state-averaged calculation data structures. It serves as a wrapper that manages
!> both energy checking arrays and state weights.
!>
!> State-averaged calculations are used for:
!> - Multi-reference calculations
!> - Excited state optimizations
!> - Conical intersection studies
!> - Transition state calculations
 module m_state_avrg
 contains

 !> Allocates all state-averaged calculation arrays.
 subroutine allocate_m_state_avrg()
      use sa_check, only: allocate_sa_check
      use sa_weights, only: allocate_sa_weights

     call allocate_sa_check()
     call allocate_sa_weights()
 end subroutine allocate_m_state_avrg

 !> Deallocates all state-averaged calculation arrays.
 subroutine deallocate_m_state_avrg()
      use sa_check, only: deallocate_sa_check
      use sa_weights, only: deallocate_sa_weights

     call deallocate_sa_check()
     call deallocate_sa_weights()
 end subroutine deallocate_m_state_avrg
 end module
