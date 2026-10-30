!> @brief Module for trial energy values used in DMC population control.
!> @author CHAMP developers
!> @date 2022
!>
!> @details This module stores trial energy values that are critical for the
!> population control algorithm in Diffusion Monte Carlo (DMC) calculations.
!> The trial energy is used in the branching factor calculation to maintain
!> a stable walker population.
!>
!> Key variables:
!> - etrial: Current trial energy estimate
!> - esigmatrial: Trial energy with statistical uncertainty
!>
module const
      use precision_kinds, only: dp
    implicit none

    !> Trial energy estimate used for DMC branching and population control
    real(dp) :: etrial

    !> Trial energy with statistical uncertainty
    real(dp) :: esigmatrial

    save
end module const

!> @brief Module for fundamental physical and mathematical constants.
!> @author CHAMP developers
!> @date 2022
!>
!> @details This module defines fundamental physical and mathematical constants
!> used throughout the quantum Monte Carlo calculations. All constants are defined
!> as parameters (compile-time constants) to ensure consistency and prevent
!> accidental modification during runtime.
!>
!> Constants provided:
!> - hb: Reduced Planck constant (ℏ) in atomic units
!> - pi: Mathematical constant π
!> - twopi: Mathematical constant 2π
!>
!> @note In atomic units (a.u.), ℏ = 1, but here it is set to 0.5 for specific
!> computational purposes. The value of π is computed using the arctangent function
!> for maximum precision.
module constants
      use precision_kinds, only: dp
    implicit none

    !> Reduced Planck's constant (ℏ) in atomic units, set to 0.5 for computational purposes
    real(dp), parameter :: hb = 0.5

    !> Mathematical constant π, computed as 4*arctan(1) for precision
    real(dp), parameter :: pi = 4.0d0*datan(1.0d0)

    !> Mathematical constant 2π, computed as 8*arctan(1) for precision
    real(dp), parameter :: twopi = 8.d0*datan(1.0d0)

end module
