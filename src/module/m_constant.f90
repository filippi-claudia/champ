
!> Module containing esigma and esigmatrial
module const
      use precision_kinds, only: dp
    implicit none

    !> esigma
    real(dp) :: etrial

    !> esigmatrial
    real(dp) :: esigmatrial

    save
end module const

!> Module containing physical constants
module constants
      use precision_kinds, only: dp
    implicit none

    !> Plank's constant
    real(dp), parameter :: hb = 0.5

    !> Value of Pi
    real(dp), parameter :: pi = 4.0d0*datan(1.0d0)

    !> Value of 2*Pi
    real(dp), parameter :: twopi = 8.d0*datan(1.0d0)

end module
