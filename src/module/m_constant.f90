module const
    !> Arguments: pi, hb, etrial, delta, deltai, fbias, nelec, imetro, ipr
      use precision_kinds, only: dp

    implicit none

    real(dp) :: delta
    real(dp) :: deltai
    real(dp) :: etrial
    real(dp) :: fbias
    integer  :: imetro
    private
    public   :: etrial, delta, deltai, fbias, imetro
    save
end module const

module const2
    !> Arguments: deltar, deltat
      use precision_kinds, only: dp

    implicit none

    real(dp) :: deltar
    real(dp) :: deltat

    private
    public :: deltar, deltat
    save
end module const2

module constants
      use precision_kinds, only: dp
  
    real(dp), parameter :: hb = 0.5
    real(dp), parameter :: pi = 4.0d0*datan(1.0d0)
    real(dp), parameter :: twopi = 8.d0*datan(1.0d0)
    
  end module