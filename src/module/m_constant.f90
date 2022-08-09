module const
    use precision_kinds, only: dp
    implicit none

    real(dp) :: etrial

    save
end module const

module constants
    use precision_kinds, only: dp
    implicit none
  
    real(dp), parameter :: hb = 0.5
    real(dp), parameter :: pi = 4.0d0*datan(1.0d0)
    real(dp), parameter :: twopi = 8.d0*datan(1.0d0)
    
end module
