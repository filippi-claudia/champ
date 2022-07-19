module constants
  use precision_kinds, only: dp

  real(dp), parameter :: hb = 0.5
  real(dp), parameter :: pi = 4.0d0*datan(1.0d0)
  real(dp), parameter :: twopi = 8.d0*datan(1.0d0)
end module
