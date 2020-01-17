module numeric_kinds
  ! named constants for 4, 2, and 1 byte integers:
  integer, parameter :: &
       i4b = selected_int_kind(9), &
       i2b = selected_int_kind(4), &
       i1b = selected_int_kind(2)
  ! single, double and quadruple precision reals:
  integer, parameter :: &
       sp = kind(1.0), &
       dp = selected_real_kind(2 * precision(1.0_sp)), &
       qp = selected_real_kind(2 * precision(1.0_dp))
end module numeric_kinds
