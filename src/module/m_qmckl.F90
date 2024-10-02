#ifdef QMCKL_FOUND
#include "qmckl_f.F90"
#endif

!> Module that contains the qmckl context
module qmckl_data
#ifdef QMCKL_FOUND

  use qmckl

  !> QMCKL context
  integer(qmckl_context) :: qmckl_ctx

  !> Logical to check if QMCKL is used
  logical                :: use_qmckl = .True.

  public :: qmckl_ctx, use_qmckl
  save

#else

  !> QMCKL context
  integer :: qmckl_ctx = 0

  !> Logical to check if QMCKL is successfully linked
  integer :: QMCKL_SUCCESS = 0

  !> Logical to check if QMCKL is used
  logical :: use_qmckl = .False.

  public :: qmckl_ctx, use_qmckl, QMCKL_SUCCESS
  save

#endif
end module

