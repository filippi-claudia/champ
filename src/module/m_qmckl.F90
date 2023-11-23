#ifdef QMCKL_FOUND
#include "qmckl_f.F90"
#endif

module qmckl_data
#ifdef QMCKL_FOUND

  use qmckl
  integer(qmckl_context) :: qmckl_ctx
  logical                :: use_qmckl = .True.

  public :: qmckl_ctx, use_qmckl
  save

#else

  integer :: qmckl_ctx = 0
  integer :: QMCKL_SUCCESS = 0
  logical :: use_qmckl = .False.

  public :: qmckl_ctx, use_qmckl, QMCKL_SUCCESS
  save

#endif
end module

