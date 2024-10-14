#ifdef QMCKL_FOUND
#include "qmckl_f.F90"
#endif

module qmckl_data
#ifdef QMCKL_FOUND

  use qmckl
  integer, parameter     :: qmckl_no_ctx_max = 2
  integer(qmckl_context) :: qmckl_ctx(qmckl_no_ctx_max)
  integer                :: qmckl_no_ctx 
  logical                :: use_qmckl = .True.

  public :: qmckl_no_ctx_max, qmckl_no_ctx, qmckl_ctx, use_qmckl
  save

#else

  integer :: qmckl_no_ctx_max = 0
  integer :: qmckl_ctx(1) = 0
  integer :: qmckl_no_ctx = 0
  integer :: QMCKL_SUCCESS = 0
  logical :: use_qmckl = .False.

  public :: qmckl_no_ctx_max, qmckl_no_ctx, qmckl_ctx, use_qmckl, QMCKL_SUCCESS
  save

#endif
end module

