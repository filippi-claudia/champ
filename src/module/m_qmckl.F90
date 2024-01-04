#ifdef QMCKL_FOUND
#include "qmckl_gpu_f.f90"

#endif

module qmckl_data
#ifdef QMCKL_FOUND

  use qmckl_gpu_f
  integer(qmckl_context_device) :: qmckl_ctx
  logical                :: use_qmckl = .True.
  !integer :: QMCKL_SUCCESS_DEVICE = 0
  !public :: qmckl_ctx, use_qmckl, QMCKL_SUCCESS_DEVICE
  public :: qmckl_ctx, use_qmckl
  save

#else

  integer :: qmckl_ctx = 0
  integer :: QMCKL_SUCCESS_DEVICE = 0
  logical :: use_qmckl = .False.

  public :: qmckl_ctx, use_qmckl, QMCKL_SUCCESS_DEVICE
  save

#endif
end module

