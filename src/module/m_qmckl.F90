#ifdef QMCKL_FOUND
#include "qmckl_f.F90"
#endif

module qmckl_data
#ifdef QMCKL_FOUND

  use qmckl
  integer(qmckl_context), public :: qmckl_ctx
  logical, public :: use_qmckl = .True.

#else

  integer, public :: qmckl_ctx = 0
  integer, public :: QMCKL_SUCCESS = 0
  logical, public :: use_qmckl = .False.

#endif
end module
