#ifdef QMCKL_FOUND
#include "qmckl_f.F90"
#endif

!> Module that contains the qmckl context
module qmckl_data
#ifdef QMCKL_FOUND

  use qmckl
  !> The maximum amount of QMCKL contexts
  integer, parameter     :: qmckl_no_ctx_max = 4
  !> QMCKL context
  integer(qmckl_context) :: qmckl_ctx(qmckl_no_ctx_max)
  !> The amount of QMCKL contexts
  integer                :: qmckl_no_ctx
  !> Logical to check if QMCKL is used
  logical                :: use_qmckl_orbitals = .True.
  logical                :: use_qmckl_jastrow = .True.
  logical                :: first_step = .True.

  public :: qmckl_no_ctx_max, qmckl_no_ctx, qmckl_ctx, use_qmckl_orbitals, use_qmckl_jastrow
  save

#else

  !> The maximum amount of QMCKL contexts
  integer, parameter :: qmckl_no_ctx_max = 0
  !> QMCKL context
  integer :: qmckl_ctx(1) = 0
  !> The amount of QMCKL contexts
  integer :: qmckl_no_ctx = 0
  !> Logical to check if QMCKL is successfully linked
  integer :: QMCKL_SUCCESS = 0
  !> Logical to check if QMCKL is used
  logical :: use_qmckl_orbitals = .False.
  logical :: use_qmckl_jastrow = .False.
  public :: qmckl_no_ctx_max, qmckl_no_ctx, qmckl_ctx, use_qmckl_orbitals, use_qmckl_jastrow, QMCKL_SUCCESS
  save

#endif
end module

