!> @brief QMCkl library interface for quantum chemistry kernels.
!> @author CHAMP developers
!> @date 2025
!>
!> @details This file conditionally includes the QMCkl Fortran interface when the
!> QMCkl library is found during compilation. QMCkl is a library providing high-performance
!> implementations of quantum chemistry kernels for quantum Monte Carlo calculations.
!>
!> The QMCkl library provides:
!> - Optimized orbital calculations
!> - Jastrow factor computations
!> - Distance and derivative calculations
!> - Platform-specific optimizations (SIMD, GPU)
!>
!> @note This file requires the QMCKL_FOUND preprocessor flag to be defined
!> during compilation to include the QMCkl interface and enable functionality.
!>
!> @see https://github.com/TREX-CoE/qmckl for more information about QMCkl.

#ifdef QMCKL_FOUND
#include "qmckl_f.F90"
#endif

!> @brief Module for QMCkl context management and configuration.
!> @author CHAMP developers
!> @date 2025
!>
!> @details This module manages QMCkl computation contexts and configuration flags.
!> When QMCkl is available, it provides multiple contexts for parallel computation
!> and flags to enable/disable specific QMCkl features (orbitals, Jastrow).
!> When QMCkl is not available, it provides stub definitions for compatibility.
module qmckl_data
#ifdef QMCKL_FOUND

  use qmckl
  
  implicit none

  !> Maximum number of QMCkl contexts (fixed at 4).
  integer, parameter     :: qmckl_no_ctx_max = 4

  !> Array of QMCkl context handles for parallel computation.
  integer(qmckl_context) :: qmckl_ctx(qmckl_no_ctx_max)

  !> Actual number of QMCkl contexts currently in use.
  integer                :: qmckl_no_ctx

  !> Flag to enable QMCkl-based orbital calculations.
  logical                :: use_qmckl_orbitals = .True.

  !> Flag to enable QMCkl-based Jastrow factor calculations.
  logical                :: use_qmckl_jastrow = .True.

  public :: qmckl_no_ctx_max, qmckl_no_ctx, qmckl_ctx, use_qmckl_orbitals, use_qmckl_jastrow
  save

#else

  implicit none

  !> Maximum number of QMCkl contexts (0 when library not available).
  integer, parameter :: qmckl_no_ctx_max = 0

  !> Dummy context array when QMCkl is not available.
  integer :: qmckl_ctx(1) = 0

  !> Number of QMCkl contexts (0 when library not available).
  integer :: qmckl_no_ctx = 0

  !> QMCkl success return code (stub when library not available).
  integer :: QMCKL_SUCCESS = 0

  !> Flag indicating QMCkl orbitals are disabled (library not available).
  logical :: use_qmckl_orbitals = .False.

  !> Flag indicating QMCkl Jastrow is disabled (library not available).
  logical :: use_qmckl_jastrow = .False.

  public :: qmckl_no_ctx_max, qmckl_no_ctx, qmckl_ctx, use_qmckl_orbitals, use_qmckl_jastrow, QMCKL_SUCCESS
  save

#endif
end module

