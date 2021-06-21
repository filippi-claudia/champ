#if defined HAVE_CONFIG_H
#  include "config.h"
#endif
!>
!!=====================================================================
!!
!! This file is part of the FDF package.
!!
!! This module provides precision for integer and reals in FDF library.
!! At this moment this module contains precision specification for:
!!
!!   a) Integer precision     (ip)
!!   b) Single Real precision (sp)
!!   c) Double Real precision (dp)
!!
!!
!! September 2007
!!
!!
!!=====================================================================

MODULE prec

!
! Precision handling
! Kind parameters
!
  integer, parameter :: ip  = selected_int_kind(9)
  integer, parameter :: i16 = selected_int_kind(16)
  integer, parameter :: sp  = selected_real_kind(6,20)
  integer, parameter :: dp  = selected_real_kind(14,100)
END MODULE prec
