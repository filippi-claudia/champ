#if defined HAVE_CONFIG_H
#  include "config.h"
#endif

#define THIS_FILE "io_fdf.F90"
!>
!!=====================================================================
!!
!! This file is part of the FDF package.
!!
!! This module implements an interface to the FORTRAN logical unit
!! system. Based on code by Richard Maine.
!!
!! Logical unit management. Units 0 to min_lun-1 are "reserved",
!! since most of the "typical" files (output, etc) use them.
!!
!! Logical units min_lun to min_max are managed by this module.
!!
!!
!! @date September 2007
!!
!!
!!=====================================================================

#define ERROR_UNIT  0
#define OUTPUT_UNIT 6

MODULE io_fdf
  USE utils
  USE prec
  USE iso_fortran_env
  implicit none

! General callable functions
  public :: io_seterr, io_setout
  public :: io_geterr, io_getout

  public :: io_assign, io_reserve
  public :: io_close, io_status


! Error and Output Units
  integer(ip), private :: stderr = ERROR_UNIT,                          &
                          stdout = OUTPUT_UNIT

! Unit control variables
  integer(ip), parameter, private :: min_lun = 10, max_lun = 99
  integer(ip), parameter, private :: nunits = max_lun-min_lun+1
  logical, private                :: lun_is_free(min_lun:max_lun) = .TRUE.


  CONTAINS

!
!   Set IO error unit
!
    SUBROUTINE io_seterr(unit)
      implicit none
!-------------------------------------------------------------- Output Variables
      integer(ip), intent(inout) :: unit
!------------------------------------------------------------------------- BEGIN
      stderr = unit
      RETURN
!--------------------------------------------------------------------------- END
    END SUBROUTINE io_seterr

!
!   Set IO output unit
!
    SUBROUTINE io_setout(unit)
      implicit none
!-------------------------------------------------------------- Output Variables
      integer(ip), intent(inout) :: unit
!------------------------------------------------------------------------- BEGIN
      stdout = unit
      RETURN
!--------------------------------------------------------------------------- END
    END SUBROUTINE io_setout

!
!   Get IO error unit
!
    SUBROUTINE io_geterr(unit)
      implicit none
!-------------------------------------------------------------- Output Variables
      integer(ip), intent(inout) :: unit
!------------------------------------------------------------------------- BEGIN
      unit = stderr
      RETURN
!--------------------------------------------------------------------------- END
    END SUBROUTINE io_geterr

!
!   Get IO output unit
!
    SUBROUTINE io_getout(unit)
      implicit none
!-------------------------------------------------------------- Output Variables
      integer(ip), intent(inout) :: unit
!------------------------------------------------------------------------- BEGIN
      unit = stdout
      RETURN
!--------------------------------------------------------------------------- END
    END SUBROUTINE io_getout


!
!   Looks for a free unit and assigns it to lun
!
    SUBROUTINE io_assign(lun)
      implicit none
!-------------------------------------------------------------- Output Variables
      integer(ip), intent(inout) :: lun

!--------------------------------------------------------------- Local Variables
      logical                    :: used, found
      integer(ip)                :: i, iostat

!------------------------------------------------------------------------- BEGIN

      i     = min_lun
      found = .FALSE.
      do while((.not. found) .and. (i .le. max_lun))
        if (lun_is_free(i)) then

          INQUIRE(unit=i, opened=used, iostat=iostat)
          if (iostat .ne. 0) used = .TRUE.
          if (.not. used) then
            lun   = i
            found = .TRUE.
          endif

          lun_is_free(i) = .FALSE.
        endif

        i = i + 1
      enddo

      if (.not. found) then
        call die('IO module: io_assign', 'No LUNs available',           &
                 THIS_FILE, __LINE__)
      endif

      RETURN
!--------------------------------------------------------------------------- END
    END SUBROUTINE io_assign

!
!   Useful to specify that one needs to use a particular unit number
!   for example, assume some legacy code expects to work with unit 15:
!   call io_reserve(15)   ! this call at the beginning of the program
!   ...
!   open(15,....)
!
    SUBROUTINE io_reserve(lun)
      implicit none
!-------------------------------------------------------------- Output Variables
      integer(ip), intent(inout) :: lun

!--------------------------------------------------------------- Local Variables
      logical                    :: used
      character(80)              :: msg
      integer(ip)                :: iostat

!------------------------------------------------------------------------- BEGIN
      INQUIRE(unit=lun, opened=used, iostat=iostat)
      if (iostat .ne. 0) used = .TRUE.
      if (used) then
        write(msg,'(a,i3,a)')                                           &
             'Cannot reserve unit',lun,'. Already connected'
        call die('IO module: io_reserve', msg, THIS_FILE, __LINE__)
      endif

      if ((lun .ge. min_lun) .and. (lun .le. max_lun))                  &
        lun_is_free(lun) = .FALSE.

      RETURN
!--------------------------------------------------------------------------- END
    END SUBROUTINE io_reserve

!
!   Use this routine instead of a simple close
!
    SUBROUTINE io_close(lun)
      implicit none
!-------------------------------------------------------------- Output Variables
      integer(ip), intent(inout) :: lun

!------------------------------------------------------------------------- BEGIN
      CLOSE(lun)
      if ((lun .ge. min_lun) .and. (lun .le. max_lun))                  &
        lun_is_free(lun) = .TRUE.

      RETURN
!--------------------------------------------------------------------------- END
    END SUBROUTINE io_close

!
!   Prints a list of the connected logical units and the names of
!   the associated files
!
    SUBROUTINE io_status()
      implicit none
!--------------------------------------------------------------- Local Variables
      logical       :: opened, named
      character(80) :: filename
      character(11) :: form
      integer(ip)   :: i, iostat

!------------------------------------------------------------------------- BEGIN
      write(stdout,'(a)') '******** io_status ********'
      do i= 0, max_lun
        INQUIRE(i, opened=opened, named=named, name=filename,           &
                form=form, iostat=iostat)
        if (iostat .eq. 0) then
          if (opened) then
            if (named) then
              write(stdout,'(i4,5x,a,5x,a)') i, form, filename
            else
              write(stdout,'(i4,5x,a,5x,a)') i, form, 'No name available'
            endif
          endif
        else
          write(stdout,'(i4,5x,a,5x,a)') i, 'IOSTAT error'
        endif
      enddo
      write(stdout,'(a)') '********           ********'

      RETURN
!--------------------------------------------------------------------------- END
    END SUBROUTINE io_status

END MODULE io_fdf
