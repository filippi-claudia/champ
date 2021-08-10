#if defined HAVE_CONFIG_H
#  include "config.h"
#endif
!>
!!=====================================================================
!!
!! This file is part of the FDF package.
!!
!! This module provides useful functions and subroutines for FDF library.
!! At this moment this module contains functions for:
!!
!!   a) String manipulation
!!   b) Warning, Die (Abort/Terminate) operations
!!
!!
!! September 2007
!!
!!
!!=====================================================================

#define ERROR_UNIT  0
#define OUTPUT_UNIT 6

MODULE utils
  USE prec
  implicit none

! String functions
  public :: leqi, leqi_strict
  public :: labeleq, packlabel
  public :: chrcap, chrlen

! Conversors between formats
  public :: s2i, s2r, arr2s, s2arr, i2s
  public :: convert_string_to_array_of_chars
  public :: convert_array_of_chars_to_string

! Warning and Terminate functions
  public :: warn, die

! Maximum size of a string
  integer(ip), parameter, public :: MAX_LENGTH = 256

  CONTAINS

!
!   Case-insensitive lexical equal-to comparison
!
    FUNCTION leqi(string1, string2)
      implicit none
!--------------------------------------------------------------- Input Variables
      character(len=*) :: string1, string2

!-------------------------------------------------------------- Output Variables
      logical          :: leqi

!--------------------------------------------------------------- Local Variables
      logical          :: completed
      character        :: char1, char2
      integer(ip)      :: i, len1, len2, lenc

!------------------------------------------------------------------------- BEGIN
      len1 = len(string1)
      len2 = len(string2)
      lenc = min(len1, len2)

      i = 1
      leqi      = .TRUE.
      completed = .FALSE.
      do while((.not. completed) .and. (i .le. lenc))
        char1 = string1(i:i)
        char2 = string2(i:i)
        call chrcap(char1, 1)
        call chrcap(char2, 1)
        if (char1 .ne. char2) then
          leqi      = .FALSE.
          completed = .TRUE.
        endif

        i = i + 1
      enddo

      if (leqi) then
        if ((len1 .gt. lenc) .and. (string1(lenc+1:len1) .ne. ' '))     &
          leqi = .FALSE.
        if ((len2 .gt. lenc) .and. (string2(lenc+1:len2) .ne. ' '))     &
          leqi = .FALSE.
      endif

      RETURN
!--------------------------------------------------------------------------- END
    END FUNCTION leqi

!
!   Examples of eq_func's for search function (Case sensitive)
!
    FUNCTION leqi_strict(str1, str2)
      implicit none
!--------------------------------------------------------------- Input Variables
      character(len=*) :: str1, str2

!-------------------------------------------------------------- Output Variables
      logical          :: leqi_strict

!------------------------------------------------------------------------- BEGIN
      leqi_strict = (str1 .eq. str2)
      RETURN
!--------------------------------------------------------------------------- END
    END FUNCTION leqi_strict

!
!   Compares s1 and s2 without regard for case, or appearance
!   of '_', '.', '-'.
!
    FUNCTION labeleq(s1, s2, logunit)
      implicit none
!--------------------------------------------------------------- Input Variables
      character(*)                    :: s1, s2
      integer(ip), optional           :: logunit

!-------------------------------------------------------------- Output Variables
      logical                         :: labeleq

!--------------------------------------------------------------- Local Variables
      character(max(len(s1),len(s2))) :: n1, n2

!------------------------------------------------------------------------- BEGIN
      call packlabel(s1, n1)
      call packlabel(s2, n2)
      labeleq = leqi(n1, n2)

      if (PRESENT(logunit) .and. labeleq .and.                          &
          (.not. leqi(s1, s2))) then
  !!        write(logunit,'(a,/,a,/,a)')                                    &
  !!           '--------- Considered equivalent:', s1, s2
      endif

      RETURN
!--------------------------------------------------------------------------- END
    END FUNCTION labeleq

!
!   Removes occurrences of '_ .-'  from s
!
    SUBROUTINE packlabel(s, n)
      implicit none
!--------------------------------------------------------------- Input Variables
      character(*) :: s

!-------------------------------------------------------------- Output Variables
      character(*) :: n

!--------------------------------------------------------------- Local Variables
      character   :: c
      integer(ip) :: i, j

      logical     :: is_sep
      is_sep(i) = (i.eq.95) .or. (i.eq.46) .or. (i.eq.45)

!------------------------------------------------------------------------- BEGIN
      n = ' '
      j = 0
      do i= 1, len(s)
        c = s(i:i)
        if (.not. is_sep(ichar(c))) then
          j = j+1
          n(j:j) = c
        endif
      enddo

      RETURN
!--------------------------------------------------------------------------- END
    END SUBROUTINE packlabel

!
!   CHRCAP accepts a STRING of NCHAR characters and replaces
!   any lowercase letters by uppercase ones.
!
    SUBROUTINE chrcap(string, nchar)
      implicit none
!--------------------------------------------------------------- Input Variables
      integer(ip)  :: nchar

!-------------------------------------------------------------- Output Variables
      character(*) :: string

!--------------------------------------------------------------- Local Variables
      integer(ip)  :: i, itemp, ncopy

!------------------------------------------------------------------------- BEGIN
      if (nchar .le. 0) then
        ncopy = LEN(string)
      else
        ncopy = nchar
      endif

      do i= 1, ncopy
        if (LGE(string(i:i),'a') .and. LLE(string(i:i),'z')) then
          itemp = ICHAR(string(i:i)) + ICHAR('A') - ICHAR('a')
          string(i:i) = CHAR(itemp)
        endif
      enddo

      RETURN
!--------------------------------------------------------------------------- END
    END SUBROUTINE chrcap

!
!   CHRLEN accepts a STRING of NCHAR characters and returns LCHAR,
!   the length of the string up to the last NONBLANK, NONNULL.
!
    SUBROUTINE chrlen(string, nchar, lchar)
      implicit none
!--------------------------------------------------------------- Input Variables
      character(*) :: string
      integer(ip)  :: nchar

!-------------------------------------------------------------- Output Variables
      integer(ip)  :: lchar

!------------------------------------------------------------------------- BEGIN
      lchar = nchar
      if (lchar .le. 0) lchar = LEN(string)

      do while(((string(lchar:lchar) .eq. ' ') .or. (string(lchar:lchar) &
               .eq. CHAR(0))) .and. (lchar .gt. 0))
        lchar = lchar - 1
      enddo

      RETURN
!--------------------------------------------------------------------------- END
    END SUBROUTINE chrlen

!
!   String to Integer translator
!
    FUNCTION s2i(string)
      implicit none
!--------------------------------------------------------------- Input Variables
      character(*), intent(in) :: string

!-------------------------------------------------------------- Output Variables
!      integer(ip)              :: s2i
      integer(i16)             :: s2i

!--------------------------------------------------------------- Local Variables
      integer(ip)              :: ierr

!------------------------------------------------------------------------- BEGIN
      read(string, fmt=*, iostat=ierr) s2i
      if (ierr .ne. 0) then
        call die('UTILS module: s2i', 'Integer conversion error',       &
                 'utils.F90', __LINE__, ERROR_UNIT)
      endif

      RETURN
!--------------------------------------------------------------------------- END
    END FUNCTION s2i

!
!   String to Real translator
!
    FUNCTION s2r(string)
      implicit none
!--------------------------------------------------------------- Input Variables
      character(*), intent(in) :: string

!-------------------------------------------------------------- Output Variables
      real(dp)                 :: s2r

!--------------------------------------------------------------- Local Variables
      integer(ip)              :: ierr

!------------------------------------------------------------------------- BEGIN
      read(string, fmt=*, iostat=ierr) s2r
      if (ierr .ne. 0) then
        call die('UTILS module: s2r', 'Real conversion error',          &
                 'utils.F90', __LINE__, ERROR_UNIT)
      endif

      RETURN
!--------------------------------------------------------------------------- END
    END FUNCTION s2r

!
!   Converts from Array of Characters to String
!
    FUNCTION arr2s(string_arr, string_size)
      implicit none
!--------------------------------------------------------------- Input Variables
      character                 :: string_arr(*)
      integer(ip)               :: string_size

!-------------------------------------------------------------- Output Variables
      character(len=MAX_LENGTH) :: arr2s

!--------------------------------------------------------------- Local Variables
      character(len=MAX_LENGTH) :: str
      character                 :: str_arr(MAX_LENGTH)
      equivalence                  (str, str_arr)

!------------------------------------------------------------------------- BEGIN
      str = ' '
      str_arr(1:string_size) = string_arr(1:string_size)
      arr2s = TRIM(str)
      RETURN
!--------------------------------------------------------------------------- END
    END FUNCTION arr2s

!
!   Arbitrary length version of s2arr, but with matching sizes
!
    subroutine convert_string_to_array_of_chars(str,arr)
    character(len=*), intent(in) :: str
    character, dimension(:), intent(out) :: arr

    integer :: n, i

    n = len(str)
    if (size(arr) /= n) call die("convert_str_to_arr","Size mismatch")
    do i = 1, n
       arr(i) = str(i:i)
    enddo
  end subroutine convert_string_to_array_of_chars
!
!   Arbitrary length version of arr2s, but with matching sizes
!
    subroutine convert_array_of_chars_to_string(arr,str)
    character, dimension(:), intent(in) :: arr
    character(len=*), intent(out) :: str

    integer :: n, i

    n = size(arr)
    if (len(str) /= n) call die("convert_arr_to_str","Size mismatch")
    do i = 1, n
       str(i:i) = arr(i)
    enddo
  end subroutine convert_array_of_chars_to_string
!
!   Converts from String to Array of Characters
!
    FUNCTION s2arr(string)
      implicit none
!--------------------------------------------------------------- Input Variables
      character(len=*)          :: string

!-------------------------------------------------------------- Output Variables
      character                 :: s2arr(MAX_LENGTH)

!--------------------------------------------------------------- Local Variables
      character(len=MAX_LENGTH) :: str
      character                 :: str_arr(MAX_LENGTH)
      equivalence                  (str, str_arr)

!------------------------------------------------------------------------- BEGIN
      str = ' '
      str(1:LEN_TRIM(string)) = string(1:LEN_TRIM(string))
      s2arr = str_arr
      RETURN
!--------------------------------------------------------------------------- END
    END FUNCTION s2arr

!
!   Converts an integer number to string
!
    FUNCTION i2s(num)
      implicit none
!--------------------------------------------------------------- Input Variables
      integer(ip)  :: num

!-------------------------------------------------------------- Output Variables
      character(5) :: i2s

!--------------------------------------------------------------- Local Variables
      integer(ip)  :: i, ntmp, zero
      character    :: cc(5)

!------------------------------------------------------------------------- BEGIN
      if (num > 99999 .OR. num < 0) then
        call die('UTILS module: i2s', 'Number is out of range',         &
                 'utils.F90', __LINE__, ERROR_UNIT)
      endif

      zero = ICHAR('0')  ! 48 is the ascii code of zero
      ntmp = num
      do i= 5, 1, -1
        cc(i) = CHAR(zero + MOD(ntmp,10))
        ntmp = ntmp/10
      enddo
      i2s = cc(1)//cc(2)//cc(3)//cc(4)//cc(5)

      RETURN
!--------------------------------------------------------------------------- END
    END function i2s

!
!   Warning routine
!
    SUBROUTINE warn(string)
      implicit none
!--------------------------------------------------------------- Input Variables
      character(len=*) :: string

!------------------------------------------------------------------------- BEGIN
      write(OUTPUT_UNIT,'(a,a)') '*** WARNING: ', TRIM(string)

      RETURN
!--------------------------------------------------------------------------- END
    END SUBROUTINE warn

!
!   Die routine (Abort/Terminate program)
!   MPI-awareness (MPI_Abort)
!
    SUBROUTINE die(routine, msg, file, line, unit, rc, cline)
#if defined(CLUSTER) || defined(BLOCKING)
      !use mpi_siesta
      use mpi
#endif
      implicit none
!--------------------------------------------------------------- Input Variables
      character(len=*), intent(in)           :: routine, msg
      character(len=*), intent(in), optional :: file, cline
      integer(ip), intent(in), optional      :: line, unit, rc

!--------------------------------------------------------------- Local Variables
      integer(ip)                            :: die_unit
#if defined(CLUSTER) || defined(BLOCKING)
      integer(ip)                            :: ierr, rc2
#endif

!------------------------------------------------------------------------- BEGIN
      if (PRESENT(unit)) then
        die_unit = unit
      else
        die_unit = ERROR_UNIT
      endif

      write(die_unit,'(a)') '*************************************************************'
      write(die_unit,'(a)') 'ERROR'
      write(die_unit,'(a)') ' '
      write(die_unit,'(3a)') TRIM(routine), ': ', TRIM(msg)
      write(die_unit,'(a)') ' '
      if (PRESENT(cline)) write(die_unit,'(5x,2a)') 'Input line: ', trim(cline)
      if (PRESENT(file)) write(die_unit,'(5x,2a)') 'File: ', trim(file)
      if (PRESENT(line)) write(die_unit,'(5x,a,i5)') 'Line: ', line
      write(die_unit,'(a)') '*************************************************************'

      if (die_unit .ne. ERROR_UNIT) then
        write(die_unit,'(a)') 'Stopping Program'
      endif

#if defined(CLUSTER) || defined(BLOCKING)
      if (.not. PRESENT(rc)) then
        rc2 = 0
      else
        rc2 = rc
      endif

      call MPI_Abort(MPI_COMM_WORLD, rc2, ierr)
#else
      STOP 'Stopping Program'
#endif
!--------------------------------------------------------------------------- END
    END SUBROUTINE die

END MODULE utils
