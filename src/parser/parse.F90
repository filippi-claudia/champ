#if defined HAVE_CONFIG_H
#  include "config.h"
#endif

#define THIS_FILE "parse.F90"
!>
!!=====================================================================
!!
!! This file is part of the FDF package.
!!
!! This module provides a simple yet powerful way to analyze the information
!! in a string (such as an input line).
!!
!! Routine, 'digest' takes as input a string 'line' and returns a pointer
!! to a derived type 'parsed_line':
!!
!!   Parsed line info (ntokens, token info and identification)
!!   Note that the token characters are stored in a single "line",
!!   and addressed using the starting and ending points.
!!   This avoids the use of dynamic memory without loss of functionality.
!!
!!  type, public :: parsed_line
!!    integer(ip)               :: ntokens
!!    character(len=MAX_LENGTH) :: line
!!    integer(ip)               :: first(MAX_NTOKENS)
!!    integer(ip)               :: last(MAX_NTOKENS)
!!    character(len=1)          :: id(MAX_NTOKENS)
!!  end type parsed_line
!!
!! which holds a list of tokens and token tags (id). The
!! parsing (split string into tokens) is done by a helper routine
!! 'parses' which currently behaves according to the FDF standard.
!! Each token is classified by helper routine 'morphol' and a token
!! id is assigned in the following way:
!!
!! * Tokens that can be read as real numbers are assigned to class
!! 'values' and given a token id 'v'. These are further classified as
!! 'integers' (id 'i') or 'reals' (id 'r').
!! * There are two list classes:
!!    'a' == integer list
!!    'c' == real list
!!    'e' == real or integer list
!! * All other tokens are tagged as 'names' (id 'n').
!!
!! The recommended usage follows the outline:
!!
!!     use parse
!!     character(len=?) line
!!     type(parsed_line), pointer :: p
!!     ...
!!     p=>digest(line)
!!     (extract information from p)
!!     call destroy(p)
!!
!! Note the pointer assignment and the explicit call to a destroyer
!! routine that frees the storage associated to p.
!!
!! The information is extracted by module procedures that fall into three
!! classes:
!!
!! a) Enquiry functions: 'search' and 'match'
!!
!! *  'search' determines whether a token in 'line' matches the given
!!    string, optionally returning an index. The search is
!!    case-insensitive by default, but this can be changed by supplying
!!    an extra procedure argument 'eq_func' with interface:
!!
!!       interface
!!         function eq_func(s1,s2)
!!           logical eq_func
!!           character(len=*), intent(in) :: s1,s2
!!         end function eq_func
!!       end interface
!!
!!    We have two different implementations of 'search' function,
!!    through a wrapper (function overload):
!!
!!       interface search
!!         module procedure search_fun
!!         module procedure search_sub
!!       end interface
!!
!!     1. %FUNCTION search_fun(string, pline_fun, after, eq_func)
!!       New search implementation. 'search' function returns
!!       the index of the token that matches with the string or
!!       -1 if not found. Leaves 'pline_fun' structure pointing
!!       to the token in the FDF structure.
!!
!!     2. %FUNCTION search_sub(pline_sub, string, ind, after, eq_func)
!!       This is the old prototype for backward compatibility.
!!       Returns .TRUE. if the string is found in the parsed line
!!       else .FALSE. Moreover can return the index of the token
!!       in the line if 'ind' is specified.
!!
!!    Example:  if (search('Mary', p) .ne. -1) ...
!!    will return the index of the first token that matches
!!    "Mary", or -1 if not found.
!!
!!    This function can take an optional keyword 'after=' (see below).
!!
!! *  'substring_search' does not match whole tokens, but substrings in
!!    tokens. And it uses the *case sensitive* Fortran 'index' function.
!
!! *  'match' is probably the most powerful routine in the module. It
!!    checks whether the token morphology of 'line' conforms to the
!!    sequence of characters specified. For example,
!!
!!    if (match(p,'nii')) ...
!!
!!    returns .TRUE. if 'line' contains at least three tokens and they are
!!    a 'name' and two 'integers'.
!!    Apart from the 'primitive' one-character ids, there is the
!!    possibility of using 'compound' virtual ids for generalized matchings:
!!
!!    - A 'v' ('value') is matched by both an 'integer' and a 'real'.
!!    - A 'j' is matched by both an 'integer' and a 'name'.
!!    - A 's' is matched by an 'integer', a 'real', and a 'name'.
!!    - A 'x' is matched by any kind of token.
!!    - A 'a' is matched by a list with integers
!!    - A 'c' is matched by a list with reals
!!    - A 'e' is matched by a list with integers or reals
!!    - A 'd' is reserved for future dictionaries...
!
!!    This function can take an optional keyword 'after=' (see below).
!!
!! b) Number functions: ntokens ('n|i|r|b|e|l|a'), nnames ('n'), nreals ('r'),
!!                      nintegers ('i'), nvalues ('i|r'), nblocks ('b'),
!!                      nendblocks ('e'), nlabels ('l'), nlists('a|c'),
!!                      nintegerlists ('a'), nreallists('c')
!!
!!    These functions return the number of tokens of each kind in 'line':
!!
!!    number_of_energies = nreals(p)
!!
!!    These functions can take an optional keyword 'after=' (see below).
!!
!! c) Extraction functions: tokens ('n|i|r|b|e|l|a|c'), names ('n'), reals ('r'),
!!                          characters,
!!                          integers ('i'), values ('i|r'), blocks ('b'),
!!                          endblocks ('e'), labels ('l'),
!!                          integerlists('a') <- a subroutine
!!                          reallists('c') <- a subroutine
!!                          valuelists('a|c') <- a subroutine
!!
!!    These functions return a piece of data which corresponds to a token
!!    of the specified kind with sequence number matching the index
!!    provided. For example,
!!
!!    nlevels = integers(p,2)
!!
!!    assigns to variable 'nlevels' the second integer in 'line'.
!!    Execution stops in the assignment cannot be made. The user should
!!    call the corresponding 'number' routine to make sure there are
!!    enough tokens of the given kind.
!!
!!    Function 'characters' returns a string of characters spanning
!!    several tokens (with the original whitespace)
!!
!!    These functions can take an optional keyword 'after=' (see below).
!!
!!
!! By default, the routines in the module perform any indexing from the
!! beginning of 'line', in such a way that the first token is assigned the
!! index 1. It is possible to specify a given token as 'origin' by using
!! the 'after=' optional keyword. For example:
!!
!!     if (search(p, 'P', ind=jp)) then            # Old implementation
!!       if (match(p, 'i', after=jp) npol = integers(p, 1, after=jp)
!!     endif
!!
!! first checks whether 'P' is found in 'line'. If so, 'match' is used to
!! check whether it is followed by at least an 'integer'. If so, its
!! valued is assigned to variable 'npol'.
!!
!! If the 'after=' optional keyword is used in routine 'search', the
!! returned index is absolute, not relative. For example, to get the
!! real number coming right after the first 'Q' which appears to the
!! right of the 'P' found above:
!!
!!     if (search(p, 'Q', ind=jq, after=jp)) then  # Old implementation
!!       if (match(p, 'r', after=jq) energy = reals(p, 1, after=jq)
!!     endif
!!
!! @authors Alberto Garcia, 1995-2007, original implementation
!! @authors Raul de la Cruz, September 2007
!! @authors Alberto Garcia, July 2008
!! @remarks Modification made to libfdf by
!! @authors Ravindra Shinde (r.l.shinde@utwente.nl)
!! @date (2021)
!!========================================================================

#define ERROR_UNIT  0

MODULE parse
  USE utils
  USE prec
  implicit none


! Serialization functions
  public :: serialize_pline, recreate_pline

! Internal functions: build parsed line and morphology
  private :: create
  private :: parses, morphol

! Digest, match and search
  public :: digest, destroy
  public :: match, search, substring_search

! Routines to get number and items
  public :: nintegers, nreals, nvalues, nnames
  public :: nblocks, nendblocks, nlabels, ntokens
  public :: integers, reals, values, names
  public :: blocks, endblocks, labels, tokens, characters
  public :: nlists, nintegerlists, nreallists
  public :: integerlists, reallists, valuelists


! Change morphology
  public :: setmorphol

! Integer|Real check routines
  private :: is_integer, is_value

! Debugging config routines
  public :: setdebug, setlog

! Internal constants
  logical, private                :: parse_debug = .FALSE.
  integer(ip), private            :: parse_log   = ERROR_UNIT
  integer(ip), parameter, private :: MAX_NTOKENS = 100

! Length of string encoding plines
  integer, parameter, public :: SERIALIZED_LENGTH =  MAX_LENGTH + 4 + 10*MAX_NTOKENS

!   Parsed line info (ntokens, token info and identification)
!   Note that the token characters are stored in a single "line",
!   and addressed using the starting and ending points.
!   This avoids the use of dynamic memory without loss of functionality.

  type, public :: parsed_line
    integer(ip)               :: ntokens
    character(len=MAX_LENGTH) :: line
    integer(ip)               :: first(MAX_NTOKENS)
    integer(ip)               :: last(MAX_NTOKENS)
    character(len=1)          :: id(MAX_NTOKENS)
  end type parsed_line

! Search wrapper (return index as function or subroutine)
  interface search
    module procedure search_fun
    module procedure search_sub
  end interface

  CONTAINS

!
!   Creates parsed_line structure
!
    SUBROUTINE create(pline)
      implicit none
!------------------------------------------------ Output Variables
      type(parsed_line), pointer :: pline

!------------------------------------------------- Local Variables
      integer(ip)                :: ierr

!----------------------------------------------------------- BEGIN
      if (ASSOCIATED(pline)) call destroy(pline)
      ALLOCATE(pline, stat=ierr)
      if (ierr .ne. 0) then
        call die('PARSE module: create', 'Error allocating pline',      &
                 THIS_FILE, __LINE__, rc=ierr)
      endif
!------------------------------------------------------------- END
    END SUBROUTINE create

!
!   Frees parsed_line structure
!
    SUBROUTINE destroy(pline)
      implicit none
!------------------------------------------------ Output Variables
      type(parsed_line), pointer :: pline

!----------------------------------------------------------- BEGIN
      if (ASSOCIATED(pline)) then
        DEALLOCATE(pline)
        NULLIFY(pline)
      endif
!------------------------------------------------------------- END
    END SUBROUTINE destroy

!
!   Return the number of items of a certain class among the tokens.
!
    FUNCTION nitems(class, pline, after)
      implicit none
!------------------------------------------------- Input Variables
      character                         :: class
      integer(ip), intent(in), optional :: after
      type(parsed_line), pointer        :: pline

!------------------------------------------------ Output Variables
      integer(ip)                       :: nitems

!------------------------------------------------- Local Variables
      character(80)                     :: msg
      integer(ip)                       :: i, starting_pos

!----------------------------------------------------------- BEGIN
      if (PRESENT(after)) then
        if (after .lt. 0) then
          write(msg,*) 'Wrong starting position when processing class: ', &
                       class
          call die('PARSE module: nitems', msg, THIS_FILE, __LINE__)
        endif
        starting_pos = after
      else
        starting_pos = 0
      endif

      nitems = 0
      do i= starting_pos+1, pline%ntokens
        if (leqi(pline%id(i), class)) nitems = nitems + 1
      enddo

      RETURN
!------------------------------------------------------------- END
    END FUNCTION nitems

!
!   Return the number of integers in the tokens.
!
    FUNCTION nintegers(pline, after)
      implicit none
!------------------------------------------------- Input Variables
      integer(ip), intent(in), optional :: after
      type(parsed_line), pointer        :: pline

!------------------------------------------------ Output Variables
      integer(ip)                       :: nintegers

!----------------------------------------------------------- BEGIN
      nintegers = nitems('i', pline, after)

      RETURN
!------------------------------------------------------------- END
    END FUNCTION nintegers

!
!   Return the number of reals in the tokens.
!
    FUNCTION nreals(pline, after)
      implicit none
!------------------------------------------------- Input Variables
      integer(ip), intent(in), optional :: after
      type(parsed_line), pointer        :: pline

!------------------------------------------------ Output Variables
      integer(ip)                       :: nreals

!----------------------------------------------------------- BEGIN
      nreals = nitems('r', pline, after)

      RETURN
!------------------------------------------------------------- END
    END FUNCTION nreals

!
!   Return the number of values in the tokens.
!
    FUNCTION nvalues(pline, after)
      implicit none
!------------------------------------------------- Input Variables
      integer(ip), intent(in), optional :: after
      type(parsed_line), pointer        :: pline

!------------------------------------------------ Output Variables
      integer(ip)                       :: nvalues

!----------------------------------------------------------- BEGIN
      nvalues = nitems('i', pline, after) + nitems('r', pline, after)

      RETURN
!------------------------------------------------------------- END
    END FUNCTION nvalues

!
!   Return the number of lists in the tokens.
!
    FUNCTION nlists(pline, after)
      implicit none
!------------------------------------------------- Input Variables
      integer(ip), intent(in), optional :: after
      type(parsed_line), pointer        :: pline

!------------------------------------------------ Output Variables
      integer(ip)                       :: nlists

!----------------------------------------------------------- BEGIN
      nlists = nitems('a', pline, after) + nitems('c', pline, after)

      RETURN
!------------------------------------------------------------- END
    END FUNCTION nlists

!
!   Return the number of integer lists in the tokens.
!
    FUNCTION nintegerlists(pline, after)
      implicit none
!------------------------------------------------- Input Variables
      integer(ip), intent(in), optional :: after
      type(parsed_line), pointer        :: pline

!------------------------------------------------ Output Variables
      integer(ip)                       :: nintegerlists

!----------------------------------------------------------- BEGIN
      nintegerlists = nitems('a', pline, after)

      RETURN
!------------------------------------------------------------- END
    END FUNCTION nintegerlists

!
!   Return the number of real lists in the tokens.
!
    FUNCTION nreallists(pline, after)
      implicit none
!------------------------------------------------- Input Variables
      integer(ip), intent(in), optional :: after
      type(parsed_line), pointer        :: pline

!------------------------------------------------ Output Variables
      integer(ip)                       :: nreallists

!----------------------------------------------------------- BEGIN
      nreallists = nitems('c', pline, after)

      RETURN
!------------------------------------------------------------- END
    END FUNCTION nreallists

!
!   Return the number of names in the tokens.
!
    FUNCTION nnames(pline, after)
      implicit none
!------------------------------------------------- Input Variables
      integer(ip), intent(in), optional :: after
      type(parsed_line), pointer        :: pline

!------------------------------------------------ Output Variables
      integer(ip)                       :: nnames

!----------------------------------------------------------- BEGIN
      nnames = nitems('n', pline, after)

      RETURN
!------------------------------------------------------------- END
    END FUNCTION nnames

!
!   Return the number of blocks in the tokens.
!
    FUNCTION nblocks(pline, after)
      implicit none
!------------------------------------------------- Input Variables
      integer(ip), intent(in), optional :: after
      type(parsed_line), pointer        :: pline

!------------------------------------------------ Output Variables
      integer(ip)                       :: nblocks

!----------------------------------------------------------- BEGIN
      nblocks = nitems('b', pline, after)

      RETURN
!------------------------------------------------------------- END
    END FUNCTION nblocks

!
!   Return the number of endblocks in the tokens.
!
    FUNCTION nendblocks(pline, after)
      implicit none
!------------------------------------------------- Input Variables
      integer(ip), intent(in), optional :: after
      type(parsed_line), pointer        :: pline

!------------------------------------------------ Output Variables
      integer(ip)                       :: nendblocks

!----------------------------------------------------------- BEGIN
      nendblocks = nitems('e', pline, after)

      RETURN
!------------------------------------------------------------- END
    END FUNCTION nendblocks

!
!   Return the number of labels in the tokens.
!
    FUNCTION nlabels(pline, after)
      implicit none
!------------------------------------------------- Input Variables
      integer(ip), intent(in), optional :: after
      type(parsed_line), pointer        :: pline

!------------------------------------------------ Output Variables
      integer(ip)                       :: nlabels

!----------------------------------------------------------- BEGIN
      nlabels = nitems('l', pline, after)

      RETURN
!------------------------------------------------------------- END
    END FUNCTION nlabels

!
!   Return the number of tokens.
!
    FUNCTION ntokens(pline, after)
      implicit none
!------------------------------------------------- Input Variables
      integer(ip), intent(in), optional :: after
      type(parsed_line), pointer        :: pline

!------------------------------------------------ Output Variables
      integer(ip)                       :: ntokens

!------------------------------------------------- Local Variables
      integer(ip)                       :: starting_pos

!----------------------------------------------------------- BEGIN
      if (PRESENT(after)) then
        if (after .lt. 0) then
          call die('PARSE module: ntokens', 'Wrong starting position',  &
                   THIS_FILE, __LINE__)
        endif
        starting_pos = after
      else
        starting_pos = 0
      endif

      ntokens = pline%ntokens - starting_pos
      if (ntokens .lt. 0) then
        call die('PARSE module: ntokens', 'Wrong starting position',    &
                 THIS_FILE, __LINE__)
      endif

      RETURN
!------------------------------------------------------------- END
    END FUNCTION ntokens

!
!   Return a given integer token, specifying it by its sequence
!   number. It is also possible to make the sequence start after
!   a given token number in the line.
!
    FUNCTION integers(pline, ind, after)
      implicit none
!------------------------------------------------- Input Variables
      integer(ip), intent(in)           :: ind
      integer(ip), intent(in), optional :: after
      type(parsed_line), pointer        :: pline

!------------------------------------------------ Output Variables
      integer(ip)                       :: integers

!------------------------------------------------- Local Variables
      logical                           :: found
      integer(ip)                       :: i, j, starting_pos

!----------------------------------------------------------- BEGIN
      if (PRESENT(after)) then
        if (after .lt. 0) then
          call die('PARSE module: integers', 'Wrong starting position', &
                   THIS_FILE, __LINE__, cline=characters(pline,1,-1))
        endif
        starting_pos = after
      else
        starting_pos = 0
      endif

      i = starting_pos+1
      j = 0
      found = .FALSE.
      do while((.not. found) .and. (i .le. pline%ntokens))
        if (leqi(pline%id(i), 'i')) j = j + 1
        if (j .eq. ind) then
          integers = s2i(tokens(pline,i))
          found = .TRUE.
        endif
        i = i + 1
      enddo

      if (.not. found) then
        call die('PARSE module: integers', 'Not enough integers in line', &
                 THIS_FILE, __LINE__,cline=characters(pline,1,-1))
      endif

      RETURN
!------------------------------------------------------------- END
    END FUNCTION integers

!
!   Return a given real token, specifying it by its sequence
!   number. It is also possible to make the sequence start after
!   a given token number in the line.
!
    FUNCTION reals(pline, ind, after)
      implicit none
!------------------------------------------------- Input Variables
      integer(ip), intent(in)           :: ind
      integer(ip), intent(in), optional :: after
      type(parsed_line), pointer        :: pline

!------------------------------------------------ Output Variables
      real(dp)                          :: reals

!------------------------------------------------- Local Variables
      logical                           :: found
      integer(ip)                       :: i, j, starting_pos

!----------------------------------------------------------- BEGIN
      if (PRESENT(after)) then
        if (after .lt. 0) then
          call die('PARSE module: reals', 'Wrong starting position',    &
                   THIS_FILE, __LINE__,cline=characters(pline,1,-1))
        endif
        starting_pos = after
      else
        starting_pos = 0
      endif

      i = starting_pos+1
      j = 0
      found = .FALSE.
      do while((.not. found) .and. (i .le. pline%ntokens))
        if (leqi(pline%id(i), 'r')) j = j + 1
        if (j .eq. ind) then
          reals = s2r(tokens(pline,i))
          found = .TRUE.
        endif
        i = i + 1
      enddo

      if (.not. found) then
        call die('PARSE module: reals', 'Not enough reals in line',     &
                 THIS_FILE, __LINE__,cline=characters(pline,1,-1))
      endif

      RETURN
!------------------------------------------------------------- END
    END FUNCTION reals

!
!   Return a given list token, specifying it by its sequence
!   number. It is also possible to make the sequence start after
!   a given token number in the line.
!
    SUBROUTINE reallists(pline, ind, nv, list, after)
      implicit none
!------------------------------------------------- Input Variables
      integer(ip), intent(in)           :: ind
      integer(ip), intent(in), optional :: after
      type(parsed_line), pointer        :: pline

!------------------------------------------------ Output Variables
      integer(ip)                       :: nv
      real(dp)                          :: list(nv)

!------------------------------------------------- Local Variables
      logical                           :: found
      integer(ip)                       :: i, j, starting_pos

      character(len=MAX_LENGTH)         :: llist, sep
      type(parsed_line), pointer        :: lpline
      integer(ip)                       :: iR
      real(dp)                          :: lR, uR, sR
      integer(ip)                       :: ri
      integer(ip)                       :: ti, lprev, li
      logical                           :: count, is_del

!----------------------------------------------------------- BEGIN
      if (PRESENT(after)) then
        if (after .lt. 0) then
          call die('PARSE module: reallists', 'Wrong starting position',    &
                   THIS_FILE, __LINE__)
        endif
        starting_pos = after
      else
        starting_pos = 0
      endif

      i = starting_pos+1
      j = 0
      found = .FALSE.
      count = .FALSE.
      do while((.not. found) .and. (i .le. pline%ntokens))
        if (leqi(pline%id(i), 'c')) j = j + 1
        if (j .eq. ind) then

           found = .TRUE.

           ! Parse token list
           llist = tokens(pline,i)
           ! The list does have the markers attached (remove them)
           li = len_trim(llist)-1
           llist = trim(llist(2:li))
           lpline => digest(llist)

           ! We now have converted the list into a
           ! parseable line

           ! Does the user request length?
           count = nv <= 0
           li = 0 ! counter for the number of items in the list
           ti = 1 ! the current token iterator
           is_del = .false.
           do while ( ti < lpline%ntokens )

              ! First we need to check whether we have a list delimiter next
              if (leqi(lpline%id(ti+1),'n')) then
                 sep = names(lpline,1,after=ti)
                 is_del = leqi(sep,',') ! apparently ',' is not a token
              end if

              if (leqi(lpline%id(ti+1),'n').and. .not. is_del) then
                 ! We have a range

                 if ( lpline%ntokens <= ti + 1 ) then
                    call die('PARSE module: reallists', 'Missing end range', &
                         THIS_FILE, __LINE__)
                 end if
                 if ( .not.scan(lpline%id(ti),'ir')>0 .or. &
                      .not.scan(lpline%id(ti+2),'ir')>0 ) then
                    call die('PARSE module: reallists', 'Range is not well-defined', &
                         THIS_FILE, __LINE__)
                 end if

                 ! grab the seperator
                 sep = names(lpline,1,after=ti)
                 if ( leqi(sep,'to') .or. leqi(sep,':') .or. &
                      leqi(sep,'--') .or. leqi(sep,'---') ) then

                    ! Sort the range
                    lR = values(lpline,1,after=ti-1)
                    uR = values(lpline,1,after=ti+1)
                    sR = 1._dp

                    ! Figure out if we have a step in the range
                    if ( ti + 3 < lpline%ntokens ) then
                       if ( leqi(lpline%id(ti+3),'n') .and. &
                          scan(lpline%id(ti+4),'ir')>0 ) then
                          sep = names(lpline,1,after=ti+2)
                          if ( leqi(sep,'step') ) then
                             sR = values(lpline,1,after=ti+3)
                             ! step after the 'step <val>'
                             ti = ti + 2
                          end if
                       end if
                    end if

                    ! Correct sign of stepper
                    if ( lR <= uR ) sR = abs(sR)
                    if ( uR <  lR ) sR = -abs(sR)
                    if ( sR == 0._dp ) call die('PARSE module: reallists', &
                         'Stepping a list cannot be stepped by 0', &
                         THIS_FILE, __LINE__ )
                    ! By adding 0.01 % we should capture a large
                    ! percentage of ill-defined ranges
                    !    lR = 1. ; uR = 1.9999 ; sR = 0.5
                    do iR = 0, int( (uR - lR) / sR + sR * 0.0001_dp )
                       call add_exit(count,li,nv,lR + sR * iR)
                    end do

                    ! jump across the range
                    ti = ti + 2
                 else
                    call die('PARSE module: reallists', 'Unknown token in list', &
                         THIS_FILE, __LINE__)
                 end if

              elseif (scan(lpline%id(ti),'ir')>0) then

                 call add_exit(count,li,nv,values(lpline,1,after=ti-1))
              end if

              ti = ti + 1
              if ( is_del ) then
                 ti = ti + 1
                 is_del = .false.
              end if

           end do

           ! Read last element (or the only element if one is given)
           if ( ti == lpline%ntokens ) then
              if (leqi(lpline%id(ti),'v')) then
                 call add_exit(count,li,nv,values(lpline,1,after=ti-1))
              end if
           end if

           ! Clean-up parsed list-line
           call destroy(lpline)

           if ( count ) then
             ! User explicitly asked for acount
             nv = li
           else if ( nv /= li ) then
             ! Update the number of elements returned
             nv = li
           end if

        endif
        i = i + 1
      enddo

      if (.not. found) then
        call die('PARSE module: reallists', 'Not enough lists in line', &
                 THIS_FILE, __LINE__)
      end if

      RETURN
!------------------------------------------------------------- END

    contains

      subroutine add_exit(is_counting,idx,nv,val)
        logical, intent(in) :: is_counting
        integer, intent(inout) :: idx
        integer, intent(in) :: nv
        real(dp), intent(in) :: val
        idx = idx + 1
        if ( .not. is_counting ) then
           if ( idx > nv ) then
              found = .false.
           else
              list(idx) = val
           end if
        end if
      end subroutine add_exit

    END SUBROUTINE reallists


!
!   Return a given list token, specifying it by its sequence
!   number. It is also possible to make the sequence start after
!   a given token number in the line.
!
    SUBROUTINE valuelists(pline, ind, nv, list, after)
      implicit none
!------------------------------------------------- Input Variables
      integer(ip), intent(in)           :: ind
      integer(ip), intent(in), optional :: after
      type(parsed_line), pointer        :: pline

!------------------------------------------------ Output Variables
      integer(ip)                       :: nv
      real(dp)                          :: list(nv)

!------------------------------------------------- Local Variables
      logical                           :: found
      integer(ip)                       :: i, j, starting_pos

      character(len=MAX_LENGTH)         :: llist, sep
      type(parsed_line), pointer        :: lpline
      integer(ip)                       :: iR
      real(dp)                          :: lR, uR, sR
      integer(ip)                       :: ri
      integer(ip)                       :: ti, lprev, li
      logical                           :: count, is_del

!----------------------------------------------------------- BEGIN
      if (PRESENT(after)) then
        if (after .lt. 0) then
          call die('PARSE module: valuelists', 'Wrong starting position',    &
                   THIS_FILE, __LINE__)
        endif
        starting_pos = after
      else
        starting_pos = 0
      endif

      i = starting_pos+1
      j = 0
      found = .FALSE.
      count = .FALSE.
      do while((.not. found) .and. (i .le. pline%ntokens))
        if (scan(pline%id(i), 'ac')>0) j = j + 1
        if (j .eq. ind) then

           found = .TRUE.

           ! Parse token list
           llist = tokens(pline,i)
           ! The list does have the markers attached (remove them)
           li = len_trim(llist)-1
           llist = trim(llist(2:li))
           lpline => digest(llist)

           ! We now have converted the list into a
           ! parseable line

           ! Does the user request length?
           count = nv <= 0
           li = 0 ! counter for the number of items in the list
           ti = 1 ! the current token iterator
           is_del = .false.
           do while ( ti < lpline%ntokens )

              ! First we need to check whether we have a list delimiter next
              if (leqi(lpline%id(ti+1),'n')) then
                 sep = names(lpline,1,after=ti)
                 is_del = leqi(sep,',') ! apparently ',' is not a token
              end if

              if (leqi(lpline%id(ti+1),'n').and. .not. is_del) then
                 ! We have a range

                 if ( lpline%ntokens <= ti + 1 ) then
                    call die('PARSE module: valuelists', 'Missing end range', &
                         THIS_FILE, __LINE__)
                 end if
                 if ( .not.scan(lpline%id(ti),'ir')>0 .or. &
                      .not.scan(lpline%id(ti+2),'ir')>0 ) then
                    call die('PARSE module: valuelists', 'Range is not well-defined', &
                         THIS_FILE, __LINE__)
                 end if

                 ! grab the seperator
                 sep = names(lpline,1,after=ti)
                 if ( leqi(sep,'to') .or. leqi(sep,':') .or. &
                      leqi(sep,'--') .or. leqi(sep,'---') ) then

                    ! Sort the range
                    lR = values(lpline,1,after=ti-1)
                    uR = values(lpline,1,after=ti+1)
                    sR = 1._dp

                    ! Figure out if we have a step in the range
                    if ( ti + 3 < lpline%ntokens ) then
                       if ( leqi(lpline%id(ti+3),'n') .and. &
                          scan(lpline%id(ti+4),'ir')>0 ) then
                          sep = names(lpline,1,after=ti+2)
                          if ( leqi(sep,'step') ) then
                             sR = values(lpline,1,after=ti+3)
                             ! step after the 'step <val>'
                             ti = ti + 2
                          end if
                       end if
                    end if

                    ! Correct sign of stepper
                    if ( lR <= uR ) sR = abs(sR)
                    if ( uR <  lR ) sR = -abs(sR)
                    if ( sR == 0._dp ) call die('PARSE module: valuelists', &
                         'Stepping a list cannot be stepped by 0', &
                         THIS_FILE, __LINE__ )
                    ! By adding 0.01 % we should capture a large
                    ! percentage of ill-defined ranges
                    !    lR = 1. ; uR = 1.9999 ; sR = 0.5
                    do iR = 0, int( (uR - lR) / sR + sR * 0.0001_dp )
                       call add_exit(count,li,nv,lR + sR * iR)
                    end do

                    ! jump across the range
                    ti = ti + 2
                 else
                    call die('PARSE module: valuelists', 'Unknown token in list', &
                         THIS_FILE, __LINE__)
                 end if

              elseif (scan(lpline%id(ti),'ir')>0) then

                 call add_exit(count,li,nv,values(lpline,1,after=ti-1))
              end if

              ti = ti + 1
              if ( is_del ) then
                 ti = ti + 1
                 is_del = .false.
              end if

           end do

           ! Read last element (or the only element if one is given)
           if ( ti == lpline%ntokens ) then
              if (scan(lpline%id(ti),'ir')>0) then
                 call add_exit(count,li,nv,values(lpline,1,after=ti-1))
              end if
           end if

           ! Clean-up parsed list-line
           call destroy(lpline)

           if ( count ) then
             ! User explicitly asked for acount
             nv = li
           else if ( nv /= li ) then
             ! Update the number of elements returned
             nv = li
           end if

        endif
        i = i + 1
      enddo

      if (.not. found) then
        call die('PARSE module: valuelists', 'Not enough lists in line', &
                 THIS_FILE, __LINE__)
      end if

      RETURN
!------------------------------------------------------------- END

    contains

      subroutine add_exit(is_counting,idx,nv,val)
        logical, intent(in) :: is_counting
        integer, intent(inout) :: idx
        integer, intent(in) :: nv
        real(dp), intent(in) :: val
        idx = idx + 1
        if ( .not. is_counting ) then
           if ( idx > nv ) then
              found = .false.
           else
              list(idx) = val
           end if
        end if
      end subroutine add_exit

    END SUBROUTINE valuelists


!
!   Return a given integer list token, specifying it by its sequence
!   number. It is also possible to make the sequence start after
!   a given token number in the line.
!
    SUBROUTINE integerlists(pline, ind, ni, list, after)
      implicit none
!------------------------------------------------- Input Variables
      integer(ip), intent(in)           :: ind
      integer(ip), intent(in), optional :: after
      type(parsed_line), pointer        :: pline

!------------------------------------------------ Output Variables
      integer(ip)                       :: ni, list(ni)

!------------------------------------------------- Local Variables
      logical                           :: found
      integer(ip)                       :: i, j, starting_pos

      character(len=MAX_LENGTH)         :: llist, sep
      type(parsed_line), pointer        :: lpline
      integer(ip)                       :: iR, lR, uR, sR
      integer(ip)                       :: ti, lprev, li
      logical                           :: count, is_del

!----------------------------------------------------------- BEGIN
      if (PRESENT(after)) then
        if (after .lt. 0) then
          call die('PARSE module: lists', 'Wrong starting position',    &
               THIS_FILE, __LINE__)
        endif
        starting_pos = after
      else
        starting_pos = 0
      endif

      i = starting_pos+1
      j = 0
      found = .FALSE.
      count = .FALSE.
      do while((.not. found) .and. (i .le. pline%ntokens))
        if (leqi(pline%id(i), 'a')) j = j + 1
        if (j .eq. ind) then

           found = .TRUE.

           ! Parse token list
           llist = tokens(pline,i)
           ! The list does have the markers attached (remove them)
           li = len_trim(llist)-1
           llist = trim(llist(2:li))
           lpline => digest(llist)

           ! We now have converted the list into a
           ! parseable line

           ! Does the user request length?
           count = ni <= 0
           li = 0 ! counter for the number of items in the list
           ti = 1 ! the current token iterator
           is_del = .false.
           do while ( ti < lpline%ntokens )

              ! First we need to check whether we have a list delimiter next
              if (leqi(lpline%id(ti+1),'n')) then
                 sep = names(lpline,1,after=ti)
                 is_del = leqi(sep,',') ! apparently ',' is not a token
              end if

              if (leqi(lpline%id(ti+1),'n').and. .not. is_del) then
                 ! We have a range

                 if ( lpline%ntokens <= ti + 1 ) then
                   call die('PARSE module: lists', 'Missing end range', &
                        THIS_FILE, __LINE__)
                 end if
                 if ( .not.leqi(lpline%id(ti),'i') .or. &
                      .not.leqi(lpline%id(ti+2),'i') ) then
                   call die('PARSE module: lists', 'Range is not well-defined', &
                        THIS_FILE, __LINE__)
                 end if

                 ! grab the seperator
                 sep = names(lpline,1,after=ti)
                 if ( leqi(sep,'to') .or. &
                      leqi(sep,':') .or. &
                      leqi(sep,'--') .or. &
                      leqi(sep,'---') ) then

                    ! Sort the range
                    lR = integers(lpline,1,after=ti-1)
                    uR = integers(lpline,1,after=ti+1)
                    sR = 1

                    ! Figure out if we have a step in the range
                    if ( ti + 3 < lpline%ntokens ) then
                       if ( leqi(lpline%id(ti+3),'n') .and. &
                          leqi(lpline%id(ti+4),'i') ) then
                          sep = names(lpline,1,after=ti+2)
                          if ( leqi(sep,'step') ) then
                             sR = integers(lpline,1,after=ti+3)
                             ! step after the 'step <val>'
                             ti = ti + 2
                          end if
                       end if
                    end if

                    ! Correct sign of stepper
                    if ( lR <= uR ) sR = abs(sR)
                    if ( uR <  lR ) sR = -abs(sR)
                    if ( sR == 0 ) call die('PARSE module: lists', &
                         'Stepping a list cannot be stepped by 0', &
                         THIS_FILE, __LINE__ )
                    do iR = lR , uR, sR
                       call add_exit(count,li,ni,iR)
                    end do

                    ! jump across the range
                    ti = ti + 2
                 else
                   call die('PARSE module: lists', 'Unknown token in list', &
                        THIS_FILE, __LINE__)
                 end if

              elseif (leqi(lpline%id(ti),'i')) then

                 call add_exit(count,li,ni,integers(lpline,1,after=ti-1))

              end if

              ti = ti + 1
              if ( is_del ) then
                 ti = ti + 1
                 is_del = .false.
              end if

           end do

           ! Read last element (or the only element if one is given)
           if ( ti == lpline%ntokens ) then
              if (leqi(lpline%id(ti),'i')) then
                 call add_exit(count,li,ni,integers(lpline,1,after=ti-1))
              end if
           end if

           ! Clean-up parsed list-line
           call destroy(lpline)

           if ( count ) ni = li

        endif
        i = i + 1
      enddo

      if (.not. found) then
        call die('PARSE module: lists', 'Not enough lists in line', &
            THIS_FILE, __LINE__)
      end if

      RETURN
!------------------------------------------------------------- END

    contains

      subroutine add_exit(is_counting,idx,ni,val)
        logical, intent(in) :: is_counting
        integer, intent(inout) :: idx
        integer, intent(in) :: ni, val
        idx = idx + 1
        if ( is_counting ) return
        if ( idx > ni ) then
          found = .false.
        else
          list(idx) = val
        end if
      end subroutine add_exit

    END SUBROUTINE integerlists

!
!   Return a given [integer|real] token, specifying it by its sequence
!   number. It is also possible to make the sequence start after
!   a given token number in the line.
!
    FUNCTION values(pline, ind, after)
      implicit none
!------------------------------------------------- Input Variables
      integer(ip), intent(in)           :: ind
      integer(ip), intent(in), optional :: after
      type(parsed_line), pointer        :: pline

!------------------------------------------------ Output Variables
      real(dp)                          :: values

!------------------------------------------------- Local Variables
      logical                           :: found
      integer(ip)                       :: i, j, starting_pos

!----------------------------------------------------------- BEGIN
      if (PRESENT(after)) then
        if (after .lt. 0) then
          call die('PARSE module: values', 'Wrong starting position',   &
                   THIS_FILE, __LINE__,cline=characters(pline,1,-1))
        endif
        starting_pos = after
      else
        starting_pos = 0
      endif

      i = starting_pos+1
      j = 0
      found = .FALSE.
      do while((.not. found) .and. (i .le. pline%ntokens))
        if ((leqi(pline%id(i), 'i')) .or. (leqi(pline%id(i), 'r')))     &
          j = j + 1
        if (j .eq. ind) then
          values = s2r(tokens(pline,i))
          found = .TRUE.
        endif
        i = i + 1
      enddo

      if (.not. found) then
        call die('PARSE module: values', 'Not enough values in line',   &
                 THIS_FILE, __LINE__,cline=characters(pline,1,-1))
      endif

      RETURN
!------------------------------------------------------------- END
    END FUNCTION values

!
!   Return a given name token, specifying it by its sequence
!   number. It is also possible to make the sequence start after
!   a given token number in the line.
!
    FUNCTION names(pline, ind, after)
      implicit none
!------------------------------------------------- Input Variables
      integer(ip), intent(in)           :: ind
      integer(ip), intent(in), optional :: after
      type(parsed_line), pointer        :: pline

!------------------------------------------------ Output Variables
      character(len=MAX_LENGTH)         :: names

!------------------------------------------------- Local Variables
      logical                           :: found
      integer(ip)                       :: i, j, starting_pos

!----------------------------------------------------------- BEGIN
      if (PRESENT(after)) then
        if (after .lt. 0) then
          call die('PARSE module: names', 'Wrong starting position',    &
                   THIS_FILE, __LINE__,cline=characters(pline,1,-1))
        endif
        starting_pos = after
      else
        starting_pos = 0
      endif

      i = starting_pos+1
      j = 0
      found = .FALSE.
      do while((.not. found) .and. (i .le. pline%ntokens))
        if (leqi(pline%id(i), 'n')) j = j + 1
        if (j .eq. ind) then
          names = trim(tokens(pline,i))
          found = .TRUE.
        endif
        i = i + 1
      enddo

      if (.not. found) then
        call die('PARSE module: names', 'Not enough names in line', &
                 THIS_FILE, __LINE__,cline=characters(pline,1,-1))
      endif

      RETURN
!------------------------------------------------------------- END
    END FUNCTION names

!
!   Return a given block label if it is found, else returns ''
!   Syntax must be: '%block Label' (bl) as stored in fdf structure
!
    FUNCTION blocks(pline)
      implicit none
!------------------------------------------------- Input Variables
      type(parsed_line), pointer    :: pline

!------------------------------------------------ Output Variables
      character(len=MAX_LENGTH)     :: blocks

!----------------------------------------------------------- BEGIN
      if (match(pline, 'bl')) then
        blocks = tokens(pline, 2)
      else
        blocks = ' '
      endif

      RETURN
!------------------------------------------------------------- END
    END FUNCTION blocks



!     FUNCTION modules(pline)
!       implicit none
! !------------------------------------------------- Input Variables
!       type(parsed_line), pointer    :: pline

! !------------------------------------------------ Output Variables
!       character(len=MAX_LENGTH)     :: modulenames

! !----------------------------------------------------------- BEGIN
!       if (match(pline, 'bl')) then
!         modulenames = trim(pline%line)(8:)
!       else
!         modulenames = ' '
!       endif

!       RETURN
! !------------------------------------------------------------- END
!     END FUNCTION modules



!
!   Return a given endblock label if it is found, else returns ''
!   Syntax must be: '%endblock Label' (el) as stored in fdf structure
!
    FUNCTION endblocks(pline)
      implicit none
!------------------------------------------------- Input Variables
      type(parsed_line), pointer    :: pline

!------------------------------------------------ Output Variables
      character(len=MAX_LENGTH)     :: endblocks

!----------------------------------------------------------- BEGIN
      if (match(pline, 'el')) then
        endblocks = tokens(pline, 2)
      else
        endblocks = ' '
      endif

      RETURN
!------------------------------------------------------------- END
    END FUNCTION endblocks

!
!   Return a given label name if it is found, else returns ''
!   Syntax must be: 'Label Value' (li|lr|ln|l) as stored in fdf structure
!
    FUNCTION labels(pline)
      implicit none
!------------------------------------------------- Input Variables
      type(parsed_line), pointer :: pline

!------------------------------------------------ Output Variables
      character(len=MAX_LENGTH)  :: labels

!----------------------------------------------------------- BEGIN
      if (match(pline, 'l')) then
        labels = tokens(pline, 1)
      else
        labels = ' '
      endif

      RETURN
!------------------------------------------------------------- END
    END FUNCTION labels

!
!   Return a given token as character, specifying it by its sequence
!   number. It is also possible to make the sequence start after
!   a given token number in the line.
!
    FUNCTION tokens(pline, ind, after)
      implicit none
!------------------------------------------------- Input Variables
      integer(ip), intent(in)           :: ind
      integer(ip), intent(in), optional :: after
      type(parsed_line), pointer        :: pline

!------------------------------------------------ Output Variables
      character(len=MAX_LENGTH)         :: tokens

!------------------------------------------------- Local Variables
      integer(ip)                       :: starting_pos, loc

!----------------------------------------------------------- BEGIN
      if (PRESENT(after)) then
        if ((after .lt. 0) .or. (after .ge. pline%ntokens))             &
          call die('PARSE module: tokens', 'Wrong starting position',   &
                   THIS_FILE, __LINE__,cline=characters(pline,1,-1))
        starting_pos = after
      else
        starting_pos = 0
      endif

      if (starting_pos+ind .gt. pline%ntokens)                          &
        call die('PARSE module: tokens', 'Wrong starting position',     &
                 THIS_FILE, __LINE__,cline=characters(pline,1,-1))

      loc = starting_pos+ind
      tokens = pline%line(pline%first(loc):pline%last(loc))

      RETURN
!------------------------------------------------------------- END
    END FUNCTION tokens

!
!   Return a given token as character, specifying it by its sequence
!   number. It is also possible to make the sequence start after
!   a given token number in the line.
!
!     FUNCTION modulenames(pline, ind, after)
!       implicit none
! !------------------------------------------------- Input Variables
!       integer(ip), intent(in)           :: ind
!       integer(ip), intent(in), optional :: after
!       type(parsed_line), pointer        :: pline

! !------------------------------------------------ Output Variables
!       character(len=MAX_LENGTH)         :: tokens

! !------------------------------------------------- Local Variables
!       integer(ip)                       :: starting_pos, loc

! !----------------------------------------------------------- BEGIN
!       if (PRESENT(after)) then
!         if ((after .lt. 0) .or. (after .ge. pline%ntokens)) &
!           call die('PARSE module: tokens', 'Wrong starting position', &
!                    THIS_FILE, __LINE__,cline=characters(pline,1,-1))
!         starting_pos = after
!       else
!         starting_pos = 0
!       endif

!       if (starting_pos+ind .gt. pline%ntokens) &
!         call die('PARSE module: tokens', 'Wrong starting position', &
!                  THIS_FILE, __LINE__,cline=characters(pline,1,-1))

!       loc = starting_pos+ind
!       modulenames = pline%line(pline%first(loc):pline%last(loc))

!       RETURN
! !------------------------------------------------------------- END
!     END FUNCTION modulenames

!
!   Return a piece of the input line, given specifying it by the sequence
!   numbers of the initial and final tokens.
!   A negative final index means that it is counted from the end, e.g.
!   ind_final=-1 refers to the last token.
!   It is also possible to make the sequence start after
!   a given token number in the line.
!
    FUNCTION characters(pline, ind_init, ind_final, after)
      implicit none
!------------------------------------------------- Input Variables
      integer(ip), intent(in)           :: ind_init
      integer(ip), intent(in)           :: ind_final
      integer(ip), intent(in), optional :: after
      type(parsed_line), pointer        :: pline

!------------------------------------------------ Output Variables
      character(len=MAX_LENGTH)         :: characters

!------------------------------------------------- Local Variables
      integer(ip)                       :: starting_pos
      integer(ip)                       :: loc_init, loc_final

!----------------------------------------------------------- BEGIN
      if (PRESENT(after)) then
        if ((after .lt. 0) .or. (after .ge. pline%ntokens))             &
          call die('PARSE module: tokens', 'Wrong starting position',   &
                   THIS_FILE, __LINE__)
        starting_pos = after
      else
        starting_pos = 0
      endif

      loc_init = starting_pos+ind_init
      if (ind_final < 0 ) then
         loc_final = pline%ntokens + ind_final + 1
      else
         loc_final = starting_pos+ind_final
      endif

      if (      (loc_init .lt. 0)                   &
           .OR. (loc_init .gt. pline%ntokens)       &
           .OR. (loc_final .lt. 0)                  &
           .OR. (loc_final .gt. pline%ntokens)      &
           .OR. (loc_final .lt. loc_init)           &
         )  then
         call die('PARSE module: characters', 'Wrong limits',     &
                 THIS_FILE, __LINE__)
      endif

      characters = pline%line(pline%first(loc_init):pline%last(loc_final))

      RETURN
!------------------------------------------------------------- END
    END FUNCTION characters

!
!   Main processing function. Digest a character line array
!   building a digested parsed_line structure.
!
    FUNCTION digest(line) result(pline)
      implicit none
!------------------------------------------------- Input Variables
      character(len=*), intent(in) :: line

!------------------------------------------------ Output Variables
      type(parsed_line), pointer   :: pline

!------------------------------------------------- Local Variables
      character                    :: token_id(MAX_NTOKENS)
      integer(ip)                  :: i, ntokens
      integer(ip)                  :: first(MAX_NTOKENS), last(MAX_NTOKENS)

!----------------------------------------------------------- BEGIN
!     Parse line, and get morphology
      call parses(ntokens, line, first, last)
      call morphol(ntokens, line, first, last, token_id)

!     Build parsed_line structure
      NULLIFY(pline)
      call create(pline)
      pline%ntokens = ntokens

      if (ntokens .gt. MAX_NTOKENS) then
         call die('PARSE module: digest', 'Too many tokens', &
                   THIS_FILE, __LINE__, rc=666)
      endif

      pline%line        = line
      do i= 1, ntokens
         pline%first(i) = first(i)
         pline%last(i)  = last(i)
         pline%id(i)    = token_id(i)
      enddo

      RETURN
!------------------------------------------------------------- END
    END FUNCTION digest

!
!   Parses a character line, filling ntokens (# of tokens)
!   first and last (beginning and ending of each one token)
!
    SUBROUTINE parses(ntokens, line, first, last)
      implicit none
!------------------------------------------------- Input Variables
      character(len=*)             :: line

!------------------------------------------------ Output Variables
      integer(ip)                  :: ntokens
      integer(ip)                  :: first(MAX_NTOKENS), last(MAX_NTOKENS)

!------------------------------------------------- Local Variables
      logical                      :: intoken, instring, completed
      logical                      :: inlist
      integer(ip)                  :: i, c, stringdel, length

!     Character statement functions
      logical :: is_digit, is_upper, is_lower, is_alpha,                &
                 is_alnum, is_extra, is_tokch
      logical :: is_comment, is_delstr, is_dellist, is_special

      is_digit(i) = (i .ge. 48) .and. (i .le. 57)
      is_upper(i) = (i .ge. 65) .and. (i .le. 90)
      is_lower(i) = (i .ge. 97) .and. (i .le. 122)
      is_alpha(i) = is_upper(i) .or. is_lower(i)
      is_alnum(i) = is_digit(i) .or. is_alpha(i)

!     Extra characters allowed in tokens:  $ % * + & - . / @ ^ _ | ~
      is_extra(i) = ((i .ge. 36) .and. (i .le. 38))                     &
                     .or. (i .eq. 42) .or. (i .eq. 43) .or. (i .eq. 45) &
                     .or. (i .eq. 46) .or. (i .eq. 47) .or. (i .eq. 64) &
                     .or. (i .eq. 94) .or. (i .eq. 95) .or. (i .eq. 124)&
                     .or. (i .eq. 126) .or. (i .eq. 58)

      is_tokch(i) = is_alnum(i) .or. is_extra(i)

!     Comments are signaled by:  !  #  ;
      is_comment(i) = (i .eq. 33) .or. (i .eq. 35) .or. (i .eq. 59)

!     String delimiters: "  '  `
      is_delstr(i)  = (i .eq. 34) .or. (i .eq. 39) .or. (i .eq. 96)

!     List delimiters: [ ]
      is_dellist(i)  = (i .eq. 91) .or. (i .eq. 93)

!     Dictionary delimiters: { }
!      is_deldict(i)  = (i .eq. 123) .or. (i .eq. 125)

!     Special characters which are tokens by themselves: <
      is_special(i) = (i .eq. 60)

!----------------------------------------------------------- BEGIN
      ntokens = 0

      intoken  = .FALSE.
      instring = .FALSE.
      inlist   = .FALSE.
      stringdel = 0

      ! Trim space at the end (not from the left)
      length = len_trim(line)

      i = 1
      completed = .FALSE.
      do while((i <= length) .and. (.not. completed))
        c = ichar(line(i:i))

!       Possible comment...
        if (is_comment(c)) then
          if (instring.or.inlist) then
            last(ntokens) = i
          else
            completed = .TRUE.
          endif

!       Character allowed in a token...
        elseif (is_tokch(c)) then
          if (.not. intoken) then
            intoken = .TRUE.
            ntokens = ntokens + 1
            first(ntokens) = i
          endif
          last(ntokens) = i

!       Character that forms a token by itself...
        elseif (is_special(c)) then
          if (.not. instring .and. .not. inlist) then
            ntokens = ntokens + 1
            first(ntokens) = i
            intoken = .FALSE.
          endif
          last(ntokens) = i

!      List delimiter... We only allow single lists, not nested lists
       elseif (is_dellist(c)) then
          if (.not. instring .and. .not. inlist) then
             inlist  = .TRUE.
             intoken = .TRUE.
             ntokens = ntokens + 1
             first(ntokens) = i
          elseif (inlist) then
             ! end list (skip last token)
             intoken = .FALSE.
             inlist  = .FALSE.
             last(ntokens) = i
          else
             last(ntokens) = i
          end if

!       String delimiter... make sure it is the right one before closing.
!       If we are currently in a token, the delimiter is appended to it.
        elseif (is_delstr(c)) then
          if (instring) then
            if (c .eq. stringdel) then
              instring = .FALSE.
              intoken  = .FALSE.
              stringdel = 0
            else
              last(ntokens) = i
            endif
          else
            if (intoken) then
              last(ntokens) = i
            else
              instring = .TRUE.
              intoken  = .TRUE.
              stringdel = c
              ntokens = ntokens + 1
              first(ntokens) = i + 1
              last(ntokens)  = i + 1
            endif
          endif

!       Token delimiter...
        else
          if (instring.or.inlist) then
            last(ntokens) = i
          else
            if (intoken) intoken = .FALSE.
          endif
        endif

        i = i + 1

        ! Check whether the parsing is correctly handled
        if ( i > MAX_LENGTH ) then
          ! Because we will limit search to the len_trim length,
          ! then this should only be found when the line has "content" too long.
          ! Note that this will *never* be executed if a comment is too
          ! long because it is checked as the first requirement and then
          ! completes parsing the line.
          call die('PARSE module: parses', 'Too long line (132 char): ' // &
              trim(line), THIS_FILE, __LINE__)
        end if

      enddo

      if (parse_debug) then
        write(parse_log,*) 'PARSER:', ntokens, 'token(s)'
        do i= 1, ntokens
          write(parse_log,*) '   Token:', '|',line(first(i):last(i)),'|'
        enddo
        write(parse_log,*) ' '
      endif
!------------------------------------------------------------- END
    END SUBROUTINE parses

!
!   Classifies the tokens according to their morphology
!
    SUBROUTINE morphol(ntokens, line, first, last, token_id)
      implicit none
!------------------------------------------------- Input Variables
      character(len=*)          :: line
      integer(ip)               :: ntokens
      integer(ip)               :: first(MAX_NTOKENS), last(MAX_NTOKENS)

!------------------------------------------------ Output Variables
      character                 :: token_id(MAX_NTOKENS)

!------------------------------------------------- Local Variables
      character(len=MAX_LENGTH) :: token, msg
      integer(ip)               :: i, j, ierr
      real(dp)                  :: real_value

!----------------------------------------------------------- BEGIN
      do i= 1, ntokens
        token = line(first(i):last(i))
        j = last(i) - first(i) + 1
        if ( ichar(token(1:1)) .eq. 91 .and. &
            ichar(token(j:j)) .eq. 93 ) then
          ! if the token starts with [ and ends with ], it will be a list
          ! We do a simple check for the list type.
          ! Since we are only dealing with integer/real lists
          ! we can simply check for a . (comma separation) which
          ! will enable an easy distinguishment between integers and reals.
          if ( index(token(1:j), '.') > 0 ) then
            token_id(i) = 'c'
          else
            token_id(i) = 'a'
          end if

!        else if ( ichar(token(1:1)) .eq. 123 .and. &
!             ichar(token(j:j)) .eq. 125 ) then
! if the token starts with { and ends with }, it will be a dictionary
!           token_id(i) = 'd'
        elseif (is_value(token)) then

!         This read also serves to double check the token for
!         real meaning (for example, ".d0" should give an error)
          read(token, fmt=*, iostat=ierr) real_value
          if (ierr .ne. 0) then
            write(msg,'(a,i3,1x,a,/,a)') 'Error in numeric conversion ', &
                 'at token number', i, ' in line ''', TRIM(line), ''''
            call die('PARSE module: morphol', msg, THIS_FILE, __LINE__)
          endif

          if (is_integer(token)) then
            token_id(i) = 'i'
          else
            token_id(i) = 'r'
          endif
        else
          token_id(i) = 'n'
        endif
      enddo

      if (parse_debug) then
        write(parse_log,*) 'MORPHOL:', ntokens, 'token(s)'
        do i= 1, ntokens
          write(parse_log,*) '   Token:', '|', token_id(i), '|'
        enddo
        write(parse_log,*) ' '
      endif
!------------------------------------------------------------- END
    END SUBROUTINE morphol

!
!   Set the morphology of a specific token in a parsed line
!
    SUBROUTINE setmorphol(ntoken, token_id, pline)
      implicit none
!------------------------------------------------- Input Variables
      character                  :: token_id
      integer(ip)                :: ntoken

!------------------------------------------------ Output Variables
      type(parsed_line), pointer :: pline

!------------------------------------------------- Local Variables
      character(len=MAX_LENGTH)  :: msg

!----------------------------------------------------------- BEGIN

!     Check if token_id is a valid morphology id
!     'a' -> List (integers)
!     'c' -> List (reals)
!     'l' -> Label
!     'b' -> BeginBlock
!     'e' -> EndBlock
!     'i' -> Integer
!     'r' -> Real
!     'n' -> Name
      if ((token_id .ne. 'a') .and. (token_id .ne. 'c') .and. &
          (token_id .ne. 'l') .and. (token_id .ne. 'b') .and. &
          (token_id .ne. 'e') .and. (token_id .ne. 'i') .and. &
          (token_id .ne. 'r') .and. (token_id .ne. 'n')) then
        write(msg,*) 'Morphology id = ''', token_id,                    &
                     ''' not valid for token = ''', tokens(pline, ntoken), ''''
        call die('PARSE module: setmorphol', msg, THIS_FILE, __LINE__)
      endif

      pline%id(ntoken) = token_id
!------------------------------------------------------------- END
    END SUBROUTINE setmorphol

!
!   Search a string along a parsed line tokens. If found, it returns
!   the index in the list of the token that matches with the string.
!   Otherwise it returns -1.
!
    FUNCTION search_fun(string, pline_fun, after, eq_func)
      implicit none
!------------------------------------------------- Input Variables
      character(len=*)                  :: string
      integer(ip), intent(in), optional :: after
      type(parsed_line), pointer        :: pline_fun
      optional                          :: eq_func

      interface
        function eq_func(s1, s2)
          logical                       :: eq_func
          character(len=*), intent(in)  :: s1, s2
        end function eq_func
      end interface

!------------------------------------------------ Output Variables
      integer(ip)                       :: search_fun

!------------------------------------------------- Local Variables
      integer(ip)                       :: i, starting_pos

!----------------------------------------------------------- BEGIN
      if (PRESENT(after)) then
        if (after .lt. 0) then
          call die('PARSE module: search_fun', 'Wrong starting position', &
                   THIS_FILE, __LINE__)
        endif
        starting_pos = after
      else
        starting_pos = 0
      endif

      search_fun = -1
      if (.not. ASSOCIATED(pline_fun)) then
        call die('PARSE module: search_fun', 'parsed_line not associated', &
                 THIS_FILE, __LINE__)
      endif

!     The default comparison routine is 'leqi' (case-insensitive)
      if (PRESENT(eq_func)) then
        i = starting_pos+1
        do while((search_fun .eq. -1) .and. (i .le. pline_fun%ntokens))
          if (eq_func(string, tokens(pline_fun, i))) search_fun = i
          i = i + 1
        enddo
      else
        i = starting_pos+1
        do while((search_fun .eq. -1) .and. (i .le. pline_fun%ntokens))
          if (leqi(string, tokens(pline_fun, i))) search_fun = i
          i = i + 1
        enddo
      endif

      RETURN
!------------------------------------------------------------- END
    END FUNCTION search_fun

!
!   Search a string along a parsed line tokens. If found, leaves
!   in 'ind' the index token in the list that matches with the string
!   and it returns .TRUE. Otherwise it returns .FALSE. and -1 in 'ind'
!
    FUNCTION search_sub(pline_sub, string, ind, after, eq_func)
      implicit none
!------------------------------------------------- Input Variables
      character(len=*)                   :: string
      integer(ip), intent(in), optional  :: after
      type(parsed_line), pointer         :: pline_sub
      optional                           :: eq_func

      interface
        function eq_func(s1, s2)
          logical                        :: eq_func
          character(len=*), intent(in)   :: s1, s2
        end function eq_func
      end interface

!------------------------------------------------ Output Variables
      logical                            :: search_sub
      integer(ip), intent(out), optional :: ind

!------------------------------------------------- Local Variables
      integer(ip)                        :: i, starting_pos

!----------------------------------------------------------- BEGIN
      if (PRESENT(after)) then
        if (after .lt. 0) then
          call die('PARSE module: search_sub', 'Wrong starting position', &
                   THIS_FILE, __LINE__)
        endif
        starting_pos = after
      else
        starting_pos = 0
      endif

      if (PRESENT(ind)) ind = -1
      search_sub = .FALSE.
      if (.not. ASSOCIATED(pline_sub)) then
        call die('PARSE module: search_sub', 'parsed_line not associated', &
                 THIS_FILE, __LINE__)
      endif

!     The default comparison routine is 'leqi' (case-insensitive)
      if (PRESENT(eq_func)) then
        i = starting_pos+1
        do while((.not. search_sub) .and. (i .le. pline_sub%ntokens))
          if (eq_func(string, tokens(pline_sub, i))) then
            if (PRESENT(ind)) ind = i
            search_sub = .TRUE.
          endif
          i = i + 1
        enddo
      else
        i = starting_pos+1
        do while((.not. search_sub) .and. (i .le. pline_sub%ntokens))
          if (leqi(string, tokens(pline_sub, i))) then
            if (PRESENT(ind)) ind = i
            search_sub = .TRUE.
          endif
          i = i + 1
        enddo
      endif

      RETURN
!------------------------------------------------------------- END
    END FUNCTION search_sub
!
!   Search a sub-string along a parsed line tokens. If found, leaves
!   in 'ind' (if present) the index token in the list that has the
!   string as a substring and it returns .TRUE. Otherwise it returns
!   .FALSE. and -1 in 'ind'
!
    FUNCTION substring_search(pline_sub, string, ind, after)
      implicit none
!------------------------------------------------- Input Variables
      character(len=*)                   :: string
      integer(ip), intent(in), optional  :: after
      type(parsed_line), pointer         :: pline_sub

!------------------------------------------------ Output Variables
      logical                            :: substring_search
      integer(ip), intent(out), optional :: ind

!------------------------------------------------- Local Variables
      integer(ip)                        :: i, starting_pos

!----------------------------------------------------------- BEGIN
      if (PRESENT(after)) then
        if (after .lt. 0) then
          call die('PARSE module: substring_search', &
                   'Wrong starting position', &
                   THIS_FILE, __LINE__)
        endif
        starting_pos = after
      else
        starting_pos = 0
      endif

      if (PRESENT(ind)) ind = -1
      substring_search = .FALSE.
      if (.not. ASSOCIATED(pline_sub)) then
        call die('PARSE module: substring_search', &
                 'parsed_line not associated', &
                 THIS_FILE, __LINE__)
      endif

!     NOTE that the we use the case-sensitive Fortran 'index' function
      i = starting_pos+1
      do while((.not. substring_search) .and. (i .le. pline_sub%ntokens))
         if (index(tokens(pline_sub, i),string) > 0) then
            if (PRESENT(ind)) ind = i
            substring_search = .TRUE.
         endif
         i = i + 1
      enddo

      RETURN
!------------------------------------------------------------- END
    END FUNCTION substring_search

!
!   Checks whether the morphology of the line or part of it
!   matches the 'signature' string str.
!   If 'after' is present, try to match the 'signature' after
!   that number of tokens.
!
    FUNCTION match(pline, str, after)
      implicit none
!------------------------------------------------- Input Variables
      character(*), intent(in)          :: str
      integer(ip), intent(in), optional :: after
      type(parsed_line), pointer        :: pline

!------------------------------------------------ Output Variables
      logical                           :: match

!------------------------------------------------- Local Variables
      character                         :: c, id
      integer(ip)                       :: i, nids, shift

!----------------------------------------------------------- BEGIN
      if (PRESENT(after)) then
        if (after .lt. 0) then
          call die('PARSE module: match', 'Wrong starting position',    &
                   THIS_FILE, __LINE__,cline=characters(pline,1,-1))
        endif
        shift = after
      else
        shift = 0
      endif

      nids = LEN_TRIM(str)
      if (pline%ntokens - shift .lt. nids) then
        match = .FALSE.
      else
        i = 1
        match = .TRUE.
        do while (match .and. (i .le. nids))
          c  = str(i:i)
          id = pline%id(shift+i)

          if (.not. leqi(c,id)) then

          !  x: matches anything
            if (leqi(c,'x')) then
               ! do nothing -- match stays .true.

          !  v: integer or real
            else if (leqi(c,'v')) then
               if (.not.(leqi(id,'i') .or. leqi(id,'r'))) then
                  match = .false.
               endif

          !  s: integer or real or name (symbol)
            else if (leqi(c,'s')) then
               if (.not.(leqi(id,'i') .or.   &
                         leqi(id,'r') .or.   &
                         leqi(id,'n')  )) then
                  match = .false.
               endif

          !  a: array (integer list)
            else if (leqi(c,'a')) then
               if (.not.(leqi(id,'a')) ) then
                  match = .false.
               endif

          !  c: array (real list)
            else if (leqi(c,'c')) then
               if (.not.(leqi(id,'c')) ) then
                  match = .false.
               endif

          !  e: array (list)
            else if (leqi(c,'e')) then
               if (.not.(leqi(id,'a') .or. &
                   leqi(id,'c') )) then
                  match = .false.
               endif

          !  j: integer or name (integer-symbol)
            else if (leqi(c,'j')) then
               if (.not.(leqi(id,'i') .or. leqi(id,'n'))) then
                  match = .false.
               endif

          !  cannot find a match
            else
               match = .false.
            endif

          endif

          i = i + 1
        enddo
      endif

      RETURN
!------------------------------------------------------------- END
    END FUNCTION match

!
!   Checks if the string has a valid integer format
!
    FUNCTION is_integer(string)
      implicit none
!------------------------------------------------- Input Variables
      character(len=*) :: string

!------------------------------------------------ Output Variables
      logical          :: is_integer

!------------------------------------------------- Local Variables
      character        :: c
      integer(ip)      :: i, length

      logical          :: is_digit, is_sign

      is_digit(c) = ((ichar(c) .ge. 48) .and. (ichar(c) .le. 57))
      is_sign(c)  = ((c .eq. '+') .or. (c .eq. '-'))

!----------------------------------------------------------- BEGIN
      length = LEN_TRIM(string)
      if (length .gt. 0) then
        c = string(1:1)
        if ((is_digit(c)) .or. (is_sign(c))) then
          i = 2
          is_integer = .TRUE.
          do while (is_integer .and. (i .le. length))
            c = string(i:i)
            if (.not. (is_digit(c))) then
              is_integer = .FALSE.
            endif
            i = i + 1
          enddo
        else
          is_integer = .FALSE.
        endif
      else
        is_integer = .FALSE.
      endif

      RETURN
!------------------------------------------------------------- END
    END FUNCTION is_integer

!
!   Checks if the string has a valid value format [real|integer]
!
    FUNCTION is_value(string)
      implicit none
!------------------------------------------------- Input Variables
      character(len=*) :: string

!------------------------------------------------ Output Variables
      logical          :: is_value

!------------------------------------------------- Local Variables
      character        :: c
      logical          :: dotsok
      integer(ip)      :: i, length, exp_mark

      logical          :: is_digit, is_sign, is_dot, is_expmark

      is_digit(c)   = ((ichar(c) .ge. 48) .and. (ichar(c) .le. 57))
      is_sign(c)    = ((c .eq. '+') .or. (c .eq. '-'))
      is_dot(c)     = ((c .eq. '.') .and. dotsok)
      is_expmark(c) = ((c .eq. 'e') .or. (c .eq. 'E') .or.              &
                       (c .eq. 'd') .or. (c .eq. 'D'))

!----------------------------------------------------------- BEGIN
      length = LEN_TRIM(string)

      is_value = .FALSE.
      dotsok   = .TRUE.

!     Find the starting point of a possible exponent
      exp_mark = length+1
      do i= 1, length
        c = string(i:i)
        if (is_expmark(c)) exp_mark = i
      enddo
      if (exp_mark .eq. length) return    ! Form: XXXXXd

      c = string(1:1)
      if (.not. (is_digit(c) .or. is_sign(c))) then
        if (is_dot(c)) then
          dotsok = .FALSE.
        else
          return
        endif
      endif

      do i= 2, exp_mark-1
        c = string(i:i)
        if (.not. (is_digit(c))) then
          if (is_dot(c)) then
            dotsok = .FALSE.
          else
            return
          endif
        endif
      enddo

!     Is the exponent an integer?
      if (exp_mark .lt. length) then
        if (.not. is_integer(string(exp_mark+1:length))) return
      endif

!     Here we could do some extra checks to see if the string still makes
!     sense... For example, "." and ".d0" pass the above tests but are not
!     readable as numbers. I believe this should be reported by the
!     conversion routine, to warn the user of a mis-typed number, instead
!     of reporting it as a string and break havoc somewhere else.

!     This cases should not be accepted, since
!     now we are scanning the whole input file blindly

      if (length == 1) then
         if (is_sign(string(1:1))) return  ! Remove '+' and '-'
         if (string(1:1) == "."  ) return  ! Remove '.'
      endif

      is_value = .TRUE.

      RETURN
!------------------------------------------------------------- END
    END FUNCTION is_value

!
!   Set debugging level for parses/morphol routines
!
    SUBROUTINE setdebug(level)
      implicit none
!------------------------------------------------- Input Variables
      integer(ip) :: level

!----------------------------------------------------------- BEGIN
      parse_debug = (level .eq. 1)
      RETURN
!------------------------------------------------------------- END
    END SUBROUTINE setdebug

!
!   Set log unit for parses/morphol routines
!
    SUBROUTINE setlog(unit)
      implicit none
!------------------------------------------------- Input Variables
      integer(ip) :: unit

!----------------------------------------------------------- BEGIN
      parse_log = unit
      RETURN
!------------------------------------------------------------- END
    END SUBROUTINE setlog

!
    subroutine serialize_pline(pline,string,length)
    type(parsed_line)   :: pline
    character(len=*), intent(out) :: string
    integer, intent(out) :: length

    integer :: pos, i
    character(len=10) buffer

    length = SERIALIZED_LENGTH
    if (len(string) < length) then
       call die('PARSE module: serialize_pline', &
            "String too short", &
            THIS_FILE, __LINE__)
    endif

    string = ""
    string(1:MAX_LENGTH) = pline%line
    pos = MAX_LENGTH
    write(string(pos+1:pos+4),"(i4)") pline%ntokens
    pos = pos + 4

    do i=1,pline%ntokens
       write(buffer,"(1x,a1,2i4)") pline%id(i), pline%first(i), pline%last(i)
       string(pos+1:pos+10) = buffer
       pos = pos + 10
    enddo

  end subroutine serialize_pline

    subroutine recreate_pline(pline,string)
    type(parsed_line), pointer   :: pline
    character(len=*), intent(in) :: string

    integer :: pos, i

    if (len(string) < SERIALIZED_LENGTH)  then
       call die('PARSE module: recreate_pline', &
            "String too short", &
            THIS_FILE, __LINE__,cline=characters(pline,1,-1))
    endif

    pline%line = string(1:MAX_LENGTH)
    pos = MAX_LENGTH
    read(string(pos+1:pos+4),"(i4)") pline%ntokens
    pos = pos + 4
    do i=1,pline%ntokens
       read(string(pos+1:pos+10),"(1x,a1,2i4)") pline%id(i), pline%first(i), pline%last(i)
       pos = pos + 10
    enddo

  end subroutine recreate_pline

END MODULE parse
