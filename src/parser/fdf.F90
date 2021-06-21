#if defined HAVE_CONFIG_H
#  include "config.h"
#endif

#define THIS_FILE "fdf.F90"
!>
!!=====================================================================
!!
!! This file is part of the FDF package.
!!
!! This module implements an extended Fortran 90/95 interface
!! to the Flexible Data Format library of A. Garcia and J.M. Soler,
!! originally written in Fortran 77.
!!
!! FEATURES:
!!
!! a) Block pointers.
!!
!! Block content can be flexibly handled by means of a pointer
!! to a derived type 'block_fdf'. Typical usage:
!!
!!     use fdf
!!     type(block_fdf)            :: bfdf
!!     type(parsed_line), pointer :: pline
!!
!!     if (fdf_block('SomeBlock', bfdf)) then
!!       do while(fdf_bline(bfdf, pline))
!!         (process line 'integers|reals|values|names ...')
!!       enddo
!!       call fdf_bclose(bfdf)
!!     endif
!!
!! The subroutine 'fdf_block' returns in 'bfdf' a structure used
!! to read the contents of the block.
!!
!! Routine fdf_bline returns in 'pline' the next non-blank parsed
!! line, non-comment line from the block, unless there are no more
!! lines, in which case it returns .FALSE. and 'pline' is undefined.
!!
!! Routine fdf_bclose runs the remaining lines in the block and ensures
!! the log may be used as input in subsequent entries.
!!
!! Routine 'backspace' moves the internal pointer of 'block_fdf'
!! structure to the previous line returned.
!!
!! Routine 'rewind' moves the internal pointer of 'block_fdf' structure
!! to the beginning of the block.
!!
!! b) Generic interface to scalar routines.
!!
!! The generic function 'fdf_get' can be used instead of any of the
!! scalar routines. The specific names are also accepted.
!!
!! c) Architecture support: this FDF implementation supports the following
!!    architectures environments.
!!
!!    1) Thread-safe: The new implementation is thread-safe and will support
!!       calling it from several OMP-threads executing in the same node.
!!
!!       The implementation is as follows: fdf_init and fdf_shutdown are
!!       SINGLE/CRITICAL sections that only one thread must execute.
!!       On the other hand 'get'/'test' routines in FDF library are
!!       thread-safe because each thread keeps its relative information
!!       about the search/query that the caller program requests.
!!
!!    2) MPI-aware: For MPI executions, FDF library renames output/debugging
!!       or log files to prevent overlaps of these files in a
!!       shared/parallel filesystem.
!!
!!       This option is enabled with _MPI_ macro, and superseded by
!!       CLUSTER and BLOCKING macros.
!!
!!    3) Cluster filesystem: It is able to read the input FDF file from a
!!       non-shared filesystem and broadcast the information to the rest
!!       of the nodes in the execution.
!!
!!       The implementation is as follows: the first node in the MPI rank
!!       that is the owner of the input FDF file, process the file and
!!       send it to the rest of the nodes in the MPI communicator.
!!
!!       This option is enabled with CLUSTER macro. This option cannot be
!!       used if BLOCKING macro is enabled.
!!
!!    4) Blocking reading: For huge executions (> 1.000 nodes) this option
!!       is useful for shared/parallel filesystems where reading the same
!!       file by several nodes could be a problem (collapsing system).
!!
!!       The implementation is as follows: reading phase is done in blocking
!!       pattern of size BLOCKSIZE nodes (configurable at compile time).
!!       This means that the number of steps needed is STEPS = #NODES/BLOCKSIZE.
!!
!!       This option is enabled with BLOCKING macro. This option cannot be
!!       used if CLUSTER macro is enabled.
!!
!! @authors Alberto Garcia, 1996-2007
!! @authors Raul de la Cruz (BSC), September 2007
!! @remarks Modifications to the original libfdf made by:
!! @authors Ravindra Shinde (r.l.shinde@utwente.nl)
!! @date (2021)
!========================================================================

!! AG: If running under MPI (flagged by the MPI symbol, choose the CLUSTER
!!     mode of operation.

#ifdef MPI
# define _MPI_
# define CLUSTER
# undef MPI
#endif

#ifndef _MPI_
# if defined(CLUSTER) || defined(BLOCKING)
#   define _MPI_
# endif
#endif

MODULE fdf
  USE io_fdf

  USE parse, only: parsed_line
  USE parse, only: nintegers, nreals
  USE parse, only: nvalues, nnames, ntokens
  USE parse, only: integers, reals
  USE parse, only: values, names, tokens, characters
  USE parse, only: match
  USE parse, only: digest, blocks, endblocks, labels
  USE parse, only: destroy, setdebug, setlog, setmorphol
  USE parse, only: nlists, nintegerlists, nreallists
  USE parse, only: integerlists, reallists, valuelists

  USE parse, only: search
  USE parse, only: fdf_bsearch => search
  USE parse, only: fdf_substring_search => substring_search

  USE parse, only: serialize_pline, recreate_pline
  USE parse, only: SERIALIZED_LENGTH

  USE utils
  USE prec
  implicit none

! User callable routines in FDF library

! Start, stop and check parallelism in FDF system
  public :: fdf_init, fdf_shutdown, fdf_parallel

! Reading label functions
  public :: fdf_get
  public :: fdf_integer, fdf_single, fdf_double
  public :: fdf_string, fdf_boolean
  public :: fdf_physical, fdf_convfac
  public :: fdf_load_filename         ! Ravindra

  ! Lists
  public :: fdf_islist, fdf_islinteger, fdf_islreal
  public :: fdf_list, fdf_linteger, fdf_ldouble

! Returns the string associated with a mark line
  public :: fdf_getline

! Test if label is defined
  public :: fdf_defined, fdf_isphysical, fdf_isblock
  public :: fdf_load_defined

! Allow to overwrite things in the FDF
  public :: fdf_overwrite, fdf_removelabel, fdf_addline

! Test if a label is used in obsolete or a deprecated state
  public :: fdf_deprecated, fdf_obsolete

! %block reading (processing each line)
  public :: fdf_block, fdf_block_linecount
  public :: fdf_bline, fdf_bbackspace, fdf_brewind, fdf_bclose
  public :: fdf_bnintegers, fdf_bnreals, fdf_bnvalues, fdf_bnnames, fdf_bntokens
  public :: fdf_bintegers, fdf_breals, fdf_bvalues, fdf_bnames, fdf_btokens
  public :: fdf_bboolean, fdf_bphysical
  public :: fdf_bnlists, fdf_bnilists, fdf_bnrlists, fdf_bnvlists
  public :: fdf_bilists, fdf_brlists, fdf_bvlists

! Match, search over blocks, and destroy block structure
  public :: fdf_bmatch, fdf_bsearch, fdf_substring_search

  public :: fdf_setoutput, fdf_setdebug

! Private functions, non-callable

! Main functions to build FDF structure (called in fdf_init)
  private :: fdf_initdata, fdf_addtoken, fdf_readline
  private :: fdf_read, fdf_readlabel, fdf_searchlabel
  private :: fdf_read_xyz
  private :: fdf_open, fdf_close

! Input/Output configuration
  private :: fdf_input
!  private :: fdf_set_output_file

! Destroy dynamic list of FDF structure (called in fdf_shutdown)
  private :: fdf_destroy, fdf_destroy_dl

! MPI init/finalize functions
#ifdef _MPI_
  private :: fdf_mpi_init, fdf_mpi_finalize
#endif

! Reading functions for CLUSTER and BLOCKING configuration
#ifdef CLUSTER
  private :: setup_fdf_cluster, broadcast_fdf_struct
#endif
#ifdef BLOCKING
  private :: fdf_readblocking
#endif

! Debugging functions, level and prints debugging info
  public :: fdf_printfdf

! Finds a label in the FDF herarchy
  private :: fdf_locate
  private :: fdf_load_locate

! Dump function (for blocks)
  private :: fdf_dump

! Wrappers functions for block access, search, matching,
! number and elements in the block (call to parse module)
  interface fdf_bnintegers
    module procedure nintegers
  end interface

  interface fdf_bnlists
    module procedure nlists
  end interface

  interface fdf_bnilists
    module procedure nintegerlists
  end interface

  interface fdf_bnrlists
    module procedure nreallists
  end interface

  interface fdf_bnvlists
    module procedure nlists
  end interface

  interface fdf_bnreals
    module procedure nreals
  end interface

  interface fdf_bnvalues
    module procedure nvalues
  end interface

  interface fdf_bnnames
    module procedure nnames
  end interface

  interface fdf_bntokens
    module procedure ntokens
  end interface

  interface fdf_bintegers
    module procedure integers
  end interface

  interface fdf_bilists
    module procedure integerlists
  end interface

  interface fdf_brlists
    module procedure reallists
  end interface

  interface fdf_bvlists
    module procedure valuelists
  end interface

  interface fdf_breals
    module procedure reals
  end interface

  interface fdf_bvalues
    module procedure values
  end interface

  interface fdf_bnames
    module procedure names
  end interface

  interface fdf_btokens
    module procedure tokens
  end interface

  interface fdf_bmatch
    module procedure match
  end interface

! fdf_get wrapper for label functions
  interface fdf_get
    module procedure fdf_integer
    module procedure fdf_single
    module procedure fdf_double
    module procedure fdf_boolean
    module procedure fdf_string
    module procedure fdf_physical
  end interface

  ! fdf_list wrapper for integer/real list functions
  interface fdf_list
    module procedure fdf_linteger
    module procedure fdf_ldouble
  end interface


! Unit numbers for input, output, error notification, and
! debugging output (the latter active if fdf_debug is true)
  logical, private                :: fdf_debug   = .FALSE.,             &
                                     fdf_debug2  = .FALSE.,             &
                                     fdf_started = .FALSE.,             &
                                     fdf_output  = .FALSE.

  integer(ip), parameter, private :: maxdepth   = 7
  integer(ip), parameter, private :: maxFileNameLength = 300
  integer(ip), private            :: ndepth
  integer(ip), private            :: fdf_in(maxdepth)
  integer(ip), private            :: fdf_out, fdf_err, fdf_log

! MPI variables (id, number of MPI tasks)
#ifdef _MPI_
  logical, private                :: mpiflag
  integer(ip), private            :: rank, ntasks
#endif

! Structure for searching inside fdf blocks
  type, public :: block_fdf
    character(len=MAX_LENGTH) :: label
    type(line_dlist), pointer :: mark => null()
    character(len=MAX_LENGTH) :: modulename = ""
  end type block_fdf

! Dynamic list for parsed_line structures
  type, public :: line_dlist
    character(len=MAX_LENGTH)  :: str
    type(parsed_line), pointer :: pline => null()
    type(line_dlist), pointer  :: next => null()
    type(line_dlist), pointer  :: prev => null()
  end type line_dlist

! FDF data structure (first and last lines)
  type, private :: fdf_file
    integer(ip)               :: nlines
    type(line_dlist), pointer :: first => null()
    type(line_dlist), pointer :: last => null()
  end type fdf_file

! Input FDF file
!  type(fdf_file), private :: file_in
! original
  type(fdf_file), pointer, private :: file_in => null()

! Export the following to enable serialization by clients of the library

  public :: serialize_fdf_struct
  public :: recreate_fdf_struct
  public :: fdf_set_started

! Define by default all the others inherit module entities as privated
! avoiding redefinitions of entities in several module files with same name
  public :: parsed_line   ! Structure for searching inside fdf blocks
  public :: leqi          ! For legacy support (old codes)
  private


CONTAINS

!
!   Initialization for fdf.
!
      SUBROUTINE fdf_init( fileInput, fileOutput, unitInput )
      implicit none
!------------------------------------------------------------- Input Variables
      character(len=*),optional,intent(in):: fileInput, fileOutput
      integer,         optional,intent(in):: unitInput

#ifndef FDF_DEBUG
!------------------------------------------------------------- Local Variables
      integer(ip)  :: debug_level, output_level
#endif
      character(len=256) :: filedebug
      character(len=maxFileNameLength):: filein, fileout

!----------------------------------------------------------------------- BEGIN
!$OMP SINGLE
      ! Prevent the user from opening two head files
      if (fdf_started) then
        call die('FDF module: fdf_init', 'Head file already set',       &
                 THIS_FILE, __LINE__, fdf_err)
      endif

#ifdef _MPI_
      call fdf_mpi_init()
#endif
      call fdf_initdata()

      call io_geterr(fdf_err)

      ! Set in/out file names, if fileInput and fileOutput are not present
      call set_file_names( filein, fileout, &
                           fileInput, fileOutput, unitInput )
      filedebug = trim(fileout) // ".debug"

#ifdef FDF_DEBUG
      ! To monitor the parsing and the build-up of the
      ! fdf data structures in all nodes
      call fdf_setdebug(2,filedebug)
      call fdf_setoutput(2,fileout)   ! All nodes print output
#endif

!!      call fdf_set_output_file(fileout)

      call fdf_input(filein)

      fdf_started = .TRUE.

#ifndef FDF_DEBUG
      ! Flags within the fdf file itself.

      ! At this point only the final fdf data structure will be shown,
      ! for level >= 2
      debug_level = fdf_get('fdf-debug', 0)
      call fdf_setdebug(debug_level,filedebug)

      ! The default is to have output only in the master node
      output_level = fdf_get('fdf-output', 1)
      call fdf_setoutput(output_level,fileout)
#endif

      if (debug_level >= 2) call fdf_printfdf()

#ifdef _MPI_
      call fdf_mpi_finalize()
#endif
!$OMP END SINGLE
      RETURN
!------------------------------------------------------------------------- END
      END SUBROUTINE fdf_init


      SUBROUTINE set_file_names( fileIn, fileOut, &
                                 optFileIn, optFileOut, unitIn )
      ! If present, copies input arguments optFileIn/Out to fileIn/Out.
      ! If absent, generates In/Out file names. If unitIn is present, and it is
      ! a named file, returns it as fileIn. If not, it copies input to a new
      ! file and returns its name. If .not.present(unitIn) => unitIn=5.
      ! If optFileIn is present, unitIn is ignored.
      implicit none
      character(len=*),intent(out):: &
        fileIn,    &! Name of file to be used as input
        fileOut     ! Name of file to be used as output
      character(len=*),optional,intent(in):: &
        optFileIn, &! Optional argument with input file name
        optFileOut  ! Optional argument with output file name
      integer,optional,intent(in):: &
        unitIn      ! Optional input file unit (not used if present(optFileIn))

      integer:: count, ierr, iostat, iu, iuIn
      logical:: named, opened
      character(len=MAX_LENGTH*2) line
      character(len=maxFileNameLength) fileName

!------------------------------------------------------------------------- BEGIN
#ifdef _MPI_
      if (rank==0) then
#endif

      ! Find a job-specific number
      call system_clock( count )
      count = mod(count,100000)

      ! Set output file name
      if (present(optFileOut)) then
        if (len(trim(optFileOut)) > len(fileOut)) &
          call die('FDF module: set_file_names', &
                   'Parameter maxFileNameLength too small.' // &
                   'Terminating.', THIS_FILE, __LINE__, fdf_err, rc=ierr)
        fileOut = optFileOut
      else                  ! set a job-specific file name
        write(fileOut,'(a,i5.5,a)') 'fdf_',count,'.log'
      endif

      ! Set input file
      if (present(optFileIn)) then     ! just copy the file name
        if (len(trim(optFileIn)) > len(fileIn)) &
          call die('FDF module: set_file_names', &
                   'Parameter maxFileNameLength too small.' // &
                   'Terminating.', THIS_FILE, __LINE__, fdf_err, rc=ierr)
        fileIn = optFileIn
      else                             ! find or set a file name

        ! Find input file unit
        if (present(unitIn)) then      ! use given unit (possibly 5)
          iuIn = unitIn
        else                           ! assume standard input
          iuIn = 5
        endif

        ! Find file name associated with given unit
        if (iuIn==5) then              ! no valid file name
           fileName = ' '
        else                           ! check if this is a named file
          inquire(unit=iuIn,opened=opened)
          if (opened) then
            inquire(unit=iuIn,named=named)
            if (named) then            ! inquire file name
              inquire(unit=iuIn,name=fileName)
            else                       ! no valid file name
              fileName = ' '
            endif ! (named)
          else
            call die('FDF module: set_file_names', 'Input unit not opened.' // &
                     'Terminating.', THIS_FILE, __LINE__, fdf_err, rc=ierr)
          endif ! (opened)
        endif ! (iuIn==5)

        ! Set input file name, possibly after copying input to it
        if (fileName==' ') then                       ! not a valid file
          write(fileIn,'(a,i5.5,a)') &
            'INPUT_TMP_',count,'.fdf'                 ! new file's name
          call io_assign(iu)                          ! new file's unit
          open(iu,file=trim(fileIn),form='formatted') ! open new file
          do
            read(iuIn,iostat=iostat,fmt='(a)') line   ! read line from old unit
            if (iostat/=0 ) exit
            write(iu,'(a)') trim(line)                ! write line to new file
          enddo
          call io_close(iu)                           ! close new file
        else                                          ! valid file
          fileIn = fileName
        endif ! (fileName=='stdin')

      endif ! (present(optFileIn))

#ifdef _MPI_
endif ! (rank==0)
#endif
!--------------------------------------------------------------------------- END
      END SUBROUTINE set_file_names

!
!   Initialize MPI subsystem if the application calling/using FDF
!   library is not running with MPI enabled.
!
#ifdef _MPI_
      SUBROUTINE fdf_mpi_init()
        !
        use mpi
        implicit none
  !--------------------------------------------------------------- Local Variables
        integer(ip) :: ierr

  !------------------------------------------------------------------------- BEGIN
        call MPI_Initialized(mpiflag, ierr)
        if (ierr .ne. MPI_SUCCESS) then
          call die('FDF module: fdf_mpi_init', 'Error initializing MPI system.' // &
                   'Terminating.', THIS_FILE, __LINE__, fdf_err, rc=ierr)
        endif

        if (.not. mpiflag) then
          call MPI_Init(ierr)
          if (ierr .ne. MPI_SUCCESS) then
            call die('FDF module: fdf_mpi_init', 'Error initializing MPI system.' // &
                     'Terminating.', THIS_FILE, __LINE__, fdf_err, rc=ierr)
          endif
        endif

        call MPI_Comm_Rank(MPI_COMM_WORLD, rank, ierr)
        if (ierr .ne. MPI_SUCCESS) then
          call die('FDF module: fdf_mpi_init', 'Error getting MPI comm rank.' // &
                   'Terminating.', THIS_FILE, __LINE__, fdf_err, rc= ierr)
        endif

        call MPI_Comm_Size(MPI_COMM_WORLD, ntasks, ierr)
        if (ierr .ne. MPI_SUCCESS) then
          call die('FDF module: fdf_mpi_init', 'Error getting MPI comm size.' // &
                   'Terminating.', THIS_FILE, __LINE__, fdf_err, rc=ierr)
        endif
        RETURN
  !--------------------------------------------------------------------------- END
      END SUBROUTINE fdf_mpi_init
#endif

  !
  !   Finalize MPI subsystem if the application calling/using FDF
  !   library is not running with MPI enabled.
  !
#ifdef _MPI_
      SUBROUTINE fdf_mpi_finalize()
        !
        use mpi
        implicit none
  !--------------------------------------------------------------- Local Variables
        integer(ip)   :: ierr

  !------------------------------------------------------------------------- BEGIN
        if (.not. mpiflag) then
          call MPI_Finalize(ierr)
          if (ierr .ne. MPI_SUCCESS) then
            call die('FDF module: fdf_mpi_finalize', 'Error finalizing MPI system.' // &
                     'Terminating.', THIS_FILE, __LINE__, fdf_err, rc=ierr)
          endif
        endif

        RETURN
  !--------------------------------------------------------------------------- END
      END SUBROUTINE fdf_mpi_finalize
#endif

!
!   Reads the input file depending on the configuration of the system:
!   Shared Filesystem, Cluster Filesystem or Blocking input access
!
    SUBROUTINE fdf_input(filein)
      implicit none
!--------------------------------------------------------------- Input Variables
      character(*)  :: filein

!------------------------------------------------------------------------- BEGIN
#ifdef CLUSTER
!!      call fdf_readcluster(filein)
  call setup_fdf_cluster(filein)
#elif defined(BLOCKING)
  call fdf_readblocking(filein)
#else
  call fdf_read_custom(filein)  ! fdf_read(filein) was the original
#endif

      if (fdf_output) write(fdf_out,'(a,a,a,i3)') '#FDF module: Opened ', filein,   &
                                  ' for input. Unit:', fdf_in(1)

      RETURN
!--------------------------------------------------------------------------- END
    END SUBROUTINE fdf_input

!
!   Reading code for Cluster Architecture, where the filesystem
!   is not shared along all the nodes in the system.
!
#ifdef CLUSTER
  SUBROUTINE fdf_readcluster(filein)
    !
    use mpi
    implicit none
!--------------------------------------------------------------- Input Variables
    character(*)  :: filein

!--------------------------------------------------------------- Local Variables
    character(80)  :: msg
    character(256) :: fileinTmp
    integer(ip)    :: ierr, texist_send, texist_recv
#ifdef SOPHISTICATED_SEARCH
    logical        :: file_exist
#endif
!------------------------------------------------------------------------- BEGIN
!     Tests if the running node has the input file:
!       If found: texist_send = rank
!       Else    : texist_send = error_code (ntasks + 1)

#ifdef SOPHISTICATED_SEARCH
    INQUIRE(file=filein, exist=file_exist)
    if (file_exist) then
      texist_send = rank
    else
      texist_send = ntasks + 1
    endif

    call MPI_AllReduce(texist_send, texist_recv, 1, MPI_INTEGER,      &
                       MPI_MIN, MPI_COMM_WORLD, ierr)
    if (ierr .ne. MPI_SUCCESS) then
      call die('FDF module: fdf_readcluster', 'Error in MPI_AllReduce (task_exist).' //  &
             'Terminating.', THIS_FILE, __LINE__, fdf_err, rc=ierr)
    endif
#else
!
!     Simplify: Assume node 0 has the file
!
      texist_recv = 0
#endif

!     The node owner of the input file send the data to the other ones
!     if none node has an input file with such name abort the application
    if (texist_recv .eq. ntasks + 1) then
      write(msg,*) 'No node found in the cluster with ',              &
                 ' input file ', filein,'. Terminating.'
      call die('FDF module: fdf_readcluster', msg,                    &
             THIS_FILE, __LINE__, fdf_err, rc=ierr)
    else
      if (texist_recv .eq. rank) then
        call fdf_read(filein)
        call fdf_sendInput()
!!Debug          call MPI_Barrier( MPI_COMM_WORLD, ierr )
        if (fdf_output) write(fdf_out,*) '#FDF module: Node', rank, 'reading/sending', &
                         ' input file ', filein
      else
        call fdf_recvInput(texist_recv, filein, fileinTmp)
!!Debug          call MPI_Barrier( MPI_COMM_WORLD, ierr )
        call fdf_read(fileinTmp)
        if (fdf_output) write(fdf_out,*) '#FDF module: Node', rank, 'receiving input', &
                       ' file from', texist_recv, 'to ', TRIM(fileinTmp)
      endif
    endif
    RETURN
!--------------------------------------------------------------------------- END
  END SUBROUTINE fdf_readcluster
!
!
  SUBROUTINE setup_fdf_cluster(filein)
!
!     A more efficient alternative to fdf_sendInput/fdf_recvInput
!     that avoids the creation of scratch files by non-root nodes.
!
!     Alberto Garcia, April 2011

    !
    use mpi
    implicit none
!--------------------------------------------------------------- Input Variables
    character(*), intent(in)  :: filein   ! File name

!--------------------------------------------------------------- Local Variables
    character(80)  :: msg
    character(256) :: fileinTmp
    integer(ip)    :: ierr, texist_send, reading_node
#ifdef SOPHISTICATED_SEARCH
    logical        :: file_exist
#endif
!------------------------------------------------------------------------- BEGIN
!     Tests if the running node has the input file:
!       If found: texist_send = rank
!       Else    : texist_send = error_code (ntasks + 1)

#ifdef SOPHISTICATED_SEARCH
    INQUIRE(file=filein, exist=file_exist)
    if (file_exist) then
      texist_send = rank
    else
      texist_send = ntasks + 1
    endif

    call MPI_AllReduce(texist_send, reading_node, 1, MPI_INTEGER,      &
                     MPI_MIN, MPI_COMM_WORLD, ierr)
    if (ierr .ne. MPI_SUCCESS) then
      call die('FDF module: fdf_readcluster', 'Error in MPI_AllReduce (task_exist).' //  &
             'Terminating.', THIS_FILE, __LINE__, fdf_err, rc=ierr)
    endif
#else
!
!     Simplify: Assume node 0 has the file
!
    reading_node = 0
#endif

!     If no node has an input file with such name abort the application
    if (reading_node .eq. ntasks + 1) then
      write(msg,*) 'No node found in the cluster with ',              &
          ' input file ', filein,'. Terminating.'
      call die('FDF module: fdf_readcluster', msg,                    &
          THIS_FILE, __LINE__, fdf_err, rc=ierr)
    endif

!     The root node reads and digests the input file
!     and sends the serialized fdf_file data structure to the other nodes

    if (rank == reading_node) then
      call fdf_read(filein)
    endif
    call broadcast_fdf_struct(reading_node)

    RETURN
!--------------------------------------------------------------------------- END
  END SUBROUTINE setup_fdf_cluster

!
!   Broadcast complete fdf structure
!
  SUBROUTINE broadcast_fdf_struct(reading_node)
    !
    use mpi
    implicit none

    integer, intent(in)       :: reading_node         ! Node which contains the struct

!--------------------------------------------------------------- Local Variables
    character, pointer        :: bufferFDF(:) => null()
    integer(ip)               :: i, j, k, ierr, nlines

!------------------------------------------------------------------------- BEGIN

    if (rank == reading_node) then
      nlines = file_in%nlines
    endif

    call MPI_Bcast(nlines, 1,                                 &
                  MPI_INTEGER, reading_node, MPI_COMM_WORLD, ierr)
    if (ierr .ne. MPI_SUCCESS) then
      call die('FDF module: broadcast_fdf', 'Error Broadcasting nlines.' // &
              'Terminating.', THIS_FILE, __LINE__, fdf_err, rc=ierr)
    endif

    ALLOCATE(bufferFDF(nlines*SERIALIZED_LENGTH), stat=ierr)
    if (ierr .ne. 0) then
      call die('FDF module: broadcast_fdf', 'Error allocating bufferFDF', &
              THIS_FILE, __LINE__, fdf_err, rc=ierr)
    endif

    if (rank == reading_node) then
      call serialize_fdf_struct(bufferFDF)
    endif

    call MPI_Bcast(bufferFDF, size(bufferFDF),              &
                  MPI_CHARACTER, reading_node, MPI_COMM_WORLD, ierr)
    if (ierr .ne. MPI_SUCCESS) then
      call die('FDF module: broadcast_fdf', 'Error Broadcasting bufferFDF.' // &
              'Terminating.', THIS_FILE, __LINE__, fdf_err, rc=ierr)
    endif
    if (rank /= reading_node) then
      call recreate_fdf_struct(nlines,bufferFDF)
    endif

    DEALLOCATE(bufferFDF)

    RETURN
!--------------------------------------------------------------------------- END
  END SUBROUTINE broadcast_fdf_struct
#endif
!
!   Reading code for Blocking read. The reading of the input file is
!   splitted in several steps of size BLOCKSIZE. The number of steps
!   needed to read the file in all the nodes depends on the total nodes,
!   being ntasks/BLOCKINGSIZE.
!
!   With this method we can avoid a colapse of a global parallel filesystem
!   if the execution of the application runs over a huge amount of nodes.
!
#ifdef BLOCKING
# ifndef BLOCKSIZE
#   define BLOCKSIZE 2
# endif
  SUBROUTINE fdf_readblocking(filein)
    !
    use mpi
    implicit none
  !--------------------------------------------------------------- Input Variables
    character(*) :: filein

  !--------------------------------------------------------------- Local Variables
    integer(ip)  :: i, ierr

  !------------------------------------------------------------------------- BEGIN

    do i= 0, ntasks-1, BLOCKSIZE
      if ((rank .ge. i) .and. (rank .le. i+BLOCKSIZE-1)) then
        call fdf_read(filein)
        if (fdf_output) write(fdf_out,*) '#FDF module: Task', rank, 'reading input', &
                        ' file ', filein, ' in step', (i/BLOCKSIZE)+1
      endif

      call MPI_Barrier(MPI_COMM_WORLD, ierr)
      if (ierr .ne. MPI_SUCCESS) then
        call die('FDF module: fdf_readblocking', 'Error in MPI_Barrier (fdf_read).' //    &
                'Terminating.', THIS_FILE, __LINE__, fdf_err, rc=ierr)
      endif
    enddo

    RETURN
  !--------------------------------------------------------------------------- END
  END SUBROUTINE fdf_readblocking
#endif

!


    RECURSIVE SUBROUTINE fdf_read_custom(filein, blocklabel)
      implicit none
!--------------------------------------------------------------- Input Variables
      character(*)               :: filein
      character(*), optional     :: blocklabel

!--------------------------------------------------------------- Local Variables
      logical                    :: dump
      logical, allocatable       :: found(:)
      character(80)              :: msg
      character(len=MAX_LENGTH)  :: label, inc_file, modulename
      character(len=MAX_LENGTH*2):: line
      integer(ip)                :: i, ierr, ntok, ind_less, nlstart, nlend, counter
      type(parsed_line), pointer :: pline

!------------------------------------------------------------------------- BEGIN
!     Open reading input file
      call fdf_open(filein)
      counter = 0
!     Read each input data line
      if (PRESENT(blocklabel)) then
        label = blocklabel
      else
        label = ' '
      endif
      do while (fdf_readline(line))

!       Check if valid data (tokens, non-blank)
        pline => digest(line)
        ntok = ntokens(pline)
        if (ntok .ne. 0) then

!         Find different special cases in the input files
!         (%block, %endblock, %include, Label1 Label2 ... < Filename)

!         %block directive
          ind_less = search('<', pline)
          if (search('%block', pline) .eq. 1) then
!            print*, "debug::library::  inside block construct "
!           No label found in %block directive
            if (ntok .eq. 1) then
              write(msg,*) '%block label not found in ', TRIM(filein)
              call die('FDF module: fdf_read', msg,                     &
                       THIS_FILE, __LINE__, fdf_err)
            endif

!           %block Label < Filename [ %dump ]
            if (ind_less .eq. 3) then

              if (ntok .ge. 4) then
!               Test if %dump is present
                if (search('%dump', pline) .eq. 5) then
                  dump = .TRUE.
                else
                  dump = .FALSE.
                endif

!               Add begin, body and end sections of block
                label = tokens(pline, 2)
                inc_file  = tokens(pline, 4)
                call destroy(pline)
                line = '%block ' // trim(label) // ' ' // trim(inc_file)
                pline => digest(line)
                call setmorphol(1, 'b', pline)
                call setmorphol(2, 'l', pline)
                call fdf_addtoken(line, pline)
                nullify(pline) ! it is stored in line

                nlstart = file_in%nlines
!                print*, "debug:: nlstart ", nlstart

                call fdf_read(inc_file, label)

!               Warn if block 'label' is empty
                if ((nlstart - file_in%nlines) .eq. 0) then
                  write(msg,*) 'FDF module: fdf_read: block ',          &
                               TRIM(label), ' is empty...'
                  call warn(msg)
                endif

                line = '%endblock ' // label
                pline => digest(line)
                call setmorphol(1, 'e', pline)
                call setmorphol(2, 'l', pline)
                call fdf_addtoken(line, pline)
                nullify(pline) ! it is stored in line

!               Dump included file to fileout
                if (dump) call fdf_dump(label)
                label = ' '

!             Filename not found in %block directive
              else
                write(msg,*) '%block filename not found in ', TRIM(filein)
                call die('FDF module: fdf_read', msg,                   &
                         THIS_FILE, __LINE__, fdf_err)
              endif

!           %block Label
            elseif (ind_less .eq. -1) then
              label = tokens(pline, 2)
              call setmorphol(1, 'b', pline)
              call setmorphol(2, 'l', pline)
              call fdf_addtoken(line, pline)
              nullify(pline) ! it is stored in line
              nlstart = file_in%nlines

!           Bad format in %block directive
            else
              write(msg,*) 'Bad ''<'' %block format in ', TRIM(filein)
              call die('FDF module: fdf_read', msg,                     &
                       THIS_FILE, __LINE__, fdf_err)
            endif

!         %endblock directive
          elseif (search('%endblock', pline) .eq. 1) then
!            print*, "debug::library::  inside endblock construct "
!           Check if %block exists before %endblock
            if (label .eq. ' ') then
              write(msg,*) 'Bad %endblock found in ', TRIM(filein)
              call die('FDF module: fdf_read', msg,                     &
                       THIS_FILE, __LINE__, fdf_err)
            else
!             Warn if block 'label' is empty
              if ((nlstart - file_in%nlines) .eq. 0) then
                write(msg,*) 'FDF module: fdf_read: block ',            &
                             TRIM(label), ' is empty...'
                call warn(msg)
              endif

              call destroy(pline)
              line = '%endblock ' // label
              pline => digest(line)
              call setmorphol(1, 'e', pline)
              call setmorphol(2, 'l', pline)
              call fdf_addtoken(line, pline)
              nullify(pline) ! it is stored in line
              label = ' '
            endif

!         %module construction
          elseif (search('%module', pline) .eq. 1) then
            counter = counter + 1
            pline => digest(line)
!            print*, "debug::library::  inside module construct ", pline%line
            call setmorphol(1, 'b', pline)
            call setmorphol(2, 'l', pline)
            call fdf_addtoken(line, pline)
            nullify(pline) ! it is stored in line
            nlstart = file_in%nlines
            modulename = trim(line(8:))
            write(msg,'(A,1x,i1,1x,A)') "Module #", counter , trim(line(8:))
            nullify(pline) ! it is stored in line
!           No label found in %module directive
            if (ntok .eq. 1) then
              write(msg,*) '%module label not found in ', TRIM(filein)
              call die('FDF module: fdf_read', msg,                     &
                        THIS_FILE, __LINE__, fdf_err)
            endif

! !           Add begin, body and end sections of module
!           %module Label
!              call destroy(pline)
              line = '%module ' // label
              !nlstart = file_in%nlines
!              print*, "debug:: beginning line no in the input file ", nlstart

!             structure created

!          %endmodule directive
          elseif (search('%endmodule', pline) .eq. 1) then
!            print*, "debug::library::  inside endmodule construct "
!           Check if %module exists before %endmodule
!              call destroy(pline)
              line = '%endmodule ' // label
              nlend = file_in%nlines
!              print*, "debug:: ending line no in the input file ", nlend


! custom added part ends here

!         %include Filename directive
          elseif (search('%include', pline) .eq. 1) then
!           Check if include filename is specified
            if (ntok .eq. 1) then
              write(msg,*) 'Filename on %include not found in ', TRIM(filein)
              call die('FDF module: fdf_read', msg,                     &
                       THIS_FILE, __LINE__, fdf_err)
            else
              inc_file = tokens(pline, 2)
              call fdf_read(inc_file)
            endif

            ! Clean pline (we simply insert the next file)
            call destroy(pline)


!         Label1 Label2 ... < Filename directive
          elseif (ind_less .ne. -1) then
!           Check if '<' is in a valid position
            if (ind_less .eq. 1) then
              write(msg,*) 'Bad ''<'' found in ', TRIM(filein)
              call die('FDF module: fdf_read', msg,                     &
                       THIS_FILE, __LINE__, fdf_err)

!           Check if '<' filename is specified
            elseif (ind_less .eq. ntok) then
              write(msg,*) 'Filename not found after ''<'' in ', TRIM(filein)
              call die('FDF module: fdf_read', msg,                     &
                       THIS_FILE, __LINE__, fdf_err)

            else
!             Search label(s) in Filename
              inc_file = tokens(pline, ind_less+1)
              ALLOCATE(found(ind_less-1), stat=ierr)
              if (ierr .ne. 0) then
                call die('FDF module: fdf_read', 'Error allocating found', &
                         THIS_FILE, __LINE__, fdf_err, rc=ierr)
              endif


!             If label(s) not found in such Filename throw an error


              found = .FALSE.
              if (.not. fdf_readlabel(ind_less-1, pline,                &
                                      inc_file, found)) then
                 i = 1
                 do while ((i .le. ind_less-1) .and. (found(i)))
                    i = i + 1
                 enddo
                 label = tokens(pline, i)
                 write(msg,*) 'Label ', TRIM(label),                     &
                             ' not found in ', TRIM(inc_file)
                 call die('FDF module: fdf_read', msg,                   &
                         THIS_FILE, __LINE__, fdf_err)
              endif

              call destroy(pline)
              DEALLOCATE(found)
            endif


!         Add remaining kind of tokens to dynamic list as labels
          else
            if (label .eq. ' ') call setmorphol(1, 'l', pline)
            call fdf_addtoken(line, pline)
            nullify(pline) ! it is stored in line
          endif
        else
!         Destroy parsed_line structure if no elements
          call destroy(pline)
        endif
      enddo

!     Close one level of input file
      if ((.not. PRESENT(blocklabel)) .and. (label .ne. ' ')) then
        write(msg,*) '%endblock ', TRIM(label),                         &
                     ' not found in ', TRIM(filein)
        call die('FDF module: fdf_read', msg, THIS_FILE, __LINE__, fdf_err)
      endif
      call fdf_close()

      RETURN
!--------------------------------------------------------------------------- END
    END SUBROUTINE fdf_read_custom






!   Read an input file (and include files) and builds memory
!   structure that will contain the data and will help in searching
!
    RECURSIVE SUBROUTINE fdf_read(filein, blocklabel)
      implicit none
!--------------------------------------------------------------- Input Variables
      character(*)               :: filein
      character(*), optional     :: blocklabel

!--------------------------------------------------------------- Local Variables
      logical                    :: dump
      logical, allocatable       :: found(:)
      character(80)              :: msg
      character(len=MAX_LENGTH)  :: label, inc_file
      character(len=MAX_LENGTH*2):: line
      integer(ip)                :: i, ierr, ntok, ind_less, nlstart
      type(parsed_line), pointer :: pline

!------------------------------------------------------------------------- BEGIN
!     Open reading input file
      call fdf_open(filein)

!     Read each input data line
      if (PRESENT(blocklabel)) then
        label = blocklabel
      else
        label = ' '
      endif
      do while (fdf_readline(line))

!       Check if valid data (tokens, non-blank)
        pline => digest(line)
        ntok = ntokens(pline)
        if (ntok .ne. 0) then

!         Find different special cases in the input files
!         (%block, %endblock, %include, Label1 Label2 ... < Filename)

!         %block directive
          ind_less = search('<', pline)
          if (search('%block', pline) .eq. 1) then

!           No label found in %block directive
            if (ntok .eq. 1) then
              write(msg,*) '%block label not found in ', TRIM(filein)
              call die('FDF module: fdf_read', msg,                     &
                       THIS_FILE, __LINE__, fdf_err)
            endif

!           %block Label < Filename [ %dump ]
            if (ind_less .eq. 3) then

              if (ntok .ge. 4) then
!               Test if %dump is present
                if (search('%dump', pline) .eq. 5) then
                  dump = .TRUE.
                else
                  dump = .FALSE.
                endif

!               Add begin, body and end sections of block
                label = tokens(pline, 2)
                inc_file  = tokens(pline, 4)
                call destroy(pline)
                line = '%block ' // label
                pline => digest(line)
                call setmorphol(1, 'b', pline)
                call setmorphol(2, 'l', pline)
                call fdf_addtoken(line, pline)
                nullify(pline) ! it is stored in line

                nlstart = file_in%nlines
                call fdf_read(inc_file, label)

!               Warn if block 'label' is empty
                if ((nlstart - file_in%nlines) .eq. 0) then
                  write(msg,*) 'FDF module: fdf_read: block ',          &
                               TRIM(label), ' is empty...'
                  call warn(msg)
                endif

                line = '%endblock ' // label
                pline => digest(line)
                call setmorphol(1, 'e', pline)
                call setmorphol(2, 'l', pline)
                call fdf_addtoken(line, pline)
                nullify(pline) ! it is stored in line

!               Dump included file to fileout
                if (dump) call fdf_dump(label)
                label = ' '

!             Filename not found in %block directive
              else
                write(msg,*) '%block filename not found in ', TRIM(filein)
                call die('FDF module: fdf_read', msg,                   &
                         THIS_FILE, __LINE__, fdf_err)
              endif

!           %block Label
            elseif (ind_less .eq. -1) then
              label = tokens(pline, 2)
              call setmorphol(1, 'b', pline)
              call setmorphol(2, 'l', pline)
              call fdf_addtoken(line, pline)
              nullify(pline) ! it is stored in line
              nlstart = file_in%nlines

!           Bad format in %block directive
            else
              write(msg,*) 'Bad ''<'' %block format in ', TRIM(filein)
              call die('FDF module: fdf_read', msg,                     &
                       THIS_FILE, __LINE__, fdf_err)
            endif

!         %endblock directive
          elseif (search('%endblock', pline) .eq. 1) then
!           Check if %block exists before %endblock
            if (label .eq. ' ') then
              write(msg,*) 'Bad %endblock found in ', TRIM(filein)
              call die('FDF module: fdf_read', msg,                     &
                       THIS_FILE, __LINE__, fdf_err)
            else
!             Warn if block 'label' is empty
              if ((nlstart - file_in%nlines) .eq. 0) then
                write(msg,*) 'FDF module: fdf_read: block ',            &
                             TRIM(label), ' is empty...'
                call warn(msg)
              endif

              call destroy(pline)
              line = '%endblock ' // label
              pline => digest(line)
              call setmorphol(1, 'e', pline)
              call setmorphol(2, 'l', pline)
              call fdf_addtoken(line, pline)
              nullify(pline) ! it is stored in line
              label = ' '
            endif

!         %include Filename directive
          elseif (search('%include', pline) .eq. 1) then
!           Check if include filename is specified
            if (ntok .eq. 1) then
              write(msg,*) 'Filename on %include not found in ', TRIM(filein)
              call die('FDF module: fdf_read', msg,                     &
                       THIS_FILE, __LINE__, fdf_err)
            else
              inc_file = tokens(pline, 2)
              call fdf_read(inc_file)
            endif

            ! Clean pline (we simply insert the next file)
            call destroy(pline)

!         Label1 Label2 ... < Filename directive
          elseif (ind_less .ne. -1) then
!           Check if '<' is in a valid position
            if (ind_less .eq. 1) then
              write(msg,*) 'Bad ''<'' found in ', TRIM(filein)
              call die('FDF module: fdf_read', msg,                     &
                       THIS_FILE, __LINE__, fdf_err)

!           Check if '<' filename is specified
            elseif (ind_less .eq. ntok) then
              write(msg,*) 'Filename not found after ''<'' in ', TRIM(filein)
              call die('FDF module: fdf_read', msg,                     &
                       THIS_FILE, __LINE__, fdf_err)
            else
!             Search label(s) in Filename
              inc_file = tokens(pline, ind_less+1)
              ALLOCATE(found(ind_less-1), stat=ierr)
              if (ierr .ne. 0) then
                call die('FDF module: fdf_read', 'Error allocating found', &
                         THIS_FILE, __LINE__, fdf_err, rc=ierr)
              endif

!             If label(s) not found in such Filename throw an error
              found = .FALSE.
              if (.not. fdf_readlabel(ind_less-1, pline,                &
                                      inc_file, found)) then
                 i = 1
                 do while ((i .le. ind_less-1) .and. (found(i)))
                    i = i + 1
                 enddo
                 label = tokens(pline, i)
                 write(msg,*) 'Label ', TRIM(label),                     &
                             ' not found in ', TRIM(inc_file)
                 call die('FDF module: fdf_read', msg,                   &
                         THIS_FILE, __LINE__, fdf_err)
              endif

              call destroy(pline)
              DEALLOCATE(found)
            endif

!         Add remaining kind of tokens to dynamic list as labels
          else
            if (label .eq. ' ') call setmorphol(1, 'l', pline)
            call fdf_addtoken(line, pline)
            nullify(pline) ! it is stored in line
          endif
        else
!         Destroy parsed_line structure if no elements
          call destroy(pline)
        endif
      enddo

!     Close one level of input file
      if ((.not. PRESENT(blocklabel)) .and. (label .ne. ' ')) then
        write(msg,*) '%endblock ', TRIM(label),                         &
                     ' not found in ', TRIM(filein)
        call die('FDF module: fdf_read', msg, THIS_FILE, __LINE__, fdf_err)
      endif
      call fdf_close()

      RETURN
!--------------------------------------------------------------------------- END
    END SUBROUTINE fdf_read


!   Read an input xyz coordinate file and builds memory
!   structure that will contain the data and will help in searching
!
    RECURSIVE SUBROUTINE fdf_read_xyz(filein, blocklabel)
      implicit none
!--------------------------------------------------------------- Input Variables
      character(*)               :: filein
      character(*), optional     :: blocklabel

!--------------------------------------------------------------- Local Variables
      logical                    :: dump
      logical, allocatable       :: found(:)
      character(80)              :: msg
      character(len=MAX_LENGTH)  :: label, inc_file
      character(len=MAX_LENGTH*2):: line
      integer(ip)                :: i, ierr, ntok, ind_less, nlstart
      type(parsed_line), pointer :: pline

!------------------------------------------------------------------------- BEGIN
!     Open reading input file
      call fdf_open(filein)

!     Read each input data line
      if (PRESENT(blocklabel)) then
        label = blocklabel
      else
        label = ' '
      endif
      do while (fdf_readline(line))

!       Check if valid data (tokens, non-blank)
        pline => digest(line)
        ntok = ntokens(pline)
        if (ntok .ne. 0) then

!         Find different special cases in the input files
!         (%block, %endblock, %include, Label1 Label2 ... < Filename)

!         %molecule directive
          ind_less = search('<', pline)
          if (search('%molecule', pline) .eq. 1) then

!           No label found in %block directive
            if (ntok .eq. 1) then
              write(msg,*) '%molecule label not found in ', TRIM(filein)
              call die('FDF module: fdf_read_xyz', msg,                     &
                       THIS_FILE, __LINE__, fdf_err)
            endif

!           %molecule Label < Filename [ %dump ]
            if (ind_less .eq. 3) then

              if (ntok .ge. 4) then
!               Test if %dump is present
                if (search('%dump', pline) .eq. 5) then
                  dump = .TRUE.
                else
                  dump = .FALSE.
                endif

!               Add begin, body and end sections of molecule block
                label = tokens(pline, 2)
                inc_file  = tokens(pline, 4)
                call destroy(pline)
                line = '%molecule ' // label
                write(*,*) "debug printing the line", line
                pline => digest(line)
                call setmorphol(1, 'b', pline)
                call setmorphol(2, 'l', pline)
                call fdf_addtoken(line, pline)
                nullify(pline) ! it is stored in line

                nlstart = file_in%nlines
                call fdf_read(inc_file, label)

!               Warn if block 'label' is empty
                if ((nlstart - file_in%nlines) .eq. 0) then
                  write(msg,*) 'FDF module: fdf_read_xyz: block ',          &
                               TRIM(label), ' is empty...'
                  call warn(msg)
                endif

                line = '%endmolecule ' // label
                pline => digest(line)
                call setmorphol(1, 'e', pline)
                call setmorphol(2, 'l', pline)
                call fdf_addtoken(line, pline)
                nullify(pline) ! it is stored in line

!               Dump included file to fileout
                if (dump) call fdf_dump(label)
                label = ' '

!             Filename not found in %block directive
              else
                write(msg,*) '%molecule filename not found in ', TRIM(filein)
                call die('FDF module: fdf_read_xyz', msg,                   &
                         THIS_FILE, __LINE__, fdf_err)
              endif

!           %molecule Label
            elseif (ind_less .eq. -1) then
              label = tokens(pline, 2)
              call setmorphol(1, 'b', pline)
              call setmorphol(2, 'l', pline)
              call fdf_addtoken(line, pline)
              nullify(pline) ! it is stored in line
              nlstart = file_in%nlines

!           Bad format in %block directive
            else
              write(msg,*) 'Bad ''<'' %molecule format in ', TRIM(filein)
              call die('FDF module: fdf_read_xyz', msg,                     &
                       THIS_FILE, __LINE__, fdf_err)
            endif

!         %endblock directive
          elseif (search('%endmolecule', pline) .eq. 1) then
!           Check if %block exists before %endblock
            if (label .eq. ' ') then
              write(msg,*) 'Bad %endmolecule found in ', TRIM(filein)
              call die('FDF module: fdf_read_xyz', msg,                     &
                       THIS_FILE, __LINE__, fdf_err)
            else
!             Warn if block 'label' is empty
              if ((nlstart - file_in%nlines) .eq. 0) then
                write(msg,*) 'FDF module: fdf_read_xyz: block ',            &
                             TRIM(label), ' is empty...'
                call warn(msg)
              endif

              call destroy(pline)
              line = '%endbmolecule ' // label
              pline => digest(line)
              call setmorphol(1, 'e', pline)
              call setmorphol(2, 'l', pline)
              call fdf_addtoken(line, pline)
              nullify(pline) ! it is stored in line
              label = ' '
            endif

!         %include Filename directive
          elseif (search('%include', pline) .eq. 1) then
!           Check if include filename is specified
            if (ntok .eq. 1) then
              write(msg,*) 'Filename on %include not found in ', TRIM(filein)
              call die('FDF module: fdf_read_xyz', msg,                     &
                       THIS_FILE, __LINE__, fdf_err)
            else
              inc_file = tokens(pline, 2)
              call fdf_read(inc_file)
            endif

            ! Clean pline (we simply insert the next file)
            call destroy(pline)

!         Label1 Label2 ... < Filename directive
          elseif (ind_less .ne. -1) then
!           Check if '<' is in a valid position
            if (ind_less .eq. 1) then
              write(msg,*) 'Bad ''<'' found in ', TRIM(filein)
              call die('FDF module: fdf_read_xyz', msg,                     &
                       THIS_FILE, __LINE__, fdf_err)

!           Check if '<' filename is specified
            elseif (ind_less .eq. ntok) then
              write(msg,*) 'Filename not found after ''<'' in ', TRIM(filein)
              call die('FDF module: fdf_read_xyz', msg,                     &
                       THIS_FILE, __LINE__, fdf_err)
            else
!             Search label(s) in Filename
              inc_file = tokens(pline, ind_less+1)
              ALLOCATE(found(ind_less-1), stat=ierr)
              if (ierr .ne. 0) then
                call die('FDF module: fdf_read_xyz', 'Error allocating found', &
                         THIS_FILE, __LINE__, fdf_err, rc=ierr)
              endif

!             If label(s) not found in such Filename throw an error
              found = .FALSE.
              if (.not. fdf_readlabel(ind_less-1, pline,                &
                                      inc_file, found)) then
                 i = 1
                 do while ((i .le. ind_less-1) .and. (found(i)))
                    i = i + 1
                 enddo
                 label = tokens(pline, i)
                 write(msg,*) 'Label ', TRIM(label),                     &
                             ' not found in ', TRIM(inc_file)
                 call die('FDF module: fdf_read_xyz', msg,                   &
                         THIS_FILE, __LINE__, fdf_err)
              endif

              call destroy(pline)
              DEALLOCATE(found)
            endif

!         Add remaining kind of tokens to dynamic list as labels
          else
            if (label .eq. ' ') call setmorphol(1, 'l', pline)
            call fdf_addtoken(line, pline)
            nullify(pline) ! it is stored in line
          endif
        else
!         Destroy parsed_line structure if no elements
          call destroy(pline)
        endif
      enddo

!     Close one level of input file
      if ((.not. PRESENT(blocklabel)) .and. (label .ne. ' ')) then
        write(msg,*) '%endmolecule ', TRIM(label),                         &
                     ' not found in ', TRIM(filein)
        call die('FDF module: fdf_read_xyz', msg, THIS_FILE, __LINE__, fdf_err)
      endif
      call fdf_close()

      RETURN
!--------------------------------------------------------------------------- END
    END SUBROUTINE fdf_read_xyz



!
!   Read an input file (and include files) searching labels to
!   include them in memory structure that will contain the data
!
    RECURSIVE FUNCTION fdf_readlabel(nelem, plabel, filein, found) result(readlabel)
      implicit none
!--------------------------------------------------------------- Input Variables
      character(*)               :: filein
      integer(ip)                :: nelem
      type(parsed_line), pointer :: plabel

!-------------------------------------------------------------- Output Variables
      logical                    :: readlabel
      logical                    :: found(nelem)

!--------------------------------------------------------------- Local Variables
      logical                    :: dump, found_elem
      logical, pointer           :: found_loc(:)
      character(80)              :: msg
      character(len=MAX_LENGTH*2):: line
      character(len=MAX_LENGTH)  :: label, inc_file
      integer(ip)                :: i, ierr, ntok, ind_less, nlstart
      integer(ip)                :: elem, nelem_loc
      integer(ip), pointer       :: found_index(:)
      type(parsed_line), pointer :: pline

!------------------------------------------------------------------------- BEGIN
!     Open input file with labels
      call fdf_open(filein)

!     While not reach to end of file and found all labels
      do while (fdf_readline(line) .and. (.not. ALL(found)))

!       Check if valid data (tokens, non-blank)
        pline => digest(line)
        ntok = ntokens(pline)
        if (ntok .ne. 0) then

!         Find different special cases in the input files
!         (%block, %endblock, %include, Label1 Label2 ... < Filename)

!         %block directive
          ind_less = search('<', pline)
          if (search('%block', pline) .eq. 1) then

!           No label found in %block directive
            if (ntok .eq. 1) then
              write(msg,*) '%block label not found in ', TRIM(filein)
              call die('FDF module: fdf_readlabel', msg,                &
                       THIS_FILE, __LINE__, fdf_err)
            endif

!           %block Label < Filename [ %dump ]
            if (ind_less .eq. 3) then

              if (ntok .ge. 4) then
!               Test if %dump is present
                if (search('%dump', pline) .eq. 5) then
                  dump = .TRUE.
                else
                  dump = .FALSE.
                endif

                label = tokens(pline, 2)
                elem  = fdf_searchlabel(found, nelem, label, plabel)

                inc_file = tokens(pline, 4)
                call destroy(pline)

!               If match with any label add [begin, body, end] of block
                if (elem .ne. -1) then
                  line = '%block ' // label
                  pline => digest(line)
                  call setmorphol(1, 'b', pline)
                  call setmorphol(2, 'l', pline)
                  call fdf_addtoken(line, pline)

                  nlstart = file_in%nlines
                  call fdf_read(inc_file, label)

!                 Warn if block 'label' is empty
                  if ((nlstart - file_in%nlines) .eq. 0) then
                    write(msg,*) 'FDF module: fdf_readlabel: block ',   &
                                 TRIM(label), ' is empty...'
                    call warn(msg)
                  endif

                  line = '%endblock ' // label
                  pline => digest(line)
                  call setmorphol(1, 'e', pline)
                  call setmorphol(2, 'l', pline)
                  call fdf_addtoken(line, pline)

!                 Dump included file to fileout
                  if (dump) call fdf_dump(label)

                  found(elem) = .TRUE.
                  label = ' '
                endif

!             Filename not found in %block directive
              else
                write(msg,*) 'Filename on %block not found in ', TRIM(filein)
                call die('FDF module: fdf_readlabel', msg,              &
                         THIS_FILE, __LINE__, fdf_err)
              endif

!           %block Label
            elseif (ind_less .eq. -1) then
              label = tokens(pline, 2)
              elem  = fdf_searchlabel(found, nelem, label, plabel)
              found_elem = .TRUE.

!             If match with any label add [begin,body,end] of block
              if (elem .ne. -1) then
                call setmorphol(1, 'b', pline)
                call setmorphol(2, 'l', pline)
                call fdf_addtoken(line, pline)
                nlstart = file_in%nlines

                found_elem = .FALSE.
                do while (fdf_readline(line) .and. (.not. found_elem))
                  pline => digest(line)
                  if (ntokens(pline) .ne. 0) then
                    if (search('%endblock', pline) .eq. 1) then
!                     Warn if block 'label' is empty
                      if ((nlstart - file_in%nlines) .eq. 0) then
                        write(msg,*) 'FDF module: fdf_readlabel: block ', &
                                     TRIM(label), ' is empty...'
                        call warn(msg)
                      endif

                      call destroy(pline)
                      line = '%endblock ' // label
                      pline => digest(line)
                      call setmorphol(1, 'e', pline)
                      call setmorphol(2, 'l', pline)
                      label = ' '

                      found_elem  = .TRUE.
                      found(elem) = .TRUE.
                    endif
                    call fdf_addtoken(line, pline)
                  endif
                enddo

!             Move to the end of the block
              else
                call destroy(pline)

                found_elem = .FALSE.
                do while (fdf_readline(line) .and. (.not. found_elem))
                  pline => digest(line)
                  if (search('%endblock', pline) .eq. 1) then
                    label = ' '
                    found_elem = .TRUE.
                  endif
                  call destroy(pline)
                enddo
              endif

!             Error due to %endblock not found
              if (.not. found_elem) then
                write(msg,*) '%endblock ', TRIM(label),                 &
                             ' not found in ', TRIM(filein)
                call die('FDF module: fdf_readlabel', msg,              &
                         THIS_FILE, __LINE__, fdf_err)
              endif

!           Bad format in %block directive
            else
              write(msg,*) 'Bad ''<'' %block format in ', TRIM(filein)
              call die('FDF module: fdf_readlabel', msg,                &
                       THIS_FILE, __LINE__, fdf_err)
            endif

!         %endblock directive
          elseif (search('%endblock', pline) .eq. 1) then
!           Bad if %endblock exists before %block
            write(msg,*) 'Bad %endblock found in ', TRIM(filein)
            call die('FDF module: fdf_readlabel', msg,                  &
                     THIS_FILE, __LINE__, fdf_err)

!         %include Filename directive
          elseif (search('%include', pline) .eq. 1) then
!           Check if include filename is specified
            if (ntok .eq. 1) then
              write(msg,*) 'Filename on %include not found in ', TRIM(filein)
              call die('FDF module: fdf_readlabel', msg,                &
                       THIS_FILE, __LINE__, fdf_err)
            else
              inc_file = tokens(pline, 2)
              call destroy(pline)
              readlabel = fdf_readlabel(nelem, plabel, inc_file, found)
            endif

!         Label1 Label2 ... < Filename directive
          elseif (ind_less .ne. -1) then
!           Check if '<' is in a valid position
            if (ind_less .eq. 1) then
              write(msg,*) 'Bad ''<'' found in ', TRIM(filein)
              call die('FDF module: fdf_readlabel', msg,                &
                       THIS_FILE, __LINE__, fdf_err)

!           Check if '<' filename is specified
            elseif (ind_less .eq. ntok) then
              write(msg,*) 'Filename not found after ''<'' in ', TRIM(filein)
              call die('FDF module: fdf_readlabel', msg,                &
                       THIS_FILE, __LINE__, fdf_err)
            else
!             Search label(s) in Filename
              line = ' '
              nelem_loc = 0
              ALLOCATE(found_index(ind_less-1), stat=ierr)
              if (ierr .ne. 0) then
                call die('FDF module: fdf_readlabel', 'Error allocating found_index', &
                         THIS_FILE, __LINE__, fdf_err, rc=ierr)
              endif
              do i= 1, ind_less-1
                label = tokens(pline, i)
                elem = fdf_searchlabel(found, nelem, label, plabel)
                if (elem .ne. -1) then
                  line = TRIM(line) // ' ' // TRIM(label)
                  nelem_loc = nelem_loc + 1
                  found_index(nelem_loc) = elem
                endif
              enddo

!             Process Filename if any label found
              if (nelem_loc .ge. 1) then
                inc_file = tokens(pline, ind_less+1)
                call destroy(pline)

                ALLOCATE(found_loc(nelem_loc), stat=ierr)
                if (ierr .ne. 0) then
                  call die('FDF module: fdf_readlabel', 'Error allocating found_loc', &
                           THIS_FILE, __LINE__, fdf_err, rc=ierr)
                endif

                found_loc = .FALSE.

!               If label(s) not found in such Filename throw an error
                pline => digest(line)
                if (.not. fdf_readlabel(nelem_loc, pline,               &
                                        inc_file, found_loc)) then
                  i = 1
                  do while ((i .le. nelem_loc) .and. (found_loc(i)))
                    i = i + 1
                  enddo
                  label = tokens(pline, i)
                  write(msg,*) 'Label ', TRIM(label), ' not found in ', TRIM(inc_file)
                  call die('FDF module: fdf_readlabel', msg,            &
                           THIS_FILE, __LINE__, fdf_err)
                else
!                 Merge results if all labels found
                  do i= 1, nelem_loc
                    found(found_index(i)) = found_loc(i)
                  enddo
                endif

                DEALLOCATE(found_index)
              endif

              DEALLOCATE(found_loc)
              call destroy(pline)
            endif

!         Label [ Value ] directive
          else
            elem = fdf_searchlabel(found, nelem, tokens(pline, 1), plabel)

!           If match with any label add it
            if (elem .ne. -1) then
              call setmorphol(1, 'l', pline)
              call fdf_addtoken(line, pline)
              found(elem) = .TRUE.
            else
!             Destroy parsed_line structure if no label found
              call destroy(pline)
            endif
          endif

        else
!         Destroy parsed_line structure if no label found
          call destroy(pline)
        endif
      enddo

!     Close input file with labels
      call fdf_close()

      readlabel = ALL(found)
      RETURN
!--------------------------------------------------------------------------- END
    END FUNCTION fdf_readlabel

!
!   Search a label in a set of not 'found' tokens given by plabel and nelem.
!   Returns the index on plabel of the token that matches with label.
!
    FUNCTION fdf_searchlabel(found, nelem, label, plabel)
      implicit none
!--------------------------------------------------------------- Input Variables
      integer(ip)                :: nelem
      logical                    :: found(nelem)
      character(*)               :: label
      type(parsed_line), pointer :: plabel

!-------------------------------------------------------------- Output Variables
      integer(ip)                :: fdf_searchlabel

!--------------------------------------------------------------- Local Variables
      logical                    :: found_elem
      integer(ip)                :: i

!------------------------------------------------------------------------- BEGIN
      i = 1
      found_elem = .FALSE.
      fdf_searchlabel = -1
      do while ((i .le. nelem) .and. (.not. found_elem))

        if (.not. found(i)) then
          if (labeleq(label, tokens(plabel, i))) then
            found_elem      = .TRUE.
            fdf_searchlabel = i
          endif
        endif
        i = i + 1
      enddo

      RETURN
!--------------------------------------------------------------------------- END
    END FUNCTION fdf_searchlabel

!
!   Dumps the content of a block
!
    SUBROUTINE fdf_dump(label)
      implicit none
!--------------------------------------------------------------- Input Variables
      character(*)               :: label

!--------------------------------------------------------------- Local Variables
      character(80)              :: msg
      type(block_fdf)            :: bfdf
      type(parsed_line), pointer :: pline

!------------------------------------------------------------------------- BEGIN
      fdf_started = .TRUE.

      if (.not. fdf_block(label, bfdf)) then
        write(msg,*) 'block ', label, 'to dump not found'
        call die('FDF module: fdf_dump', msg, THIS_FILE, __LINE__, fdf_err)
      endif

!     fdf_bline prints each block line in fdf_out
      do while(fdf_bline(bfdf, pline))
      enddo

      fdf_started = .FALSE.

      RETURN
!--------------------------------------------------------------------------- END
    END SUBROUTINE fdf_dump

!
!   Init FDF file structure
!
    SUBROUTINE fdf_initdata()
      implicit none
!--------------------------------------------------------------- Local Variables
      integer(ip) :: ierr

!------------------------------------------------------------------------- BEGIN
      ndepth = 0

      ALLOCATE(file_in, stat=ierr)
      if (ierr .ne. 0) then
        call die('FDF module: fdf_initdata', 'Error allocating file_in', &
                 THIS_FILE, __LINE__, fdf_err, rc=ierr)
      endif

      file_in%nlines = 0
      NULLIFY(file_in%first)
      NULLIFY(file_in%last)

      RETURN
!--------------------------------------------------------------------------- END
    END SUBROUTINE fdf_initdata

!
!   Add a line individually to the dynamic list of parsed lines
!   This can not include block's and is restricted to key values
!
    SUBROUTINE fdf_addline(line)
      implicit none
!--------------------------------------------------------------- Input Variables
      character(len=MAX_LENGTH)  :: line

!--------------------------------------------------------------- Local Variables
      integer(ip)                :: ntok
      type(parsed_line), pointer :: pline

!------------------------------------------------------------------------- BEGIN

!     Check if valid data (tokens, non-blank)
      pline => digest(line)

      call setmorphol(1, 'l', pline)
      call fdf_addtoken(line, pline)

      if (fdf_debug2) then
         write(fdf_log,*) '***FDF_ADDLINE********************************'
         write(fdf_log,*) 'Line:', TRIM(line)
         write(fdf_log,*) '**********************************************'
      endif

    END SUBROUTINE fdf_addline

!
!   Remove a line from the dynamic list of parsed lines
!   This can not include block's and is restricted to key values
!
    SUBROUTINE fdf_removelabel(label)
      implicit none
!--------------------------------------------------------------- Input Variables
      character(len=MAX_LENGTH)  :: label

!--------------------------------------------------------------- Local Variables
      type(line_dlist), pointer  :: mark

!------------------------------------------------------------------------- BEGIN

      do while ( fdf_locate(label,mark) )

         if (fdf_debug2) then
            write(fdf_log,*) '***FDF_REMOVELABEL*******************************'
            write(fdf_log,*) 'Line:', TRIM(mark%str)
            write(fdf_log,*) 'Label:', trim(label)
            write(fdf_log,*) '**********************************************'
         endif

         ! To circumvent the first/last line in the fdf-file
         ! we have to check for the existence of the
         ! first/last mark being the one removed.
         ! That special case *must* correct the first/last
         ! tokens.
         if ( associated(mark,target=file_in%first) ) then
            file_in%first => mark%next
         end if
         if ( associated(mark,target=file_in%last) ) then
            file_in%last => mark%prev
         end if

         ! Remove the label from the dynamic list
         call destroy(mark%pline)
         if ( associated(mark%prev) ) then
            mark%prev%next => mark%next
         end if
         if ( associated(mark%next) ) then
            mark%next%prev => mark%prev
         end if
         DEALLOCATE(mark)

         NULLIFY(mark)
      end do

    END SUBROUTINE fdf_removelabel

!
!   Overwrite label line in dynamic list of parsed lines
!
    SUBROUTINE fdf_overwrite(line)
!--------------------------------------------------------------- Input Variables
      character(len=MAX_LENGTH)   :: line

!--------------------------------------------------------------- Local Variables
      type(parsed_line), pointer  :: pline
      character(len=MAX_LENGTH)   :: label

      integer :: ierr

      pline => digest(line)
      if ( search('%block', pline) == 1 .or. &
          search('%endblock', pline) == 1 ) then

        ! We do not allow this in a single line
        call die('FDF module: fdf_overwrite', 'Error overwriting block (not implemented)',   &
            THIS_FILE, __LINE__, fdf_err, rc=ierr)

      else if ( search('%include', pline) == 1 ) then

        ! We do not allow this in a single line
        call die('FDF module: fdf_overwrite', 'Error overwriting flags from input file (not implemented)',   &
            THIS_FILE, __LINE__, fdf_err, rc=ierr)

      else if ( search('<', pline) /= -1 ) then

        ! We do not allow this in a single line
        call die('FDF module: fdf_overwrite', 'Error piping in overwriting (not implemented)',   &
            THIS_FILE, __LINE__, fdf_err, rc=ierr)

      else

      label = tokens(pline,1)
      call setmorphol(1, 'l', pline)
      call fdf_removelabel(label)

        ! Add token to the list of fdf-flags
        ! Since we add it directly we shouldn't destroy the pline
      call fdf_addtoken(line, pline)
      if ( fdf_debug ) then
        write(fdf_log,'(2a)') '---> Overwriting token: ', trim(label)
      end if

    end if

    END SUBROUTINE fdf_overwrite

!
!   Add a token to the dynamic list of parsed lines
!
    SUBROUTINE fdf_addtoken(line, pline)
      implicit none
!--------------------------------------------------------------- Input Variables
      character(len=MAX_LENGTH)  :: line
      type(parsed_line), pointer :: pline

!--------------------------------------------------------------- Local Variables
      integer(ip)                :: i, ierr
      type(line_dlist), pointer  :: mark

!------------------------------------------------------------------------- BEGIN
      ALLOCATE(mark, stat=ierr)
      if (ierr .ne. 0) then
        call die('FDF module: fdf_addtoken', 'Error allocating mark',   &
                 THIS_FILE, __LINE__, fdf_err, rc=ierr)
      endif

      mark%str   =  line
      mark%pline => pline
      NULLIFY(mark%next)

      if (ASSOCIATED(file_in%first)) then
        mark%prev         => file_in%last
        file_in%last%next => mark
      else
        NULLIFY(mark%prev)
        file_in%first => mark
      endif

      file_in%last => mark
      file_in%nlines = file_in%nlines + 1

      if (fdf_debug2) then
        write(fdf_log,*) '***FDF_ADDTOKEN*******************************'
        write(fdf_log,*) 'Line:', TRIM(mark%str)
        write(fdf_log,*) 'Ntokens:', mark%pline%ntokens
        do i= 1, mark%pline%ntokens
          write(fdf_log,*) '  Token:', trim(tokens(pline,i)), &
                           ' (', mark%pline%id(i), ')'
        enddo
        write(fdf_log,*) '**********************************************'
      endif

      RETURN
!--------------------------------------------------------------------------- END
    END SUBROUTINE fdf_addtoken

!
!   Opens a file for FDF processing.
!
    SUBROUTINE fdf_open(filename)
      implicit none
!-------------------------------------------------------------- Output Variables
      character(*)  :: filename

!--------------------------------------------------------------- Local Variables
      logical       :: file_exists
      character(80) :: msg
      integer(ip)   :: lun

!------------------------------------------------------------------------- BEGIN
      ndepth = ndepth + 1
      if (ndepth .gt. maxdepth) then
        call die('FDF module: fdf_open', 'Too many nested fdf files...', &
                 THIS_FILE, __LINE__, fdf_err)
      endif

      if (leqi(filename, 'stdin')) then
        lun = INPUT_UNIT
        if (fdf_debug) write(fdf_log,'(a,i1,a)')                        &
          '---> Reading from standard input [DEPTH:', ndepth,'] '
      else
        call io_assign(lun)

        INQUIRE(file=filename, exist=file_exists)
        if (file_exists) then
          open(unit=lun, file=filename, status='old', form='formatted')
          REWIND(lun)
          if (fdf_debug) write(fdf_log,'(a,i1,a,a)')                    &
            '---> Opened [DEPTH:', ndepth,'] ', TRIM(filename)
        else
          write(msg,'(a,a)') 'Cannot open ', TRIM(filename)
          call die('FDF module: fdf_open', msg, THIS_FILE, __LINE__, fdf_err)
        endif
      endif

      fdf_in(ndepth) = lun
      REWIND(fdf_in(ndepth))

      RETURN
!--------------------------------------------------------------------------- END
    END SUBROUTINE fdf_open

!
!   Closes currently opened fdf file
!
    SUBROUTINE fdf_close()
      implicit none
!------------------------------------------------------------------------- BEGIN
      if (ndepth .ge. 1) then
        call io_close(fdf_in(ndepth))
        if (fdf_debug)                                                  &
          write(fdf_log,'(a,i1,a)') '---> Closed [DEPTH:', ndepth,']'
        ndepth = ndepth - 1
      endif

      RETURN
!--------------------------------------------------------------------------- END
    END SUBROUTINE fdf_close

!
!   Set output file for FDF subsystem
!
    SUBROUTINE fdf_set_output_file(fileout)
      implicit none
!----------------------------------------------------- Input Variables
      character(len=*), intent(in)   :: fileout

!----------------------------------------------------- Local Variables
      character(256) :: fileouttmp
!----------------------------------------------------- BEGIN
      call io_assign(fdf_out)

#ifdef _MPI_
      if (rank /= 0) then
         fileouttmp =  trim(fileout) // "." // i2s(rank)
      else
         fileouttmp = fileout
      endif
#else
      fileouttmp = fileout
#endif

#ifdef FDF_DEBUG
      !
      !     If debugging, all the nodes use named log files
      !
      open( unit=fdf_out, file=TRIM(fileouttmp), form='formatted', &
           access='sequential', status='replace' )
#else
      !
      !     Only the master node opens a named log file in the current dir.
      !     Non-master nodes use a scratch file.
      !     These log files tend to be quite small, so there
      !     should not be problems such as filling up filesystems.
      !     ... your mileage might vary. This is a grey area.
      !     Some compilers allow the user to specify where scratch
      !     files go.
      !
#ifdef _MPI_
      if (rank /= 0) then
         open( unit=fdf_out, form='formatted', &
              access='sequential', status='scratch' )
      else
         open( unit=fdf_out, file=TRIM(fileouttmp), form='formatted', &
              access='sequential', status='replace' )
      end if
#else
      open( unit=fdf_out, file=TRIM(fileouttmp), form='formatted', &
           access='sequential', status='replace' )
#endif
#endif

      RETURN
!--------------------------------------------------------------------------- END
    END SUBROUTINE fdf_set_output_file

!
!   Frees and shutdown FDF system
!
    SUBROUTINE fdf_shutdown()
      implicit none
!------------------------------------------------------------------------- BEGIN
!$OMP SINGLE
      if (fdf_started) then
        call fdf_destroy(file_in)
        fdf_started = .FALSE.

        call io_close(fdf_out)
      endif
!$OMP END SINGLE

      RETURN
!--------------------------------------------------------------------------- END
    END SUBROUTINE fdf_shutdown

!
!   Destroy the fdf_file structure for the input file
!
    SUBROUTINE fdf_destroy(fdfp)
      implicit none
!-------------------------------------------------------------- Output Variables
      type(fdf_file) :: fdfp

!------------------------------------------------------------------------- BEGIN
      if (ASSOCIATED(fdfp%first)) call fdf_destroy_dl(fdfp%first)

      RETURN
!--------------------------------------------------------------------------- END
    END SUBROUTINE fdf_destroy

!
!   Destroy recursively the dynamic list of parsed lines
!
    RECURSIVE SUBROUTINE fdf_destroy_dl(dlp)
      implicit none
!-------------------------------------------------------------- Output Variables
      type(line_dlist), pointer :: dlp

      !! Use for tail recursion later:      type(line_dlist), pointer :: pnext

      !------------------------------------------------------------------------- BEGIN
      ! This is NOT tail-recursive!!
      if (ASSOCIATED(dlp%next)) call fdf_destroy_dl(dlp%next)
      call destroy(dlp%pline)
      DEALLOCATE(dlp)

      RETURN
!--------------------------------------------------------------------------- END
    END SUBROUTINE fdf_destroy_dl

!
!   Tells when FDF library has been compiled with parallel
!   execution support. Depends on _MPI_ macro.
!
    FUNCTION fdf_parallel()
      implicit none
!-------------------------------------------------------------- Output Variables
      logical      :: fdf_parallel

!------------------------------------------------------------------------- BEGIN
#ifdef _MPI_
      fdf_parallel = .TRUE.
#else
      fdf_parallel = .FALSE.
#endif

      RETURN
!--------------------------------------------------------------------------- END
    END FUNCTION fdf_parallel

!
!   Read a line of the 'ndepth' input file, returning .TRUE. if
!   there are more lines to read from input file, .FALSE. otherwise.
!
    FUNCTION fdf_readline(line)
      implicit none
!-------------------------------------------------------------- Output Variables
      logical      :: fdf_readline
      character(*) :: line

!--------------------------------------------------------------- Local Variables
      integer(ip)  :: stat

!------------------------------------------------------------------------- BEGIN
      read(fdf_in(ndepth), '(a)', iostat=stat) line

      if (stat .eq. 0) then
        fdf_readline = .TRUE.
        if (fdf_debug2) write(fdf_log, '(a,a76)') 'fdf_readline > ', line
      else
        fdf_readline = .FALSE.
      endif

      RETURN
!--------------------------------------------------------------------------- END
    END FUNCTION fdf_readline

!
!   Returns the string line associated with the mark pointer
!   in the FDF herarchy. If mark is not associated returns ''.
!
    FUNCTION fdf_getline(mark)
      implicit none
!--------------------------------------------------------------- Input Variables
      type(line_dlist), pointer :: mark

!-------------------------------------------------------------- Output Variables
      character(len=MAX_LENGTH) :: fdf_getline

!------------------------------------------------------------------------- BEGIN
      if (ASSOCIATED(mark)) then
        fdf_getline = mark%str
      else
        fdf_getline = ' '
      endif

      RETURN
!--------------------------------------------------------------------------- END
    END FUNCTION fdf_getline

!
!   Prints all the fdf_file structure of the input file(s)
!
    SUBROUTINE fdf_printfdf()
      implicit none
!--------------------------------------------------------------- Local Variables
      integer(ip)               :: i, ntokens
      character*1               :: id
      type(line_dlist), pointer :: dlp
      character(len=MAX_LENGTH) :: tok

!------------------------------------------------------------------------- BEGIN
      dlp => file_in%first

      write(fdf_log,*) '*** FDF Memory Structure Summary: ************'
      do while (ASSOCIATED(dlp))
        ntokens = dlp%pline%ntokens
        write(fdf_log,*) 'Line:', TRIM(dlp%str)
        write(fdf_log,*) 'Ntokens:', ntokens
        do i= 1, ntokens
          tok = tokens(dlp%pline,i)
          id  = dlp%pline%id(i)
          write(fdf_log,*) '  Token:', trim(tok), '(', dlp%pline%id(i), ')'
        enddo
        dlp => dlp%next
      enddo
      write(fdf_log,*) '**********************************************'

      RETURN
!--------------------------------------------------------------------------- END
    END SUBROUTINE fdf_printfdf

!
!   Returns an integer associated with label 'label', or the default
!   value if label is not found in the fdf file.
!   Optionally can return a pointer to the line found.
!
    FUNCTION fdf_integer(label, default, line)
      implicit none
!--------------------------------------------------------------- Input Variables
      character(*)                        :: label
      integer(ip)                         :: default

!-------------------------------------------------------------- Output Variables
      integer(ip)                         :: fdf_integer
      type(line_dlist), pointer, optional :: line

!--------------------------------------------------------------- Local Variables
      character(80)                       :: msg
      type(line_dlist), pointer           :: mark

!------------------------------------------------------------------------- BEGIN
!     Prevents using FDF routines without initialize
      if (.not. fdf_started) then
        call die('FDF module: fdf_integer', 'FDF subsystem not initialized', &
                 THIS_FILE, __LINE__, fdf_err)
      endif

      if (fdf_locate(label, mark)) then
        if (.not. match(mark%pline, 'li')) then
          write(msg,*) 'no integer value for ', label
          call die('FDF module: fdf_integer', msg, THIS_FILE, __LINE__, fdf_err)
        endif

        fdf_integer = integers(mark%pline, 1, 1)
        if (fdf_output) write(fdf_out,'(a,t30,i0)') label, fdf_integer
      else
        fdf_integer = default
        if (fdf_output) write(fdf_out,'(a,t30,i0,t60,a)') label, default, '# default value'
      endif

      if (PRESENT(line)) line = mark

      RETURN
!--------------------------------------------------------------------------- END
    END FUNCTION fdf_integer

!
!   Returns true or false whether or not the label 'label' is
!   a value with units or not.
!   I.e. it returns true if the line has the form lvn, if not found, or not lvn,
!   it returns false.
!
    FUNCTION fdf_isphysical(label)
      implicit none
!--------------------------------------------------------------- Input Variables
      character(*)                        :: label

!-------------------------------------------------------------- Output Variables
      logical                             :: fdf_isphysical

!--------------------------------------------------------------- Local Variables
      type(line_dlist), pointer           :: mark

!------------------------------------------------------------------------- BEGIN
!     Prevents using FDF routines without initialize
      if (.not. fdf_started) then
         call die('FDF module: fdf_isphysical', 'FDF subsystem not initialized', &
                 THIS_FILE, __LINE__, fdf_err)
      endif

      if (fdf_locate(label, mark)) then
         fdf_isphysical = match(mark%pline, 'lvn')
      else
         fdf_isphysical = .false.
      endif
      if (fdf_output) write(fdf_out,'(a,t30,l10)') "#:physical? " // label, fdf_isphysical

      RETURN
!--------------------------------------------------------------------------- END
    END FUNCTION fdf_isphysical

!
!   Returns true or false whether or not the label 'label' is
!   a list or not, you cannot get the line out from this routine
!
    FUNCTION fdf_islist(label)
      implicit none
!--------------------------------------------------------------- Input Variables
      character(*)                        :: label

!-------------------------------------------------------------- Output Variables
      logical                             :: fdf_islist

!--------------------------------------------------------------- Local Variables
      type(line_dlist), pointer           :: mark

!------------------------------------------------------------------------- BEGIN
!     Prevents using FDF routines without initialize
      if (.not. fdf_started) then
         call die('FDF module: fdf_islist', 'FDF subsystem not initialized', &
                 THIS_FILE, __LINE__, fdf_err)
      endif

      if (fdf_locate(label, mark)) then
        ! if it is a list:
        fdf_islist = match(mark%pline, 'le')
      else
         fdf_islist = .false.
      endif
      if (fdf_output) write(fdf_out,'(a,t30,l10)') "#:list? " // label, fdf_islist

      RETURN
!--------------------------------------------------------------------------- END
    END FUNCTION fdf_islist

    FUNCTION fdf_islinteger(label)
      implicit none
!--------------------------------------------------------------- Input Variables
      character(*)                        :: label

!-------------------------------------------------------------- Output Variables
      logical                             :: fdf_islinteger

!--------------------------------------------------------------- Local Variables
      type(line_dlist), pointer           :: mark

!------------------------------------------------------------------------- BEGIN
!     Prevents using FDF routines without initialize
      if (.not. fdf_started) then
         call die('FDF module: fdf_islinteger', 'FDF subsystem not initialized', &
                 THIS_FILE, __LINE__, fdf_err)
      endif

      if (fdf_locate(label, mark)) then
        ! if it is an integer list:
        fdf_islinteger = match(mark%pline, 'la')
      else
         fdf_islinteger = .false.
      endif
      if (fdf_output) write(fdf_out,'(a,t30,l10)') "#:linteger? " // label, &
          fdf_islinteger

      RETURN
!--------------------------------------------------------------------------- END
    END FUNCTION fdf_islinteger

    FUNCTION fdf_islreal(label)
      implicit none
!--------------------------------------------------------------- Input Variables
      character(*)                        :: label

!-------------------------------------------------------------- Output Variables
      logical                             :: fdf_islreal

!--------------------------------------------------------------- Local Variables
      type(line_dlist), pointer           :: mark

!------------------------------------------------------------------------- BEGIN
!     Prevents using FDF routines without initialize
      if (.not. fdf_started) then
         call die('FDF module: fdf_islreal', 'FDF subsystem not initialized', &
                 THIS_FILE, __LINE__, fdf_err)
      endif

      if (fdf_locate(label, mark)) then
        ! if it is a reallist:
        fdf_islreal = match(mark%pline, 'lc')
      else
         fdf_islreal = .false.
      endif
      if (fdf_output) write(fdf_out,'(a,t30,l10)') "#:lreal? " // label, &
          fdf_islreal

      RETURN
!--------------------------------------------------------------------------- END
    END FUNCTION fdf_islreal

!
!   Returns a list with label 'label', or the default
!   value if label is not found in the fdf file.
!
    SUBROUTINE fdf_linteger(label,ni,list,line)
      implicit none
!--------------------------------------------------------------- Input Variables
      character(*)                        :: label
      integer(ip)                         :: ni

!-------------------------------------------------------------- Output Variables
      integer(ip)                         :: list(ni)
      type(line_dlist), pointer, optional :: line

!--------------------------------------------------------------- Local Variables
      character(80)                       :: msg
      type(line_dlist), pointer           :: mark
      integer(ip)                         :: lni, llist(1)

!------------------------------------------------------------------------- BEGIN
!     Prevents using FDF routines without initialize
      if (.not. fdf_started) then
         call die('FDF module: fdf_linteger', 'FDF subsystem not initialized', &
              THIS_FILE, __LINE__, fdf_err)
      endif

      if (fdf_locate(label, mark)) then
         if (.not. match(mark%pline, 'la')) then
            write(msg,*) 'no list value for ', label
            call die('FDF module: fdf_linteger', msg, THIS_FILE, __LINE__, fdf_err)
         endif

         ! Retrieve length of list
         lni = -1
         call integerlists(mark%pline,1,lni,llist)
         if ( ni <= 0 ) then
            ! the user has requested size...
            ni = lni
         else
            ! the list is not long enough
            if ( ni < lni ) then
              write(msg, '(2a,2(a,i0))')'List ', trim(label), &
                  ' container too small: ', ni, ' versus ', lni
              call die('FDF module: fdf_linteger', trim(msg), &
                  THIS_FILE, __LINE__, fdf_err)
            end if
            call integerlists(mark%pline,1,ni,list)
         end if

         ! find a way to write out the list anyway
         if (fdf_output) write(fdf_out,'(a,t30,i10)') label, lni
      else
         write(msg,*) 'no list value for ', label
         call die('FDF module: fdf_linteger', msg, THIS_FILE, __LINE__, fdf_err)
      endif

      if (PRESENT(line)) line = mark

      RETURN
!--------------------------------------------------------------------------- END
    END SUBROUTINE fdf_linteger

!
!   Returns a list with label 'label', or the default
!   value if label is not found in the fdf file.
!
    SUBROUTINE fdf_ldouble(label,nv,list,line)
      implicit none
!--------------------------------------------------------------- Input Variables
      character(*)                        :: label
      integer(ip)                         :: nv

!-------------------------------------------------------------- Output Variables
      real(dp)                            :: list(nv)
      type(line_dlist), pointer, optional :: line

!--------------------------------------------------------------- Local Variables
      character(80)                       :: msg
      type(line_dlist), pointer           :: mark
      integer(ip)                         :: lnv
      real(dp)                            :: llist(1)

!------------------------------------------------------------------------- BEGIN
!     Prevents using FDF routines without initialize
      if (.not. fdf_started) then
         call die('FDF module: fdf_ldouble', 'FDF subsystem not initialized', &
              THIS_FILE, __LINE__, fdf_err)
      endif

      if (fdf_locate(label, mark)) then
         if (.not. match(mark%pline, 'le')) then
            write(msg,*) 'no list value for ', label
            call die('FDF module: fdf_ldouble', msg, THIS_FILE, __LINE__, fdf_err)
         endif

         ! Retrieve length of list
         lnv = -1
         call valuelists(mark%pline,1,lnv,llist)
         if ( nv <= 0 ) then
            ! the user has requested size...
            nv = lnv
         else
            ! the list is not long enough
            if ( nv < lnv ) then
              write(msg, '(2a,2(a,i0))')'List ', trim(label), &
                  ' container too small: ', nv, ' versus ', lnv
              call die('FDF module: fdf_ldouble', trim(msg), &
                  THIS_FILE, __LINE__, fdf_err)
            end if
            call valuelists(mark%pline,1,nv,list)
         end if

         if (fdf_output) write(fdf_out,'(a,t30,i10)') label, lnv
      else
         write(msg,*) 'no list value for ', label
         call die('FDF module: fdf_ldouble', msg, THIS_FILE, __LINE__, fdf_err)
      endif

      if (PRESENT(line)) line = mark

      RETURN
!--------------------------------------------------------------------------- END
    END SUBROUTINE fdf_ldouble

!
!   Returns a string associated with label 'label', or the default
!   string if label is not found in the fdf file.
!   Optionally can return a pointer to the line found.
!
    FUNCTION fdf_string(label, default, line)
      implicit none
!--------------------------------------------------------------- Input Variables
      character(*)                        :: label
      character(*)                        :: default

!-------------------------------------------------------------- Output Variables
      character(80)                       :: fdf_string
      type(line_dlist), pointer, optional :: line

!--------------------------------------------------------------- Local Variables
      type(line_dlist), pointer           :: mark

!------------------------------------------------------------------------- BEGIN
!     Prevents using FDF routines without initialize
      if (.not. fdf_started) then
        call die('FDF module: fdf_string', 'FDF subsystem not initialized', &
                 THIS_FILE, __LINE__, fdf_err)
      endif

      if (fdf_locate(label, mark)) then
         if (ntokens(mark%pline) < 2) then
            fdf_string = ""
            if (fdf_output) write(fdf_out,'(a,t30,a)') label, &
             "#  *** Set to empty string *** "
         else
            ! Get all the characters spanning the space from the second to
            ! the last token
            fdf_string = characters(mark%pline, ind_init=2, ind_final=-1)
            if (fdf_output) write(fdf_out,'(a,t30,a)') label, fdf_string
         endif
      else
        fdf_string = default
        if (fdf_output) write(fdf_out,'(a,t30,a,t60,a)') label, default, '# default value'
      endif

      if (PRESENT(line)) line = mark

      RETURN
!--------------------------------------------------------------------------- END
    END FUNCTION fdf_string

!
!   Returns the name of the file appearing after the label , or the default
!   string if label is not found in the fdf file.
!   Optionally can return a pointer to the line found.
!
    FUNCTION fdf_load_filename(label, default, line)
      implicit none
!--------------------------------------------------------------- Input Variables
      character(*)                        :: label
      character(*)                        :: default

!-------------------------------------------------------------- Output Variables
      character(80)                       :: fdf_load_filename
      type(line_dlist), pointer, optional :: line

!--------------------------------------------------------------- Local Variables
      type(line_dlist), pointer           :: mark

!------------------------------------------------------------------------- BEGIN
!     Prevents using FDF routines without initialize
      if (.not. fdf_started) then
        call die('FDF module: fdf_load_filename', 'FDF subsystem not initialized', &
                 THIS_FILE, __LINE__, fdf_err)
      endif

      if (fdf_load_locate(label, mark)) then
        if (tokens(mark%pline, 1) == "load" ) then
          if (ntokens(mark%pline) < 3) then
              fdf_load_filename = ""
              if (fdf_output) write(fdf_out,'(a,t30,a)') label, &
              "#  *** Set to empty string *** "
          else
              ! Get all the characters spanning the space from the second to
              ! the last token
              fdf_load_filename = trim(characters(mark%pline, ind_init=3, ind_final=-1))
              if (fdf_output) write(fdf_out,'(a,t30,a)') label, fdf_load_filename
          endif
        else
          call die('FDF module: fdf_load_filename', 'Incorrect load statement', THIS_FILE, __LINE__, fdf_err)
        endif
      else
        fdf_load_filename = default
        if (fdf_output) write(fdf_out,'(a,t30,a,t60,a)') label, default, '# default value'
      endif

      if (PRESENT(line)) line = mark

      RETURN
!--------------------------------------------------------------------------- END
    END FUNCTION fdf_load_filename


!   Returns the name of the module appearing after the %module syntax , or the default
!   module name if label is not found in the fdf file.
!   Optionally can return a pointer to the line found.
!
    FUNCTION fdf_locate_module(label, default, line)
      implicit none
!--------------------------------------------------------------- Input Variables
      character(*)                        :: label
      character(*)                        :: default

!-------------------------------------------------------------- Output Variables
      character(80)                       :: fdf_locate_module
      type(line_dlist), pointer, optional :: line

!--------------------------------------------------------------- Local Variables
      type(line_dlist), pointer           :: mark

!------------------------------------------------------------------------- BEGIN
!     Prevents using FDF routines without initialize
      if (.not. fdf_started) then
        call die('FDF module: fdf_locate_module', 'FDF subsystem not initialized', &
                 THIS_FILE, __LINE__, fdf_err)
      endif

      if (fdf_load_locate(label, mark)) then
        if (tokens(mark%pline, 1) == "load" ) then
          if (ntokens(mark%pline) < 3) then
              fdf_locate_module = ""
              if (fdf_output) write(fdf_out,'(a,t30,a)') label, &
              "#  *** Set to empty string *** "
          else
              ! Get all the characters spanning the space from the second to
              ! the last token
              fdf_locate_module = trim(characters(mark%pline, ind_init=3, ind_final=-1))
              if (fdf_output) write(fdf_out,'(a,t30,a)') label, fdf_locate_module
          endif
        else
          call die('FDF module: fdf_locate_module', 'Incorrect load statement', THIS_FILE, __LINE__, fdf_err)
        endif
      else
        fdf_locate_module = default
        if (fdf_output) write(fdf_out,'(a,t30,a,t60,a)') label, default, '# default value'
      endif

      if (PRESENT(line)) line = mark

      RETURN
!--------------------------------------------------------------------------- END
    END FUNCTION fdf_locate_module




!
!   Returns true if label 'label' appears by itself or in the form
!   label {yes,true,.true.,t,y,1} (case insensitive).
!
!   Returns false if label 'label' appears in the form
!   label {no,false,.false.,f,n,0} (case insensitive).
!
!   If label is not found in the fdf file, fdf_boolean returns the
!   LOGICAL variable default.
!
!   Optionally can return a pointer to the line found.
!
    FUNCTION fdf_boolean(label, default, line)
      implicit none
!--------------------------------------------------------------- Input Variables
      character(*)                        :: label
      logical                             :: default

!-------------------------------------------------------------- Output Variables
      logical                             :: fdf_boolean
      type(line_dlist), pointer, optional :: line

!--------------------------------------------------------------- Local Variables
      character(80)                       :: msg, valstr
      integer                             :: valint
      type(line_dlist), pointer           :: mark


!------------------------------------------------------------------------- BEGIN
!     Prevents using FDF routines without initialize
      if (.not. fdf_started) then
        call die('FDF module: fdf_boolean', 'FDF subsystem not initialized', &
                 THIS_FILE, __LINE__, fdf_err)
      endif

      if (fdf_locate(label, mark)) then

!       If the label appears by itself, we interpret it as .true.
        if (ntokens(mark%pline) .ne. 1) then

!         Look for second word
          if (nintegers(mark%pline) == 0) then
            valstr = names(mark%pline, 1, 1)

            if (is_true(valstr)) then
              fdf_boolean = .TRUE.
              if (fdf_output) write(fdf_out,'(a,t30,l10)') label, fdf_boolean

            elseif (is_false(valstr)) then
              fdf_boolean = .FALSE.
              if (fdf_output) write(fdf_out,'(a,t30,l10)') label, fdf_boolean

            else
              write(msg,*) 'unexpected logical value ', label, ' = ', valstr
              call die('FDF module: fdf_boolean', msg,                    &
                      THIS_FILE, __LINE__, fdf_err)
            endif
          elseif (nintegers(mark%pline) == 1) then
            valint = integers(mark%pline, 1, 1)

            if ( valint == 1 ) then
              fdf_boolean = .TRUE.
              if (fdf_output) write(fdf_out,'(a,t30,l10)') label, fdf_boolean

            elseif (valint == 0) then
              fdf_boolean = .FALSE.
              if (fdf_output) write(fdf_out,'(a,t30,l10)') label, fdf_boolean

            else
              write(msg,*) 'unexpected logical value ', label, ' = ', valint
              call die('FDF module: fdf_boolean', msg,                    &
                      THIS_FILE, __LINE__, fdf_err)
            endif
          else
            fdf_boolean = .TRUE.
            if (fdf_output) write(fdf_out,'(a,t30,l10,5x,a)') label, fdf_boolean,          &
                                           '# label by itself'
          endif
        endif
      else
        fdf_boolean = default
        if (fdf_output) write(fdf_out,'(a,t30,l10,t60,a)') label, default, '# default value'
      endif

      if (PRESENT(line)) line = mark

      RETURN

      CONTAINS

      logical function is_true(valstr)  result(a)
      character(len=*), intent(in) :: valstr
      a = leqi(valstr, 'yes')    .or. leqi(valstr, 'true') .or. &
          leqi(valstr, '.true.') .or. leqi(valstr, 't')    .or. &
          leqi(valstr, 'y')
      end function is_true

      logical function is_false(valstr)  result(a)
      character(len=*), intent(in) :: valstr
      a = leqi(valstr, 'no')      .or. leqi(valstr, 'false') .or. &
          leqi(valstr, '.false.') .or. leqi(valstr, 'f')     .or. &
          leqi(valstr, 'n')
      end function is_false

!--------------------------------------------------------------------------- END
    END FUNCTION fdf_boolean

!
!   Returns true if label 'label' appears by itself or in the form
!   label {yes,true,.true.,t,y} (case insensitive).
!
!   Returns false if label 'label' appears in the form
!   label {no,false,.false.,f,n} (case insensitive).
!
!   If label is not found in the fdf file, fdf_boolean returns the
!   LOGICAL variable default.
!
!   Optionally can return a pointer to the line found.
!
    FUNCTION fdf_bboolean(pline,ind,after)
      implicit none
!--------------------------------------------------------------- Input Variables
      integer(ip), intent(in)           :: ind
      integer(ip), intent(in), optional :: after
      type(parsed_line), pointer        :: pline

!-------------------------------------------------------------- Output Variables
      logical                             :: fdf_bboolean

!--------------------------------------------------------------- Local Variables
      character(80)                       :: msg, valstr
      type(line_dlist), pointer           :: mark


!------------------------------------------------------------------------- BEGIN
!     Prevents using FDF routines without initialize
      if (.not. fdf_started) then
        call die('FDF module: fdf_bboolean', 'FDF subsystem not initialized', &
                 THIS_FILE, __LINE__, fdf_err)
      endif

      if (ind <= nnames(pline,after=after)) then

        valstr = names(pline,ind,after=after)

        if (is_true(valstr)) then
           fdf_bboolean = .TRUE.
           if (fdf_output) write(fdf_out,'(a,t30,l10)') valstr, fdf_bboolean

        elseif (is_false(valstr)) then
           fdf_bboolean = .FALSE.
           if (fdf_output) write(fdf_out,'(a,t30,l10)') valstr, fdf_bboolean

        else
           write(msg,*) 'unexpected logical value ', valstr
           call die('FDF module: fdf_bboolean', msg,                    &
                THIS_FILE, __LINE__, fdf_err)
        endif
      else
         fdf_bboolean = .TRUE.
         if (fdf_output) write(fdf_out,'(l10,5x,a)') fdf_bboolean, &
                                       '# block designation by itself'
      endif

      RETURN

    CONTAINS

      logical function is_true(valstr)  result(a)
      character(len=*), intent(in) :: valstr
      a = leqi(valstr, 'yes')    .or. leqi(valstr, 'true') .or. &
          leqi(valstr, '.true.') .or. leqi(valstr, 't')    .or. &
          leqi(valstr, 'y')
      end function is_true

      logical function is_false(valstr)  result(a)
      character(len=*), intent(in) :: valstr
      a = leqi(valstr, 'no')      .or. leqi(valstr, 'false') .or. &
          leqi(valstr, '.false.') .or. leqi(valstr, 'f')     .or. &
          leqi(valstr, 'n')
      end function is_false

!--------------------------------------------------------------------------- END
    END FUNCTION fdf_bboolean

!
!   Returns a single precision value associated with label 'label',
!   or the default value if label is not found in the fdf file.
!   Optionally can return a pointer to the line found.
!   Note that integers on the line are also accepted
!
    FUNCTION fdf_single(label, default, line)
      implicit none
!--------------------------------------------------------------- Input Variables
      character(*)                        :: label
      real(sp)                            :: default

!-------------------------------------------------------------- Output Variables
      real(sp)                            :: fdf_single
      type(line_dlist), pointer, optional :: line

!--------------------------------------------------------------- Local Variables
      character(80)                       :: msg
      type(line_dlist), pointer           :: mark

!------------------------------------------------------------------------- BEGIN
!     Prevents using FDF routines without initialize
      if (.not. fdf_started) then
        call die('FDF module: fdf_single', 'FDF subsystem not initialized', &
                 THIS_FILE, __LINE__, fdf_err)
      endif

      if (fdf_locate(label, mark)) then
        if (.not. match(mark%pline, 'lv')) then
          write(msg,*) 'no real value for ', label
          call die('FDF module: fdf_single', msg, THIS_FILE, __LINE__,  fdf_err)
        endif
        fdf_single = values(mark%pline, 1, 1)
        if (fdf_output) write(fdf_out,'(a,t30,g20.10)') label, fdf_single
      else
        fdf_single = default
        if (fdf_output) write(fdf_out,'(a,t30,g20.10,t60,a)') label, default, '# default value'
      endif

      if (PRESENT(line)) line = mark

      RETURN
!--------------------------------------------------------------------------- END
    END FUNCTION fdf_single

!
!   Returns a double precision value associated with label 'label',
!   or the default value if label is not found in the fdf file.
!   Optionally can return a pointer to the line found.
!   Note that integers on the line are also accepted
!
    FUNCTION fdf_double(label, default, line)
      implicit none
!--------------------------------------------------------------- Input Variables
      character(*)                        :: label
      real(dp)                            :: default

!-------------------------------------------------------------- Output Variables
      real(dp)                            :: fdf_double
      type(line_dlist), pointer, optional :: line

!--------------------------------------------------------------- Local Variables
      character(80)                       :: msg
      type(line_dlist), pointer           :: mark

!------------------------------------------------------------------------- BEGIN
!     Prevents using FDF routines without initialize
      if (.not. fdf_started) then
        call die('FDF module: fdf_double', 'FDF subsystem not initialized', &
                 THIS_FILE, __LINE__, fdf_err)
      endif

      if (fdf_locate(label, mark)) then
        if (.not. match(mark%pline, 'lv')) then
          write(msg,*) 'no real value for ', label
          call die('FDF module: fdf_double', msg, THIS_FILE, __LINE__, fdf_err)
        endif
        fdf_double = values(mark%pline, 1, 1)
        if (fdf_output) write(fdf_out,'(a,t30,g20.10)') label, fdf_double
      else
        fdf_double = default
        if (fdf_output) write(fdf_out,'(a,t30,g20.10,t60,a)') label, default, '# default value'
      endif

      if (PRESENT(line)) line = mark

      RETURN
!--------------------------------------------------------------------------- END
    END FUNCTION fdf_double

!
!   Returns a double precision value associated with label 'label',
!   or the default value if label is not found in the fdf file.
!   Converts the units to defunit.
!   Optionally can return a pointer to the line found.
!
    FUNCTION fdf_physical(label, default, defunit, line)
      implicit none
!--------------------------------------------------------------- Input Variables
      character(*)                        :: label, defunit
      real(dp)                            :: default

!-------------------------------------------------------------- Output Variables
      real(dp)                            :: fdf_physical
      type(line_dlist), pointer, optional :: line

!--------------------------------------------------------------- Local Variables
      character(10)                       :: unitstr
      character(80)                       :: msg
      real(dp)                            :: value
      type(line_dlist), pointer           :: mark

!------------------------------------------------------------------------- BEGIN
!     Prevents using FDF routines without initialize
      if (.not. fdf_started) then
        call die('FDF module: fdf_physical', 'FDF subsystem not initialized', &
                 THIS_FILE, __LINE__, fdf_err)
      endif

!     Label found
      if (fdf_locate(label, mark)) then
        if (.not. match(mark%pline, 'lv')) then
          write(msg,*) 'no real value for ', label
          call die('FDF module: fdf_physical', msg, THIS_FILE, __LINE__, fdf_err)
        endif

!       Label with value
        value = values(mark%pline, 1, 1)
        fdf_physical = value

!       Look for unit
        if (.not. match(mark%pline, 'lvn')) then
          write(msg,*) 'no unit specified for ', label
          call die('FDF module: fdf_physical', msg, THIS_FILE, __LINE__, fdf_err)
        endif

        unitstr = names(mark%pline, 1, 2)
        if (.not. leqi(unitstr, defunit))                               &
          fdf_physical = value * fdf_convfac(unitstr, defunit)

        if (fdf_output) write(fdf_out,'(a,t30,g20.10,1x,a10)') label, fdf_physical, defunit
        if (fdf_output) write(fdf_out,'(a,a,t30,g20.10,1x,a10)')                         &
             '# above item originally: ', label, value, unitstr
      else
        fdf_physical = default
        if (fdf_output) write(fdf_out,'(a,t30,g20.10,1x,a,t60,a)')                        &
             label, default, defunit, '# default value'
      endif

      if (PRESENT(line)) line = mark

      RETURN
!--------------------------------------------------------------------------- END
    END FUNCTION fdf_physical

!
!   Returns a double precision value from a block-line after a certain input value
!   or the default value if label is not found in the fdf file.
!   Converts the units to defunit.
!
    FUNCTION fdf_bphysical(pline, default, defunit, after)
      implicit none
!--------------------------------------------------------------- Input Variables
      type(parsed_line), pointer        :: pline
      real(dp)                          :: default
      character(*)                      :: defunit
      integer(ip), intent(in), optional :: after

!-------------------------------------------------------------- Output Variables
      real(dp)                            :: fdf_bphysical

!--------------------------------------------------------------- Local Variables
      character(10)                       :: unitstr
      character(80)                       :: msg
      real(dp)                            :: value
      type(line_dlist), pointer           :: mark

!------------------------------------------------------------------------- BEGIN
!     Prevents using FDF routines without initialize
      if (.not. fdf_started) then
         call die('FDF module: fdf_bphysical', 'FDF subsystem not initialized', &
              THIS_FILE, __LINE__, fdf_err)
      endif

      if (.not. match(pline, 'vn', after)) then
         write(msg,*) 'no real value for line: '//pline%line
         call die('FDF module: fdf_bphysical', msg, THIS_FILE, &
              __LINE__, fdf_err)
      endif

      ! get value in block-line
      value = values(pline, 1, after)

      ! get unit in block-line
      unitstr = names(pline, 1, after)
      if ( leqi(unitstr, defunit) ) then
         fdf_bphysical = value
      else
         fdf_bphysical = value * fdf_convfac(unitstr, defunit)
      end if

      if ( fdf_output ) then
         if ( present(after) ) then
            write(fdf_out,'(5x,g20.10,1x,a10,1x,i0)') fdf_bphysical, &
                 defunit, after
         else
            write(fdf_out,'(5x,g20.10,1x,a10)') fdf_bphysical, defunit
         end if
         write(fdf_out,'(a,a,t30,g20.10,1x,a10)') &
              '# above item on line: ', pline%line
      end if

!--------------------------------------------------------------------------- END
    END FUNCTION fdf_bphysical

!
!   Returns conversion factor between a subset of physical units
!   Written by j.m.soler. dec'96.
!   Modified by alberto garcia, jan'97.
!   Added more units by Nick Papior, Aug'17.
!
    FUNCTION fdf_convfac(from, to)
      implicit none
!--------------------------------------------------------------- Input Variables
      character(*)           :: from, to

!-------------------------------------------------------------- Output Variables
      real(dp)               :: fdf_convfac

!--------------------------------------------------------------- Local Variables
      character(80)          :: msg
      integer(ip)            :: iu, ifrom, ito

!------------------------------------------------------------------------- BEGIN

!
!     We allow case variations in the units. this could be dangerous
!     (meV --> MeV!!) in real life, but not in this restricted
!     field.

      ! README BEFORE ADDING UNITS:
      !
      ! Units should be added through the small Python code:
      !   fdf_units.py
      ! Add the appropriate unit in the designated unit-type and
      ! run the python script. It will then create (to std-out)
      ! a drop-in replacement for the following lines.

      integer(ip), parameter :: nu = 81
      character(8) :: dimm(nu)
      character(10) :: name(nu)
      real(dp) :: unit(nu)
      data (dimm(iu), name(iu), unit(iu), iu=1, 3) / &
          'mass    ', 'g         ', 1.d-3, &
          'mass    ', 'kg        ', 1.d0, &
          'mass    ', 'amu       ', 1.66054d-27 /

      data (dimm(iu), name(iu), unit(iu), iu=4, 9) / &
          'length  ', 'm         ', 1.d0, &
          'length  ', 'cm        ', 1.d-2, &
          'length  ', 'nm        ', 1.d-9, &
          'length  ', 'pm        ', 1.d-12, &
          'length  ', 'ang       ', 1.d-10, &
          'length  ', 'bohr      ', 0.529177d-10 /

      data (dimm(iu), name(iu), unit(iu), iu=10, 19) / &
          'energy  ', 'j         ', 1.d0, &
          'energy  ', 'kj        ', 1.d3, &
          'energy  ', 'erg       ', 1.d-7, &
          'energy  ', 'mev       ', 1.60219d-22, &
          'energy  ', 'ev        ', 1.60219d-19, &
          'energy  ', 'mry       ', 2.17991d-21, &
          'energy  ', 'ry        ', 2.17991d-18, &
          'energy  ', 'mha       ', 4.35982d-21, &
          'energy  ', 'mhartree  ', 4.35982d-21, &
          'energy  ', 'ha        ', 4.35982d-18 /
      data (dimm(iu), name(iu), unit(iu), iu=20, 29) / &
          'energy  ', 'hartree   ', 4.35982d-18, &
          'energy  ', 'k         ', 1.38066d-23, &
          'energy  ', 'kelvin    ', 1.38066d-23, &
          'energy  ', 'kcal/mol  ', 6.94780d-21, &
          'energy  ', 'kj/mol    ', 1.6606d-21, &
          'energy  ', 'hz        ', 6.6262d-34, &
          'energy  ', 'thz       ', 6.6262d-22, &
          'energy  ', 'cm-1      ', 1.986d-23, &
          'energy  ', 'cm^-1     ', 1.986d-23, &
          'energy  ', 'cm**-1    ', 1.986d-23 /

      data (dimm(iu), name(iu), unit(iu), iu=30, 39) / &
          'time    ', 's         ', 1.d0, &
          'time    ', 'au        ', 2.418884326505e-17, &
          'time    ', 'ns        ', 1.d-9, &
          'time    ', 'ps        ', 1.d-12, &
          'time    ', 'fs        ', 1.d-15, &
          'time    ', 'min       ', 60.d0, &
          'time    ', 'mins      ', 60.d0, &
          'time    ', 'hour      ', 3600.d0, &
          'time    ', 'hours     ', 3600.d0, &
          'time    ', 'day       ', 86400.d0 /
      data (dimm(iu), name(iu), unit(iu), iu=40, 40) / &
          'time    ', 'days      ', 86400.d0 /

      data (dimm(iu), name(iu), unit(iu), iu=41, 44) / &
          'force   ', 'n         ', 1.d0, &
          'force   ', 'ev/ang    ', 1.60219d-9, &
          'force   ', 'ry/bohr   ', 4.11943d-8, &
          'force   ', 'ha/bohr   ', 2.059715d-08 /

      data (dimm(iu), name(iu), unit(iu), iu=45, 54) / &
          'pressure', 'pa        ', 1.d0, &
          'pressure', 'gpa       ', 1.d9, &
          'pressure', 'atm       ', 1.01325d5, &
          'pressure', 'bar       ', 1.d5, &
          'pressure', 'mbar      ', 1.d11, &
          'pressure', 'ev/ang**3 ', 1.60219d11, &
          'pressure', 'ev/ang^3  ', 1.60219d11, &
          'pressure', 'ry/bohr**3', 1.47108d13, &
          'pressure', 'ry/bohr^3 ', 1.47108d13, &
          'pressure', 'ha/bohr^3 ', 2.94216d13 /
      data (dimm(iu), name(iu), unit(iu), iu=55, 55) / &
          'pressure', 'ha/bohr**3', 2.94216d13 /

      data (dimm(iu), name(iu), unit(iu), iu=56, 57) / &
          'charge  ', 'c         ', 1.d0, &
          'charge  ', 'e         ', 1.602177d-19 /

      data (dimm(iu), name(iu), unit(iu), iu=58, 62) / &
          'dipole  ', 'c*m       ', 1.d0, &
          'dipole  ', 'd         ', 3.33564d-30, &
          'dipole  ', 'debye     ', 3.33564d-30, &
          'dipole  ', 'e*bohr    ', 8.47835d-30, &
          'dipole  ', 'e*ang     ', 1.602177d-29 /

      data (dimm(iu), name(iu), unit(iu), iu=63, 64) / &
          'mominert', 'kg*m**2   ', 1.d0, &
          'mominert', 'ry*fs**2  ', 2.17991d-48 /

      data (dimm(iu), name(iu), unit(iu), iu=65, 71) / &
          'efield  ', 'v/m       ', 1.d0, &
          'efield  ', 'v/nm      ', 1.d9, &
          'efield  ', 'v/ang     ', 1.d10, &
          'efield  ', 'v/bohr    ', 1.8897268d10, &
          'efield  ', 'ry/bohr/e ', 2.5711273d11, &
          'efield  ', 'ha/bohr/e ', 5.1422546d11, &
          'efield  ', 'har/bohr/e', 5.1422546d11 /

      data (dimm(iu), name(iu), unit(iu), iu=72, 73) / &
          'angle   ', 'deg       ', 1.d0, &
          'angle   ', 'rad       ', 5.72957795d1 /

      data (dimm(iu), name(iu), unit(iu), iu=74, 81) / &
          'torque  ', 'mev/deg   ', 1.0d-3, &
          'torque  ', 'mev/rad   ', 1.745533d-5, &
          'torque  ', 'ev/deg    ', 1.0d0, &
          'torque  ', 'ev/rad    ', 1.745533d-2, &
          'torque  ', 'mry/deg   ', 13.6058d-3, &
          'torque  ', 'mry/rad   ', 0.237466d-3, &
          'torque  ', 'ry/deg    ', 13.6058d0, &
          'torque  ', 'ry/rad    ', 0.237466d0 /

!
      ifrom = 0
      ito   = 0
      do iu= 1, nu
        if (leqi(name(iu), from)) ifrom = iu
        if (leqi(name(iu), to))   ito   = iu
      end do

      if (ifrom .eq. 0) then
        write(msg,*) 'unknown unit = ', from
        call die('FDF module: fdf_convfac', msg, THIS_FILE, __LINE__, fdf_err)
      endif

      if (ito .eq. 0) then
        write(msg,*) 'unknown unit = ', to
        call die('FDF module: fdf_convfac', msg, THIS_FILE, __LINE__, fdf_err)
      endif

      if (leqi(dimm(ifrom), dimm(ito))) then
        fdf_convfac = unit(ifrom) / unit(ito)
      else
        write(msg,*) 'unit''s physical dimensions don''t match: ',      &
                     from, ', ', to
        call die('FDF module: fdf_convfac', msg, THIS_FILE, __LINE__, fdf_err)
      endif
!--------------------------------------------------------------------------- END
    END FUNCTION fdf_convfac

!
!   Searches for label in the fdf hierarchy. If it appears the function
!   returns .TRUE. and leaves mark pointer positioned at the line.
!   Otherwise, it returns .FALSE. and mark points to NULL.
!
    FUNCTION fdf_locate(label, mark)
      implicit none
!--------------------------------------------------------------- Input Variables
      character(*)              :: label

!-------------------------------------------------------------- Output Variables
      logical                   :: fdf_locate
      type(line_dlist), pointer :: mark

!--------------------------------------------------------------- Local Variables
      character(80)             :: strlabel

!------------------------------------------------------------------------- BEGIN
      fdf_locate = .FALSE.

!      if (fdf_donothing) return

      mark => file_in%first
      do while ((.not. fdf_locate) .and. (ASSOCIATED(mark)))

        if (match(mark%pline, 'l')) then
          strlabel = labels(mark%pline)

          if (labeleq(strlabel, label, fdf_log)) then
            fdf_locate = .TRUE.
          else
            mark => mark%next
          endif
        else
          mark => mark%next
        endif
      enddo

      if (.not. fdf_locate) NULLIFY(mark)

      RETURN
!--------------------------------------------------------------------------- END
    END FUNCTION fdf_locate


!
!   Searches for label in the fdf hierarchy. If it appears the function
!   returns .TRUE. and leaves mark pointer positioned at the line.
!   Otherwise, it returns .FALSE. and mark points to NULL.
!
    FUNCTION fdf_load_locate(label, mark)
      implicit none
!--------------------------------------------------------------- Input Variables
      character(*)              :: label

!-------------------------------------------------------------- Output Variables
      logical                   :: fdf_load_locate
      type(line_dlist), pointer :: mark

!--------------------------------------------------------------- Local Variables
      character(80)             :: strlabel

!------------------------------------------------------------------------- BEGIN
      fdf_load_locate = .FALSE.

!      if (fdf_donothing) return

      mark => file_in%first
      do while ((.not. fdf_load_locate) .and. (ASSOCIATED(mark)))

        if (match(mark%pline, 'l')) then
          strlabel = tokens(mark%pline, 2)

          if (labeleq(strlabel, label, fdf_log)) then
            fdf_load_locate = .TRUE.
          else
            mark => mark%next
          endif
        else
          mark => mark%next
        endif
      enddo

      if (.not. fdf_load_locate) NULLIFY(mark)

      RETURN
!--------------------------------------------------------------------------- END
    END FUNCTION fdf_load_locate



!
!   Returns true or false whether or not the label 'label' is
!   a block.
!   I.e. it returns true if the line has the form bl, if not found, or not bl
!   it returns false.
!
    FUNCTION fdf_isblock(label)
      implicit none
!--------------------------------------------------------------- Input Variables
      character(*)                        :: label

!-------------------------------------------------------------- Output Variables
      logical                             :: fdf_isblock

!--------------------------------------------------------------- Local Variables
      type(line_dlist), pointer :: mark
      character(80) :: strlabel

!------------------------------------------------------------------------- BEGIN
!     Prevents using FDF routines without initialize
      if (.not. fdf_started) then
         call die('FDF module: fdf_isblock', 'FDF subsystem not initialized', &
                 THIS_FILE, __LINE__, fdf_err)
      endif

      fdf_isblock = .false.

      mark => file_in%first
      do while ( associated(mark) )

!!$        if ( match(mark%pline, 'l') ) then
!!$          strlabel = labels(mark%pline)
!!$
!!$          if ( labeleq(strlabel, label, fdf_log) ) then
!!$            ! fdf has first-encounter acceptance.
!!$            ! I.e. for an input
!!$            !   Label_Name 1
!!$            !   %block Label_Name
!!$            !     1
!!$            !   %endblock Label_Name
!!$            ! the former will be accepted first.
!!$            exit
!!$          end if

        if ( match(mark%pline, 'bl') ) then
          strlabel = blocks(mark%pline)

          if ( labeleq(strlabel, label, fdf_log) ) then
            fdf_isblock = .true.
            exit
          end if
        end if

        mark => mark%next
      end do

      if (fdf_output) write(fdf_out,'(a,t30,l10)') "Is block " // trim(label) // ' present?',  fdf_isblock

      RETURN
!--------------------------------------------------------------------------- END
    END FUNCTION fdf_isblock

!
!   Searches for block label in the fdf hierarchy. If it appears returns
!   .TRUE. and leaves block mark pointer positioned at the first line.
!   Otherwise, it returns .FALSE. and block mark points to NULL.
!
    FUNCTION fdf_block(label, bfdf)
      implicit none
!--------------------------------------------------------------- Input Variables
      character(*)    :: label

!-------------------------------------------------------------- Output Variables
      logical         :: fdf_block
      type(block_fdf) :: bfdf

!--------------------------------------------------------------- Local Variables
      character(80)   :: strlabel

!------------------------------------------------------------------------- BEGIN
!     Prevents using FDF routines without initialize
      if (.not. fdf_started) then
        call die('FDF module: fdf_block', 'FDF subsystem not initialized', &
                 THIS_FILE, __LINE__, fdf_err)
      endif

      fdf_block = .FALSE.

      bfdf%mark => file_in%first
      do while ((.not. fdf_block) .and. (ASSOCIATED(bfdf%mark)))

        if (match(bfdf%mark%pline, 'bl')) then
          strlabel = blocks(bfdf%mark%pline)

          if (labeleq(strlabel, label, fdf_log)) then
            fdf_block = .TRUE.
            bfdf%label = label

            if (fdf_output) write(fdf_out,'(a,a)') '%block ', TRIM(label)
          endif
        endif

        bfdf%mark => bfdf%mark%next
      enddo

      if (.not. fdf_block) NULLIFY(bfdf%mark)

      RETURN
!--------------------------------------------------------------------------- END
    END FUNCTION fdf_block

!
!   Get successive parsed lines from block returning
!   .TRUE. while more lines exist in the block bfdf.
!
    FUNCTION fdf_bline(bfdf, pline)
      implicit none
!--------------------------------------------------------------- Input Variables
      type(block_fdf)            :: bfdf

!-------------------------------------------------------------- Output Variables
      logical                    :: fdf_bline
      type(parsed_line), pointer :: pline

!--------------------------------------------------------------- Local Variables
      character(80)              :: strlabel

!------------------------------------------------------------------------- BEGIN
!     Prevents using FDF routines without initialize
      if (.not. fdf_started) then
        call die('FDF module: fdf_bline', 'FDF subsystem not initialized', &
                 THIS_FILE, __LINE__, fdf_err)
      endif

      if (.not. ASSOCIATED(bfdf%mark)) then
        call die('FDF module: fdf_bline', 'block_fdf structure not initialized', &
                 THIS_FILE, __LINE__, fdf_err)
      endif

      fdf_bline = .TRUE.

!     If we are in the head of the block move to the content
      if (match(bfdf%mark%pline, 'bl')) then
        strlabel = blocks(bfdf%mark%pline)

        if (labeleq(strlabel, bfdf%label, fdf_log)) then
          bfdf%mark => bfdf%mark%next

          if (fdf_output) write(fdf_out,'(a,a)') '%block ', TRIM(bfdf%label)
        endif
      endif

      if (match(bfdf%mark%pline, 'el')) then
        strlabel = endblocks(bfdf%mark%pline)

        if (labeleq(strlabel, bfdf%label, fdf_log)) then
          fdf_bline = .FALSE.
          NULLIFY(pline)

          if (fdf_output) write(fdf_out,'(a,a)') '%endblock ', TRIM(bfdf%label)
        endif
      endif

      if (fdf_bline) then
        if (fdf_output) write(fdf_out,'(1x,a)') TRIM(bfdf%mark%str)

        pline     => bfdf%mark%pline
        bfdf%mark => bfdf%mark%next
      endif

      RETURN
!--------------------------------------------------------------------------- END
    END FUNCTION fdf_bline



!
!   Backspace to the previous physical line in the block
!   returning .TRUE. while more lines exist in the block bfdf.
!
    FUNCTION fdf_bbackspace(bfdf,pline)
      implicit none
!--------------------------------------------------------------- Input Variables
      type(block_fdf)            :: bfdf

!-------------------------------------------------------------- Output Variables
      logical                    :: fdf_bbackspace
      type(parsed_line), pointer, optional :: pline

!--------------------------------------------------------------- Local Variables
      character(80)              :: strlabel

!------------------------------------------------------------------------- BEGIN
!     Prevents using FDF routines without initialize
      if (.not. fdf_started) then
        call die('FDF module: fdf_bbackspace', 'FDF subsystem not initialized', &
                 THIS_FILE, __LINE__, fdf_err)
      endif

      if (.not. ASSOCIATED(bfdf%mark)) then
        call die('FDF module: fdf_bbackspace', 'block_fdf structure not initialized', &
                 THIS_FILE, __LINE__, fdf_err)
      endif

      fdf_bbackspace = .TRUE.

!     If we are in the bottom of the block move to the content

      if (match(bfdf%mark%pline, 'el')) then

        strlabel = endblocks(bfdf%mark%pline)

        if (labeleq(strlabel, bfdf%label, fdf_log)) then
          bfdf%mark => bfdf%mark%prev

          if (fdf_output) write(fdf_out,'(1x,a)') "#:(Backspace to) " // "|" //  &
                                TRIM(bfdf%mark%str) // "|"
        endif

!     If we are at the head we cannot backspace

      else if (match(bfdf%mark%pline, 'bl')) then
        strlabel = blocks(bfdf%mark%pline)

        if (labeleq(strlabel, bfdf%label, fdf_log)) then
          fdf_bbackspace = .FALSE.
          if (fdf_output) write(fdf_out,'(1x,a)') "#:(Cannot backspace) " // "|" //  &
                                TRIM(bfdf%mark%str) // "|"
        endif

      else

        bfdf%mark => bfdf%mark%prev
        if (fdf_output) write(fdf_out,'(1x,a)') "#:(Backspace to) " // "|" //  &
                                TRIM(bfdf%mark%str) // "|"
      endif

      if ( present(pline) ) pline => bfdf%mark%pline

      RETURN
!--------------------------------------------------------------------------- END
    END FUNCTION fdf_bbackspace

!
!   Moves the pointer of the working line (bfdf%mark)
!   to the beginning of the block 'label' structure.
!
    SUBROUTINE fdf_brewind(bfdf)
      implicit none
!-------------------------------------------------------------- Output Variables
      type(block_fdf) :: bfdf

!--------------------------------------------------------------- Local Variables
      character(80)   :: msg

!------------------------------------------------------------------------- BEGIN
!     Prevents using FDF routines without initialize
      if (.not. fdf_started) then
        call die('FDF module: fdf_brewind', 'FDF subsystem not initialized', &
                 THIS_FILE, __LINE__, fdf_err)
      endif

      if (.not. ASSOCIATED(bfdf%mark)) then
        call die('FDF module: fdf_brewind', 'block_fdf structure not initialized', &
                 THIS_FILE, __LINE__, fdf_err)
      endif

      if (.not. fdf_block(bfdf%label, bfdf)) then
        write(msg,*) 'Block ', bfdf%label, ' not found in FDF structure'
        call die('FDF module: fdf_brewind', msg, &
                 THIS_FILE, __LINE__, fdf_err)
      endif

      RETURN
!--------------------------------------------------------------------------- END
    END SUBROUTINE fdf_brewind

!
!   Closes the opened block by looping the remaining lines of the working line.
!   This is only needed to complete the fdf-*.log files output for direct
!   usage later. It does nothing internally.
!
    SUBROUTINE fdf_bclose(bfdf)
      implicit none
!-------------------------------------------------------------- Output Variables
      type(block_fdf) :: bfdf

!--------------------------------------------------------------- Local Variables
      type(parsed_line), pointer :: pline
      integer(ip) :: i
      character(80) :: msg

!------------------------------------------------------------------------- BEGIN
!     Prevents using FDF routines without initialize
      if (.not. fdf_started) then
        call die('FDF module: fdf_bclose', 'FDF subsystem not initialized', &
                 THIS_FILE, __LINE__, fdf_err)
      endif

      ! Quick return (no need for errors)
      if ( .not. associated(bfdf%mark) ) return

      ! This should hopefully discourage compilers to optimize the loop away...
      i = 0
      do while ( fdf_bline(bfdf, pline) )
        i = i + fdf_bnvalues(pline)
      end do
      write(msg,'(a,i10)') 'Block ', i

      RETURN
!--------------------------------------------------------------------------- END
    END SUBROUTINE fdf_bclose


!
!   Count number of lines with an optional specification.
!   I.e. this will read through the block and return the number of lines in the
!   block which matches the morphology (morph)
!   This may be used to easily digest number of non-empty lines in the block.
!   Note that a match on the morphology only compares the number of ID's in
!   morph. I.e. a line with 'vvvil' will match 'vvvi'.
!
    FUNCTION fdf_block_linecount(label, morph)
      implicit none
!--------------------------------------------------------------- Input Variables
      character(len=*) :: label
      character(len=*), optional :: morph
!-------------------------------------------------------------- Output Variables
      integer(ip) :: fdf_block_linecount

!--------------------------------------------------------------- Local Variables
      type(block_fdf) :: bfdf
      type(parsed_line), pointer :: pline
      logical :: orig_fdf_output

!------------------------------------------------------------------------- BEGIN
!     Prevents using FDF routines without initialize
      if (.not. fdf_started) then
        call die('FDF module: fdf_block_linecount', 'FDF subsystem not initialized', &
                 THIS_FILE, __LINE__, fdf_err)
      endif

      ! Store the fdf_output variable (suppress writing to log)
      orig_fdf_output = fdf_output
      fdf_output = .false.

      ! Find the block and search for morhp
      fdf_block_linecount = 0
      if ( fdf_block(label, bfdf) ) then

        do while ( fdf_bline(bfdf, pline) )
          if ( present(morph) ) then
            if ( fdf_bmatch(pline, morph) ) then
              fdf_block_linecount = fdf_block_linecount + 1
            end if
          else
            fdf_block_linecount = fdf_block_linecount + 1
          end if
        end do

      end if

      ! Restore output
      fdf_output = orig_fdf_output

      if ( fdf_output ) then
        if ( present(morph) ) then
          write(fdf_out,'(3a,3x,i0)') "#:block-line-count? ", &
              trim(label), ' ('//trim(morph)//')', fdf_block_linecount
        else
          write(fdf_out,'(2a,3x,i0)') "#:block-line-count? ", &
              trim(label), fdf_block_linecount
        end if
      end if

      RETURN
!--------------------------------------------------------------------------- END
    END FUNCTION fdf_block_linecount

!
!   Check if label is defined
!
    logical FUNCTION fdf_defined(label)
      implicit none
!--------------------------------------------------------------- Input Variables
      character(*)              :: label

!--------------------------------------------------------------- Local Variables
      type(line_dlist), pointer :: mark

!--------------------------------------------------------------------- BEGIN
      ! First, check whether a single label exists:
      fdf_defined = fdf_locate(label, mark)
      if (.not. fdf_defined) then
         ! Check whether there is a block with that label
         fdf_defined = fdf_isblock(label)
      endif
      if ( fdf_output ) then
        write(fdf_out,'(a,t30,l10)') 'Is defined? ' // label, fdf_defined
      endif

      RETURN
!----------------------------------------------------------------------- END
    END FUNCTION fdf_defined

!   Check if label is defined
!
    logical FUNCTION fdf_load_defined(label)
      implicit none
!--------------------------------------------------------------- Input Variables
      character(*)              :: label

!--------------------------------------------------------------- Local Variables
      type(line_dlist), pointer :: mark

!--------------------------------------------------------------------- BEGIN
      ! First, check whether a single label exists:
      fdf_load_defined = fdf_load_locate(label, mark)
      if (.not. fdf_load_defined) then
         ! Check whether there is a block with that label
         fdf_load_defined = fdf_isblock(label)
      endif
      if ( fdf_output ) then
        write(fdf_out,'(a,t30,l10)') 'Is  load ' //  trim(label) //  ' defined? ' , fdf_load_defined
      endif

      RETURN
!----------------------------------------------------------------------- END
    END FUNCTION fdf_load_defined





!
!   Output levels:
!   level <= 0: nothing
!   level  = 1: standard
!
    SUBROUTINE fdf_setoutput(level,fileout_in)
      implicit none
!------------------------------------------------------------- Input Variables
      integer(ip)                  :: level
      character(len=*), intent(in) :: fileout_in


      character(len=256) :: fileout

      fileout = fileout_in
      if (level .le. 0) then
        if (fdf_output) then
          call io_close(fdf_out)
          fdf_output = .FALSE.
        endif
      else
        if (.not. fdf_output) then
#ifdef _MPI_
          if (rank /= 0) then
             if (level .ge. 2) then
                fileout =  trim(fileout) // "." // i2s(rank)
                call io_assign(fdf_out)
                open(fdf_out, file=fileout, form='formatted',               &
                     status='unknown')
                REWIND(fdf_out)
                fdf_output = .TRUE.
             endif
          else
             call io_assign(fdf_out)
             open(fdf_out, file=fileout, form='formatted',               &
                  status='unknown')
             REWIND(fdf_out)
             fdf_output = .TRUE.
          endif
#else
          call io_assign(fdf_out)
          open(fdf_out, file=fileout, form='formatted',               &
               status='unknown')
          REWIND(fdf_out)
          fdf_output = .TRUE.
#endif
        endif
      endif
!----------------------------------------------------------------------- END
    END SUBROUTINE fdf_setoutput

!
!   Debugging levels:
!   level <= 0: nothing
!   level  = 1: standard
!   level >= 2: exhaustive

    SUBROUTINE fdf_setdebug(level,filedebug)
      implicit none
!------------------------------------------------------------- Input Variables
      integer(ip)      :: level
      character(len=*) :: filedebug

!----------------------------------------------------------------------- BEGIN
      if (level .le. 0) then
        if (fdf_debug) then
          call io_close(fdf_log)
          fdf_debug = .FALSE.
        endif
      else
        if (.not. fdf_debug) then
          call io_assign(fdf_log)

          print*, "filedebug master before ", filedebug
#ifdef _MPI_
          if (rank /= 0) then
             filedebug =  trim(filedebug) // "." // i2s(rank)
             print*, "filedebug ", filedebug
          endif
#endif
          open(fdf_log, file=filedebug, form='formatted',               &
               status='unknown')
          REWIND(fdf_log)
          fdf_debug = .TRUE.

!         Set logging/debugging info for PARSE module also
          call setlog(fdf_log)
          call setdebug(1)
        endif
      endif

      fdf_debug2 = (level .ge. 2)

      RETURN
!----------------------------------------------------------------------- END
    END SUBROUTINE fdf_setdebug

!
!   For handling deprecated labels.
!   Also there is an optional "newlabel" if it has been changed into
!   a new label.
!
    subroutine fdf_deprecated(label,newlabel)
      implicit none
!--------------------------------------------------------------- Input Variables
      character(*)           :: label
      character(*)           :: newlabel

!------------------------------------------------------------------------- BEGIN
      if ( fdf_defined(label) ) then
         if (fdf_output) write(fdf_out,'(a)') "#**Warning: FDF symbol '"//trim(label)// &
              "' is deprecated."
         if ( fdf_defined(newlabel) ) then
            if (fdf_output) write(fdf_out,'(a)') "#           FDF symbol '"//trim(newlabel)// &
                 "' will be used instead."
         else
            if (fdf_output) write(fdf_out,'(a)') "#           FDF symbol '"//trim(newlabel)// &
                 "' replaces '"//trim(label)//"'."
         end if
      end if

!--------------------------------------------------------------------------- END
    end subroutine fdf_deprecated

!
!   For handling obsoleted labels.
!
    subroutine fdf_obsolete(label)
      implicit none
!--------------------------------------------------------------- Input Variables
      character(*)           :: label

!------------------------------------------------------------------------- BEGIN
      if ( fdf_defined(label) ) then
         if (fdf_output) write(fdf_out,'(a)') "#**Warning: FDF symbol '"//trim(label)// &
              "' is obsolete."
      end if

!--------------------------------------------------------------------------- END
    end subroutine fdf_obsolete

!===================== Serialization utilities for clients

    subroutine serialize_fdf_struct(buffer)
    character, pointer        :: buffer(:)
!    character(len=1), intent(inout), allocatable   :: buffer(:)

    character(len=SERIALIZED_LENGTH)  bufline
    type(line_dlist), pointer :: mark
    integer(ip) :: i, length, init, final

!    integer :: nchars ! total size of serialized content
!
!    if (allocated(buffer)) deallocate(buffer)
!    nchars = file_in%nlines * SERIALIZED_LENGTH
!    allocate(buffer(nchars))

    mark => file_in%first
    do i= 1, file_in%nlines
       call serialize_pline(mark%pline,bufline,length)
       init  = (i-1)*SERIALIZED_LENGTH+1
       final = (i)*SERIALIZED_LENGTH
       call convert_string_to_array_of_chars(bufline,buffer(init:final))
       mark => mark%next
    enddo
  end subroutine serialize_fdf_struct

    subroutine recreate_fdf_struct(nlines,bufferFDF)
    character(len=1), intent(in)    :: bufferFDF(:)
    integer, intent(in)       :: nlines

    character(len=SERIALIZED_LENGTH)  bufline
    type(parsed_line), pointer    :: pline
    integer(ip) :: i, init, final

!    nlines = size(bufferFDF) / SERIALIZED_LENGTH

    do i= 1, nlines
       init  = (i-1)*SERIALIZED_LENGTH+1
       final = (i)*SERIALIZED_LENGTH
       call convert_array_of_chars_to_string(bufferFDF(init:final),bufline)
       allocate(pline)
       call recreate_pline(pline,bufline)
       call fdf_addtoken(pline%line,pline)
    enddo

  end subroutine recreate_fdf_struct

!=================================================================
!--------- Obsolete routines ----------------
!
!   The owner node of the input file sends data file to the
!   other processes in the MPI communicator.
!
#ifdef CLUSTER
  SUBROUTINE fdf_sendInput()

    use mpi
    implicit none
!--------------------------------------------------------------- Local Variables
    character, pointer        :: bufferFDF(:) => null()
    integer(ip)               :: i, j, k, ierr
    type(line_dlist), pointer :: mark

!------------------------------------------------------------------------- BEGIN
    ALLOCATE(bufferFDF(file_in%nlines*MAX_LENGTH), stat=ierr)
    if (ierr .ne. 0) then
      call die('FDF module: fdf_sendInput', 'Error allocating bufferFDF', &
               THIS_FILE, __LINE__, fdf_err, rc=ierr)
    endif

    mark => file_in%first
    do i= 1, file_in%nlines*MAX_LENGTH, MAX_LENGTH
      bufferFDF(i:i+MAX_LENGTH-1) = s2arr(mark%str)
      mark => mark%next
    enddo
    print *, " printing from fdf_sendInput cluster mode "
    call MPI_Bcast(file_in%nlines, 1,                                 &
                   MPI_INTEGER, rank, MPI_COMM_WORLD, ierr)

    if (ierr .ne. MPI_SUCCESS) then
      call die('FDF module: fdf_sendInput', 'Error Broadcasting nlines.' // &
               'Terminating.', THIS_FILE, __LINE__, fdf_err, rc=ierr)
    endif

    call MPI_Bcast(bufferFDF, file_in%nlines*MAX_LENGTH,              &
                   MPI_CHARACTER, rank, MPI_COMM_WORLD, ierr)
    if (ierr .ne. MPI_SUCCESS) then
      call die('FDF module: fdf_sendInput', 'Error Broadcasting bufferFDF.' // &
               'Terminating.', THIS_FILE, __LINE__, fdf_err, rc=ierr)
    endif

    DEALLOCATE(bufferFDF)

    RETURN
!--------------------------------------------------------------------------- END
  END SUBROUTINE fdf_sendInput
#endif

!
!   Remaining nodes recieve the input data file from the owner
!   of the FDF file inside the cluster.
!
#ifdef CLUSTER
  SUBROUTINE fdf_recvInput(root, filein, fileinTmp)
    !use mpi_siesta
    use mpi
    implicit none
!--------------------------------------------------------------- Input Variables
    character(*)       :: filein, fileinTmp
    integer(ip)        :: root

!--------------------------------------------------------------- Local Variables
    integer(ip)        :: i, j, lun, ierr, nlines
    character, pointer :: bufferFDF(:) => null()
    character(len=10)  :: fmt
!------------------------------------------------------------------------- BEGIN
    call MPI_Bcast(nlines, 1,                                         &
                   MPI_INTEGER, root, MPI_COMM_WORLD, ierr)
    if (ierr .ne. MPI_SUCCESS) then
      call die('FDF module: fdf_recvInput', 'Error Broadcasting nlines.' // &
               'Terminating.', THIS_FILE, __LINE__, fdf_err, rc=ierr)
    endif

    ALLOCATE(bufferFDF(nlines*MAX_LENGTH), stat=ierr)
    if (ierr .ne. 0) then
      call die('FDF module: fdf_recvInput', 'Error allocating bufferFDF', &
               THIS_FILE, __LINE__, fdf_err, rc=ierr)
    endif

    call MPI_Bcast(bufferFDF, nlines*MAX_LENGTH,                      &
                   MPI_CHARACTER, root, MPI_COMM_WORLD, ierr)
    if (ierr .ne. MPI_SUCCESS) then
      call die('FDF module: fdf_recvInput', 'Error Broadcasting bufferFDF.' // &
               'Terminating.', THIS_FILE, __LINE__, fdf_err, rc=ierr)
    endif

    call io_assign(lun)
    fileinTmp = TRIM(filein) // '.' // i2s(rank)
    open(unit=lun, file=fileinTmp, form='formatted',                  &
         status='unknown')

    if (MAX_LENGTH.lt.10) then
      write(fmt,"(a1,I1,a2)") "(", MAX_LENGTH, "a)"
    else if (MAX_LENGTH.lt.100) then
      write(fmt,"(a1,I2,a2)") "(", MAX_LENGTH, "a)"
    else if (MAX_LENGTH.lt.1000) then
      write(fmt,"(a1,I3,a2)") "(", MAX_LENGTH, "a)"
    else
      write(fmt,"(a1,I4,a2)") "(", MAX_LENGTH, "a)"
    endif

    do i= 1, nlines*MAX_LENGTH, MAX_LENGTH
      write(lun,fmt) (bufferFDF(j),j=i,i+MAX_LENGTH-1)
    enddo
    call io_close(lun)

    DEALLOCATE(bufferFDF)

    RETURN
!--------------------------------------------------------------------------- END
  END SUBROUTINE fdf_recvInput
#endif

      ! To enable client-side setting
    SUBROUTINE fdf_set_started(status)
      logical, intent(in) :: status

      fdf_started = status
    end SUBROUTINE fdf_set_started

END MODULE fdf
