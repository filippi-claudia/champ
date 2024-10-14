!> Module to control the flow of the program
module control

    implicit none

    !> mode of the program
    character(len=12) :: mode

    !> ipr flag; decides the verbosity of the output
    integer  :: ipr

    private
    public :: mode, ipr, init_control_mode
    save
    contains
    subroutine init_control_mode(str_mode)
        implicit None
        character(12), intent(IN) :: str_mode
        mode = str_mode
    end subroutine init_control_mode

end module control

!> Module to control the VMC calculation mode
!> @var vmc_idump: Whether to store the dump files at the end of the VMC calculation
!> @var vmc_irstar: Whether to restart the VMC calculation
!> @var vmc_isite: Whether to read the configuration from a file or generate it
!> @var vmc_nconf: Number of configurations to generate
!> @var vmc_nblk: Number of blocks to run
!> @var vmc_nblk_max: Maximum number of blocks to run
!> @var vmc_nblkeq: Number of blocks to equilibrate
!> @var vmc_nconf_new: Number of configurations to generate in the new block
!> @var vmc_nstep: Number of steps to run in a block
!> @var vmc_icharged_atom: Index of the charged atom
!> @var vmc_nblk_ci: Number of blocks to run in the CI calculation
module control_vmc

    !> Whether to store the dump files at the end of the VMC calculation
    integer :: vmc_idump

    !> Whether to restart the VMC calculation
    integer :: vmc_irstar

    !> Whether to read the configuration from a file or generate it
    integer :: vmc_isite

    !> Number of configurations to generate
    integer :: vmc_nconf

    !> Number of blocks to run
    integer :: vmc_nblk

    !> Maximum number of blocks to run
    integer :: vmc_nblk_max

    !> Number of blocks to equilibrate
    integer :: vmc_nblkeq

    !> Number of configurations to generate in the new block
    integer :: vmc_nconf_new

    !> Number of steps to run in a block
    integer :: vmc_nstep

    !> Index of the charged atom
    integer :: vmc_icharged_atom

    !> Number of blocks to run in the CI calculation
    integer :: vmc_nblk_ci

    private
    public :: vmc_idump, vmc_irstar, vmc_isite, vmc_nconf, vmc_nblk, vmc_nblk_max
    public :: vmc_nblkeq, vmc_nconf_new, vmc_nstep, vmc_icharged_atom, vmc_nblk_ci
    save
end module control_vmc


!> Module to control the DMC calculation mode
!> @var dmc_idump: Whether to store the dump files at the end of the DMC calculation
!> @var dmc_irstar: Whether to restart the DMC calculation
!> @var dmc_isite: Whether to read the configuration from a file or generate it
!> @var dmc_nconf: Number of configurations to generate
!> @var dmc_nblk: Number of blocks to run
!> @var dmc_nblkeq: Number of blocks to equilibrate
!> @var dmc_nconf_new: Number of configurations to generate in the new block
!> @var dmc_nstep: Number of steps to run in a block

module control_dmc

    !> Whether to store the dump files at the end of the DMC calculation
    integer :: dmc_idump

    !> Whether to restart the DMC calculation
    integer :: dmc_irstar

    !> Whether to read the configuration from a file or generate it
    integer :: dmc_isite

    !> Number of configurations to generate
    integer :: dmc_nconf

    !> Number of blocks to run
    integer :: dmc_nblk

    !> Number of blocks to equilibrate
    integer :: dmc_nblkeq

    !> Number of configurations to generate in the new block
    integer :: dmc_nconf_new

    !> Number of steps to run in a block
    integer :: dmc_nstep

    private
    public :: dmc_idump, dmc_irstar, dmc_isite, dmc_nconf, dmc_nblk
    public :: dmc_nblkeq, dmc_nconf_new, dmc_nstep
    save
end module control_dmc

!> Module to control the periodic calculations
!> @var iperiodic: periodic boundary conditions
!> @var ibasis: the basis set
module contrl_per

    implicit none

    !> periodic boundary conditions
    integer :: iperiodic

    !> the basis set
    integer :: ibasis

    private
    public   :: iperiodic, ibasis
    save
end module contrl_per

module contrldmc
    !> Arguments: iacc_rej, icross, icuspg, icut_br, icut_e, idiv_v, idmc, ipq, itau_eff, nfprod, rttau, tau, taueff, tautot
      use multiple_geo, only: MFORCE
      use precision_kinds, only: dp

    implicit none

    integer :: iacc_rej
    integer :: icross
    integer :: icuspg
    integer :: icut_br
    integer :: icut_e
    integer :: idiv_v
    integer :: idmc
    integer :: ipq
    integer :: itau_eff
    integer :: nfprod
    integer :: limit_wt_dmc
    real(dp) :: ibranching_c
    real(dp) :: rttau
    real(dp) :: tau
    real(dp), dimension(:), allocatable :: taueff !(MFORCE)
    real(dp) :: tautot

    private
    public :: iacc_rej, icross, icuspg, icut_br, icut_e, ibranching_c, idiv_v, idmc, ipq, itau_eff, nfprod, rttau, tau, taueff, tautot, limit_wt_dmc
    public :: allocate_contrldmc, deallocate_contrldmc
    save
contains
    subroutine allocate_contrldmc()
      use multiple_geo, only: MFORCE
        if (.not. allocated(taueff)) allocate (taueff(MFORCE))
    end subroutine allocate_contrldmc

    subroutine deallocate_contrldmc()
        if (allocated(taueff)) deallocate (taueff)
    end subroutine deallocate_contrldmc

end module contrldmc


!> Module to control the file handling
module contrl_file
#if defined(TREXIO_FOUND)
      use trexio,  only: trexio_back_end_t
#endif
    implicit none

    !> Name of the log file
    character(20) :: log_filename

    !> Name of the proc file
    character(20) :: proc_filename

    !> Name of the input file
    character(80) :: file_input

    !> Name of the output file
    character(80) :: file_output

    !> Name of the error file
    character(80) :: file_error

    !> Unit number for the input file
    integer       :: iunit

    !> Unit number for the output file
    integer       :: ounit

    !> Unit number for the error file
    integer       :: errunit

#if defined(TREXIO_FOUND)
    integer(trexio_back_end_t) :: backend
#endif

    private
    public :: log_filename, proc_filename
    public :: file_input, file_output, file_error
    public :: close_files
    public :: init_procfile, init_logfile, initialize
    public :: iunit, ounit, errunit
#if defined(TREXIO_FOUND)
    public :: backend
#endif
    save
contains

    !> Subroutine to close the files
    subroutine close_files()
      use mpiconf, only: wid
        close (5)
        close (6)
        if (wid) close (45)
    end subroutine close_files

    !> Subroutine to initialize the log file
    subroutine init_logfile()
      use mpiconf, only: wid

        !> Open the standard output and the log file only on the master
        if (wid) then
            log_filename = 'output.log'
            open (45, file=log_filename, status='unknown')
        else
            log_filename = '/dev/null'
            close (6)
            open (6, file='/dev/null')
        endif
    end subroutine init_logfile

    !> Subroutine to initialize the file handling
    !> @author Ravindra Shinde
    subroutine initialize()
        use mpiconf, only: wid      ! logical :: true only for mpirank=0

        use, intrinsic :: iso_fortran_env !, only: stdin=>input_unit, stdout=>output_unit, stderr=>error_unit


        implicit none
        character(len=40), parameter            :: VERSION = 'v2.0.0'
        character(len=80), allocatable          :: arg(:)
        integer                                 :: i, j, iostat, argcount
        character(len=10), dimension(12)        :: extensions
        character(len=100)                      :: string_format  = '(A, T40, A)'

        ! Initialize file names to empty variables to avoid garbage in the output
        file_error = ''
        file_input = ''
        file_output = ''



        ! Make sure ounit default is stdout, and errunit is stderr
        ounit = output_unit
        errunit = error_unit
        ! Get all the command line arguments


        argcount = command_argument_count()
        if ( .not. allocated(arg)) allocate(arg(12))
        arg = ""
        do i = 1, argcount
            call get_command_argument(i, arg(i))
        end do

        do i = 1,argcount
            select case (arg(i))
                case ('-v', '-ver', '--version')
                    write(*, '(a)') 'version ', VERSION
                    stop

                case ('-V', '--verbose')
                    write(output_unit, '(a)') 'verbose mode '

                case ('-d', '-debug', '--debug')
                    write(output_unit, '(a)') 'debug mode '

                case ('-p', '-prefix', '--prefix')
                    write(output_unit, '(a)') 'prefix text for generated files '

                case ('-i', '-in', '-inp', '-input', '--input')
                    if ((index(arg(i+1), ".in") /= 0) .or. (index(arg(i+1), ".inp") /= 0) .or. &
                        (index(arg(i+1), ".dat") /= 0) ) then
                        file_input = arg(i+1)
                    else
                        write(error_unit,*) "input file should have an extention .in / .inp / .dat"
                        stop
                    endif

                case ('-o', '-ou', '-out', '-output', '--output')
                    file_output = arg(i+1)
                    if (.not. wid ) then
                        file_output = '/dev/null'
                        close (6)
                        open (6, file='/dev/null')
                    else
                      open (newunit=ounit,file=file_output, iostat=iostat, action='write', status='unknown' )
                      if (iostat /= 0) error stop "error in opening output unit"
                    endif

                case ('-e', '-er', '-err', '-error', '--error')
                    file_error = arg(i+1)
                    if (.not. wid ) then
                        file_error = '/dev/null'
                    endif
                    open (newunit=errunit,file=file_error, iostat=iostat, action='write' )

                case ('-h', '--help')
                    call print_help()
                    stop

                case default
                    extensions(1) = ".in"  ; extensions(2) = ".inp" ; extensions(3) = ".dat" ; extensions(4) = ".out" ;
                    extensions(5) = ".log" ; extensions(6) = ".err" ; extensions(7) = ".e"   ; extensions(8) = "error"

                    ! default error file if not mentioned.
                    do j = 1, 8
                        if (index(arg(i+1), extensions(j)) /= 0) then
                            write(output_unit, '(2a)') 'unrecognised command-line option: ', arg(i)
                            call print_help()
                        endif
                    enddo

            end select
        enddo
        if ( allocated(arg)) deallocate(arg)


        contains
        subroutine print_help()
            write(output_unit, '(a, /)') "Command line options::"
            write(output_unit, '(a, /)') "  '-i', '-in', '-inp', '-input',  '--input'   Provide input  file i.e. '-i input.inp' "
            write(output_unit, '(a, /)') "  '-o', '-ou', '-out', '-output', '--output'  Provide output file i.e. '-o output.log' "
            write(output_unit, '(a, /)') "  '-e', '-er', '-err', '-error',  '--error'   Provide error  file i.e. '-e error' "

            write(output_unit, '(a, /)') "  '-v', '-ver', '--version'                   print version information"
            write(output_unit, '(a, /)') "  '-h', '--help'                              print usage information"
            write(output_unit, '(a, /)') "  '-V', '--verbose'                           verbose mode printing"
            write(output_unit, '(a, /)') "  '-d', '-debug', '--debug'                   run in debug mode"
            write(output_unit, '(a, /)') "  '-p', '-prefix', '--prefix'                 Attach prefix to generated filenames"
        end subroutine print_help

    end subroutine initialize


    !> Subroutine to initialize the proc file
    subroutine init_procfile()
        use mpiconf, only: idtask
        use control, only: ipr

        implicit none

        integer  :: iunit

        if (ipr .gt. 1) then
            write(proc_filename, '(a,i5.5)') "check.", idtask
            open (newunit=iunit, form='formatted', file=proc_filename)
        endif
    end subroutine init_procfile

end module contrl_file

!> Module to control the allocation and deallocation of the control modules
module m_control
contains
!> Subroutine to allocate the control modules
subroutine allocate_m_control()
    use contrldmc, only: allocate_contrldmc

    implicit none

    call allocate_contrldmc()

end subroutine allocate_m_control

!> Subroutine to deallocate the control modules
subroutine deallocate_m_control()
    use contrldmc, only: deallocate_contrldmc

    implicit none

    call deallocate_contrldmc()

end subroutine deallocate_m_control
end module
