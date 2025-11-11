!> @brief Module for global program control and execution mode.
!> @author CHAMP developers
!> @date 2022
!>
!> @details This module manages the overall execution mode of the program and
!> output verbosity level. The mode determines whether the program runs VMC,
!> DMC, optimization, or other calculation types. The verbosity level controls
!> the amount of diagnostic output.
!>
!> Key parameters:
!> - mode: Execution mode string (e.g., 'vmc', 'dmc')
!> - ipr: Output verbosity level (higher values = more output)
!>
!> @note The mode is typically set during input parsing and remains constant
!> throughout the calculation.
module control

    implicit none

    !> Execution mode of the program (e.g., 'vmc', 'dmc')
    character(len=12) :: mode

    !> Output verbosity level: -1=minimal, 0=normal, 1+=verbose with debug info
    integer  :: ipr

    private
    public :: mode, ipr, init_control_mode
    save
    contains

    !> Initializes the program execution mode.
    !>
    !> @details Sets the mode string that controls which calculation type
    !> will be executed (VMC, DMC, etc.).
    !>
    !> @param[in] str_mode String specifying the execution mode
    subroutine init_control_mode(str_mode)
        implicit None
        character(12), intent(IN) :: str_mode
        mode = str_mode
    end subroutine init_control_mode

end module control

!> @brief Module for VMC (Variational Monte Carlo) calculation control parameters.
!> @author CHAMP developers
!> @date 2022
!>
!> @details This module stores control parameters specific to VMC calculations,
!> including block structure, configuration management, restart options, and
!> output control. These parameters determine how the VMC sampling is performed.
!>
!> Key parameter groups:
!> - Block structure: nblk, nblkeq, nstep (equilibration and production)
!> - Configuration management: nconf, nconf_new, isite (walker setup)
!> - I/O control: idump, irstar (checkpoint and restart)
!> - Special modes: icharged_atom, nblk_ci (advanced calculations)
!>
!> @note These parameters are typically read from the vmc.inp input file.
module control_vmc

    !> Flag to write walker configurations to dump file at end (0=no, 1=yes)
    integer :: vmc_idump

    !> Flag to restart from previous VMC run (0=new run, 1=restart)
    integer :: vmc_irstar

    !> Configuration initialization mode (0=generate, 1=read from file)
    integer :: vmc_isite

    !> Total number of walker configurations for VMC sampling
    integer :: vmc_nconf

    !> Number of production blocks (excluding equilibration)
    integer :: vmc_nblk

    !> Maximum number of blocks allowed in optimization
    integer :: vmc_nblk_max

    !> Number of equilibration blocks (statistics discarded)
    integer :: vmc_nblkeq

    !> Number of new configurations to generate 
    integer :: vmc_nconf_new

    !> Number of Monte Carlo steps per block
    integer :: vmc_nstep

    !> Index of charged atom for special calculations (e.g., ionization)
    integer :: vmc_icharged_atom

    !> Number of blocks for CI (Configuration Interaction) coefficient optimization
    integer :: vmc_nblk_ci

    private
    public :: vmc_idump, vmc_irstar, vmc_isite, vmc_nconf, vmc_nblk, vmc_nblk_max
    public :: vmc_nblkeq, vmc_nconf_new, vmc_nstep, vmc_icharged_atom, vmc_nblk_ci
    save
end module control_vmc


!> @brief Module for DMC (Diffusion Monte Carlo) calculation control parameters.
!> @author CHAMP developers
!> @date 2022
!>
!> @details This module stores control parameters specific to DMC calculations,
!> which provide more accurate ground-state energies than VMC through fixed-node
!> diffusion quantum Monte Carlo. Parameters control walker population, branching,
!> equilibration, and output management.
!>
!> Key parameter groups:
!> - Block structure: nblk, nblkeq, nstep (equilibration and sampling)
!> - Walker population: nconf, nconf_new (dynamic population control)
!> - I/O control: idump, irstar (checkpoint and restart)
!> - Configuration source: isite (walker initialization)
!>
!> @note These parameters are typically read from the dmc.inp input file.
!> @note DMC walker populations may fluctuate due to branching/merging.
module control_dmc

    !> Flag to write walker configurations to dump file at end (0=no, 1=yes)
    integer :: dmc_idump

    !> Flag to restart from previous DMC run (0=new run, 1=restart)
    integer :: dmc_irstar

    !> Configuration initialization mode (0=from VMC, 1=read from file)
    integer :: dmc_isite

    !> Target number of walker configurations (subject to branching)
    integer :: dmc_nconf

    !> Number of production blocks for DMC statistics
    integer :: dmc_nblk

    !> Number of equilibration blocks (allows population stabilization)
    integer :: dmc_nblkeq

    !> Number of new configurations after population equilibration
    integer :: dmc_nconf_new

    !> Number of DMC steps per block
    integer :: dmc_nstep

    private
    public :: dmc_idump, dmc_irstar, dmc_isite, dmc_nconf, dmc_nblk
    public :: dmc_nblkeq, dmc_nconf_new, dmc_nstep
    save
end module control_dmc

!> @brief Module for periodic boundary condition control parameters.
!> @author CHAMP developers
!> @date 2022
!>
!> @details This module controls the use of periodic boundary conditions (PBC)
!> in simulations of extended systems (crystals, solids). The iperiodic flag
!> enables PBC treatment, while ibasis determines how basis functions handle
!> periodicity.
!>
!> Typical values:
!> - iperiodic: 0=molecular/finite system, 1=periodic system (crystal)
!> - ibasis: Basis set type appropriate for PBC
!>
!> @note Required for solid-state calculations with periodic symmetry.
module contrl_per

    implicit none

    !> Flag for periodic boundary conditions (0=molecular, 1=periodic/crystal)
    integer :: iperiodic

    !> Basis set type for periodic systems
    integer :: ibasis

    private
    public   :: iperiodic, ibasis
    save
end module contrl_per

!> @brief Module for advanced DMC algorithm control parameters and time step arrays.
!> @author CHAMP developers
!> @date 2022
!>
!> @details This module contains detailed control parameters for the DMC algorithm,
!> including time step management, branching control, acceptance/rejection tracking,
!> and force calculations. The taueff array allows different effective time steps
!> for different force components (multiple geometry optimization).
!>
!> Key parameter categories:
!> - Time steps: tau (base), taueff (force-dependent), tautot (total), rttau (sqrt(tau))
!> - Branching control: ibranching_c, icut_br, icut_e, limit_wt_dmc
!> - Algorithm flags: idmc (DMC type), ipq, itau_eff
!> - Statistics: iacc_rej (accept/reject), icross (node crossing), icuspg (cusp)
!> - Forces: nfprod (force products), idiv_v (velocity divergence)
!>
!> @note The taueff array is allocated dynamically based on MFORCE.
module contrldmc
      use multiple_geo, only: MFORCE
      use precision_kinds, only: dp

    implicit none

    !> Flag for acceptance/rejection statistics (0=off, 1=collect stats)
    integer :: iacc_rej

    !> Node crossing control flag (affects node avoidance algorithms)
    integer :: icross

    !> Cusp correction flag for Gaussian basis sets?
    integer :: icuspg

    !> Branching cutoff control (limits max weight for walker splitting)
    integer :: icut_br

    !> Energy cutoff flag (removes high-energy walkers)
    integer :: icut_e

    !> Velocity divergence control for velocity Verlet integration?
    integer :: idiv_v

    !> DMC algorithm type 
    integer :: idmc

    !> Population control method?
    integer :: ipq

    !> Effective time step flag
    integer :: itau_eff

    !> Number of force products to compute for multiple geometries
    integer :: nfprod

    !> Maximum walker weight limit to prevent divergence
    integer :: limit_wt_dmc

    !> Branching parameter C for weight w = exp(-C*(E_L - E_ref)*tau)
    real(dp) :: ibranching_c

    !> Square root of time step (for diffusion step variance)
    real(dp) :: rttau

    !> Base DMC time step (a.u., controls diffusion and bias)
    real(dp) :: tau

    !> Effective time steps for force calculations (MFORCE), allocated dynamically
    real(dp), dimension(:), allocatable :: taueff

    !> Total accumulated time step for multiple time step algorithms
    real(dp) :: tautot

    private
    public :: iacc_rej, icross, icuspg, icut_br, icut_e, ibranching_c, idiv_v, idmc, ipq, itau_eff, nfprod, rttau, tau, taueff, tautot, limit_wt_dmc
    public :: allocate_contrldmc, deallocate_contrldmc
    save
contains
    !> @brief Allocate effective time step array
    !> @details Allocates taueff(MFORCE) for storing effective time steps.
    subroutine allocate_contrldmc()
      use multiple_geo, only: MFORCE
        if (.not. allocated(taueff)) allocate (taueff(MFORCE))
    end subroutine allocate_contrldmc

    !> @brief Deallocate effective time step array.
    !> @details Frees memory for taueff array.
    subroutine deallocate_contrldmc()
        if (allocated(taueff)) deallocate (taueff)
    end subroutine deallocate_contrldmc

end module contrldmc


!> @brief Module for file handling, command-line parsing, and I/O unit management.
!> @author CHAMP developers
!> @date 2022
!>
!> @details This module manages all file I/O operations for the CHAMP program,
!> including input/output file names, unit numbers, log files, proc files, and
!> command-line argument processing. It provides subroutines for initialization,
!> command-line parsing, and proper file closure.
!>
!> Key functionality:
!> - Command-line parsing: -i (input), -o (output), -e (error)
!> - MPI-aware I/O: master writes to files, workers redirect to /dev/null
!> - Log files: output.log (master only), proc files (per-process debug)
!> - TREXIO backend support: when compiled with TREXIO library
!>
!> @note Master process (MPI rank 0) handles output; workers are silenced.
!> @note Unit 5=stdin, 6=stdout, 45=log file, ounit/errunit=custom files.
module contrl_file
#if defined(TREXIO_FOUND)
      use trexio,  only: trexio_back_end_t
#endif
    implicit none

    !> Name of the log file (output.log on master, /dev/null on workers)
    character(20) :: log_filename

    !> Name of the per-process check file (check.XXXXX where XXXXX is MPI rank)
    character(20) :: proc_filename

    !> Name of the main input file (from -i, .in/.inp/.dat extension)
    character(80) :: file_input

    !> Name of the main output file (from -o, master writes, workers use /dev/null)
    character(80) :: file_output

    !> Name of the error log file (from -e error)
    character(80) :: file_error

    !> Unit number for input file
    integer       :: iunit

    !> Unit number for output file (dynamically assigned with newunit)
    integer       :: ounit

    !> Unit number for error file (dynamically assigned with newunit)
    integer       :: errunit

#if defined(TREXIO_FOUND)
    !> TREXIO backend type for quantum chemistry data I/O
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

    !> @brief Close standard input, output, and log files.
    !> @details Closes units 5 (stdin), 6 (stdout), and 45 (log file on master).
    !> Called at program termination.
    subroutine close_files()
      use mpiconf, only: wid
        close (5)
        close (6)
        if (wid) close (45)
    end subroutine close_files

    !> @brief Initialize the log file for master process.
    !> @details Opens output.log on master (wid=true), redirects stdout to /dev/null on workers.
    !> Ensures only master process generates visible output.
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

    !> @brief Parse command-line arguments and initialize file I/O.
    !> @author Ravindra Shinde
    !> @details Parses command-line options for input/output/error files, version,
    !> help, verbose, and debug modes. Sets up file units and redirects worker output.
    !>
    !> Recognized options:
    !> - -i, --input: Input file (.in/.inp/.dat)
    !> - -o, --output: Output file (master only)
    !> - -e, --error: Error log file
    !> - -v, --version: Print version and exit
    !> - -h, --help: Print help and exit
    !> - -V, --verbose: Enable verbose output
    !> - -d, --debug: Enable debug mode
    !> - -p, --prefix: Add prefix to output files
    !>
    !> @note Workers redirect output to /dev/null to avoid cluttered output.
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


    !> @brief Initialize per-process debug file (proc file).
    !> @details Creates check.XXXXX file (XXXXX = MPI rank) for per-process debugging
    !> when verbosity level ipr > 1. Each MPI process writes to its own proc file.
    !>
    !> @note Only created when ipr > 1 (verbose debug mode).
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

!> @brief Master module for control parameter allocation and deallocation.
!> @author CHAMP developers
!> @date 2022
!>
!> @details This module provides top-level allocation and deallocation routines
!> for all control modules. It coordinates memory management for control parameters
!> by calling allocation routines from sub-modules.
!>
!> Current sub-modules managed:
!> - contrldmc: DMC control arrays (taueff)
!>
!> @note Add new allocation calls here when introducing new control modules.
module m_control
contains
!> @brief Allocate all control module arrays.
!> @details Calls allocation routines for all control sub-modules.
!> Currently allocates contrldmc arrays (taueff for force-dependent time steps).
subroutine allocate_m_control()
    use contrldmc, only: allocate_contrldmc

    implicit none

    call allocate_contrldmc()

end subroutine allocate_m_control

!> @brief Deallocate all control module arrays.
!> @details Calls deallocation routines for all control sub-modules.
!> Frees memory for contrldmc arrays.
subroutine deallocate_m_control()
    use contrldmc, only: deallocate_contrldmc

    implicit none

    call deallocate_contrldmc()

end subroutine deallocate_m_control
end module
