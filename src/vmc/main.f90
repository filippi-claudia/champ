!>------------------------------------------------------------------------------
!>        Main Program of CHAMP
!>------------------------------------------------------------------------------
!>
!> DESCRIPTION:
!> Read the input file and run either a simple sampling or optimize the
!> wave function using the method specified in the input
!>
!> URL           : https://github.com/filippi-claudia/champ
!>---------------------------------------------------------------------------
!> @author Claudia Filippi
module main_mod
    contains
subroutine initialize_main
    !> Initialize MPI
      use contrl_file, only: init_logfile,init_procfile
      use contrl_file, only: initialize
      use mpi
      use mpiconf, only: idtask,mpiconf_init,nproc
      use mpitimer, only: elapsed_time,mpi_time,time_check1
      use mpitimer, only: time_start
      use parser_mod, only: parser

    implicit none
    integer :: ierr


    call mpi_init(ierr)
    call mpi_comm_rank(MPI_COMM_WORLD, idtask, ierr)
    call mpi_comm_size(MPI_COMM_WORLD, nproc, ierr)


    time_start = mpi_time()
    time_check1 = time_start
    !> init our own mpi vars
    call mpiconf_init()

    !> Mode gets reset in metrop_mov1...but loses mpi info
    !call init_control_mode('vmc_one_mpi ')           ! commented by ravindra. Not needed

    !> Initiaize output.log file.
    call init_logfile()
    call initialize()

    ! read the input from input file and other data files
    call elapsed_time("MPI initializations : ")
    call parser()
    call elapsed_time("Parsing all the files : ")

    !> Initiaize log check.XXX files. It needs ipr flag value.
    call init_procfile()

    call elapsed_time("Before VMC driver : ")

end subroutine

subroutine finalize_main()
    !> Finalize MPI
    use mpi
    use allocation_mod, only: deallocate_vmc
    use contrl_file,    only: close_files,ounit
    use mpitimer,    only: mpi_time, time_start, time_final

    implicit none
    integer :: ierr
    time_final = mpi_time()

    write(ounit,'(a,g16.6,a)') " REAL TIME (Total) of computation ::  ", time_final - time_start, " seconds "

    call close_files()
    call mpi_finalize(ierr)
    call deallocate_vmc()
end subroutine

end module


program main
    !> Main program of CHAMP
    use main_mod, only: initialize_main, finalize_main
    use vmc_driver_mod, only: run_vmc_driver

    implicit none

    !> Initialize MPI
    call initialize_main()

    !> Run plain VMC or wave-function optimization.
    call run_vmc_driver()

    !> Finalize MPI
    call finalize_main()
end
