!------------------------------------------------------------------------------
!        Main Program of CHAMP
!------------------------------------------------------------------------------
!> @author
!> Claudia Filippi
!
! DESCRIPTION:
!> Read the input file and run either a simple sampling or optimize the
! wave function using the method specified in the input
!
! URL           : https://github.com/filippi-claudia/champ
!---------------------------------------------------------------------------

module main_mod

contains
      use allocation_mod, only: deallocate_vmc
      use contrl_file, only: close_files,init_logfile,init_procfile
      use contrl_file, only: initialize,ounit
      use control, only: init_control_mode
      use mpi
      use mpiconf, only: idtask,mpiconf_init,nproc,wid
      use mpitimer, only: elapsed_time,time,time_check1,time_final
      use mpitimer, only: time_start
      use optwf_mod, only: optwf
      use parser_mod, only: parser
      use precision_kinds, only: dp

    implicit None
    integer :: ierr


    !> Initialize MPI
    call mpi_init(ierr)
    call mpi_comm_rank(MPI_COMM_WORLD, idtask, ierr)
    call mpi_comm_size(MPI_COMM_WORLD, nproc, ierr)


    time_start = time()
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

    call MPI_BARRIER(MPI_Comm_World, ierr)
    call elapsed_time("MPI Barrier before optwf : ")

end subroutine

subroutine finalize_main()
    use mpi_f08
    use allocation_mod, only: deallocate_vmc
    use contrl_file,    only: ounit
    use mpitimer,    only: time, time_start, time_final

    implicit none
    integer :: ierr
    ! call close_files()
    time_final = time()

    write(ounit,'(a,g16.6,a)') " REAL TIME (Total) of computation ::  ", time_final - time_start, " seconds "

    call mpi_finalize(ierr)
    call deallocate_vmc()
end subroutine

end module


program main
    use main_mod, only: initialize_main, finalize_main
    use optwf_mod, only: optwf

    implicit None
    integer :: ierr

    call initialize_main()

    ! ! run the the optimization
    call optwf()

    call finalize_main()
end
