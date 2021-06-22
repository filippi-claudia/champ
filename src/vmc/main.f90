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

program main

    use mpi
    use mpiconf, only: idtask, nproc
    use mpiconf, only: mpiconf_init
    use contr3, only: init_control_mode
    use contrl_file, only: init_logfile, init_procfile, close_files, initialize
    use allocation_mod, only: deallocate_vmc
    use optwf_mod, only: optwf
    use mpiconf, only: wid      ! logical :: true only for mpirank=0

    implicit None
    integer :: ierr

    !> Initialize MPI
    call mpi_init(ierr)
    call mpi_comm_rank(MPI_COMM_WORLD, idtask, ierr)
    call mpi_comm_size(MPI_COMM_WORLD, nproc, ierr)

    !> init our own mpi vars
    call mpiconf_init()

    !> Mode gets reset in metrop_mov1...but loses mpi info
    !call init_control_mode('vmc_one_mpi ')           ! commented by ravindra. Not needed

    !> Initiaize output.log file.
    call init_logfile()
    call initialize()

    ! read the input from input file and other data files
    call parser()

    !> Initiaize log check.XXX files. It needs ipr flag value.
    call init_procfile()

    call MPI_BARRIER(MPI_Comm_World, ierr)
    ! ! run the the optimization
    call optwf()

    ! call close_files()
    call mpi_finalize(ierr)
    call deallocate_vmc()

end
