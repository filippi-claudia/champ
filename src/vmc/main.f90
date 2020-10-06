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
    use contrl_file, only: init_files, close_files

    implicit None
    integer :: ierr

    !> Initialize MPI
    call mpi_init(ierr)
    call mpi_comm_rank(MPI_COMM_WORLD, idtask, ierr)
    call mpi_comm_size(MPI_COMM_WORLD, nproc, ierr)

    !> init our own mpi vars
    call mpiconf_init()

    !> Mode gets reset in metrop_mov1...but loses mpi info
    call init_control_mode('vmc_one_mpi ')

    !> Initiaize the different log/out files
    call init_files()

    ! read the input
    call read_input()

    ! run the the optimization
    call optwf()

    call close_files()
    call mpi_finalize(ierr)

end
