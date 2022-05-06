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
    use contrl_file, only: init_logfile, init_procfile, close_files 
    use md_mass
    use md_var
    use allocation_mod, only: deallocate_vmc
    use optwf_mod, only: optwf

    implicit None
    integer :: ierr
    integer :: md_step

    !> Initialize MPI
    call mpi_init(ierr)
    call mpi_comm_rank(MPI_COMM_WORLD, idtask, ierr)
    call mpi_comm_size(MPI_COMM_WORLD, nproc, ierr)

    !> init our own mpi vars
    call mpiconf_init()

    !> Mode gets reset in metrop_mov1...but loses mpi info
    call init_control_mode('vmc_one_mpi ')

    !> Initiaize output.log file.
    call init_logfile() 

    ! read the input
    call read_input()

    !> Initiaize log check.XXX files. It needs ipr flag value.  
    call init_procfile() 

    ! run the the optimization
    call p2gtid('mdyn:md_step',md_step,0,1) !number of iterations
    write(6,*) "MD STEP", md_step
    if(md_step.ge.1) then
       call allocate_mass()
       call allocate_md()
       call md
       call deallocate_mass()
       call deallocate_md()
    else
       call optwf()
    endif
    call close_files()
    call mpi_finalize(ierr)
    call deallocate_vmc()

end
