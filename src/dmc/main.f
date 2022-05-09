      program maindmc
c Written by Claudia Filippi
      use mpiconf, only: idtask, nproc, wid
      use mpiconf, only: mpiconf_init
      use allocation_mod, only: deallocate_dmc
      use optwf_contrl, only: ioptwf
      use contr3, only: mode
      use contrl_file, only: initialize
      use mpi
      use contrl_file,    only: ounit
      use mpitimer,    only: time, time_start, time_check1, time_final
      use parser_mod,  only: parser
      use error,       only: fatal_error
      use optwf_matrix_corsamp_mod, only: optwf_matrix_corsamp
      use dmc_f_mod,   only: dmc

      implicit none

      integer :: ierr, ibranch_elec
      character*40 :: filename

      call mpi_init(ierr)
      call mpi_comm_rank(MPI_COMM_WORLD, idtask, ierr)
      call mpi_comm_size(MPI_COMM_WORLD, nproc, ierr)

      time_start = time()
      call mpiconf_init()

      call initialize()

c Open the standard output and the log file only on the master
      if(wid) then
        open(45, file='output.log', status='unknown')
      else
        close(6)
        open(6, file='/dev/null')
        open(45, file='trash.log')
      endif

      write(filename, '(''problem.'',i0)') idtask
      open(18,file=filename, status='unknown')


!     read the input from parser
      call parser()
      call MPI_BARRIER(MPI_Comm_World, ierr)


      if(mode.eq.'dmc_one_mpi2') then
        if(ioptwf.gt.0) call fatal_error('MAIN: no DMC optimization with global population')

!        call p2gtid('dmc:ibranch_elec', ibranch_elec, 0, 1)
        ! why a local variable is used to decide the following line
        if(ibranch_elec.gt.0) call fatal_error('MAIN: no DMC single-branch with global population')
      endif

      if(ioptwf.gt.0) then
       call optwf_matrix_corsamp
      else
       call dmc
      endif

      close(5)
      close(6)
      close(45)

      time_final = time()
      !write(ounit,'(a,g16.6,a)') " Total time of computation ::  ", time_final - time_start, " seconds "

      call mpi_finalize(ierr)
      call deallocate_dmc()

      stop
      end
