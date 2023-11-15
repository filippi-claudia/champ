      program maindmc
! Written by Claudia Filippi
      use allocation_mod, only: deallocate_dmc
      use contrl_file, only: initialize,ounit, close_files, init_logfile
      use control, only: mode, ipr
      use dmc_f_mod, only: dmc
      use error,   only: fatal_error
      use mpi
      use mpiconf, only: idtask,mpiconf_init,nproc,wid
      use mpitimer, only: elapsed_time,time,time_check1,time_final
      use mpitimer, only: time_start
      use optwf_control, only: ioptwf, method
      use optwf_matrix_corsamp_mod, only: optwf_matrix_corsamp
      use optwf_sr_mod, only: optwf_sr
      use parser_mod, only: parser

      implicit none

      integer :: ierr
      character(len=40) :: filename

      call mpi_init(ierr)
      call mpi_comm_rank(MPI_COMM_WORLD, idtask, ierr)
      call mpi_comm_size(MPI_COMM_WORLD, nproc, ierr)

      time_start = time()
      time_check1 = time_start
      call mpiconf_init()

      call init_logfile()
      call initialize()

!     read the input from parser
      call elapsed_time("MPI initializations : ")
      call parser()
      call elapsed_time("Parsing all the files : ")

      if(ipr.gt.5) then
            write(filename, '(''problem.'',i0)') idtask
            open(18,file=filename, status='unknown')
      endif

      if(ioptwf.gt.0) then
            if(mode.eq.'dmc_one_mpi2') call fatal_error('MAIN: no DMC optimization with global population')
         if (method .eq. 'sr_n') then
            call optwf_sr
         else
            call fatal_error('MAIN: Only SR optimization')
         endif
      else
        call dmc
      endif

      time_final = time()
      write(ounit,'(a,g16.6,a)') " REAL TIME (Total) of computation ::  ", time_final - time_start, " seconds "

      call close_files()
      call mpi_finalize(ierr)
      call deallocate_dmc()

      stop
      end
