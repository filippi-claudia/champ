      program maindmc
c Written by Claudia Filippi
      use allocation_mod, only: deallocate_dmc
      use contrl_file, only: initialize,ounit
      use control, only: mode, ipr
      use dmc_f_mod, only: dmc
      use error,   only: fatal_error
      use mpi
      use mpiconf, only: idtask,mpiconf_init,nproc,wid
      use mpitimer, only: elapsed_time,time,time_check1,time_final
      use mpitimer, only: time_start
      use optwf_control, only: ioptwf,method
      use optwf_matrix_corsamp_mod, only: optwf_matrix_corsamp
      use optwf_sr_mod, only: optwf_sr
      use parser_mod, only: parser

      implicit none

      integer :: ierr, ibranch_elec
      character*40 :: filename

      call mpi_init(ierr)
      call mpi_comm_rank(MPI_COMM_WORLD, idtask, ierr)
      call mpi_comm_size(MPI_COMM_WORLD, nproc, ierr)

      time_start = time()
      time_check1 = time_start
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

      if(ipr.gt.5) then
            write(filename, '(''problem.'',i0)') idtask
            open(18,file=filename, status='unknown')
      endif


!     read the input from parser
      call elapsed_time("MPI initializations : ")
      call parser()
      call elapsed_time("Parsing all the files : ")
      call MPI_BARRIER(MPI_Comm_World, ierr)
      call elapsed_time("MPI Barrier before main DMC : ")


      if(mode.eq.'dmc_one_mpi2') then
        if(ioptwf.gt.0) call fatal_error('MAIN: no DMC optimization with global population')

        ! why a local variable is used to decide the following line
        if(ibranch_elec.gt.0) call fatal_error('MAIN: no DMC single-branch with global population')
      endif

        if (ioptwf .gt. 0) then
            if (method .eq. 'sr_n') then
                call optwf_sr
            else
                call fatal_error('MAIN: Only SR optimization')
            endif
       else
            call dmc
      endif

      close(5)
      close(6)
      close(45)

      time_final = time()
      write(ounit,'(a,g16.6,a)') " REAL TIME (Total) of computation ::  ", time_final - time_start, " seconds "

      call mpi_finalize(ierr)
      call deallocate_dmc()

      stop
      end
