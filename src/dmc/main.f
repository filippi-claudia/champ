      program maindmc
c Written by Claudia Filippi
      use mpiconf, only: idtask, nproc, wid
      use mpiconf, only: mpiconf_init
      use allocation_mod, only: deallocate_dmc
      use optwf_contrl, only: ioptwf
      use contr3, only: mode
      use mpi

      implicit none

      integer :: ierr
      character*40 :: filename

      call mpi_init(ierr)
      call mpi_comm_rank(MPI_COMM_WORLD, idtask, ierr)
      call mpi_comm_size(MPI_COMM_WORLD, nproc, ierr)

      call mpiconf_init()


c Open the standard output and the log file only on the master
      if(wid) then
        open(45, file='output.log', status='unknown')
      else
        close(6)
        open(6, file='/dev/null')
        open(45, file='trash.log')
      endif

      if(idtask.le.9) then
        write(filename, '(''problem.'',i1)') idtask
       elseif(idtask.le.99) then
        write(filename, '(''problem.'',i2)') idtask
       elseif(idtask.le.999) then
        write(filename, '(''problem.'',i3)') idtask
       else
        call fatal_error('MAIN: idtask ge 1000')
      endif
      open(18,file=filename, status='unknown')

!      BUG:: Ravindra. Following line needs a replacement
!      call read_input

      if(ioptwf.gt.0) then
       call optwf_matrix_corsamp
      else
       call dmc
      endif

      close(5)
      close(6)
      close(45)

      call mpi_finalize(ierr)
      call deallocate_dmc()

      stop
      end
