      program maindmc
c MPI version created by Claudia Filippi
      implicit double precision (a-h,o-z)

      include 'mpi_qmc.h'
      include 'mpif.h'

      character*40 filename

      character*12 mode
      common /contr3/ mode

      logical wid
      common /mpiconf/ idtask,nproc,wid

      call mpi_init(ierr)

      call mpi_comm_rank(MPI_COMM_WORLD,idtask,ierr)
      call mpi_comm_size(MPI_COMM_WORLD,nproc,ierr)

      if(nproc.gt.nprocx) call fatal_error('MAIN: nproc > nprocx')

      mode='dmc_one_mpi1'

      wid=(idtask.eq.0)

c Open the standard output and the log file only on the master
      if(wid) then
        open(45,file='output.log',status='unknown')
      else
        close(6)
        open(6,file='/dev/null')
        open(45,file='/dev/null')
      endif

      open(1,file='filename',status='old')
      read(1,'(a40)') filename
      close(1)

cJF closes standard input(?) on 5, probably safer to avoid 5 for files
c   or anything < 10 for that matter
      close(5)
      open(5,file=filename,status='old')

      if(idtask.le.9) then
        write(filename,'(''problem.'',i1)') idtask
       elseif(idtask.le.99) then
        write(filename,'(''problem.'',i2)') idtask
       elseif(idtask.le.999) then
        write(filename,'(''problem.'',i3)') idtask
       else
        call fatal_error('MAIN: idtask ge 1000')
      endif
      open(18,file=filename,status='unknown')

      call read_input

      call p2gtid('optwf:ioptwf',ioptwf,0,1)

      if(ioptwf.gt.0) then
       call optwf
      else
       call dmc
      endif

      close(5)
      close(6)
      close(45)

c     if(wid) then
c       close(45)
c     else
c       close(45,status='delete')
c     endif

      call mpi_finalize(ierr)

      stop
      end
