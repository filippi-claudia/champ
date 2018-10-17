      program maindmc
c Written by Claudia Filippi
      implicit double precision (a-h,o-z)

      include 'mpi_qmc.h'
      include 'mpif.h'

      character*20 filename

      character*12 mode
      common /contr3/ mode

      logical wid
      common /mpiconf/ idtask,nproc,wid

      call mpi_init(ierr)

      call mpi_comm_rank(MPI_COMM_WORLD,idtask,ierr)
      call mpi_comm_size(MPI_COMM_WORLD,nproc,ierr)

      if(nproc.gt.nprocx) call fatal_error('MAIN: nproc > nprocx')

      mode='dmc_one_mpi2'

      wid=(idtask.eq.0)

c Open the standard output and the log file only on the master
      if(wid) then
        open(45,file='output.log',status='unknown')
      else
        close(6)
        open(6,file='/dev/null')
        open(45,file='trash.log')
      endif

      open(1,file='filename',status='old')
      read(1,'(a20)') filename
      close(1)

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

      call MPI_ATTR_GET(MPI_COMM_WORLD,MPI_TAG_UB,ivalue,flag,ierr)
c     write(6,*) 'In main.f from MPI_ATTR_GET',ivalue
      if(ierr.ne.0) write(6,*) 'Warning:? in main.f from MPI_ATTR_GET',ierr

      call read_input

      call p2gtid('optwf:ioptwf',ioptwf,0,1)
      if(ioptwf.gt.0) call fatal_error('MAIN: no DMC optimization with global population')
      call p2gtid('dmc:ibranch_elec',ibranch_elec,0,1)
      if(ibranch_elec.gt.0) call fatal_error('MAIN: no DMC single-branch with global population')

      call dmc

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
