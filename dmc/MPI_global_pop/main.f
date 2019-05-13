      program maindmc
c MPI version created by Claudia Filippi
      implicit double precision (a-h,o-z)

      include 'mpi_qmc.h'
      include 'mpif.h'
      include 'vmc.h'

      character*40 filename

      character*12 mode
      common /contr3/ mode

      logical wid
      common /mpiconf/ idtask,nproc,wid
      common /mpitype/ jas_type1,jas_type2

      dimension iblocklen(MELEC),idispl(MELEC)

      mode='dmc_one_mpi3'

      call mpi_init(ierr)

      call mpi_comm_rank(MPI_COMM_WORLD,idtask,ierr)
      call mpi_comm_size(MPI_COMM_WORLD,nproc,ierr)

      if(nproc.gt.nprocx) call fatal_error('MAIN: nproc > nprocx')

      wid=(idtask.eq.0)

c To minimize the length of what is passed, use jas_typ to pick certain
c elements of matrix
      do 1 i=1,nelec
        iblocklen(i)=nelec-(i-1)
   1    idispl(i)=MELEC*(i-1)+(i-1)
      call mpi_type_indexed(nelec,iblocklen,idispl,mpi_double_precision
     &,jas_type1,ierr)
      call mpi_type_commit(jas_type1,ierr)

      do 2 i=1,nelec
        iblocklen(i)=3*nelec
   2    idispl(i)=3*MELEC*(i-1)
      call mpi_type_indexed(nelec,iblocklen,idispl,mpi_double_precision
     &,jas_type2,ierr)
      call mpi_type_commit(jas_type2,ierr)

c Open the standard output and the log file only on the master
      if(wid) then
        open(45,file='output.log',status='unknown')
      else
        close(6)
        open(6,file='/dev/null')
        open(45,file='/dev/null')
      endif

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
