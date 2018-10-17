      program main
c Written by Claudia Filippi

      implicit real*8(a-h,o-z)
      character*12 mode
      character method*20

c mpif.h is system, mpi_qmc.h is ours
      include 'mpi_qmc.h'
      include 'mpif.h'

      character*40 filename

      common /const/ pi,hb,etrial,delta,deltai,fbias,nelec,imetro,ipr

      common /contr3/ mode
      logical wid
      common /mpiconf/ idtask,nproc,wid

      call mpi_init(ierr)

      call mpi_comm_rank(MPI_COMM_WORLD,idtask,ierr)
      call mpi_comm_size(MPI_COMM_WORLD,nproc,ierr)

      if(nproc.gt.nprocx) stop 'nproc>nprocx in main'

      wid=(idtask.eq.0)

c Mode gets reset in metrop_mov1... but loses mpi info
      mode='vmc_one_mpi '

c Open the standard output and the log file only on the master
      if(wid) then
        open(45,file='output.log',status='unknown')
      else
c     endif
c     if(idtask.ne.1) then
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

      call read_input

      if(ipr.gt.1) then
        if(idtask.lt.10) then
          write(filename,'(i1)') idtask
         elseif(idtask.lt.100) then
          write(filename,'(i2)') idtask
         elseif(idtask.lt.1000) then
          write(filename,'(i3)') idtask
         else
          write(filename,'(i4)') idtask
        endif
        filename='check.'//filename(1:index(filename,' ')-1)
        open(unit=88,form='formatted',file=filename)
      endif


      call p2gtid('optwf:ioptwf',ioptwf,0,1)
      call p2gtad('optwf:method',method,'linear',1)
      call p2gtad('optwf:method',method,'linear',1)
      call p2gtid('optwf:idl_flag',idl_flag,0,1)

      if(ioptwf.gt.0) then
        if(idl_flag.gt.0) then
          call dl_optwf
        else if(method.eq.'sr_n'.or.method.eq.'lin_d'.or.method.eq.'mix_n') then
          call sr_optwf
        else
         call optwf
       endif
      else
       call vmc
      endif

      close(5)
      close(6)
      close(45)

      call mpi_finalize(ierr)

      end
