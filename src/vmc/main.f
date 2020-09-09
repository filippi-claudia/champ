      program main
c Written by Claudia Filippi

      use mpi_qmc, only: NPROCX
      use mpiconf, only: idtask, nproc, wid
      use contr3, only: mode
      use const, only: pi, hb, etrial, delta, deltai, fbias, nelec, imetro, ipr

      implicit real*8(a-h,o-z)

      character method*20
      character*40 filename

c mpif.h is system, mpi_qmc.h is ours
      include 'mpif.h'

      call mpi_init(ierr)

      call mpi_comm_rank(MPI_COMM_WORLD,idtask,ierr)
      call mpi_comm_size(MPI_COMM_WORLD,nproc,ierr)

      if(nproc.gt.NPROCX) stop 'nproc>NPROCX in main'

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
      call p2gtid('optwf:idl_flag',idl_flag,0,1)
      call p2gtid('optwf:ilbfgs_flag',ilbfgs_flag,0,1)

      if(ioptwf.gt.0) then
        if(idl_flag.gt.0) then
          call optwf_dl
        elseif(ilbfgs_flag.gt.0) then
          call optwf_olbfgs
        elseif(method.eq.'sr_n') then
          call optwf_sr
        elseif(method.eq.'lin_d') then
          call optwf_lin_d
        elseif(method.eq.'mix_n') then
          call optwf_mix
        else
         call optwf_matrix_corsamp
       endif
      else
       call vmc
      endif

      close(5)
      close(6)
      close(45)

      call mpi_finalize(ierr)

      end
