      subroutine mmpol_reduce(wgsum)

      implicit real*8(a-h,o-z)

      include 'mpif.h'
      include 'dmc.h'
      include 'mmpol.h'

      character*12 mode
      common /contr3/ mode

      logical wid
      common /mpiconf/ idtask,nproc,wid

      if(immpol.eq.0) return

      if(mode.eq.'dmc_one_mpi2') then

      call mpi_reduce(dmmpol_sum,dmmpol_collect,1
     &,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)
      dmmpol_sum=dmmpol_collect

      call mpi_reduce(cmmpol_sum,cmmpol_collect,1
     &,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)
      cmmpol_sum=cmmpol_collect

      else

      dmmpol_now=dmmpol_sum/wgsum
      cmmpol_now=cmmpol_sum/wgsum

      dmmpol2_sum=dmmpol_sum*dmmpol_now
      cmmpol2_sum=cmmpol_sum*cmmpol_now

      call mpi_reduce(dmmpol_sum,dmmpol_collect,1
     &,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

      call mpi_reduce(dmmpol2_sum,dmmpol2_collect,1
     &,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

      call mpi_reduce(cmmpol_sum,cmmpol_collect,1
     &,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

      call mpi_reduce(cmmpol2_sum,cmmpol2_collect,1
     &,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)


      call mpi_barrier(MPI_COMM_WORLD,ierr)

      if(.not.wid) return

      dmmpol_cum=dmmpol_cum+dmmpol_collect
      dmmpol_cm2=dmmpol_cm2+dmmpol2_collect

      cmmpol_cum=cmmpol_cum+cmmpol_collect
      cmmpol_cm2=cmmpol_cm2+cmmpol2_collect

      endif

      return
      end

      subroutine mmpol_send(irecv,itag_s)

      implicit real*8(a-h,o-z)

      include 'dmc.h'
      include 'force.h'
      include 'mmpol.h'
      include 'mpif.h'

      common /mmpolo/ dmmpolo(MWALK),cmmpolo(MWALK),
     &         eeko1(MWALK,MCHMM),eeko2(MWALK,MCHMM),eeko3(MWALK,MCHMM)
      common /branch/ wtgen(0:MFPRD1),ff(0:MFPRD1),eold(MWALK,MFORCE),
     &pwt(MWALK,MFORCE),wthist(MWALK,0:MFORCE_WT_PRD,MFORCE),
     &wt(MWALK),eigv,eest,wdsumo,wgdsumo,fprod,nwalk

      common /mpiconf/ idtask,nproc,wid

      dimension istatus(MPI_STATUS_SIZE)

      if(immpol.eq.0) return

      itag_s=itag_s+1
      call mpi_isend(dmmpolo(nwalk),1,mpi_double_precision,irecv
     &     ,itag_s,MPI_COMM_WORLD,irequest,ierr)
      itag_s=itag_s+1
      call mpi_isend(cmmpolo(nwalk),1,mpi_double_precision,irecv
     &     ,itag_s,MPI_COMM_WORLD,irequest,ierr)

      return

      entry mmpol_recv(isend,itag_r)

      if(immpol.eq.0) return

      itag_r=itag_r+1
      call mpi_recv(dmmpolo(nwalk),1,mpi_double_precision,isend
     &     ,itag_r,MPI_COMM_WORLD,istatus,ierr)
      itag_r=itag_r+1
      call mpi_recv(cmmpolo(nwalk),1,mpi_double_precision,isend
     &     ,itag_r,MPI_COMM_WORLD,istatus,ierr)

      return
      end

