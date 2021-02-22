      subroutine pcm_reduce(wgsum)

      use dmc_mod, only: MWALK, MFPROD, MFPRD1, MPATH
      use mpiconf, only: idtask, nproc, wid, NPROCX
      use contr3, only: mode

      implicit real*8(a-h,o-z)



      include 'mpif.h'
      include 'pcm.h'


      if(ipcm.eq.0) return

      if(mode.eq.'dmc_one_mpi2') then

      call mpi_reduce(spcmsum,spcmcollect,1
     &,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)
      spcmsum=spcmcollect

      call mpi_reduce(vpcmsum,vpcmcollect,1
     &,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)
      vpcmsum=vpcmcollect

      call mpi_reduce(qopcm_sum,qopcmcollect,1
     &,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)
      qopcm_sum=qopcmcollect

      else

      spcmnow=spcmsum/wgsum
      vpcmnow=vpcmsum/wgsum
      qopcmnow=qopcm_sum/wgsum

      spcm2sum=spcmsum*spcmnow
      vpcm2sum=vpcmsum*vpcmnow
      qopcm2sum=qopcm_sum*qopcmnow

      call mpi_reduce(spcmsum,spcmcollect,1
     &,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

      call mpi_reduce(spcm2sum,spcm2collect,1
     &,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

      call mpi_reduce(vpcmsum,vpcmcollect,1
     &,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

      call mpi_reduce(vpcm2sum,vpcm2collect,1
     &,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

      call mpi_reduce(qopcm_sum,qopcmcollect,1
     &,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

      call mpi_reduce(qopcm2sum,qopcm2collect,1
     &,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

      call mpi_barrier(MPI_COMM_WORLD,ierr)

      if(.not.wid) return

      spcmcum=spcmcum+spcmcollect
      spcmcm2=spcmcm2+spcm2collect

      vpcmcum=vpcmcum+vpcmcollect
      vpcmcm2=vpcmcm2+vpcm2collect

      qopcm_cum=qopcm_cum+qopcmcollect
      qopcm_cm2=qopcm_cm2+qopcm2collect

      endif

      return
      end

      subroutine pcm_send(irecv,itag_s)

      use dmc_mod, only: MWALK, MFPROD, MFPRD1, MPATH
      use force_mod, only: MFORCE, MFORCE_WT_PRD, MWF

      implicit real*8(a-h,o-z)

      include 'mpif.h'
      include 'pcm.h'

      common /pcmo/ spcmo(MWALK),vpcmo(MWALK),qopcmo(MWALK),enfpcmo(MWALK,MCHS)
      common /branch/ wtgen(0:MFPRD1),ff(0:MFPRD1),eold(MWALK,MFORCE),
     &pwt(MWALK,MFORCE),wthist(MWALK,0:MFORCE_WT_PRD,MFORCE),
     &wt(MWALK),eigv,eest,wdsumo,wgdsumo,fprod,nwalk

      dimension istatus(MPI_STATUS_SIZE)

      if(ipcm.eq.0) return

      itag_s=itag_s+1
      call mpi_isend(spcmo(nwalk),1,mpi_double_precision,irecv
     &     ,itag_s,MPI_COMM_WORLD,irequest,ierr)
      itag_s=itag_s+1
      call mpi_isend(vpcmo(nwalk),1,mpi_double_precision,irecv
     &     ,itag_s,MPI_COMM_WORLD,irequest,ierr)
      itag_s=itag_s+1
      call mpi_isend(qopcmo(nwalk),1,mpi_double_precision,irecv
     &     ,itag_s,MPI_COMM_WORLD,irequest,ierr)

      return

      entry pcm_recv(isend,itag_r)

      if(ipcm.eq.0) return

      itag_r=itag_r+1
      call mpi_recv(spcmo(nwalk),1,mpi_double_precision,isend
     &     ,itag_r,MPI_COMM_WORLD,istatus,ierr)
      itag_r=itag_r+1
      call mpi_recv(vpcmo(nwalk),1,mpi_double_precision,isend
     &     ,itag_r,MPI_COMM_WORLD,istatus,ierr)
      itag_r=itag_r+1
      call mpi_recv(qopcmo(nwalk),1,mpi_double_precision,isend
     &     ,itag_r,MPI_COMM_WORLD,istatus,ierr)

      return
      end
