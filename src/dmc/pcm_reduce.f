      module pcm_reduce_mod
      contains
      subroutine pcm_reduce(wgsum)

      use control, only: mode
      use mpi
      use mpiconf, only: wid
      use pcm_averages, only: qopcm_cm2,qopcm_cum,qopcm_sum,spcmcm2
      use pcm_averages, only: spcmcum,spcmsum,vpcmcm2,vpcmcum,vpcmsum
      use pcm_cntrl, only: ipcm
      use precision_kinds, only: dp

      implicit none

      integer :: ierr
      real(dp) :: qopcm2collect, qopcm2sum, qopcmcollect, qopcmnow, spcm2collect
      real(dp) :: spcm2sum, spcmcollect, spcmnow, vpcm2collect
      real(dp) :: vpcm2sum, vpcmcollect, vpcmnow, wgsum

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

      use branch,  only: nwalk
      use mpi
      use pcm_cntrl, only: ipcm
      use pcmo,    only: qopcmo_dmc,spcmo_dmc,vpcmo_dmc

      implicit none

      integer :: ierr, irecv, irequest, isend, itag_r
      integer :: itag_s
      integer, dimension(MPI_STATUS_SIZE) :: istatus




      if(ipcm.eq.0) return

      itag_s=itag_s+1
      call mpi_isend(spcmo_dmc(nwalk),1,mpi_double_precision,irecv
     &     ,itag_s,MPI_COMM_WORLD,irequest,ierr)
      itag_s=itag_s+1
      call mpi_isend(vpcmo_dmc(nwalk),1,mpi_double_precision,irecv
     &     ,itag_s,MPI_COMM_WORLD,irequest,ierr)
      itag_s=itag_s+1
      call mpi_isend(qopcmo_dmc(nwalk),1,mpi_double_precision,irecv
     &     ,itag_s,MPI_COMM_WORLD,irequest,ierr)

      return

      entry pcm_recv(isend,itag_r)

      if(ipcm.eq.0) return

      itag_r=itag_r+1
      call mpi_recv(spcmo_dmc(nwalk),1,mpi_double_precision,isend
     &     ,itag_r,MPI_COMM_WORLD,istatus,ierr)
      itag_r=itag_r+1
      call mpi_recv(vpcmo_dmc(nwalk),1,mpi_double_precision,isend
     &     ,itag_r,MPI_COMM_WORLD,istatus,ierr)
      itag_r=itag_r+1
      call mpi_recv(qopcmo_dmc(nwalk),1,mpi_double_precision,isend
     &     ,itag_r,MPI_COMM_WORLD,istatus,ierr)

      return
      end
      end module
