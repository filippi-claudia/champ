      subroutine pcm_reduce(wgsum)

      implicit real*8(a-h,o-z)

      include 'mpif.h'
      include 'dmc.h'
      include 'pcm.h'

      logical wid
      common /mpiconf/ idtask,nproc,wid

      if(ipcm.eq.0) return

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

      return
      end
