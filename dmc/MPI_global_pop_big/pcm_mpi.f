      subroutine pcm_reduce

      implicit real*8(a-h,o-z)

      include 'mpif.h'
      include 'dmc.h'
      include 'pcm.h'

      logical wid
      common /mpiconf/ idtask,nproc,wid

      if(ipcm.eq.0) return

      call mpi_reduce(spcmsum,spcmcollect,1
     &,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)
      spcmsum=spcmcollect

      call mpi_reduce(vpcmsum,vpcmcollect,1
     &,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)
      vpcmsum=vpcmcollect

      call mpi_reduce(qopcm_sum,qopcmcollect,1
     &,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)
      qopcm_sum=qopcmcollect

      call mpi_barrier(MPI_COMM_WORLD,ierr)

      return
      end


      subroutine pcm_send(irecv,itag_s)

      implicit real*8(a-h,o-z)

      include 'dmc.h'
      include 'force.h'
      include 'mpif.h'
      include 'pcm.h'

      common /pcmo/ spcmo(MWALK),vpcmo(MWALK),qopcmo(MWALK),enfpcmo(MWALK,MCHS)
      common /mpiconf/ idtask,nproc,wid

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
