      subroutine mmpol_reduce

      implicit real*8(a-h,o-z)

      include 'mpif.h'
      include 'dmc.h'
      include 'mmpol.h'

      logical wid
      common /mpiconf/ idtask,nproc,wid

      if(immpol.eq.0) return

      call mpi_reduce(dmmpol_sum,dmmpol_collect,1
     &,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)
      dmmpol_sum=dmmpol_collect

      call mpi_reduce(cmmpol_sum,cmmpol_collect,1
     &,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)
      cmmpol_sum=cmmpol_collect

      call mpi_barrier(MPI_COMM_WORLD,ierr)

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

