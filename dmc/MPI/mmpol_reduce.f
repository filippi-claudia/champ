      subroutine mmpol_reduce(wgsum)

      implicit real*8(a-h,o-z)

      include 'mpif.h'
      include 'dmc.h'
      include 'mmpol.h'

      logical wid
      common /mpiconf/ idtask,nproc,wid

      if(immpol.eq.0) return

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

      return
      end
