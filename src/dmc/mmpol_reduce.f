      module mmpol_reduce_mod
      contains
      subroutine mmpol_reduce(wgsum)

      use mpiconf, only: wid
      use mmpol_cntrl, only: immpol
      use mmpol_averages, only: cmmpol_cum, cmmpol_cm2, dmmpol_sum
      use mmpol_averages, only: cmmpol_sum, dmmpol_cum, dmmpol_cm2
      use mpi

      use precision_kinds, only: dp
      use control, only: mode
      implicit none

      integer :: ierr
      real(dp) :: cmmpol2_collect, cmmpol2_sum, cmmpol_collect, cmmpol_now, dmmpol2_collect
      real(dp) :: dmmpol2_sum, dmmpol_collect, dmmpol_now, wgsum


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

      use branch, only: nwalk
      use mmpolo, only: cmmpolo_dmc, dmmpolo_dmc
      use mmpol_cntrl, only: immpol
      use mpi

      implicit none

      integer :: ierr, irecv, irequest, isend, itag_r
      integer :: itag_s
      integer, dimension(MPI_STATUS_SIZE) :: istatus



      if(immpol.eq.0) return

      itag_s=itag_s+1
      call mpi_isend(dmmpolo_dmc(nwalk),1,mpi_double_precision,irecv
     &     ,itag_s,MPI_COMM_WORLD,irequest,ierr)
      itag_s=itag_s+1
      call mpi_isend(cmmpolo_dmc(nwalk),1,mpi_double_precision,irecv
     &     ,itag_s,MPI_COMM_WORLD,irequest,ierr)

      return

      entry mmpol_recv(isend,itag_r)

      if(immpol.eq.0) return

      itag_r=itag_r+1
      call mpi_recv(dmmpolo_dmc(nwalk),1,mpi_double_precision,isend
     &     ,itag_r,MPI_COMM_WORLD,istatus,ierr)
      itag_r=itag_r+1
      call mpi_recv(cmmpolo_dmc(nwalk),1,mpi_double_precision,isend
     &     ,itag_r,MPI_COMM_WORLD,istatus,ierr)

      return
      end

      end module
