      subroutine prop_reduce(wgsum)

      use prp003, only: vprop_cm2, cc_nuc, vprop_sum, vprop_cum
      use mpiconf, only: idtask, nproc, wid, NPROCX
      use contr3, only: mode
      implicit real*8(a-h,o-z)



      include 'mpif.h'
      include 'dmc.h'
      include 'properties.h'


      dimension vp2sum(MAXPROP), vpcollect(MAXPROP), vp2collect(MAXPROP)

      if(iprop.eq.0) return

      if(mode.eq.'dmc_one_mpi2') then

       call mpi_reduce(vprop_sum,vpcollect,nprop
     &      ,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

      do 1 i=1,nprop
  1     vprop_sum(i)=vpcollect(i)

      else

      do 10 i=1,nprop
        vpnow=vprop_sum(i)/wgsum
 10     vp2sum(i)=vprop_sum(i)*vpnow

      call mpi_reduce(vprop_sum,vpcollect,nprop
     &     ,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

      call mpi_reduce(vp2sum,vp2collect,nprop
     &     ,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

      call mpi_barrier(MPI_COMM_WORLD,ierr)

      if(.not.wid) return

      do 20 i=1,nprop
       vprop_cum(i)=vprop_cum(i)+vpcollect(i)
       vprop_cm2(i)=vprop_cm2(i)+vp2collect(i)
 20   enddo

      endif

      return
      end

      subroutine prop_send(irecv,itag_s)

      use mpiconf, only: idtask, nproc, wid, NPROCX
      implicit real*8(a-h,o-z)


      include 'dmc.h'
      include 'force.h'
      include 'mpif.h'
      include 'properties.h'
      include 'prop_dmc.h'

      common /branch/ wtgen(0:MFPRD1),ff(0:MFPRD1),eold(MWALK,MFORCE),
     &pwt(MWALK,MFORCE),wthist(MWALK,0:MFORCE_WT_PRD,MFORCE),
     &wt(MWALK),eigv,eest,wdsumo,wgdsumo,fprod,nwalk

      dimension istatus(MPI_STATUS_SIZE)

      if(iprop.eq.0) return

      itag_s=itag_s+1
      call mpi_isend(vprop_old(1,nwalk),nprop,mpi_double_precision,irecv
     &     ,itag_s,MPI_COMM_WORLD,irequest,ierr)
c     itag_s=itag_s+1

      return

      entry prop_recv(isend,itag_r)

      if(iprop.eq.0) return

      itag_r=itag_r+1
      call mpi_recv(vprop_old(1,nwalk),nprop,mpi_double_precision,isend
     &     ,itag_r,MPI_COMM_WORLD,istatus,ierr)
c     itag_r=itag_r+1

      return
      end
