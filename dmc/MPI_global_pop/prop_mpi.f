      subroutine prop_reduce

      implicit real*8(a-h,o-z)

      include 'dmc.h'
      include 'properties.h'
      include 'prop_dmc.h'
      include 'mpif.h'

      if(iprop.eq.0) return

      do 10 i=1,nprop
        call mpi_reduce(vprop_sum(i),vpcollect,1
     &       ,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)
 10     vprop_sum(i)=vpcollect

      return
      end

      subroutine prop_send(irecv,itag_s)

      implicit real*8(a-h,o-z)

      include 'dmc.h'
      include 'force.h'
      include 'mpif.h'
      include 'properties.h'
      include 'prop_dmc.h'

      common /branch/ wtgen(0:MFPRD1),ff(0:MFPRD1),eold(MWALK,MFORCE),
     &pwt(MWALK,MFORCE),wthist(MWALK,0:MFORCE_WT_PRD,MFORCE),
     &wt(MWALK),eigv,eest,wdsumo,wgdsumo,fprod,nwalk
      common /mpiconf/ idtask,nproc,wid

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
