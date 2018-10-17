      subroutine prop_reduce(wgsum)

      implicit real*8(a-h,o-z)

      include 'mpif.h'
      include 'properties.h'

      dimension vpcollect(MAXPROP), vp2collect(MAXPROP)


      logical wid
      common /mpiconf/ idtask,nproc,wid


      if(iprop.eq.0) return

      do 10 i=1,nprop

       vpnow=vprop_sum(i)/wgsum
       vp2sum=vprop_sum(i)*vpnow


       call mpi_reduce(vprop_sum(i),vpcollect(i),1
     &      ,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

       call mpi_reduce(vp2sum,vp2collect(i),1
     &      ,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

 10   enddo

      call mpi_barrier(MPI_COMM_WORLD,ierr)
      if(.not.wid) return

      do 20 i=1,nprop
       vprop_cum(i)=vprop_cum(i)+vpcollect(i)
       vprop_cm2(i)=vprop_cm2(i)+vp2collect(i)
 20   enddo

      return
      end
