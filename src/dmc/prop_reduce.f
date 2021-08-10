      subroutine prop_reduce(wgsum)

      use prp000, only: iprop, nprop
      use prp003, only: vprop_sum, vprop_cum, vprop_cm2
      use mpiconf, only: wid
      use contr3, only: mode
      use properties, only: MAXPROP
      use mpi

      use precision_kinds, only: dp
      implicit none

      integer :: i, ierr
      real(dp) :: vpnow, wgsum
      real(dp), dimension(MAXPROP) :: vp2sum
      real(dp), dimension(MAXPROP) :: vpcollect
      real(dp), dimension(MAXPROP) :: vp2collect



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

      use prp000, only: iprop, nprop
      use prp002, only: vprop_old
      use branch, only: nwalk
      use mpi

      implicit none

      integer :: ierr, irecv, irequest, isend, itag_r
      integer :: itag_s
      integer, dimension(MPI_STATUS_SIZE) :: istatus




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
