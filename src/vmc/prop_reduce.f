      subroutine prop_reduce

      use properties, only: MAXPROP
      use mpiconf, only: wid
      use prp000, only: iprop, nprop
      use prp003, only: vprop_cm2, vprop_cum
      use mpi

      use precision_kinds, only: dp
      implicit none

      integer :: i, ierr

      real(dp), dimension(MAXPROP) :: vpcollect
      real(dp), dimension(MAXPROP) :: vp2collect


      if(iprop.eq.0) return

      call mpi_reduce(vprop_cum,vpcollect,nprop
     &,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

      call mpi_reduce(vprop_cm2,vp2collect,nprop
     &,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

      call mpi_barrier(MPI_COMM_WORLD,ierr)

      if(wid) then
        do 10 i=1,nprop
          vprop_cum(i)=vpcollect(i)
   10     vprop_cm2(i)=vp2collect(i)
       else
        do 20 i=1,nprop
          vprop_cum(i)=0
   20     vprop_cm2(i)=0
      endif

      return
      end
