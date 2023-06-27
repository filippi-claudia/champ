module prop_reduce_mod
contains
      subroutine prop_reduce

      use mpi
      use mpiconf, only: wid
      use precision_kinds, only: dp
      use properties, only: MAXPROP
      use prp000,  only: iprop,nprop
      use prp003,  only: vprop_cm2,vprop_cum

      implicit none

      integer :: i, ierr

      real(dp), dimension(MAXPROP) :: vpcollect
      real(dp), dimension(MAXPROP) :: vp2collect


      if(iprop.eq.0) return

      call mpi_reduce(vprop_cum,vpcollect,nprop &
      ,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

      call mpi_reduce(vprop_cm2,vp2collect,nprop &
      ,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

      call mpi_barrier(MPI_COMM_WORLD,ierr)

      if(wid) then
        do i=1,nprop
          vprop_cum(i)=vpcollect(i)
          vprop_cm2(i)=vp2collect(i)
        enddo
       else
        do i=1,nprop
          vprop_cum(i)=0
          vprop_cm2(i)=0
        enddo
      endif

      return
      end
end module
