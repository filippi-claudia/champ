      subroutine optx_jas_ci_reduce
c Written by Claudia Filippi

      use optwf_parms, only: nparmj
      use dets, only: ndet
      use mix_jas_ci, only: de_o_ci, dj_de_ci, dj_o_ci, dj_oe_ci
      use optwf_contrl, only: ioptci, ioptjas
      use optwf_parms, only: nparmj
      use ci000, only: nciterm
      use method_opt, only: method
      use mpi
      use precision_kinds, only: dp

      implicit none

      integer :: i, ierr, j

      real(dp), dimension(nparmj,ndet) :: collect


      if(ioptjas.eq.0.or.ioptci.eq.0.or.method.eq.'sr_n'.or.method.eq.'lin_d') return

      call mpi_reduce(dj_o_ci,collect,nparmj*nciterm
     &     ,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

      call mpi_bcast(collect,nparmj*nciterm
     &     ,mpi_double_precision,0,MPI_COMM_WORLD,ierr)

      do 10 i=1,nparmj
        do 10 j=1,nciterm
  10     dj_o_ci(i,j)=collect(i,j)

      call mpi_reduce(dj_de_ci,collect,nparmj*nciterm
     &     ,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

      call mpi_bcast(collect,nparmj*nciterm
     &     ,mpi_double_precision,0,MPI_COMM_WORLD,ierr)

      do 20 i=1,nparmj
        do 20 j=1,nciterm
  20     dj_de_ci(i,j)=collect(i,j)

      call mpi_reduce(de_o_ci,collect,nparmj*nciterm
     &     ,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

      call mpi_bcast(collect,nparmj*nciterm
     &     ,mpi_double_precision,0,MPI_COMM_WORLD,ierr)

      do 30 i=1,nparmj
        do 30 j=1,nciterm
  30     de_o_ci(i,j)=collect(i,j)

      call mpi_reduce(dj_oe_ci,collect,nparmj*nciterm
     &     ,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

      call mpi_bcast(collect,nparmj*nciterm
     &     ,mpi_double_precision,0,MPI_COMM_WORLD,ierr)

      do 40 i=1,nparmj
        do 40 j=1,nciterm
  40     dj_oe_ci(i,j)=collect(i,j)

      call mpi_barrier(MPI_COMM_WORLD,ierr)

      return
      end
