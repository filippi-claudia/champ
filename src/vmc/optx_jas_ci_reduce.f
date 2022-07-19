      module optx_jas_ci_reduce_mod
      contains
      subroutine optx_jas_ci_reduce
c Written by Claudia Filippi

      use optwf_parms, only: nparmj
      use dets, only: ndet
      use mix_jas_ci, only: de_o_ci, dj_de_ci, dj_o_ci, dj_oe_ci
      use optwf_control, only: ioptci, ioptjas
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

      do i=1,nparmj
        do j=1,nciterm
         dj_o_ci(i,j)=collect(i,j)
        enddo
      enddo

      call mpi_reduce(dj_de_ci,collect,nparmj*nciterm
     &     ,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

      call mpi_bcast(collect,nparmj*nciterm
     &     ,mpi_double_precision,0,MPI_COMM_WORLD,ierr)

      do i=1,nparmj
        do j=1,nciterm
         dj_de_ci(i,j)=collect(i,j)
        enddo
      enddo

      call mpi_reduce(de_o_ci,collect,nparmj*nciterm
     &     ,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

      call mpi_bcast(collect,nparmj*nciterm
     &     ,mpi_double_precision,0,MPI_COMM_WORLD,ierr)

      do i=1,nparmj
        do j=1,nciterm
         de_o_ci(i,j)=collect(i,j)
        enddo
      enddo

      call mpi_reduce(dj_oe_ci,collect,nparmj*nciterm
     &     ,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

      call mpi_bcast(collect,nparmj*nciterm
     &     ,mpi_double_precision,0,MPI_COMM_WORLD,ierr)

      do i=1,nparmj
        do j=1,nciterm
         dj_oe_ci(i,j)=collect(i,j)
        enddo
      enddo

      call mpi_barrier(MPI_COMM_WORLD,ierr)

      return
      end
      end module
