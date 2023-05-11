      module optx_orb_ci_reduce_mod
      contains
      subroutine optx_orb_ci_reduce
c Written by Claudia Filippi

      use optorb_mod, only: mxreduced
      use optci, only: mxciterm
      use optwf_control, only: ioptci, ioptorb
      use mix_orb_ci, only: ci_de_o, ci_o_ho, ci_o_o, ci_o_oe
      use ci000, only: nciterm
      use optwf_control, only: method
      use optorb_cblock, only: nreduced
      use mpi
      use optci,   only: mxciterm
      use optorb_cblock, only: nreduced
      use optorb_mod, only: mxreduced
      use optwf_control, only: ioptci,ioptorb,method
      use precision_kinds, only: dp

      implicit none

      integer :: i, ierr, j

      real(dp), dimension(mxciterm,mxreduced) :: collect


      if(ioptci.eq.0.or.ioptorb.eq.0.or.method.eq.'sr_n'.or.method.eq.'lin_d') return

      call mpi_reduce(ci_o_o,collect,mxciterm*nreduced
     &     ,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

      call mpi_bcast(collect,mxciterm*nreduced
     &     ,mpi_double_precision,0,MPI_COMM_WORLD,ierr)

      do j=1,nreduced
        do i=1,nciterm
         ci_o_o(i,j)=collect(i,j)
        enddo
      enddo

      call mpi_reduce(ci_o_oe,collect,mxciterm*nreduced
     &     ,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

      call mpi_bcast(collect,mxciterm*nreduced
     &     ,mpi_double_precision,0,MPI_COMM_WORLD,ierr)

      do j=1,nreduced
        do i=1,nciterm
         ci_o_oe(i,j)=collect(i,j)
        enddo
      enddo

      call mpi_reduce(ci_de_o,collect,mxciterm*nreduced
     &     ,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

      call mpi_bcast(collect,mxciterm*nreduced
     &     ,mpi_double_precision,0,MPI_COMM_WORLD,ierr)

      do j=1,nreduced
        do i=1,nciterm
         ci_de_o(i,j)=collect(i,j)
        enddo
      enddo

      call mpi_reduce(ci_o_ho,collect,mxciterm*nreduced
     &     ,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

      call mpi_bcast(collect,mxciterm*nreduced
     &     ,mpi_double_precision,0,MPI_COMM_WORLD,ierr)

      do j=1,nreduced
        do i=1,nciterm
         ci_o_ho(i,j)=collect(i,j)
        enddo
      enddo

      call mpi_barrier(MPI_COMM_WORLD,ierr)

      return
      end
      end module
