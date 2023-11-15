module optx_jas_orb_reduce_mod
contains
      subroutine optx_jas_orb_reduce
! Written by Claudia Filippi

      use optorb_mod, only: mxreduced
      use optwf_parms, only: nparmj
      use csfs, only: nstates
      use optwf_control, only: ioptjas, ioptorb
      use optwf_parms, only: nparmj
      use mix_jas_orb, only: de_o, dj_ho, dj_o, dj_oe
      use optwf_control, only: method
      use mpi
      use optorb_mod, only: mxreduced
      use optwf_control, only: ioptjas,ioptorb,method
      use optwf_parms, only: nparmj
      use precision_kinds, only: dp

      implicit none

      integer :: i, ierr, istate, j, nreduced
      real(dp), dimension(nparmj,mxreduced) :: collect


      if(ioptjas.eq.0.or.ioptorb.eq.0.or.method.eq.'sr_n'.or.method.eq.'lin_d') return

      do istate=1,nstates

        call mpi_reduce(dj_o(1,1,istate),collect,nparmj*nreduced &
             ,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

        call mpi_bcast(collect,nparmj*nreduced &
             ,mpi_double_precision,0,MPI_COMM_WORLD,ierr)

        do j=1,nreduced
          do i=1,nparmj
           dj_o(i,j,istate)=collect(i,j)
          enddo
        enddo

        call mpi_reduce(dj_oe(1,1,istate),collect,nparmj*nreduced &
             ,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

        call mpi_bcast(collect,nparmj*nreduced &
             ,mpi_double_precision,0,MPI_COMM_WORLD,ierr)

        do j=1,nreduced
          do i=1,nparmj
           dj_oe(i,j,istate)=collect(i,j)
          enddo
        enddo

        call mpi_reduce(de_o(1,1,istate),collect,nparmj*nreduced &
             ,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

        call mpi_bcast(collect,nparmj*nreduced &
             ,mpi_double_precision,0,MPI_COMM_WORLD,ierr)

        do j=1,nreduced
          do i=1,nparmj
           de_o(i,j,istate)=collect(i,j)
          enddo
        enddo

        call mpi_reduce(dj_ho(1,1,istate),collect,nparmj*nreduced &
             ,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

        call mpi_bcast(collect,nparmj*nreduced &
             ,mpi_double_precision,0,MPI_COMM_WORLD,ierr)

        do j=1,nreduced
          do i=1,nparmj
           dj_ho(i,j,istate)=collect(i,j)
          enddo
        enddo
      enddo


      return
      end
end module
