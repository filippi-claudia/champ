      module optjas_reduce_mod
      contains
      subroutine optjas_reduce
c Written by Claudia Filippi
      use mpi
      use optwf_parms, only: nparmj
      use csfs, only: nstates
      use mstates_mod, only: MSTATES
      use gradjerr, only: grad_jas_bcm2, grad_jas_bcum
      use optwf_control, only: ioptjas
      use optwf_parms, only: nparmj
      use gradhessj, only: d2j, d2j_e, de, de_de, de_e, dj, dj_de, dj_dj
      use gradhessj, only: dj_dj_e, dj_e, dj_e2
      use gradhessj, only: e2
      use gradjerrb, only: ngrad_jas_bcum, ngrad_jas_blocks
      use method_opt, only: method
      use precision_kinds, only: dp

      implicit none

      integer :: i, ierr, istate, j, ngrad_jas_collect

      real(dp), dimension(nparmj, MSTATES) :: collect
      real(dp), dimension(nparmj, nparmj, MSTATES) :: collect2


      if(ioptjas.eq.0.or.method.eq.'sr_n'.or.method.eq.'lin_d') return

c Note, to do: error is not collected

      call mpi_reduce(dj,collect,nparmj*nstates
     &     ,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

      call mpi_bcast(collect,nparmj*nstates
     &     ,mpi_double_precision,0,MPI_COMM_WORLD,ierr)

      do istate=1,nstates
        do i=1,nparmj
          dj(i,istate)=collect(i,istate)
        enddo
      enddo

      call mpi_reduce(dj_e,collect,nparmj*nstates
     &     ,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

      call mpi_bcast(collect,nparmj*nstates
     &     ,mpi_double_precision,0,MPI_COMM_WORLD,ierr)

      do istate=1,nstates
        do i=1,nparmj
          dj_e(i,istate)=collect(i,istate)
        enddo
      enddo

      call mpi_reduce(de,collect,nparmj*nstates
     &     ,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

      call mpi_bcast(collect,nparmj*nstates
     &     ,mpi_double_precision,0,MPI_COMM_WORLD,ierr)

      do istate=1,nstates
        do i=1,nparmj
          de(i,istate)=collect(i,istate)
        enddo
      enddo

      call mpi_reduce(de_e,collect,nparmj*nstates
     &     ,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

      call mpi_bcast(collect,nparmj*nstates
     &     ,mpi_double_precision,0,MPI_COMM_WORLD,ierr)

      do istate=1,nstates
        do i=1,nparmj
          de_e(i,istate)=collect(i,istate)
        enddo
      enddo

      call mpi_reduce(e2,collect,nparmj*nstates
     &     ,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

      call mpi_bcast(collect,nparmj*nstates
     &     ,mpi_double_precision,0,MPI_COMM_WORLD,ierr)

      do istate=1,nstates
        do i=1,nparmj
          e2(i,istate)=collect(i,istate)
        enddo
      enddo

      call mpi_reduce(dj_e2,collect,nparmj*nstates
     &     ,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

      call mpi_bcast(collect,nparmj*nstates
     &     ,mpi_double_precision,0,MPI_COMM_WORLD,ierr)

      do istate=1,nstates
        do i=1,nparmj
          dj_e2(i,istate)=collect(i,istate)
        enddo
      enddo

      call mpi_reduce(dj_de,collect2,nparmj*nparmj*nstates
     &     ,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

      call mpi_bcast(collect2,nparmj*nparmj*nstates
     &     ,mpi_double_precision,0,MPI_COMM_WORLD,ierr)

      do istate=1,nstates
        do i=1,nparmj
          do j=1,nparmj
            dj_de(i,j,istate)=collect2(i,j,istate)
          enddo
        enddo
      enddo

      call mpi_reduce(dj_dj,collect2,nparmj*nparmj*nstates
     &     ,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

      call mpi_bcast(collect2,nparmj*nparmj*nstates
     &     ,mpi_double_precision,0,MPI_COMM_WORLD,ierr)

      do istate=1,nstates
        do i=1,nparmj
          do j=1,nparmj
            dj_dj(i,j,istate)=collect2(i,j,istate)
          enddo
        enddo
      enddo

      call mpi_reduce(dj_dj_e,collect2,nparmj*nparmj*nstates
     &     ,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

      call mpi_bcast(collect2,nparmj*nparmj*nstates
     &     ,mpi_double_precision,0,MPI_COMM_WORLD,ierr)

      do istate=1,nstates
        do i=1,nparmj
          do j=1,nparmj
            dj_dj_e(i,j,istate)=collect2(i,j,istate)
          enddo
        enddo
      enddo

      call mpi_reduce(d2j,collect2,nparmj*nparmj*nstates
     &     ,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

      call mpi_bcast(collect2,nparmj*nparmj*nstates
     &     ,mpi_double_precision,0,MPI_COMM_WORLD,ierr)

      do istate=1,nstates
        do i=1,nparmj
          do j=1,nparmj
            d2j(i,j,istate)=collect2(i,j,istate)
          enddo
        enddo
      enddo

      call mpi_reduce(d2j_e,collect2,nparmj*nparmj*nstates
     &     ,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

      call mpi_bcast(collect2,nparmj*nparmj*nstates
     &     ,mpi_double_precision,0,MPI_COMM_WORLD,ierr)

      do istate=1,nstates
        do i=1,nparmj
          do j=1,nparmj
            d2j_e(i,j,istate)=collect2(i,j,istate)
          enddo
        enddo
      enddo

      call mpi_reduce(de_de,collect2,nparmj*nparmj*nstates
     &     ,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

      call mpi_bcast(collect2,nparmj*nparmj*nstates
     &     ,mpi_double_precision,0,MPI_COMM_WORLD,ierr)

      do istate=1,nstates
        do i=1,nparmj
          do j=1,nparmj
            de_de(i,j,istate)=collect2(i,j,istate)
          enddo
        enddo
      enddo

      if(ngrad_jas_blocks.gt.0) then
        call mpi_reduce(ngrad_jas_bcum,ngrad_jas_collect,1
     &     ,mpi_integer,mpi_sum,0,MPI_COMM_WORLD,ierr)

        call mpi_bcast(ngrad_jas_collect,1
     &     ,mpi_integer,0,MPI_COMM_WORLD,ierr)

        ngrad_jas_bcum=ngrad_jas_collect

        call mpi_reduce(grad_jas_bcum,collect,nparmj*nstates
     &     ,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

        call mpi_bcast(collect,nparmj*nstates
     &     ,mpi_double_precision,0,MPI_COMM_WORLD,ierr)

        do istate=1,nstates
          do i=1,nparmj
            grad_jas_bcum(i,istate)=collect(i,istate)
          enddo
        enddo

        call mpi_reduce(grad_jas_bcm2,collect,nparmj*nstates
     &     ,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

        call mpi_bcast(collect,nparmj*nstates
     &     ,mpi_double_precision,0,MPI_COMM_WORLD,ierr)

        do istate=1,nstates
          do i=1,nparmj
            grad_jas_bcm2(i,istate)=collect(i,istate)
          enddo
        enddo

      endif

c these averages should be set to zero on the slaves but optjas_reduce
c is only called at the end of run (differently than prop_reduce) and
c only the master writes to output and dumper

      call mpi_barrier(MPI_COMM_WORLD,ierr)

      return
      end
      end module
