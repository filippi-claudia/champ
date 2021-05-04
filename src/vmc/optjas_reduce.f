      subroutine optjas_reduce
c Written by Claudia Filippi
      use mpi
      use optjas, only: MPARMJ
      use csfs, only: nstates
      use mstates_mod, only: MSTATES
      use gradjerr, only: grad_jas_bcm2, grad_jas_bcum
      use optwf_contrl, only: ioptjas
      use optwf_parms, only: nparmj
      use gradhessj, only: d2j, d2j_e, de, de_de, de_e, dj, dj_de, dj_dj
      use gradhessj, only: dj_dj_e, dj_e, dj_e2
      use gradhessj, only: e2
      use gradjerrb, only: ngrad_jas_bcum, ngrad_jas_blocks
      use method_opt, only: method
      use precision_kinds, only: dp

      implicit none

      integer :: i, ierr, istate, j, ngrad_jas_collect

      real(dp), dimension(MPARMJ, MSTATES) :: collect
      real(dp), dimension(MPARMJ, MPARMJ, MSTATES) :: collect2


      if(ioptjas.eq.0.or.method.eq.'sr_n'.or.method.eq.'lin_d') return

c Note, to do: error is not collected

      call mpi_reduce(dj,collect,MPARMJ*nstates
     &     ,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

      call mpi_bcast(collect,MPARMJ*nstates
     &     ,mpi_double_precision,0,MPI_COMM_WORLD,ierr)

      do 10 istate=1,nstates
        do 10 i=1,nparmj
  10      dj(i,istate)=collect(i,istate)

      call mpi_reduce(dj_e,collect,MPARMJ*nstates
     &     ,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

      call mpi_bcast(collect,MPARMJ*nstates
     &     ,mpi_double_precision,0,MPI_COMM_WORLD,ierr)

      do 20 istate=1,nstates
        do 20 i=1,nparmj
  20      dj_e(i,istate)=collect(i,istate)

      call mpi_reduce(de,collect,MPARMJ*nstates
     &     ,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

      call mpi_bcast(collect,MPARMJ*nstates
     &     ,mpi_double_precision,0,MPI_COMM_WORLD,ierr)

      do 25 istate=1,nstates
        do 25 i=1,nparmj
  25      de(i,istate)=collect(i,istate)

      call mpi_reduce(de_e,collect,MPARMJ*nstates
     &     ,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

      call mpi_bcast(collect,MPARMJ*nstates
     &     ,mpi_double_precision,0,MPI_COMM_WORLD,ierr)

      do 26 istate=1,nstates
        do 26 i=1,nparmj
  26      de_e(i,istate)=collect(i,istate)

      call mpi_reduce(e2,collect,MPARMJ*nstates
     &     ,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

      call mpi_bcast(collect,MPARMJ*nstates
     &     ,mpi_double_precision,0,MPI_COMM_WORLD,ierr)

      do 27 istate=1,nstates
        do 27 i=1,nparmj
  27      e2(i,istate)=collect(i,istate)

      call mpi_reduce(dj_e2,collect,MPARMJ*nstates
     &     ,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

      call mpi_bcast(collect,MPARMJ*nstates
     &     ,mpi_double_precision,0,MPI_COMM_WORLD,ierr)

      do 28 istate=1,nstates
        do 28 i=1,nparmj
  28      dj_e2(i,istate)=collect(i,istate)

      call mpi_reduce(dj_de,collect2,MPARMJ*MPARMJ*nstates
     &     ,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

      call mpi_bcast(collect2,MPARMJ*MPARMJ*nstates
     &     ,mpi_double_precision,0,MPI_COMM_WORLD,ierr)

      do 30 istate=1,nstates
        do 30 i=1,nparmj
          do 30 j=1,nparmj
  30        dj_de(i,j,istate)=collect2(i,j,istate)

      call mpi_reduce(dj_dj,collect2,MPARMJ*MPARMJ*nstates
     &     ,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

      call mpi_bcast(collect2,MPARMJ*MPARMJ*nstates
     &     ,mpi_double_precision,0,MPI_COMM_WORLD,ierr)

      do 40 istate=1,nstates
        do 40 i=1,nparmj
          do 40 j=1,nparmj
  40        dj_dj(i,j,istate)=collect2(i,j,istate)

      call mpi_reduce(dj_dj_e,collect2,MPARMJ*MPARMJ*nstates
     &     ,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

      call mpi_bcast(collect2,MPARMJ*MPARMJ*nstates
     &     ,mpi_double_precision,0,MPI_COMM_WORLD,ierr)

      do 50 istate=1,nstates
        do 50 i=1,nparmj
          do 50 j=1,nparmj
  50        dj_dj_e(i,j,istate)=collect2(i,j,istate)

      call mpi_reduce(d2j,collect2,MPARMJ*MPARMJ*nstates
     &     ,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

      call mpi_bcast(collect2,MPARMJ*MPARMJ*nstates
     &     ,mpi_double_precision,0,MPI_COMM_WORLD,ierr)

      do 52 istate=1,nstates
        do 52 i=1,nparmj
          do 52 j=1,nparmj
  52        d2j(i,j,istate)=collect2(i,j,istate)

      call mpi_reduce(d2j_e,collect2,MPARMJ*MPARMJ*nstates
     &     ,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

      call mpi_bcast(collect2,MPARMJ*MPARMJ*nstates
     &     ,mpi_double_precision,0,MPI_COMM_WORLD,ierr)

      do 54 istate=1,nstates
        do 54 i=1,nparmj
          do 54 j=1,nparmj
  54        d2j_e(i,j,istate)=collect2(i,j,istate)

      call mpi_reduce(de_de,collect2,MPARMJ*MPARMJ*nstates
     &     ,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

      call mpi_bcast(collect2,MPARMJ*MPARMJ*nstates
     &     ,mpi_double_precision,0,MPI_COMM_WORLD,ierr)

      do 56 istate=1,nstates
        do 56 i=1,nparmj
          do 56 j=1,nparmj
  56        de_de(i,j,istate)=collect2(i,j,istate)

      if(ngrad_jas_blocks.gt.0) then
        call mpi_reduce(ngrad_jas_bcum,ngrad_jas_collect,1
     &     ,mpi_integer,mpi_sum,0,MPI_COMM_WORLD,ierr)

        call mpi_bcast(ngrad_jas_collect,1
     &     ,mpi_integer,0,MPI_COMM_WORLD,ierr)

        ngrad_jas_bcum=ngrad_jas_collect

        call mpi_reduce(grad_jas_bcum,collect,MPARMJ*nstates
     &     ,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

        call mpi_bcast(collect,MPARMJ*nstates
     &     ,mpi_double_precision,0,MPI_COMM_WORLD,ierr)

        do 60 istate=1,nstates
          do 60 i=1,nparmj
  60        grad_jas_bcum(i,istate)=collect(i,istate)

        call mpi_reduce(grad_jas_bcm2,collect,MPARMJ*nstates
     &     ,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

        call mpi_bcast(collect,MPARMJ*nstates
     &     ,mpi_double_precision,0,MPI_COMM_WORLD,ierr)

        do 70 istate=1,nstates
          do 70 i=1,nparmj
  70        grad_jas_bcm2(i,istate)=collect(i,istate)

      endif

c these averages should be set to zero on the slaves but optjas_reduce
c is only called at the end of run (differently than prop_reduce) and
c only the master writes to output and dumper

      call mpi_barrier(MPI_COMM_WORLD,ierr)

      return
      end
