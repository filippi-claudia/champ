      subroutine optorb_reduce

      use optorb_mod, only: nmatdim
      use csfs, only: nstates
      use optorb_cblock, only: norbterm
      use optwf_contrl, only: ioptorb
      use orb_mat_003, only: orb_o_cum
      use orb_mat_004, only: orb_oe_cum
      use orb_mat_005, only: orb_ho_cum
      use orb_mat_006, only: orb_oo_cum
      use orb_mat_007, only: orb_oho_cum
      use orb_mat_024, only: orb_f_bcm2, orb_f_bcum
      use orb_mat_030, only: orb_ecum, orb_wcum
      use method_opt, only: method
      use mpi
      use precision_kinds, only: dp

      implicit none

      integer :: i, iefpsample, ierr, isample_cmat, istate
      integer :: matdim, norb_f_bcum, norb_f_collect, nreduced
      real(dp), dimension(norbterm+nmatdim) :: collect

      if(ioptorb.eq.0.or.method.eq.'sr_n'.or.method.eq.'lin_d') return

      call mpi_reduce(norb_f_bcum,norb_f_collect,1
     &     ,mpi_integer,mpi_sum,0,MPI_COMM_WORLD,ierr)

      call mpi_bcast(norb_f_collect,1
     &     ,mpi_integer,0,MPI_COMM_WORLD,ierr)

      norb_f_bcum=norb_f_collect

      do 60 istate=1,nstates

        call mpi_reduce(orb_o_cum(1,istate),collect,norbterm
     &        ,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

        call mpi_bcast(collect,norbterm
     &        ,mpi_double_precision,0,MPI_COMM_WORLD,ierr)

        do 10 i=1,norbterm
   10     orb_o_cum(i,istate)=collect(i)

        call mpi_reduce(orb_oe_cum(1,istate),collect,norbterm
     &        ,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

        call mpi_bcast(collect,norbterm
     &        ,mpi_double_precision,0,MPI_COMM_WORLD,ierr)

        do 20 i=1,norbterm
   20     orb_oe_cum(i,istate)=collect(i)

        call mpi_reduce(orb_ho_cum(1,istate),collect,norbterm
     &        ,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

        call mpi_bcast(collect,norbterm
     &        ,mpi_double_precision,0,MPI_COMM_WORLD,ierr)

        do 30 i=1,norbterm
   30     orb_ho_cum(i,istate)=collect(i)

        call mpi_reduce(orb_f_bcum(1,istate),collect,norbterm
     &        ,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

        call mpi_bcast(collect,norbterm
     &        ,mpi_double_precision,0,MPI_COMM_WORLD,ierr)

        do 50 i=1,norbterm
   50     orb_f_bcum(i,istate)=collect(i)

        call mpi_reduce(orb_f_bcm2(1,istate),collect,norbterm
     &        ,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

        call mpi_bcast(collect,norbterm
     &        ,mpi_double_precision,0,MPI_COMM_WORLD,ierr)

        do 60 i=1,norbterm
   60     orb_f_bcm2(i,istate)=collect(i)

      if(iefpsample.ne.1) then
        call mpi_reduce(orb_wcum,collect,nstates
     &       ,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)
        call mpi_bcast(collect,nstates
     &        ,mpi_double_precision,0,MPI_COMM_WORLD,ierr)

        do 62 istate=1,nstates
   62     orb_wcum(istate)=collect(istate)

        call mpi_reduce(orb_ecum,collect,nstates
     &       ,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)
        call mpi_bcast(collect,nstates
     &        ,mpi_double_precision,0,MPI_COMM_WORLD,ierr)

        do 64 istate=1,nstates
   64     orb_ecum(istate)=collect(istate)
      endif

      if(isample_cmat.eq.0) return

      matdim=nreduced*(nreduced+1)/2

      do 70 istate=1,nstates
        call mpi_reduce(orb_oo_cum(1,istate),collect,matdim
     &        ,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

        call mpi_bcast(collect,matdim
     &        ,mpi_double_precision,0,MPI_COMM_WORLD,ierr)

        do 70 i=1,matdim
   70     orb_oo_cum(i,istate)=collect(i)

      matdim=nreduced*nreduced

      do 75 istate=1,nstates
        call mpi_reduce(orb_oho_cum(1,istate),collect,matdim
     &        ,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

        call mpi_bcast(collect,matdim
     &        ,mpi_double_precision,0,MPI_COMM_WORLD,ierr)

        do 75 i=1,matdim
   75     orb_oho_cum(i,istate)=collect(i)

c these averages should be set to zero on the slaves but optorb_reduce
c is only called at the end of run (differently than prop_reduce) and
c only the master writes to output and dumper

      call mpi_barrier(MPI_COMM_WORLD,ierr)

      return
      end
