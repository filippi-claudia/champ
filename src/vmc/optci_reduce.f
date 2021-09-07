      subroutine optci_reduce

      use precision_kinds, only: dp
      use optci, only: MXCITERM, MXCIREDUCED
      use optwf_contrl, only: ioptci
      use mstates_ctrl, only: iefficiency, nstates_psig
      use mstates2, only: effcm2, effcum
      use ci000, only: nciterm
      use ci005_blk, only: ci_o_cum
      use ci006_blk, only: ci_de_cum
      use ci008_blk, only: ci_oe_cm2, ci_oe_cum
      use ci009_blk, only: ci_oo_cm2, ci_oo_cum
      use ci010_blk, only: ci_ooe_cum
      use method_opt, only: method
      use mpi

      implicit none

      integer :: i, ierr, j, matdim
      real(dp) :: nmatdim, MXORBTERM

c     parameter(MXTMP=max(MXORBTERM,nmatdim))
c     max does not work with g77
    !   parameter(MXTMP=MXCITERM+MXCIMATDIM)
    !   dimension collect(MXTMP),collect2(MXCITERM,MXCIREDUCED)

      integer :: MXTMP
      real(dp), DIMENSION(:), ALLOCATABLE :: optci_reduce_collect
      real(dp), DIMENSION(:, :), ALLOCATABLE :: optci_reduce_collect2

      MXTMP=max(MXORBTERM,nmatdim)

      allocate(optci_reduce_collect(MXTMP))
      allocate(optci_reduce_collect2(MXCITERM,MXCIREDUCED))

      if (iefficiency.gt.0) then
        call mpi_reduce(effcum,optci_reduce_collect,nstates_psig
     &      ,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

        call mpi_bcast(optci_reduce_collect,nstates_psig
     &     ,mpi_double_precision,0,MPI_COMM_WORLD,ierr)

       do 80 j=1,nstates_psig
  80     effcum(j)=optci_reduce_collect(j)

       call mpi_reduce(effcm2,optci_reduce_collect,nstates_psig
     &      ,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

        call mpi_bcast(optci_reduce_collect,nstates_psig
     &     ,mpi_double_precision,0,MPI_COMM_WORLD,ierr)

       do 90 j=1,nstates_psig
  90     effcm2(j)=optci_reduce_collect(j)
      endif

      if(ioptci.eq.0.or.method.eq.'sr_n'.or.method.eq.'lin_d') then
        deallocate(optci_reduce_collect)
        deallocate(optci_reduce_collect2)
        return
      endif

      call mpi_reduce(ci_o_cum(1),optci_reduce_collect(1),nciterm
     &      ,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

      call mpi_bcast(optci_reduce_collect,nciterm
     &     ,mpi_double_precision,0,MPI_COMM_WORLD,ierr)

      do 10 i=1,nciterm
  10    ci_o_cum(i)=optci_reduce_collect(i)

      call mpi_reduce(ci_de_cum(1),optci_reduce_collect(1),nciterm
     &      ,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

      call mpi_bcast(optci_reduce_collect,nciterm
     &     ,mpi_double_precision,0,MPI_COMM_WORLD,ierr)

      do 20 i=1,nciterm
  20    ci_de_cum(i)=optci_reduce_collect(i)

      call mpi_reduce(ci_oe_cum,optci_reduce_collect2,MXCITERM*nciterm
     &     ,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

      call mpi_bcast(optci_reduce_collect2,MXCITERM*nciterm
     &     ,mpi_double_precision,0,MPI_COMM_WORLD,ierr)

      do 30 i=1,nciterm
       do 30 j=1,nciterm
  30    ci_oe_cum(i,j) = optci_reduce_collect2(i,j)

      call mpi_reduce(ci_oe_cm2,optci_reduce_collect2,MXCITERM*nciterm
     &     ,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

      call mpi_bcast(optci_reduce_collect2,MXCITERM*nciterm
     &     ,mpi_double_precision,0,MPI_COMM_WORLD,ierr)

      do 40 i=1,nciterm
       do 40 j=1,nciterm
  40    ci_oe_cm2(i,j) = optci_reduce_collect2(i,j)

      matdim=nciterm*(nciterm+1)/2

      call mpi_reduce(ci_oo_cum(1),optci_reduce_collect(1),matdim
     &      ,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

      call mpi_bcast(optci_reduce_collect,matdim
     &     ,mpi_double_precision,0,MPI_COMM_WORLD,ierr)

      do 50 i=1,matdim
  50    ci_oo_cum(i)=optci_reduce_collect(i)

      call mpi_reduce(ci_oo_cm2(1),optci_reduce_collect(1),matdim
     &      ,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

      call mpi_bcast(optci_reduce_collect,matdim
     &     ,mpi_double_precision,0,MPI_COMM_WORLD,ierr)

      do 60 i=1,matdim
  60    ci_oo_cm2(i)=optci_reduce_collect(i)

      call mpi_reduce(ci_ooe_cum(1),optci_reduce_collect(1),matdim
     &      ,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

      call mpi_bcast(optci_reduce_collect,matdim
     &     ,mpi_double_precision,0,MPI_COMM_WORLD,ierr)

      do 70 i=1,matdim
  70    ci_ooe_cum(i)=optci_reduce_collect(i)

      deallocate(optci_reduce_collect)
      deallocate(optci_reduce_collect2)

      call mpi_barrier(MPI_COMM_WORLD,ierr)

      return
      end
