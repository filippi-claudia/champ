      module optci_reduce_mod
      contains
      subroutine optci_reduce

      use ci000,   only: nciterm
      use ci005_blk, only: ci_o_cum
      use ci006_blk, only: ci_de_cum
      use ci008_blk, only: ci_oe_cm2,ci_oe_cum
      use ci009_blk, only: ci_oo_cm2,ci_oo_cum
      use ci010_blk, only: ci_ooe_cum
      use mpi
      use mstates2, only: effcm2,effcum
      use mstates_ctrl, only: iefficiency,nstates_psig
      use optci,   only: mxcireduced,mxciterm,ncimatdim
      use optorb_cblock, only: norbterm
      use optorb_mod, only: nmatdim
      use optwf_control, only: ioptci,method
      use precision_kinds, only: dp

      implicit none

      integer :: i, ierr, j, matdim
      ! real(dp) :: nmatdim, MXORBTERM

c     parameter(MXTMP=max(MXORBTERM,nmatdim))
c     max does not work with g77
    !   parameter(MXTMP=mxciterm+ncimatdim)
    !   dimension collect(MXTMP),collect2(mxciterm,mxcireduced)

      integer :: MXTMP
      real(dp), DIMENSION(:), ALLOCATABLE :: optci_reduce_collect
      real(dp), DIMENSION(:, :), ALLOCATABLE :: optci_reduce_collect2

      ! defined like that in reference branch
      ! https://github.com/filippi-claudia/champ/blob/847f0e5e94d77035d957158406aac47c3e27af54/src/vmc/optci_reduce.f#L16
      ! MXTMP = mxciterm + ncimatdim
      MXTMP=max(norbterm,ncimatdim)

      allocate(optci_reduce_collect(MXTMP))
      allocate(optci_reduce_collect2(mxciterm,mxcireduced))

      if (iefficiency.gt.0) then
        call mpi_reduce(effcum,optci_reduce_collect,nstates_psig
     &      ,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

        call mpi_bcast(optci_reduce_collect,nstates_psig
     &     ,mpi_double_precision,0,MPI_COMM_WORLD,ierr)

       do j=1,nstates_psig
         effcum(j)=optci_reduce_collect(j)
       enddo

       call mpi_reduce(effcm2,optci_reduce_collect,nstates_psig
     &      ,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

        call mpi_bcast(optci_reduce_collect,nstates_psig
     &     ,mpi_double_precision,0,MPI_COMM_WORLD,ierr)

       do j=1,nstates_psig
         effcm2(j)=optci_reduce_collect(j)
       enddo
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

      do i=1,nciterm
        ci_o_cum(i)=optci_reduce_collect(i)
      enddo

      call mpi_reduce(ci_de_cum(1),optci_reduce_collect(1),nciterm
     &      ,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

      call mpi_bcast(optci_reduce_collect,nciterm
     &     ,mpi_double_precision,0,MPI_COMM_WORLD,ierr)

      do i=1,nciterm
        ci_de_cum(i)=optci_reduce_collect(i)
      enddo

      call mpi_reduce(ci_oe_cum,optci_reduce_collect2,mxciterm*nciterm
     &     ,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

      call mpi_bcast(optci_reduce_collect2,mxciterm*nciterm
     &     ,mpi_double_precision,0,MPI_COMM_WORLD,ierr)

      do i=1,nciterm
       do j=1,nciterm
        ci_oe_cum(i,j) = optci_reduce_collect2(i,j)
       enddo
      enddo

      call mpi_reduce(ci_oe_cm2,optci_reduce_collect2,mxciterm*nciterm
     &     ,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

      call mpi_bcast(optci_reduce_collect2,mxciterm*nciterm
     &     ,mpi_double_precision,0,MPI_COMM_WORLD,ierr)

      do i=1,nciterm
       do j=1,nciterm
        ci_oe_cm2(i,j) = optci_reduce_collect2(i,j)
       enddo
      enddo

      matdim=nciterm*(nciterm+1)/2

      call mpi_reduce(ci_oo_cum(1),optci_reduce_collect(1),matdim
     &      ,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

      call mpi_bcast(optci_reduce_collect,matdim
     &     ,mpi_double_precision,0,MPI_COMM_WORLD,ierr)

      do i=1,matdim
        ci_oo_cum(i)=optci_reduce_collect(i)
      enddo

      call mpi_reduce(ci_oo_cm2(1),optci_reduce_collect(1),matdim
     &      ,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

      call mpi_bcast(optci_reduce_collect,matdim
     &     ,mpi_double_precision,0,MPI_COMM_WORLD,ierr)

      do i=1,matdim
        ci_oo_cm2(i)=optci_reduce_collect(i)
      enddo

      call mpi_reduce(ci_ooe_cum(1),optci_reduce_collect(1),matdim
     &      ,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

      call mpi_bcast(optci_reduce_collect,matdim
     &     ,mpi_double_precision,0,MPI_COMM_WORLD,ierr)

      do i=1,matdim
        ci_ooe_cum(i)=optci_reduce_collect(i)
      enddo

      deallocate(optci_reduce_collect)
      deallocate(optci_reduce_collect2)

      call mpi_barrier(MPI_COMM_WORLD,ierr)

      return
      end
      end module
