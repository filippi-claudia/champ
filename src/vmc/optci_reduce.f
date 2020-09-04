      subroutine optci_reduce

      use optwf_contrl, only: ioptci
      use mstates_ctrl, only: iefficiency, iguiding, nstates_psig
      use mstates2, only: effcm2, effcum
      use ci000, only: iciprt, nciprim, nciterm

      implicit real*8(a-h,o-z)






      include 'vmc.h'
      include 'optci.h'
      include 'optci_cblk.h'
      include 'mstates.h'
      include 'mpif.h'

c     parameter(MXTMP=max(MXORBTERM,MXMATDIM))
c     max does not work with g77
      parameter(MXTMP=MXCITERM+MXCIMATDIM)
      dimension collect(MXTMP),collect2(MXCITERM,MXCIREDUCED)


      if (iefficiency.gt.0) then
        call mpi_reduce(effcum,collect,nstates_psig
     &      ,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

        call mpi_bcast(collect,nstates_psig
     &     ,mpi_double_precision,0,MPI_COMM_WORLD,ierr)

       do 80 j=1,nstates_psig
  80     effcum(j)=collect(j)

       call mpi_reduce(effcm2,collect,nstates_psig
     &      ,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

        call mpi_bcast(collect,nstates_psig
     &     ,mpi_double_precision,0,MPI_COMM_WORLD,ierr)

       do 90 j=1,nstates_psig
  90     effcm2(j)=collect(j)
      endif

      if(ioptci.eq.0.or.method.eq.'sr_n'.or.method.eq.'lin_d') return

      call mpi_reduce(ci_o_cum(1),collect(1),nciterm
     &      ,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

      call mpi_bcast(collect,nciterm
     &     ,mpi_double_precision,0,MPI_COMM_WORLD,ierr)

      do 10 i=1,nciterm
  10    ci_o_cum(i)=collect(i)

      call mpi_reduce(ci_de_cum(1),collect(1),nciterm
     &      ,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

      call mpi_bcast(collect,nciterm
     &     ,mpi_double_precision,0,MPI_COMM_WORLD,ierr)

      do 20 i=1,nciterm
  20    ci_de_cum(i)=collect(i)

      call mpi_reduce(ci_oe_cum,collect2,MXCITERM*nciterm
     &     ,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

      call mpi_bcast(collect2,MXCITERM*nciterm
     &     ,mpi_double_precision,0,MPI_COMM_WORLD,ierr)

      do 30 i=1,nciterm
       do 30 j=1,nciterm
  30    ci_oe_cum(i,j) = collect2(i,j)

      call mpi_reduce(ci_oe_cm2,collect2,MXCITERM*nciterm
     &     ,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

      call mpi_bcast(collect2,MXCITERM*nciterm
     &     ,mpi_double_precision,0,MPI_COMM_WORLD,ierr)

      do 40 i=1,nciterm
       do 40 j=1,nciterm
  40    ci_oe_cm2(i,j) = collect2(i,j) 

      matdim=nciterm*(nciterm+1)/2

      call mpi_reduce(ci_oo_cum(1),collect(1),matdim
     &      ,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

      call mpi_bcast(collect,matdim
     &     ,mpi_double_precision,0,MPI_COMM_WORLD,ierr)

      do 50 i=1,matdim
  50    ci_oo_cum(i)=collect(i)

      call mpi_reduce(ci_oo_cm2(1),collect(1),matdim
     &      ,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

      call mpi_bcast(collect,matdim
     &     ,mpi_double_precision,0,MPI_COMM_WORLD,ierr)

      do 60 i=1,matdim
  60    ci_oo_cm2(i)=collect(i)

      call mpi_reduce(ci_ooe_cum(1),collect(1),matdim
     &      ,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

      call mpi_bcast(collect,matdim
     &     ,mpi_double_precision,0,MPI_COMM_WORLD,ierr)

      do 70 i=1,matdim
  70    ci_ooe_cum(i)=collect(i)

      call mpi_barrier(MPI_COMM_WORLD,ierr)

      return
      end
