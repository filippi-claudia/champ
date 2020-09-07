      subroutine optx_orb_ci_reduce
c Written by Claudia Filippi

      use vmc, only: MELEC, MORB, MBASIS, MDET, MCENT, MCTYPE, MCTYP3X
      use vmc, only: NSPLIN, nrad, MORDJ, MORDJ1, MMAT_DIM, MMAT_DIM2, MMAT_DIM20
      use vmc, only: radmax, delri
      use vmc, only: NEQSX, MTERMS
      use vmc, only: MCENT3, NCOEF, MEXCIT
      use optwf_contrl, only: ioptci, ioptorb
      use mix_orb_ci, only: ci_de_o, ci_o_ho, ci_o_o, ci_o_oe

      use ci000, only: nciterm

      use method_opt, only: method

      implicit real*8(a-h,o-z)






      include 'optorb.h'
      include 'optci.h'
      include 'mpif.h'



      dimension collect(MXCITERM,MXREDUCED)

      if(ioptci.eq.0.or.ioptorb.eq.0.or.method.eq.'sr_n'.or.method.eq.'lin_d') return

      call mpi_reduce(ci_o_o,collect,MXCITERM*nreduced
     &     ,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

      call mpi_bcast(collect,MXCITERM*nreduced
     &     ,mpi_double_precision,0,MPI_COMM_WORLD,ierr)

      do 10 j=1,nreduced
        do 10 i=1,nciterm
  10     ci_o_o(i,j)=collect(i,j)

      call mpi_reduce(ci_o_oe,collect,MXCITERM*nreduced
     &     ,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

      call mpi_bcast(collect,MXCITERM*nreduced
     &     ,mpi_double_precision,0,MPI_COMM_WORLD,ierr)

      do 20 j=1,nreduced
        do 20 i=1,nciterm
  20     ci_o_oe(i,j)=collect(i,j)

      call mpi_reduce(ci_de_o,collect,MXCITERM*nreduced
     &     ,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

      call mpi_bcast(collect,MXCITERM*nreduced
     &     ,mpi_double_precision,0,MPI_COMM_WORLD,ierr)

      do 30 j=1,nreduced
        do 30 i=1,nciterm
  30     ci_de_o(i,j)=collect(i,j)

      call mpi_reduce(ci_o_ho,collect,MXCITERM*nreduced
     &     ,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

      call mpi_bcast(collect,MXCITERM*nreduced
     &     ,mpi_double_precision,0,MPI_COMM_WORLD,ierr)

      do 40 j=1,nreduced
        do 40 i=1,nciterm
  40     ci_o_ho(i,j)=collect(i,j)

      call mpi_barrier(MPI_COMM_WORLD,ierr)

      return
      end
