      subroutine optx_orb_ci_reduce
c Written by Claudia Filippi

      use optwf_contrl, only: ioptci, ioptjas, ioptorb, nparm
      use optwf_parms, only: nparmd, nparme, nparmg, nparmj, nparml, nparms
      implicit real*8(a-h,o-z)



      include 'vmc.h'
      include 'optorb.h'
      include 'optci.h'
      include 'mpif.h'


      common /mix_orb_ci/ ci_o_o(MXCITERM,MXREDUCED),ci_o_oe(MXCITERM,MXREDUCED),
     &ci_de_o(MXCITERM,MXREDUCED),ci_o_ho(MXCITERM,MXREDUCED)

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
