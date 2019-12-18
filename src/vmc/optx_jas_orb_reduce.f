      subroutine optx_jas_orb_reduce
c Written by Claudia Filippi

      implicit real*8(a-h,o-z)

      include 'vmc.h'
      include 'force.h'
      include 'mstates.h'
      include 'optjas.h'
      include 'optorb.h'
      include 'mpif.h'

      common /optwf_contrl/ ioptjas,ioptorb,ioptci,nparm
      common /optwf_parms/ nparml,nparme,nparmd,nparms,nparmg,nparmj

      common /mix_jas_orb/ dj_o(MPARMJ,MXREDUCED,MSTATES),dj_oe(MPARMJ,MXREDUCED,MSTATES),
     &de_o(MPARMJ,MXREDUCED,MSTATES),dj_ho(MPARMJ,MXREDUCED,MSTATES)

      common /csfs/ ccsf(MDET,MSTATES,MWF),cxdet(MDET*MDETCSFX)
     &,icxdet(MDET*MDETCSFX),iadet(MDET),ibdet(MDET),ncsf,nstates

      dimension collect(MPARMJ,MXREDUCED)

      if(ioptjas.eq.0.or.ioptorb.eq.0.or.method.eq.'sr_n'.or.method.eq.'lin_d') return

      do 40 istate=1,nstates

        call mpi_reduce(dj_o(1,1,istate),collect,MPARMJ*nreduced
     &       ,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

        call mpi_bcast(collect,MPARMJ*nreduced
     &       ,mpi_double_precision,0,MPI_COMM_WORLD,ierr)

        do 10 j=1,nreduced
          do 10 i=1,nparmj
  10       dj_o(i,j,istate)=collect(i,j)

        call mpi_reduce(dj_oe(1,1,istate),collect,MPARMJ*nreduced
     &       ,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

        call mpi_bcast(collect,MPARMJ*nreduced
     &       ,mpi_double_precision,0,MPI_COMM_WORLD,ierr)

        do 20 j=1,nreduced
          do 20 i=1,nparmj
  20       dj_oe(i,j,istate)=collect(i,j)

        call mpi_reduce(de_o(1,1,istate),collect,MPARMJ*nreduced
     &       ,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

        call mpi_bcast(collect,MPARMJ*nreduced
     &       ,mpi_double_precision,0,MPI_COMM_WORLD,ierr)

        do 30 j=1,nreduced
          do 30 i=1,nparmj
  30       de_o(i,j,istate)=collect(i,j)

        call mpi_reduce(dj_ho(1,1,istate),collect,MPARMJ*nreduced
     &       ,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

        call mpi_bcast(collect,MPARMJ*nreduced
     &       ,mpi_double_precision,0,MPI_COMM_WORLD,ierr)

        do 40 j=1,nreduced
          do 40 i=1,nparmj
  40       dj_ho(i,j,istate)=collect(i,j)

      call mpi_barrier(MPI_COMM_WORLD,ierr)

      return
      end
