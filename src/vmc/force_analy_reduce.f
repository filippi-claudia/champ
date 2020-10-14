      subroutine force_analy_reduce
      use vmc_mod, only: MELEC, MORB, MBASIS, MDET, MCENT, MCTYPE, MCTYP3X
      use vmc_mod, only: NSPLIN, nrad, MORDJ, MORDJ1, MMAT_DIM, MMAT_DIM2, MMAT_DIM20
      use vmc_mod, only: radmax, delri
      use vmc_mod, only: NEQSX, MTERMS
      use vmc_mod, only: MCENT3, NCOEF, MEXCIT
      use atom, only: ncent, ncent_tot

      use mpiconf, only: wid
      use da_energy_sumcum, only: da_energy_cm2, da_energy_cum, da_psi_cum

      use force_analy, only: iforce_analy
      implicit real*8(a-h,o-z)

      include 'mpif.h'

      dimension collect(3*ncent_tot)

      if(iforce_analy.eq.0) return

      call mpi_reduce(da_energy_cum,collect,3*ncent
     &,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

      if(wid) then
        ii=0
        do 10 ic=1,ncent
          do 10 k=1,3
            ii=ii+1
  10        da_energy_cum(k,ic)=collect(ii)
       else
        do 15 ic=1,ncent
          do 15 k=1,3
  15        da_energy_cum(k,ic)=0.d0
      endif

      call mpi_reduce(da_psi_cum,collect,3*ncent
     &,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

      if(wid) then
        ii=0
        do 20 ic=1,ncent
          do 20 k=1,3
            ii=ii+1
  20        da_psi_cum(k,ic)=collect(ii)
       else
        do 25 ic=1,ncent
          do 25 k=1,3
  25        da_psi_cum(k,ic)=0.d0
      endif

      call mpi_reduce(da_energy_cm2,collect,3*ncent
     &,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

      if(wid) then
        ii=0
        do 30 ic=1,ncent
          do 30 k=1,3
            ii=ii+1
  30        da_energy_cm2(k,ic)=collect(ii)
       else
        do 35 ic=1,ncent
          do 35 k=1,3
  35        da_energy_cm2(k,ic)=0.d0
      endif

      call mpi_barrier(MPI_COMM_WORLD,ierr)

      return
      end
