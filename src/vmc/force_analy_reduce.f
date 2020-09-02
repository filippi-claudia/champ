      subroutine force_analy_reduce
      use atom, only: ncent

      use mpiconf, only: wid
      use da_energy_sumcum, only: da_energy_cm2, da_energy_cum, da_psi_cum

      use force_analy, only: iforce_analy
      implicit real*8(a-h,o-z)

      include 'mpif.h'
      include 'vmc.h'

      dimension collect(3*MCENT)

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
