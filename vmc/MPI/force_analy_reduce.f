      subroutine force_analy_reduce

      implicit real*8(a-h,o-z)
      include 'mpif.h'
      include 'vmc.h'

      logical wid
      common /mpiconf/ idtask,nproc,wid

      common /atom/ znuc(MCTYPE),cent(3,MCENT),pecent
     &,iwctype(MCENT),nctype,ncent

      common /da_energy_ave/ da_energy_sum(3,MCENT),da_psi_sum(3,MCENT),
     & da_energy_cum(3,MCENT),da_psi_cum(3,MCENT),da_energy_cm2(3,MCENT)

      common /force_analy/ iforce_analy

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
