      module force_analy_reduce_mod
      contains
      subroutine force_analy_reduce
      use system, only: ncent, ncent_tot

      use mpiconf, only: wid
      use da_energy_sumcum, only: da_energy_cm2, da_energy_cum, da_psi_cum
      use m_force_analytic, only: iforce_analy
      use mpi

      use precision_kinds, only: dp
      implicit none

      integer :: ic, ierr, ii, k

      real(dp), dimension(3*ncent_tot) :: collect



      if(iforce_analy.eq.0) return

      call mpi_reduce(da_energy_cum,collect,3*ncent
     &,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

      if(wid) then
        ii=0
        do ic=1,ncent
          do k=1,3
            ii=ii+1
            da_energy_cum(k,ic)=collect(ii)
          enddo
        enddo
       else
        do ic=1,ncent
          do k=1,3
            da_energy_cum(k,ic)=0.d0
          enddo
        enddo
      endif

      call mpi_reduce(da_psi_cum,collect,3*ncent
     &,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

      if(wid) then
        ii=0
        do ic=1,ncent
          do k=1,3
            ii=ii+1
            da_psi_cum(k,ic)=collect(ii)
          enddo
        enddo
       else
        do ic=1,ncent
          do k=1,3
            da_psi_cum(k,ic)=0.d0
          enddo
        enddo
      endif

      call mpi_reduce(da_energy_cm2,collect,3*ncent
     &,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

      if(wid) then
        ii=0
        do ic=1,ncent
          do k=1,3
            ii=ii+1
            da_energy_cm2(k,ic)=collect(ii)
          enddo
        enddo
       else
        do ic=1,ncent
          do k=1,3
            da_energy_cm2(k,ic)=0.d0
          enddo
        enddo
      endif

      call mpi_barrier(MPI_COMM_WORLD,ierr)

      return
      end
      end module
