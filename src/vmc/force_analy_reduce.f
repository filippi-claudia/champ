      subroutine force_analy_reduce
      use mstates_mod, only: MSTATES
      use vmc_mod, only: MELEC, MORB, MBASIS, MDET, MCENT, MCTYPE, MCTYP3X
      use vmc_mod, only: NSPLIN, nrad, MORDJ, MORDJ1, MMAT_DIM, MMAT_DIM2, MMAT_DIM20
      use vmc_mod, only: radmax, delri
      use vmc_mod, only: NEQSX, MTERMS
      use vmc_mod, only: MCENT3, NCOEF, MEXCIT
      use atom, only: ncent
      use mpiconf, only: wid
      use da_energy_sumcum, only: da_energy_cm2, da_energy_cum, da_psi_cum
      use force_analy, only: iforce_analy
      use csfs, only: nstates
      use mpi

      implicit real*8(a-h,o-z)

      dimension collect(3*MCENT,MSTATES)

      if(iforce_analy.eq.0) return

      do istate=1,nstates
         call mpi_reduce(da_energy_cum(1,1,istate),collect(1,istate),3*ncent
     &        ,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

         if(wid) then
            ii=0
            do ic=1,ncent
               do k=1,3
                  ii=ii+1
                  da_energy_cum(k,ic,istate)=collect(ii,istate)
               enddo
            enddo
         else
            do ic=1,ncent
               do k=1,3
                  da_energy_cum(k,ic,istate)=0.0d0
               enddo
            enddo
         endif

         call mpi_reduce(da_psi_cum(1,1,istate),collect(1,istate),3*ncent
     &        ,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

         if(wid) then
            ii=0
            do ic=1,ncent
               do k=1,3
                  ii=ii+1
                  da_psi_cum(k,ic,istate)=collect(ii,istate)
               enddo
            enddo
         else
            do ic=1,ncent
               do k=1,3
                  da_psi_cum(k,ic,istate)=0.0d0
               enddo
            enddo
         endif

         call mpi_reduce(da_energy_cm2(1,1,istate),collect(1,istate),3*ncent
     &        ,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

         if(wid) then
            ii=0
            do ic=1,ncent
               do k=1,3
                  ii=ii+1
                  da_energy_cm2(k,ic,istate)=collect(ii,istate)
               enddo
            enddo
         else
            do ic=1,ncent
               do k=1,3
                  da_energy_cm2(k,ic,istate)=0.0d0
               enddo
            enddo
         endif
      enddo

      call mpi_barrier(MPI_COMM_WORLD,ierr)

      end subroutine
