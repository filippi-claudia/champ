module force_analy_reduce_mod
      contains
      subroutine force_analy_reduce()
      use da_energy_sumcum, only: da_energy_cum,da_psi_cum,da_energy_cm2
      use m_force_analytic, only: iforce_analy, da_energy_psi_block
      use m_force_analytic, only: wcum_block, da_psi_block
      use force_pth, only: PTH
      use mpi
      use mpiconf, only: wid
      use precision_kinds, only: dp
      use system,  only: ncent,ncent_tot
      use vd_mod, only: dmc_ivd, da_branch_cum
      use control_vmc, only: vmc_nblk


      implicit none

      integer :: ic, ierr, ii, k, iph, iblk

      real(dp), dimension(3*ncent_tot*PTH) :: collect
      real(dp), dimension(vmc_nblk*3*ncent_tot*PTH) :: collect2
      real(dp), dimension(vmc_nblk) :: collect3

!     no multistate indices
      if(iforce_analy.eq.0) return

      call mpi_reduce(da_energy_cum,collect,3*ncent*PTH,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

      if(wid) then
        ii=0
        do iph=1,PTH
          do ic=1,ncent
            do k=1,3
              ii=ii+1
              da_energy_cum(k,ic,iph)=collect(ii)
            enddo
          enddo
        enddo
      else
        do iph=1,PTH
          do ic=1,ncent
            do k=1,3
              da_energy_cum(k,ic,iph)=0.d0
            enddo
          enddo
        enddo
      endif

      call mpi_reduce(da_psi_cum,collect,3*ncent*PTH,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)


      if(wid) then
        ii=0
        do iph=1,PTH
          do ic=1,ncent
            do k=1,3
              ii=ii+1
              da_psi_cum(k,ic,iph)=collect(ii)
            enddo
          enddo
        enddo
       else
        do iph=1,PTH
          do ic=1,ncent
            do k=1,3
              da_psi_cum(k,ic,iph)=0.d0
            enddo
          enddo
        enddo
      endif

      if (dmc_ivd.gt.0) then
        call mpi_reduce(da_branch_cum,collect,3*ncent*PTH,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

        if(wid) then
          ii=0
          do iph=1,PTH
            do ic=1,ncent
              do k=1,3
                ii=ii+1
                da_branch_cum(k,ic,iph)=collect(ii)
              enddo
            enddo
          enddo
        else
          do iph=1,PTH
            do ic=1,ncent
              do k=1,3
                da_branch_cum(k,ic,iph)=0.d0
              enddo
            enddo
          enddo
        endif
      endif

       call mpi_reduce(da_energy_cm2,collect,3*ncent*PTH,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

       if(wid) then
         ii=0
         do iph=1,PTH
           do ic=1,ncent
             do k=1,3
               ii=ii+1
               da_energy_cm2(k,ic,iph)=collect(ii)
             enddo
           enddo
         enddo
       else
         do iph=1,PTH
           do ic=1,ncent
             do k=1,3
               da_energy_cm2(k,ic,iph)=0.d0
             enddo
           enddo
         enddo
       endif

       call mpi_reduce(da_energy_psi_block,collect2,vmc_nblk*3*ncent*PTH,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

       if(wid) then
         ii=0
         do iph=1,PTH
           do ic=1,ncent
             do k=1,3
               do iblk=1,vmc_nblk
                 ii=ii+1
                 da_energy_psi_block(iblk,k,ic,iph)=collect2(ii)
               enddo
             enddo
           enddo
         enddo
       else
         do iph=1,PTH
           do ic=1,ncent
             do k=1,3
               do iblk=1,vmc_nblk
                 da_energy_psi_block(iblk,k,ic,iph)=0.d0
               enddo
             enddo
           enddo
         enddo
       endif

       call mpi_reduce(wcum_block,collect3,vmc_nblk,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

       if(wid) then
         ii=0
         do iblk=1,vmc_nblk
           ii=ii+1
           wcum_block(iblk)=collect3(ii)
         enddo
       else
        do iblk=1,vmc_nblk
          wcum_block(iblk)=0.d0
        enddo
       endif

       call mpi_reduce(da_psi_block,collect2,3*ncent*vmc_nblk*PTH,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

       if(wid) then
         ii=0
         do iph=1,PTH
           do ic=1,ncent
             do k=1,3
               do iblk=1,vmc_nblk
                 ii=ii+1
                 da_psi_block(iblk,k,ic,iph)=collect2(ii)
               enddo
             enddo
           enddo
         enddo
       else
         do iph=1,PTH
           do ic=1,ncent
             do k=1,3
               do iblk=1,vmc_nblk
                 da_psi_block(iblk,k,ic,iph)=0.d0
               enddo
             enddo
           enddo
         enddo
       endif

      return
      end
end module
