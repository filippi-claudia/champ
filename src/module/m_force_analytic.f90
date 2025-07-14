 module m_force_analytic
     !> Arguments: iforce_analy, iuse_zmat, alfgeo
      use precision_kinds, only: dp
     implicit none

     integer :: iforce_analy
     integer :: f_analy_err
     integer :: iuse_zmat
     real(dp) :: alfgeo
     real(dp) :: energy_block
     integer :: current_block

    ! from force_fin
     real(dp), dimension(:, :, :), allocatable :: da_energy_ave !(3,MCENT,PTH)
     real(dp), dimension(:), allocatable :: da_energy_err !(3)

     real(dp), dimension(:, :, :, :), allocatable :: da_energy_psi_block !(nblk, 3, MCENT, PTH) 
     real(dp), dimension(:, :, :, :), allocatable :: da_psi_block !(nblk, 3, MCENT, PTH) 
     real(dp), dimension(:), allocatable :: wcum_block !(nblk)
     ! from force_mat_n
     real(dp), dimension(:, :), allocatable :: force_o !(6*MCENT,mconf)

     private
     public   :: iforce_analy, iuse_zmat, alfgeo, f_analy_err
     public   :: da_energy_ave, da_energy_err
     public   :: allocate_force_fin, deallocate_force_fin
     public   :: force_o
     public   :: allocate_force_mat_n, deallocate_force_mat_n
     public   :: da_energy_psi_block, da_psi_block, wcum_block
     public   :: current_block
     save

 contains
     subroutine allocate_force_fin()
      use system,  only: ncent_tot
      use force_pth, only: PTH
      use control_vmc, only: vmc_nblk

         if (.not. allocated(da_energy_ave)) allocate (da_energy_ave(3, ncent_tot, PTH))
         if (.not. allocated(da_energy_err)) allocate (da_energy_err(3))
         if (.not. allocated(da_energy_psi_block)) allocate (da_energy_psi_block(vmc_nblk, 3, ncent_tot, PTH))
         if (.not. allocated(da_psi_block)) allocate (da_psi_block(vmc_nblk, 3, ncent_tot, PTH))
         if (.not. allocated(wcum_block)) allocate (wcum_block(vmc_nblk))
     end subroutine allocate_force_fin

     subroutine deallocate_force_fin()
         if (allocated(da_energy_err)) deallocate (da_energy_err)
         if (allocated(da_energy_ave)) deallocate (da_energy_ave)
         if (allocated(da_psi_block)) deallocate (da_psi_block)
         if (allocated(da_energy_psi_block)) deallocate (da_energy_psi_block)
         if (allocated(wcum_block)) deallocate (wcum_block)
     end subroutine deallocate_force_fin

     subroutine allocate_force_mat_n()
      use sr_mod,  only: mconf
      use system,  only: ncent_tot
         if (.not. allocated(force_o)) allocate (force_o(6*ncent_tot, mconf))
     end subroutine allocate_force_mat_n

     subroutine deallocate_force_mat_n()
         if (allocated(force_o)) deallocate (force_o)
     end subroutine deallocate_force_mat_n

 end module m_force_analytic
