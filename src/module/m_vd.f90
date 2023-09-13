module vd_mod
    use precision_kinds, only: dp
    use system, only: ncent_tot
    use dmc_mod, only: mwalk
    use multiple_geo, only: MFORCE_WT_PRD
  
    implicit none
  
    integer :: dmc_ivd
    real(dp), dimension(:, :), allocatable :: da_branch !(3, MCENT, PTH)
    real(dp), dimension(:, :), allocatable :: da_branch_sum !(3, MCENT, PTH)
    real(dp), dimension(:, :), allocatable :: da_branch_cum !(3, MCENT, PTH)
    real(dp), dimension(:, :, :), allocatable :: esnake !(3, MCENT, MWALK, PTH)
    real(dp), dimension(:, :, :), allocatable :: deriv_eold !(3, MCENT, MWALK)
    real(dp), dimension(:, :, :, :), allocatable :: ehist !(3, MCENT, MWALK, 0:MFORCE_WT_PRD, PTH)
    
    private
    public :: dmc_ivd
    public :: da_branch, da_branch_sum, da_branch_cum, esnake, deriv_eold, ehist
    public :: allocate_da_branch, deallocate_da_branch
  
contains
    subroutine allocate_da_branch()
      if (dmc_ivd.gt.0) then
         if (.not. allocated(da_branch_cum)) allocate (da_branch_cum(3, ncent_tot))
         if (.not. allocated(da_branch_sum)) allocate (da_branch_sum(3, ncent_tot))
         if (.not. allocated(da_branch)) allocate (da_branch(3, ncent_tot))
         if (.not. allocated(esnake)) allocate (esnake(3, ncent_tot, mwalk))
         if (.not. allocated(deriv_eold)) allocate (deriv_eold(3, ncent_tot, mwalk))
         if (.not. allocated(ehist)) allocate (ehist(3, ncent_tot, mwalk, 0:MFORCE_WT_PRD))
      endif
    end subroutine allocate_da_branch
  
    subroutine deallocate_da_branch()
      use system, only: ncent
      use branch, only: nwalk
      use multiple_geo, only: nwprod
      
      integer :: k, ic, iw, ip
      
      if (dmc_ivd.gt.0) then
         if (allocated(da_branch_cum)) deallocate(da_branch_cum)
         if (allocated(da_branch_sum)) deallocate(da_branch_sum)
         if (allocated(da_branch)) deallocate(da_branch)
         if (allocated(esnake)) deallocate(esnake)
         if (allocated(deriv_eold)) deallocate(deriv_eold)
         if (allocated(ehist)) deallocate(ehist)
      endif
    end subroutine deallocate_da_branch
  
end module vd_mod