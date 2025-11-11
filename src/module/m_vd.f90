!> @brief Module for velocity drifts in DMC force calculations
!> @date 2025
!> @author Emiel Slootman
!> @details This module manages arrays and data structures for computing
!> velocity drifts of the local energy with respect to nuclear
!> positions in Diffusion Monte Carlo (DMC) calculations. These derivatives
!> are essential for calculating forces on nuclei.
!>
!> The module provides:
!> - Storage for branching derivatives and their accumulation
!> - Energy derivative history tracking across walkers and time steps
!> - Memory management routines for dynamic allocation/deallocation
!>
!> @note The arrays are only allocated when dmc_ivd > 0, indicating that
!> velocity drift calculations are enabled.
module vd_mod
    use precision_kinds, only: dp
    use system, only: ncent_tot
    use dmc_mod, only: mwalk
    use multiple_geo, only: MFORCE_WT_PRD
    use force_pth, only: PTH
  
    implicit none
  
    !> Flag to enable/disable velocity drift calculations.
    integer :: dmc_ivd
    
    !> Branching derivatives for current DMC step.
    real(dp), dimension(:, :, :), allocatable :: da_branch
    
    !> Sum of branching derivatives over walkers.
    real(dp), dimension(:, :, :), allocatable :: da_branch_sum
    
    !> Cumulative branching derivatives.
    real(dp), dimension(:, :, :), allocatable :: da_branch_cum
    
    !> Energy derivatives for walker ensemble.
    real(dp), dimension(:, :, :, :), allocatable :: esnake
    
    !> Previous step energy derivatives.
    real(dp), dimension(:, :, :), allocatable :: deriv_eold
    
    !> Historical energy derivatives with time weighting.
    real(dp), dimension(:, :, :, :, :), allocatable :: ehist
    
    private
    public :: dmc_ivd
    public :: da_branch, da_branch_sum, da_branch_cum, esnake, deriv_eold, ehist
    public :: allocate_da_branch, deallocate_da_branch
  
contains

    !> Allocates memory for velocity drift arrays.
    subroutine allocate_da_branch()
      if (dmc_ivd.gt.0) then
         if (.not. allocated(da_branch_cum)) allocate (da_branch_cum(3, ncent_tot, PTH))
         if (.not. allocated(da_branch_sum)) allocate (da_branch_sum(3, ncent_tot, PTH))
         if (.not. allocated(da_branch)) allocate (da_branch(3, ncent_tot, PTH))
         if (.not. allocated(esnake)) allocate (esnake(3, ncent_tot, mwalk, PTH))
         if (.not. allocated(deriv_eold)) allocate (deriv_eold(3, ncent_tot, mwalk))
         if (.not. allocated(ehist)) allocate (ehist(3, ncent_tot, mwalk, 0:MFORCE_WT_PRD, PTH))
      endif
    end subroutine allocate_da_branch
  
    !> Deallocates memory for velocity drift arrays.
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