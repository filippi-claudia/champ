!> @brief Module for saving and storing wavefunction parameters during optimization.
!> @author CHAMP developers
!> @date 2025
!>
!> @details This module stores current and best wavefunction parameters during optimization,
!> including Jastrow parameters (a, b, c), CI/CSF coefficients, and determinant coefficients.
!> It maintains both working copies and best-found parameters for backtracking and restart purposes.
module save_mod
  use precision_kinds, only: dp
  implicit none

  !> Number of type-a Jastrow parameters.
  integer :: mparmja

  !> Number of type-b Jastrow parameters.
  integer :: mparmjb

  !> Number of type-c Jastrow parameters.
  integer :: mparmjc

  !> Saved type-a (e-n) Jastrow parameters (4th-order).
  real(dp), allocatable :: a4_save(:,:,:)

  !> Saved type-b (e-e) Jastrow parameters.
  real(dp), allocatable :: b_save(:,:,:)

  !> Saved type-c (e-e-n) Jastrow parameters.
  real(dp), allocatable :: c_save(:,:,:)

  !> Saved orbital coefficients.
  real(dp), allocatable :: coef_save(:,:,:)

  !> Saved determinant coefficients.
  real(dp), allocatable :: cdet_save(:,:)

  !> Saved CSF coefficients.
  real(dp), allocatable :: ccsf_save(:,:)

  !> Number of type-a Jastrow parameters for best wavefunction.
  integer :: mparmja_best

  !> Number of type-b Jastrow parameters for best wavefunction.
  integer :: mparmjb_best

  !> Number of type-c Jastrow parameters for best wavefunction.
  integer :: mparmjc_best

  !> Best type-a Jastrow parameters (4th-order).
  real(dp), allocatable :: a4_best(:,:,:)

  !> Best type-b Jastrow parameters.
  real(dp), allocatable :: b_best(:,:,:)

  !> Best type-c Jastrow parameters.
  real(dp), allocatable :: c_best(:,:,:)

  !> Best orbital coefficients.
  real(dp), allocatable :: coef_best(:,:,:)

  !> Best determinant coefficients.
  real(dp), allocatable :: cdet_best(:,:)

  !> Best CSF coefficients.
  real(dp), allocatable :: ccsf_best(:,:)

  !> Saved number of CI terms.
  integer :: nciterm_sav

  !> Saved number of orbital terms.
  integer :: norbterm_sav

  !> Saved number of determinant parameters.
  integer :: nparmd_sav

  !> Saved number of Jastrow parameters.
  integer :: nparmj_sav

  !> Saved number of reduced parameters.
  integer :: nreduced_sav

  !> Saved number of active orbitals.
  integer :: nadorb_sav

  !> Saved pseudopotential matrix elements for nonlocal grids from dmc/nonloc_grid.f.
  real(dp), allocatable :: t_vpsp_save(:, :, :)

  save
end module
