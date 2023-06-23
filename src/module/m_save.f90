module save_mod
  use precision_kinds, only: dp
  implicit none
  integer :: mparmja, mparmjb, mparmjc
  real(dp), allocatable :: a4_save(:,:,:)
  real(dp), allocatable :: b_save(:,:,:)
  real(dp), allocatable :: c_save(:,:,:)
  real(dp), allocatable :: coef_save(:,:,:)
  real(dp), allocatable :: cdet_save(:,:)
  real(dp), allocatable :: ccsf_save(:,:)

  integer :: mparmja_best, mparmjb_best, mparmjc_best
  real(dp), allocatable :: a4_best(:,:,:)
  real(dp), allocatable :: b_best(:,:,:)
  real(dp), allocatable :: c_best(:,:,:)
  real(dp), allocatable :: coef_best(:,:,:)
  real(dp), allocatable :: cdet_best(:,:)
  real(dp), allocatable :: ccsf_best(:,:)

  integer :: nciterm_sav, norbterm_sav, nparmd_sav, nparmj_sav, nreduced_sav, nadorb_sav

  ! from dmc/nonloc_grid.f
  real(dp), allocatable :: t_vpsp_save(:, :, :)

  save
end module
