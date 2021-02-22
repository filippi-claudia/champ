!> Subroutines to allocate/deallocate memory of specific DMC-related
!> variables.
!> The allocation will occurr only in via read_input.f
subroutine allocate_dmc()
  use age, only: allocate_iage
  use contrldmc, only: allocate_contrldmc 
  use config, only: allocate_config_dmc 
  use estsum, only: allocate_estsum_dmc 
  use estcum, only: allocate_estcum_dmc
  use est2cm, only: allocate_est2cm_dmc
  use derivest, only: allocate_derivest
  use branch, only: allocate_branch
  use c_averages, only: allocate_c_averages
  !> Allocate dmc-related arrays:
  call allocate_iage()
  call allocate_contrldmc()
  call allocate_config_dmc()
  call allocate_estsum_dmc()
  call allocate_estcum_dmc()
  call allocate_est2cm_dmc()
  call allocate_derivest()
  call allocate_branch()
  call allocate_c_averages()
end subroutine allocate_dmc

subroutine deallocate_dmc()
  use age, only: deallocate_iage
  use contrldmc, only: deallocate_contrldmc 
  use config, only: deallocate_config_dmc 
  use estsum, only: deallocate_estsum_dmc 
  use estcum, only: deallocate_estcum_dmc
  use est2cm, only: deallocate_est2cm_dmc 
  use derivest, only: deallocate_derivest
  use branch, only: deallocate_branch
  use c_averages, only: deallocate_c_averages
  !> Deallocate dmc-related arrays:
  call deallocate_iage()
  call deallocate_contrldmc()
  call deallocate_config_dmc()
  call deallocate_estsum_dmc()
  call deallocate_estcum_dmc()
  call deallocate_est2cm_dmc() 
  call deallocate_derivest()
  call deallocate_branch()
  call deallocate_c_averages()
end subroutine deallocate_dmc