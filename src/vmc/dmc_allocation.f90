!> Subroutines to allocate/deallocate memory in DMC-related
!> variables:
subroutine allocate_m_dmc()
  use age, only: allocate_iage
  use contrldmc, only: allocate_contrldmc 
  use config, only: allocate_config_dmc 
  use estsum, only: allocate_estsum_dmc 
  use estcum, only: allocate_estcum_dmc
  use est2cm, only: allocate_est2cm_dmc 
  !> Allocate dmc-related arrays:
  call allocate_iage()
  call allocate_contrldmc()
  call allocate_config_dmc()
  call allocate_estsum_dmc()
  call allocate_estcum_dmc()
  call allocate_est2cm_dmc() 
end subroutine allocate_m_dmc

subroutine deallocate_m_dmc()
  use age, only: deallocate_iage
  use contrldmc, only: deallocate_contrldmc 
  use config, only: deallocate_config_dmc 
  use estsum, only: deallocate_estsum_dmc 
  use estcum, only: deallocate_estcum_dmc
  use est2cm, only: deallocate_est2cm_dmc 
  !> Deallocate dmc-related arrays:
  call deallocate_iage()
  call deallocate_contrldmc()
  call deallocate_config_dmc()
  call deallocate_estsum_dmc()
  call deallocate_estcum_dmc()
  call deallocate_est2cm_dmc() 
end subroutine deallocate_m_dmc