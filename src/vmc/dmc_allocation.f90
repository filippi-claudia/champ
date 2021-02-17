!> Subroutines to allocate/deallocate memory in DMC-related
!> variables:
subroutine allocate_m_dmc()
  use age, only: allocate_iage
  use contrldmc, only: allocate_contrldmc 
  use config, only: allocate_config_dmc 
  use estsum, only: allocate_estsum_dmc 
  !> Allocate all the dmc_mod arrays:
  call allocate_iage()
  call allocate_contrldmc()
  call allocate_config_dmc()
  call  allocate_estsum_dmc()
end subroutine allocate_m_dmc

subroutine deallocate_m_dmc()
  use age, only: deallocate_iage
  use contrldmc, only: deallocate_contrldmc 
  use config, only: deallocate_config_dmc 
  use estsum, only: deallocate_estsum_dmc 
  !> Deallocate all the dmc_mod arrays:
  call deallocate_iage()
  call deallocate_contrldmc()
  call deallocate_config_dmc()
  call deallocate_estsum_dmc()
end subroutine deallocate_m_dmc