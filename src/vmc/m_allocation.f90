module allocation_mod
!> Contains subroutines to allocate/deallocate VMC and DMC calculations.

implicit none
public
contains
  !> Subroutines to allocate/deallocate memory of specific VMC-related
  !> variables.
  !> The allocation will occurr only in via read_input.f
  subroutine allocate_vmc()
  
        call allocate_m_common()
        call allocate_m_basis
        call allocate_m_control
        call allocate_m_deriv
        call allocate_m_efield
        call allocate_m_estimators
        call allocate_m_ewald
        call allocate_m_force
        call allocate_m_gradhess
        call allocate_m_grdnt
        call allocate_m_grid
        call allocate_m_jastrow
        call allocate_m_mixderiv
        call allocate_m_mmpol
        call allocate_m_mstates
        call allocate_m_optci
        call allocate_m_optorb
        call allocate_m_optwf
        call allocate_m_pcm
        call allocate_m_prop
        call allocate_m_pseudo
        call allocate_m_sampling
        call allocate_m_sr
        call allocate_m_state_avrg
  
  end subroutine allocate_vmc
  
  subroutine deallocate_vmc()
  
        call deallocate_m_common
        call deallocate_m_basis
        call deallocate_m_control
        call deallocate_m_deriv
        call deallocate_m_efield
        call deallocate_m_estimators
        call deallocate_m_ewald
        call deallocate_m_force
        call deallocate_m_gradhess
        call deallocate_m_grdnt
        call deallocate_m_grid
        call deallocate_m_jastrow
        call deallocate_m_mixderiv
        call deallocate_m_mmpol
        call deallocate_m_mstates
        call deallocate_m_optci
        call deallocate_m_optorb
        call deallocate_m_optwf
        call deallocate_m_pcm
        call deallocate_m_prop
        call deallocate_m_pseudo
        call deallocate_m_sampling
        call deallocate_m_sr
        call deallocate_m_state_avrg
  
  end subroutine deallocate_vmc
  !> Subroutines to allocate/deallocate memory of specific DMC-related
  !> variables.
  !> Allocation will occurr only via read_input.f.
  !> Deallocation will ocurr before ending of dmc/main.f.
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
    use c_averages_index, only: allocate_c_averages_index
    use jacobsave, only: allocate_jacobsave
    use velratio, only: allocate_velratio

    implicit none
  
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
    call allocate_c_averages_index()
    call allocate_jacobsave()
    call allocate_velratio()
  
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
    use c_averages_index, only: deallocate_c_averages_index
    use jacobsave, only: deallocate_jacobsave
    use velratio, only: deallocate_velratio

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
    call deallocate_c_averages_index()
    call deallocate_jacobsave()
    call deallocate_velratio()
  
  end subroutine deallocate_dmc
end module allocation_mod
