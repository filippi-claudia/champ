!> Contains subroutines to allocate/deallocate VMC and DMC calculations.
!>@param allocate_vmc Allocates memory for VMC calculations.
!>@param deallocate_vmc Deallocates memory for VMC calculations.
!>@param allocate_dmc Allocates memory for DMC calculations.
!>@param deallocate_dmc Deallocates memory for DMC calculations.
module allocation_mod
      use jastrow, only: allocate_m_jastrow,deallocate_m_jastrow
      use m_basis, only: allocate_m_basis,deallocate_m_basis
      use m_common, only: allocate_m_common,deallocate_m_common
      use m_control, only: allocate_m_control,deallocate_m_control
      use m_deriv, only: allocate_m_deriv,deallocate_m_deriv
      use m_efield, only: allocate_m_efield,deallocate_m_efield
      use m_estimators, only: allocate_m_estimators
      use m_estimators, only: deallocate_m_estimators
      use m_ewald, only: deallocate_m_ewald
      use m_force, only: allocate_m_force,deallocate_m_force
      use m_gradhess, only: allocate_m_gradhess,deallocate_m_gradhess
      use m_grdnt, only: allocate_m_grdnt,deallocate_m_grdnt
      use m_grid,  only: allocate_m_grid,deallocate_m_grid
      use m_mixderiv, only: allocate_m_mixderiv,deallocate_m_mixderiv
      use m_mmpol, only: allocate_m_mmpol,deallocate_m_mmpol
      use m_mstates, only: allocate_m_mstates,deallocate_m_mstates
      use m_optci, only: allocate_m_optci,deallocate_m_optci
      use m_optorb, only: allocate_m_optorb,deallocate_m_optorb
      use m_optwf, only: allocate_m_optwf,deallocate_m_optwf
      use m_pcm,   only: allocate_m_pcm,deallocate_m_pcm
      use m_prop,  only: allocate_m_prop,deallocate_m_prop
      use m_pseudo, only: allocate_m_pseudo,deallocate_m_pseudo
      use m_sampling, only: allocate_m_sampling,deallocate_m_sampling
      use m_sr,    only: allocate_m_sr,deallocate_m_sr
      use m_state_avrg, only: allocate_m_state_avrg
      use m_state_avrg, only: deallocate_m_state_avrg

implicit none
public
contains
  !> Subroutines to allocate memory of specific VMC-related variables.
  !>@param allocate_m_common()
  !>@param allocate_m_basis
  !>@param allocate_m_control
  !>@param allocate_m_deriv
  !>@param allocate_m_efield
  !>@param allocate_m_estimators
  !>@param allocate_m_force
  !>@param allocate_m_gradhess
  !>@param allocate_m_grdnt
  !>@param allocate_m_grid
  !>@param allocate_m_jastrow
  !>@param allocate_m_mixderiv
  !>@param allocate_m_mmpol
  !>@param allocate_m_mstates
  !>@param allocate_m_optci
  !>@param allocate_m_optorb
  !>@param allocate_m_optwf
  !>@param allocate_m_pcm
  !>@param allocate_m_prop
  !>@param allocate_m_pseudo
  !>@param allocate_m_sampling
  !>@param allocate_m_sr
  !>@param allocate_m_state_avrg
  subroutine allocate_vmc()
    call allocate_m_common()
    call allocate_m_basis
    call allocate_m_control
    call allocate_m_deriv
    call allocate_m_efield
    call allocate_m_estimators
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

  !> Subroutines to deallocate memory of specific VMC-related variables.
  !>@param deallocate_m_common
  !>@param deallocate_m_basis
  !>@param deallocate_m_control
  !>@param deallocate_m_deriv
  !>@param deallocate_m_efield
  !>@param deallocate_m_estimators
  !>@param deallocate_m_ewald
  !>@param deallocate_m_force
  !>@param deallocate_m_gradhess
  !>@param deallocate_m_grdnt
  !>@param deallocate_m_grid
  !>@param deallocate_m_jastrow
  !>@param deallocate_m_mixderiv
  !>@param deallocate_m_mmpol
  !>@param deallocate_m_mstates
  !>@param deallocate_m_optci
  !>@param deallocate_m_optorb
  !>@param deallocate_m_optwf
  !>@param deallocate_m_pcm
  !>@param deallocate_m_prop
  !>@param deallocate_m_pseudo
  !>@param deallocate_m_sampling
  !>@param deallocate_m_sr
  !>@param deallocate_m_state_avrg
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

  !> Subroutines to allocate memory of specific DMC-related variables.
  !>@param allocate_iage()
  !>@param allocate_contrldmc()
  !>@param allocate_config_dmc()
  !>@param allocate_estsum_dmc()
  !>@param allocate_estcum_dmc()
  !>@param allocate_est2cm_dmc()
  !>@param allocate_derivest()
  !>@param allocate_branch()
  !>@param allocate_c_averages()
  !>@param allocate_c_averages_index()
  !>@param allocate_jacobsave()
  !>@param allocate_velratio()
  !>@param allocate_da_branch()
  !>@param allocate_pathak()
  subroutine allocate_dmc()
    use age,     only: allocate_iage
    use branch,  only: allocate_branch
    use c_averages, only: allocate_c_averages
    use c_averages_index, only: allocate_c_averages_index
    use config,  only: allocate_config_dmc
    use contrldmc, only: allocate_contrldmc
    use derivest, only: allocate_derivest
    use est2cm,  only: allocate_est2cm_dmc
    use estcum,  only: allocate_estcum_dmc
    use estsum,  only: allocate_estsum_dmc
    use fragments, only: allocate_fragments, nfrag
    use jacobsave, only: allocate_jacobsave
    use velratio, only: allocate_velratio
    use vd_mod, only: allocate_da_branch
    use pathak_mod, only: allocate_pathak

    implicit none

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
    call allocate_da_branch()
    call allocate_pathak()
    call allocate_fragments()
  end subroutine allocate_dmc

  !> Subroutines to deallocate memory of specific DMC-related variables.
  !> Deallocation will ocurr before ending of dmc/main.f90.
  !>@param deallocate_iage()
  !>@param deallocate_contrldmc()
  !>@param deallocate_config_dmc()
  !>@param deallocate_estsum_dmc()
  !>@param deallocate_estcum_dmc()
  !>@param deallocate_est2cm_dmc()
  !>@param deallocate_derivest()
  !>@param deallocate_branch()
  !>@param deallocate_c_averages()
  !>@param deallocate_c_averages_index()
  !>@param deallocate_jacobsave()
  !>@param deallocate_velratio()
  !>@param deallocate_da_branch()
  !>@param deallocate_pathak()
  subroutine deallocate_dmc()
    use age,     only: deallocate_iage
    use branch,  only: deallocate_branch
    use c_averages, only: deallocate_c_averages
    use c_averages_index, only: deallocate_c_averages_index
    use config,  only: deallocate_config_dmc
    use contrldmc, only: deallocate_contrldmc
    use derivest, only: deallocate_derivest
    use est2cm,  only: deallocate_est2cm_dmc
    use estcum,  only: deallocate_estcum_dmc
    use estsum,  only: deallocate_estsum_dmc
    use fragments,  only: deallocate_fragments
    use jacobsave, only: deallocate_jacobsave
    use velratio, only: deallocate_velratio
    use vd_mod, only: deallocate_da_branch
    use pathak_mod, only: deallocate_pathak

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
    call deallocate_da_branch()
    call deallocate_pathak()
    call deallocate_fragments()
  end subroutine deallocate_dmc
end module allocation_mod
