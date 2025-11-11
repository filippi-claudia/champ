!> @brief Module for centralized memory allocation and deallocation.
!> @author CHAMP developers
!> @date 2022
!>
!> @details This module provides unified interfaces for allocating and deallocating
!> memory for all data structures used in VMC and DMC calculations. It orchestrates
!> the memory management across multiple modules, ensuring proper initialization order
!> and complete cleanup.
!>
!> The module contains:
!> - allocate_vmc(): Allocates all VMC-related arrays
!> - deallocate_vmc(): Deallocates all VMC-related arrays
!> - allocate_dmc(): Allocates all DMC-specific arrays
!> - deallocate_dmc(): Deallocates all DMC-specific arrays
!>
!> @note Memory allocation should be done after all input parameters are read and
!> before any calculation begins. Deallocation should be done before program termination
!> or when switching between different calculation modes.
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
  !> Allocates memory for all VMC calculation arrays and data structures.
  !>
  !> @details This subroutine systematically allocates memory by calling allocation
  !> routines from 23 modules in the proper initialization order. It performs the following tasks:
  !> - allocate_m_common(): Allocates common arrays (electron-electron/nucleus distances, orbital values, Slater matrices)
  !> - allocate_m_basis: Allocates basis set coefficient and exponent arrays for atomic orbitals
  !> - allocate_m_control: Allocates control parameter arrays for run settings
  !> - allocate_m_deriv: Allocates derivative arrays needed for wavefunction optimization
  !> - allocate_m_efield: Allocates electric field parameter arrays for external field calculations
  !> - allocate_m_estimators: Allocates energy and property estimator accumulator arrays
  !> - allocate_m_force: Allocates force calculation arrays for atomic forces
  !> - allocate_m_gradhess: Allocates gradient and Hessian matrices for optimization
  !> - allocate_m_grdnt: Allocates gradient arrays for energy and wavefunction derivatives
  !> - allocate_m_grid: Allocates spatial grid arrays for density and orbital visualization
  !> - allocate_m_jastrow: Allocates Jastrow factor coefficient and scaling parameter arrays
  !> - allocate_m_mixderiv: Allocates mixed derivative term arrays for coupled optimizations
  !> - allocate_m_mmpol: Allocates molecular mechanics polarization data structures
  !> - allocate_m_mstates: Allocates multiple electronic state arrays for excited states
  !> - allocate_m_optci: Allocates CI coefficient optimization arrays
  !> - allocate_m_optorb: Allocates orbital coefficient optimization arrays
  !> - allocate_m_optwf: Allocates wavefunction optimization control and history arrays
  !> - allocate_m_pcm: Allocates PCM (Polarizable Continuum Model) cavity and surface charge arrays
  !> - allocate_m_prop: Allocates property calculation arrays (dipole, quadrupole, etc.)
  !> - allocate_m_pseudo: Allocates pseudopotential radial grid and coefficient arrays
  !> - allocate_m_sampling: Allocates walker configuration, random number, and statistics arrays
  !> - allocate_m_sr: Allocates Stochastic Reconfiguration overlap and Hamiltonian matrix arrays
  !> - allocate_m_state_avrg: Allocates state-averaged calculation weight and energy arrays
  !>
  !> @note Called from vmc/main.f90 after reading input files and before VMC equilibration.
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

  !> Deallocates memory for all VMC calculation arrays and data structures.
  !>
  !> @details This subroutine systematically frees memory by calling deallocation
  !> routines from 23 modules. It performs the following tasks:
  !> - deallocate_m_common: Frees common arrays (distances, orbitals, Slater matrices)
  !> - deallocate_m_basis: Releases basis set coefficient and exponent arrays
  !> - deallocate_m_control: Frees control parameter arrays
  !> - deallocate_m_deriv: Releases derivative arrays for optimization
  !> - deallocate_m_efield: Frees electric field parameter arrays
  !> - deallocate_m_estimators: Releases energy and property estimator arrays
  !> - deallocate_m_ewald: Frees Ewald summation arrays (k-vectors, structure factors)
  !> - deallocate_m_force: Releases force calculation arrays
  !> - deallocate_m_gradhess: Frees gradient and Hessian matrices
  !> - deallocate_m_grdnt: Releases gradient arrays for various quantities
  !> - deallocate_m_grid: Frees spatial grid arrays
  !> - deallocate_m_jastrow: Releases Jastrow factor coefficient and scaling arrays
  !> - deallocate_m_mixderiv: Frees mixed derivative term arrays
  !> - deallocate_m_mmpol: Releases molecular mechanics polarization data
  !> - deallocate_m_mstates: Frees multiple electronic state arrays
  !> - deallocate_m_optci: Releases CI coefficient optimization arrays
  !> - deallocate_m_optorb: Frees orbital optimization arrays
  !> - deallocate_m_optwf: Releases wavefunction optimization data structures
  !> - deallocate_m_pcm: Frees PCM (Polarizable Continuum Model) cavity and charge arrays
  !> - deallocate_m_prop: Releases property calculation arrays
  !> - deallocate_m_pseudo: Frees pseudopotential radial grid and coefficient arrays
  !> - deallocate_m_sampling: Releases walker configuration and statistics arrays
  !> - deallocate_m_sr: Frees Stochastic Reconfiguration matrix arrays
  !> - deallocate_m_state_avrg: Releases state-averaged calculation weight arrays
  !>
  !> @note Called at the end of VMC run in vmc/main.f90 or before switching to DMC mode.
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

  !> Allocates memory for DMC (Diffusion Monte Carlo) calculation arrays and data structures.
  !>
  !> @details This subroutine prepares DMC-specific memory by calling 15 allocation
  !> routines. It must be called after allocate_vmc() and performs the following tasks:
  !> - allocate_iage(): Allocates walker age arrays for population control tracking
  !> - allocate_contrldmc(): Allocates DMC control parameters (time step, target population)
  !> - allocate_config_dmc(): Allocates DMC walker coordinate and weight arrays
  !> - allocate_estsum_dmc(): Allocates block-averaged DMC estimator arrays
  !> - allocate_estcum_dmc(): Allocates cumulative DMC estimator accumulators
  !> - allocate_est2cm_dmc(): Allocates DMC second moment estimators for error analysis
  !> - allocate_derivest(): Allocates derivative estimator arrays for force calculations
  !> - allocate_branch(): Allocates branching factor and population control arrays
  !> - allocate_c_averages(): Allocates core estimator average arrays (energy, variance)
  !> - allocate_c_averages_index(): Allocates indexing arrays for core estimator access
  !> - allocate_jacobsave(): Allocates Jacobian determinant storage for T-moves
  !> - allocate_velratio(): Allocates velocity ratio arrays for T-move acceptance
  !> - allocate_da_branch(): Allocates variational derivative arrays with branching weights
  !> - allocate_pathak(): Allocates Pathak T-move correction arrays
  !> - allocate_fragments(): Allocates molecular fragment tracking arrays (if nfrag > 0)
  !>
  !> @note Called from dmc/main.f90 after allocate_vmc() and before DMC equilibration.
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

  !> Deallocates memory for all DMC calculation arrays and data structures.
  !>
  !> @details This subroutine systematically frees DMC-specific memory by calling
  !> 15 deallocation routines. It performs the following cleanup tasks:
  !> - deallocate_iage(): Frees walker age tracking arrays
  !> - deallocate_contrldmc(): Releases DMC control parameter arrays
  !> - deallocate_config_dmc(): Frees DMC walker configuration and weight arrays
  !> - deallocate_estsum_dmc(): Releases block-averaged estimator arrays
  !> - deallocate_estcum_dmc(): Frees cumulative estimator accumulator arrays
  !> - deallocate_est2cm_dmc(): Releases second moment estimator arrays
  !> - deallocate_derivest(): Frees derivative estimator arrays used for forces
  !> - deallocate_branch(): Releases branching factor and population control arrays
  !> - deallocate_c_averages(): Frees core estimator average arrays
  !> - deallocate_c_averages_index(): Releases core estimator indexing arrays
  !> - deallocate_jacobsave(): Frees Jacobian determinant storage arrays
  !> - deallocate_velratio(): Releases velocity ratio arrays for T-moves
  !> - deallocate_da_branch(): Frees variational derivative branching arrays
  !> - deallocate_pathak(): Releases Pathak T-move correction arrays
  !> - deallocate_fragments(): Frees molecular fragment arrays
  !>
  !> @note Called near the end of dmc/main.f90 to clean up DMC arrays before
  !> program termination or before switching back to VMC mode.
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
