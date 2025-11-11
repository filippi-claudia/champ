!> @brief Module for atomic derivative energy and wavefunction accumulators.
!> @author CHAMP developers
!> @date 2020
!>
!> @details This module stores statistical accumulators for energy and wavefunction
!> derivatives with respect to atomic positions. These quantities are used for
!> calculating atomic forces in VMC and DMC. The accumulators include sums,
!> cumulative averages, and second moments for error estimation.
!>
!> Key arrays (da = derivative with respect to atom):
!> - da_energy_*: Energy derivative accumulators
!> - da_psi_*: Wavefunction derivative accumulators  
!> - Dimensions: (3 Cartesian components, ncent_tot centers, PTH paths)
!>
!> @note Used in force calculations and geometry optimization.
module da_energy_sumcum
      use precision_kinds, only: dp

    implicit none

    !> Second moment of energy derivatives for error analysis, dimension (3, ncent_tot, PTH)
    real(dp), dimension(:, :, :), allocatable :: da_energy_cm2

    !> Cumulative average of energy derivatives, dimension (3, ncent_tot, PTH)
    real(dp), dimension(:, :, :), allocatable :: da_energy_cum

    !> Block sum of energy derivatives, dimension (3, ncent_tot, PTH)
    real(dp), dimension(:, :, :), allocatable :: da_energy_sum

    !> Cumulative average of wavefunction derivatives, dimension (3, ncent_tot, PTH)
    real(dp), dimension(:, :, :), allocatable :: da_psi_cum

    !> Block sum of wavefunction derivatives, dimension (3, ncent_tot, PTH)
    real(dp), dimension(:, :, :), allocatable :: da_psi_sum

    !> Block sum of energy times wavefunction derivatives, dimension (3, ncent_tot, PTH)
    real(dp), dimension(:, :, :), allocatable :: da_energy_psi_sum

    private
    public :: da_energy_cm2, da_energy_cum, da_energy_sum
    public :: da_energy_psi_sum, da_psi_cum, da_psi_sum
    public :: allocate_da_energy_sumcum, deallocate_da_energy_sumcum
    save
contains
    !> Allocates memory for atomic derivative energy accumulators.
    !>
    !> @details Allocates arrays for storing statistical accumulators of energy
    !> and wavefunction derivatives with respect to atomic positions. All arrays
    !> have dimensions (3, ncent_tot, PTH) for 3 Cartesian components, all centers,
    !> and multiple force paths.
    !>
    !> @note Called during force calculation initialization.
    subroutine allocate_da_energy_sumcum()
      use system,  only: ncent_tot
      use force_pth, only: PTH
        if (.not. allocated(da_energy_cm2)) allocate (da_energy_cm2(3, ncent_tot, PTH))
        if (.not. allocated(da_energy_cum)) allocate (da_energy_cum(3, ncent_tot, PTH))
        if (.not. allocated(da_energy_sum)) allocate (da_energy_sum(3, ncent_tot, PTH))
        if (.not. allocated(da_psi_cum)) allocate (da_psi_cum(3, ncent_tot, PTH))
        if (.not. allocated(da_psi_sum)) allocate (da_psi_sum(3, ncent_tot, PTH))
        if (.not. allocated(da_energy_psi_sum)) allocate (da_energy_psi_sum(3, ncent_tot, PTH))
    end subroutine allocate_da_energy_sumcum

    !> Deallocates memory for atomic derivative energy accumulators.
    subroutine deallocate_da_energy_sumcum()
        if (allocated(da_psi_sum)) deallocate(da_psi_sum)
        if (allocated(da_psi_cum)) deallocate(da_psi_cum)
        if (allocated(da_energy_sum)) deallocate(da_energy_sum)
        if (allocated(da_energy_cum)) deallocate(da_energy_cum)
        if (allocated(da_energy_cm2)) deallocate(da_energy_cm2)
        if (allocated(da_energy_psi_sum)) deallocate(da_energy_psi_sum)
    end subroutine deallocate_da_energy_sumcum

end module da_energy_sumcum

!> @brief Module for Jastrow factor derivatives with respect to atomic positions.
!> @author CHAMP developers
!> @date 2020
!>
!> @details This module stores derivatives of the Jastrow correlation factor
!> with respect to atomic positions. 
!>
!> Key arrays:
!> - da_d2j: Laplacian of Jastrow derivatives
!> - da_j: Jastrow factor derivatives
!> - da_vj: Gradient of Jastrow derivatives
!>
module da_jastrow
      use precision_kinds, only: dp

    implicit none

    !> Laplacian of Jastrow derivatives with respect to atomic positions, dimension (3, ncent_tot)
    real(dp), dimension(:, :), allocatable :: da_d2j

    !> Jastrow factor derivatives for electron pairs, dimension (3, nelec, nelec, ncent_tot)
    real(dp), dimension(:, :, :, :), allocatable :: da_j

    !> Gradient of Jastrow derivatives, dimension (3, 3, nelec, ncent_tot)
    real(dp), dimension(:, :, :, :), allocatable :: da_vj

    private
    public   ::  da_d2j, da_j, da_vj
    public :: allocate_da_jastrow, deallocate_da_jastrow
    save
contains
    !> Allocates memory for Jastrow factor derivative arrays.
    !>
    !> @details Allocates arrays for storing Jastrow factor derivatives with
    !> respect to atomic positions, needed for force calculations.
    !>
    subroutine allocate_da_jastrow()
      use system,  only: ncent_tot,nelec
        if (.not. allocated(da_d2j)) allocate (da_d2j(3, ncent_tot))
        if (.not. allocated(da_j)) allocate (da_j(3, nelec, nelec, ncent_tot))
        if (.not. allocated(da_vj)) allocate (da_vj(3, 3, nelec, ncent_tot))
    end subroutine allocate_da_jastrow

    !> Deallocates memory for Jastrow factor derivative arrays.
    subroutine deallocate_da_jastrow
        if (allocated(da_vj)) deallocate(da_vj)
        if (allocated(da_j)) deallocate(da_j)
        if (allocated(da_d2j)) deallocate(da_d2j)
    end subroutine deallocate_da_jastrow

end module da_jastrow

!> @brief Module for orbital value derivatives with respect to atomic positions.
!> @author CHAMP developers
!> @date 2022
!>
!> @details This module stores derivatives of molecular orbital values with
!> respect to atomic positions. 
!>
!> Key arrays:
!> - da_orb: Orbital value derivatives
!> - da_dorb: Gradient of orbital derivatives
!> - da_d2orb: Laplacian of orbital derivatives
!>
!> @note Critical for Pulay force corrections in basis-set-dependent calculations.
module da_orbval
      use precision_kinds, only: dp
      use vmc_mod, only: norb_tot

    implicit none

    !> Laplacian of orbital derivatives, dimension (3, nelec, norb_tot, ncent_tot)
    real(dp), dimension(:, :, :, :), allocatable :: da_d2orb

    !> Gradient of orbital derivatives, dimension (3, 3, nelec, norb_tot, ncent_tot)
    real(dp), dimension(:, :, :, :, :), allocatable :: da_dorb

    !> Orbital value derivatives with respect to atomic positions, dimension (3, nelec, norb_tot, ncent_tot)
    real(dp), dimension(:, :, :, :), allocatable :: da_orb

    private
    public   ::  da_d2orb, da_dorb, da_orb
    public :: allocate_da_orbval, deallocate_da_orbval
    save
contains
    !> Allocates memory for orbital derivative arrays.
    !>
    !> @details Allocates arrays for storing orbital value derivatives and their
    !> spatial derivatives with respect to atomic positions.
    !>
    !> @note Called during force calculation initialization for Pulay forces.
    subroutine allocate_da_orbval()
      use system, only: ncent, ncent_tot, nelec
      use vmc_mod, only: norb_tot
        if (.not. allocated(da_d2orb)) allocate (da_d2orb(3, nelec, norb_tot, ncent_tot))
        if (.not. allocated(da_dorb)) allocate (da_dorb(3, 3, nelec, norb_tot, ncent_tot))
        if (.not. allocated(da_orb)) allocate (da_orb(3, nelec, norb_tot, ncent_tot))
    end subroutine allocate_da_orbval

    !> Deallocates memory for orbital derivative arrays.
    subroutine deallocate_da_orbval()
        if (allocated(da_orb)) deallocate(da_orb)
        if (allocated(da_dorb)) deallocate(da_dorb)
        if (allocated(da_d2orb)) deallocate(da_d2orb)
    end subroutine deallocate_da_orbval

end module da_orbval

!> @brief Module for pseudopotential derivatives with respect to atomic positions.
!> @author CHAMP developers
!> @date 2022
!>
!> @details This module stores derivatives of pseudopotential energies and
!> matrix elements with respect to atomic positions. These are needed for
!> calculating forces in pseudopotential calculations, including both local
!> and non-local contributions.
!>
!> Key arrays:
!> - da_pecent: Local pseudopotential energy derivatives per center
!> - da_vps: Pseudopotential matrix element derivatives
!> - da_nonloc: Non-local pseudopotential contribution derivatives
!> - da_pe_en: Total pseudopotential energy derivatives
!>
!> @note Used in pseudopotential force calculations for valence-only QMC.
module da_pseudo

    use precision_kinds, only: dp

    implicit none

    !> Local pseudopotential energy derivatives per center, dimension (3, ncent_tot)
    real(dp), dimension(:, :), allocatable :: da_pecent

    !> Pseudopotential matrix element derivatives, dimension (3, nelec, ncent_tot, MPS_L)
    real(dp), dimension(:, :, :, :), allocatable :: da_vps

    !> Non-local pseudopotential contribution derivatives, dimension (3, ncent_tot)
    real(dp), dimension(:, :), allocatable :: da_nonloc

    !> Total pseudopotential energy derivatives per center, dimension (3, ncent_tot)
    real(dp), dimension(:, :), allocatable :: da_pe_en

    private
    public   :: da_pecent, da_vps, da_nonloc, da_pe_en
    public :: allocate_da_pseudo, deallocate_da_pseudo
    save
contains
    !> Allocates memory for pseudopotential derivative arrays.
    !>
    !> @details Allocates arrays for storing pseudopotential energy and matrix
    !> element derivatives. Initializes da_nonloc to zero.
    !>
    !> @note Called during pseudopotential force calculation initialization.
    subroutine allocate_da_pseudo()
      use pseudo_mod, only: MPS_L
      use system, only: ncent_tot, nelec
        if (.not. allocated(da_pecent)) allocate (da_pecent(3, ncent_tot))
        if (.not. allocated(da_vps)) allocate (da_vps(3, nelec, ncent_tot, MPS_L))
        if (.not. allocated(da_nonloc)) allocate (da_nonloc(3, ncent_tot))
        if (.not. allocated(da_pe_en)) allocate (da_pe_en(3, ncent_tot))

        da_nonloc = 0.0D0

    end subroutine allocate_da_pseudo

    !> Deallocates memory for pseudopotential derivative arrays.
    subroutine deallocate_da_pseudo()
        if (allocated(da_nonloc)) deallocate(da_nonloc)
        if (allocated(da_vps)) deallocate(da_vps)
        if (allocated(da_pecent)) deallocate(da_pecent)
        if (allocated(da_pe_en)) deallocate(da_pe_en)
    end subroutine deallocate_da_pseudo

end module da_pseudo

!> @brief Module for current step energy and wavefunction derivatives.
!> @author CHAMP developers
!> @date 2022
!>
!> @details This module stores the instantaneous (current step) values of energy
!> and wavefunction derivatives with respect to atomic positions. These are used
!> for accumulating force statistics during Monte Carlo sampling.
!>
!> @note Updated at each MC step during force calculations.
module da_energy_now
      use precision_kinds, only: dp

    implicit none

    !> Current step energy derivatives with respect to atomic positions, dimension (3, ncent_tot)
    real(dp), dimension(:, :), allocatable :: da_energy

    !> Current step wavefunction derivatives with respect to atomic positions, dimension (3, ncent_tot)
    real(dp), dimension(:, :), allocatable :: da_psi

    private
    public   ::  da_energy, da_psi
    public :: allocate_da_energy_now, deallocate_da_energy_now
    save
contains
    !> Allocates memory for current step derivative arrays.
    !>
    !> @details Allocates arrays for storing instantaneous energy and wavefunction
    !> derivatives at the current Monte Carlo step.
    subroutine allocate_da_energy_now()
      use system, only: ncent_tot
        if (.not. allocated(da_energy)) allocate (da_energy(3, ncent_tot))
        if (.not. allocated(da_psi)) allocate (da_psi(3, ncent_tot))
    end subroutine allocate_da_energy_now

    !> Deallocates memory for current step derivative arrays.
    subroutine deallocate_da_energy_now()
        if (allocated(da_psi)) deallocate(da_psi)
        if (allocated(da_energy)) deallocate(da_energy)
    end subroutine deallocate_da_energy_now

end module da_energy_now

!> @brief Module for energy derivatives with respect to Jastrow parameters.
!> @author CHAMP developers
!> @date 2022
!>
!> @details This module stores energy derivatives with respect to Jastrow factor
!> parameters. These derivatives are essential for wavefunction optimization using
!> gradient-based methods.
!>
!> @note Used in Jastrow parameter optimization (VMC optimization).
module deloc_dj_m
    use mstates_mod, only: MSTATES
    use optwf_parms, only: nparmj
    use precision_kinds, only: dp

    implicit none

    !> Energy derivatives with respect to Jastrow parameters, dimension (nparmj, MSTATES)
    real(dp), dimension(:, :), allocatable :: denergy

    private
    public :: denergy
    public :: allocate_deloc_dj_m, deallocate_deloc_dj_m
    save
contains
    !> Allocates memory for Jastrow parameter derivative arrays.
    !>
    !> @details Allocates arrays for storing energy derivatives with respect to
    !> Jastrow parameters for wavefunction optimization.
    subroutine allocate_deloc_dj_m()
      use mstates_mod, only: MSTATES
      use optwf_parms, only: nparmj
        if (.not. allocated(denergy)) allocate (denergy(nparmj, MSTATES))
    end subroutine allocate_deloc_dj_m

    !> Deallocates memory for Jastrow parameter derivative arrays.
    subroutine deallocate_deloc_dj_m()
        if (allocated(denergy)) deallocate(denergy)
    end subroutine deallocate_deloc_dj_m

end module deloc_dj_m

!> @brief Module for per-determinant energy derivatives.
!> @author CHAMP developers
!> @date 2022
!>
!> @details This module stores energy derivatives calculated separately for each
!> determinant in a multi-determinant wavefunction. Used in optimizations that
!> require determinant-resolved derivative information.
!>
!> @note Important for multi-determinant wavefunction optimizations.
module denergy_det_m
      use precision_kinds, only: dp

    implicit none

    !> Energy derivatives per determinant and block, dimension (ndet, 2, nbjx)
    real(dp), dimension(:, :, :), allocatable :: denergy_det

    private
    public :: denergy_det
    public :: allocate_denergy_det_m, deallocate_denergy_det_m
    save
contains
    !> Allocates memory for per-determinant energy derivative arrays.
    !>
    !> @details Allocates arrays for storing energy derivatives separately for
    !> each determinant in multi-determinant calculations.
    subroutine allocate_denergy_det_m()
      use slater, only: ndet
      use vmc_mod, only: nbjx
        if (.not. allocated(denergy_det)) allocate (denergy_det(ndet, 2, nbjx))
    end subroutine allocate_denergy_det_m

    !> Deallocates memory for per-determinant energy derivative arrays.
    subroutine deallocate_denergy_det_m()
        if (allocated(denergy_det)) deallocate(denergy_det)
    end subroutine deallocate_denergy_det_m

end module denergy_det_m

!> @brief Module for Jastrow factor parameter derivatives.
!> @author CHAMP developers
!> @date 2022
!>
!> @details This module stores derivatives of the Jastrow factor with respect to
!> its variational parameters. These include the Jastrow value derivatives (g),
!> Laplacian derivatives (d2g), and two-electron derivatives (go). Essential for
!> computing gradients in Jastrow parameter optimization.
!>
!> @note Core module for Jastrow factor optimization in VMC.
module derivjas
    use precision_kinds, only: dp

    implicit none

    !> Laplacian of Jastrow factor with respect to parameters, dimension (nparmj, nwftypejas)
    real(dp), dimension(:, :), allocatable :: d2g

    !> Gradient of Jastrow derivatives with respect to parameters, dimension (3, nelec, nparmj, nwftypejas)
    real(dp), dimension(:, :, :, :), allocatable :: g

    !> Two-electron Jastrow derivatives with respect to parameters, dimension (nelec, nelec, nparmj, nwftypejas)
    real(dp), dimension(:, :, :, :), allocatable :: go

    !> Jastrow factor value derivatives with respect to parameters, dimension (nparmj, nwftypejas)
    real(dp), dimension(:, :), allocatable :: gvalue

    private
    public   :: d2g, g, go, gvalue
    public :: allocate_derivjas, deallocate_derivjas
    save
contains
    !> Allocates memory for Jastrow parameter derivative arrays.
    !>
    !> @details Allocates arrays for storing Jastrow factor derivatives with
    !> respect to variational parameters for optimization.
    subroutine allocate_derivjas()
      use optwf_parms, only: nparmj
      use system, only: nelec
      use vmc_mod, only: nwftypejas
        if (.not. allocated(d2g)) allocate (d2g(nparmj, nwftypejas))
        if (.not. allocated(g)) allocate (g(3, nelec, nparmj, nwftypejas))
        if (.not. allocated(go)) allocate (go(nelec, nelec, nparmj, nwftypejas))
        if (.not. allocated(gvalue)) allocate (gvalue(nparmj,nwftypejas))
    end subroutine allocate_derivjas

    !> Deallocates memory for Jastrow parameter derivative arrays.
    subroutine deallocate_derivjas()
        if (allocated(gvalue)) deallocate(gvalue)
        if (allocated(go)) deallocate(go)
        if (allocated(g)) deallocate(g)
        if (allocated(d2g)) deallocate(d2g)
    end subroutine deallocate_derivjas

end module derivjas

!> @brief Module for orbital derivative workspace indices.
!> @author CHAMP developers
!> @date 2022
!>
!> @details This module stores workspace indices for orbital derivatives in
!> multi-determinant calculations. Used to track which orbitals need derivative
!> evaluation for each determinant.
!>
!> @note Used internally for efficient derivative calculations.
module dorb_m

    implicit none

    !> Orbital derivative workspace indices per electron and determinant, dimension (nelec, ndet)
    integer, dimension(:, :), allocatable :: iworbd

    private
    public :: iworbd
    public :: allocate_dorb_m, deallocate_dorb_m
    save

contains

    !> Allocates memory for orbital derivative index arrays.
    !>
    !> @details Allocates workspace index array for orbital derivatives, initialized to zero.
    subroutine allocate_dorb_m()
      use slater, only: ndet
      use system, only: nelec
        if (.not. allocated(iworbd)) allocate (iworbd(nelec, ndet), source=0)
    end subroutine allocate_dorb_m

    !> Deallocates memory for orbital derivative index arrays.
    subroutine deallocate_dorb_m()
        if (allocated(iworbd)) deallocate(iworbd)
    end subroutine deallocate_dorb_m

end module dorb_m

!> @brief Module for Jastrow parameter scaling factors.
!> @author CHAMP developers
!> @date 2022
!>
!> @details This module stores scaling factors for Jastrow parameters.
!> These factors (d1d2 and d2d2) are used in the parameterization of correlation
!> factors, particularly for handling different center types (a) and spin types (b).
!>
module ijasnonlin
      use precision_kinds, only: dp

    implicit none

    !> First-order to second-order scaling factors for center types, dimension (nctype_tot)
    real(dp), dimension(:), allocatable :: d1d2a

    !> First-order to second-order scaling factors for spin types, dimension (2)
    real(dp), dimension(:), allocatable :: d1d2b

    !> Second-order to second-order scaling factors for center types, dimension (nctype_tot)
    real(dp), dimension(:), allocatable :: d2d2a

    !> Second-order to second-order scaling factors for spin types, dimension (2)
    real(dp), dimension(:), allocatable :: d2d2b

    private
    public :: d1d2a, d1d2b, d2d2a, d2d2b
    public :: allocate_ijasnonlin, deallocate_ijasnonlin
    save
contains
    !> Allocates memory for Jastrow scaling factor arrays.
    !>
    !> @details Allocates arrays for parameter scaling in Jastrow factors.
    subroutine allocate_ijasnonlin()
      use system, only: nctype_tot
        if (.not. allocated(d1d2a)) allocate (d1d2a(nctype_tot))
        if (.not. allocated(d1d2b)) allocate (d1d2b(2))
        if (.not. allocated(d2d2a)) allocate (d2d2a(nctype_tot))
        if (.not. allocated(d2d2b)) allocate (d2d2b(2))
    end subroutine allocate_ijasnonlin

    !> Deallocates memory for non-linear Jastrow scaling factor arrays.
    subroutine deallocate_ijasnonlin()
        if (allocated(d2d2b)) deallocate(d2d2b)
        if (allocated(d2d2a)) deallocate(d2d2a)
        if (allocated(d1d2b)) deallocate(d1d2b)
        if (allocated(d1d2a)) deallocate(d1d2a)
    end subroutine deallocate_ijasnonlin

end module ijasnonlin

!> @brief Module for DMC derivative estimators and accumulators.
!> @author CHAMP developers
!> @date 2022
!>
!> @details This module stores statistical accumulators for force derivatives
!> in Diffusion Monte Carlo (DMC) calculations. These include block sums,
!> cumulative averages, second moments, and total averages across multiple paths.
!>
!> Key arrays:
!> - derivsum: Block sums of force derivatives
!> - derivcum: Cumulative averages of force derivatives
!> - derivcm2: Second moments for error estimation
!> - derivtotave: Total averaged forces
!>
module derivest
   use precision_kinds, only: dp
   use system,  only: ncent_tot
   use force_pth, only: PTH

   implicit none

    !> Second moment of force derivatives for error analysis, dimension (3, ncent_tot, PTH)
    real(dp), dimension(:,:,:), allocatable :: derivcm2

    !> Cumulative average of force derivatives, dimension (3, 3, ncent_tot, PTH)
    real(dp), dimension(:,:,:,:), allocatable :: derivcum

    !> Block sum of force derivatives, dimension (3, 3, ncent_tot, PTH)
    real(dp), dimension(:,:,:,:), allocatable :: derivsum

    !> Total averaged forces per center and path, dimension (3, ncent_tot, PTH)
    real(dp), dimension(:,:,:), allocatable :: derivtotave

    private
    public :: derivcm2, derivcum, derivsum, derivtotave
    public :: allocate_derivest, deallocate_derivest
    save

contains
    !> Allocates memory for DMC derivative estimator arrays.
    !>
    !> @details Allocates arrays for DMC force derivative accumulators across
    !> multiple force paths for statistical sampling.
    subroutine allocate_derivest()
        if (.not. allocated(derivcm2)) allocate(derivcm2(3,ncent_tot,PTH))
        if (.not. allocated(derivcum)) allocate(derivcum(3,3,ncent_tot, PTH))
        if (.not. allocated(derivsum)) allocate(derivsum(3,3,ncent_tot, PTH))
        if (.not. allocated(derivtotave)) allocate(derivtotave(3,ncent_tot, PTH))
    end subroutine allocate_derivest

    !> Deallocates memory for DMC derivative estimator arrays.
    subroutine deallocate_derivest
        if (allocated(derivcm2)) deallocate(derivcm2)
        if (allocated(derivcum)) deallocate(derivcum)
        if (allocated(derivsum)) deallocate(derivsum)
        if (allocated(derivtotave)) deallocate(derivtotave)
    end subroutine deallocate_derivest
 end module derivest

!> @brief Master module for coordinated derivative memory management.
!> @author CHAMP developers
!> @date 2022
!>
!> @details This module provides unified allocation and deallocation interfaces
!> for all derivative-related modules. It orchestrates memory management for:
!> - Atomic position derivatives (forces)
!> - Jastrow parameter derivatives (optimization)
!> - Orbital derivatives (Pulay forces)
!> - Pseudopotential derivatives
!> - Energy and wavefunction derivative accumulators
!>
!> @note Primary interface for all derivative calculation memory management.
module m_deriv
contains
!> Allocates memory for all derivative-related modules.
!>
!> @details This subroutine calls allocation routines from 10 derivative modules:
!> - allocate_da_energy_sumcum(): Atomic force accumulators
!> - allocate_da_jastrow(): Jastrow factor atomic derivatives
!> - allocate_da_orbval(): Orbital value atomic derivatives
!> - allocate_da_pseudo(): Pseudopotential atomic derivatives
!> - allocate_da_energy_now(): Current step derivatives
!> - allocate_deloc_dj_m(): Jastrow parameter derivatives
!> - allocate_denergy_det_m(): Per-determinant derivatives
!> - allocate_derivjas(): Jastrow parameter derivative components
!> - allocate_dorb_m(): Orbital derivative workspace
!> - allocate_ijasnonlin(): Non-linear Jastrow scaling factors
!>
!> @note Called from allocate_vmc() when derivatives/forces are needed.
subroutine allocate_m_deriv()
    use da_energy_now, only: allocate_da_energy_now
    use da_energy_sumcum, only: allocate_da_energy_sumcum
    use da_jastrow, only: allocate_da_jastrow
    use da_orbval, only: allocate_da_orbval
    use da_pseudo, only: allocate_da_pseudo
    use deloc_dj_m, only: allocate_deloc_dj_m
    use denergy_det_m, only: allocate_denergy_det_m
    use derivjas, only: allocate_derivjas
    use dorb_m, only: allocate_dorb_m
    use ijasnonlin, only: allocate_ijasnonlin

    implicit none

    call allocate_da_energy_sumcum()
    call allocate_da_jastrow()
    call allocate_da_orbval()
    call allocate_da_pseudo()
    call allocate_da_energy_now()
    call allocate_deloc_dj_m()
    call allocate_denergy_det_m()
    call allocate_derivjas()
    call allocate_dorb_m()
    call allocate_ijasnonlin()
end subroutine allocate_m_deriv

!> Deallocates memory for all derivative-related modules.
!>
!> @details This subroutine calls deallocation routines from all 10 derivative
!> modules to free memory used for force calculations and wavefunction optimization.
!>
!> @note Called from deallocate_vmc() at program termination.
subroutine deallocate_m_deriv()
    use da_energy_now, only: deallocate_da_energy_now
    use da_energy_sumcum, only: deallocate_da_energy_sumcum
    use da_jastrow, only: deallocate_da_jastrow
    use da_orbval, only: deallocate_da_orbval
    use da_pseudo, only: deallocate_da_pseudo
    use deloc_dj_m, only: deallocate_deloc_dj_m
    use denergy_det_m, only: deallocate_denergy_det_m
    use derivjas, only: deallocate_derivjas
    use dorb_m, only: deallocate_dorb_m
    use ijasnonlin, only: deallocate_ijasnonlin

    implicit none

    call deallocate_da_energy_sumcum()
    call deallocate_da_jastrow()
    call deallocate_da_orbval()
    call deallocate_da_pseudo()
    call deallocate_da_energy_now()
    call deallocate_deloc_dj_m()
    call deallocate_denergy_det_m()
    call deallocate_derivjas()
    call deallocate_dorb_m()
    call deallocate_ijasnonlin()
end subroutine deallocate_m_deriv
end module
