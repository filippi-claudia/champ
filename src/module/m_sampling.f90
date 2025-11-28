!> @brief Module for electron configuration and Monte Carlo sampling data.
!> @author CHAMP developers
!> @date 2020
!>
!> @details This module stores current and proposed electron configurations during
!> Monte Carlo sampling, including positions, velocities, wavefunctions, and energies.
!> It maintains both old (current) and new (proposed) states for Metropolis acceptance testing.
module config
      use mstates_mod, only: MSTATES
      use multiple_geo, only: MFORCE
      use precision_kinds, only: dp

    implicit none

    !> Time step for each electron.
    real(dp), dimension(:), allocatable :: delttn

    !> New (proposed) local energies.
    real(dp), dimension(:), allocatable :: enew

    !> Old (current) local energies for all states.
    real(dp), dimension(:, :), allocatable :: eold

    !> Index of nearest center for new configuration.
    integer, dimension(:), allocatable :: nearestn

    !> Index of nearest center for old configuration.
    integer, dimension(:), allocatable :: nearesto

    !> New electron probability.
    real(dp) :: pen

    !> Old electron probabilities for all states.
    real(dp), dimension(:), allocatable :: peo

    !> New wavefunction squared values.
    real(dp), dimension(:), allocatable :: psi2n

    !> Old wavefunction squared values for all states.
    real(dp), dimension(:, :), allocatable :: psi2o

    !> Old determinant part of wavefunction.
    real(dp), dimension(:), allocatable :: psido

    !> Old Jastrow part of wavefunction.
    real(dp), dimension(:), allocatable :: psijo

    !> Minimum distance to centers for new configuration.
    real(dp), dimension(:), allocatable :: rminn

    !> Minimum distance to centers (new-old mix).
    real(dp), dimension(:), allocatable :: rminno

    !> Minimum distance to centers for old configuration.
    real(dp), dimension(:), allocatable :: rmino

    !> Minimum distance to centers (old-new mix).
    real(dp), dimension(:), allocatable :: rminon

    !> Vector to nearest center for new configuration.
    real(dp), dimension(:, :), allocatable :: rvminn

    !> Vector to nearest center (new-old mix).
    real(dp), dimension(:, :), allocatable :: rvminno

    !> Vector to nearest center for old configuration.
    real(dp), dimension(:, :), allocatable :: rvmino

    !> Vector to nearest center (old-new mix).
    real(dp), dimension(:, :), allocatable :: rvminon

    !> New transition Jastrow factor.
    real(dp) :: tjfn

    !> Old transition Jastrow factors for all states.
    real(dp), dimension(:), allocatable :: tjfo

    !> Old-old transition Jastrow factor.
    real(dp) :: tjfoo

    !> New electron velocities.
    real(dp), dimension(:, :), allocatable :: vnew

    !> Old electron velocities.
    real(dp), dimension(:, :), allocatable :: vold

    !> New electron positions.
    real(dp), dimension(:, :), allocatable :: xnew

    !> Old electron positions.
    real(dp), dimension(:, :), allocatable :: xold

    !> Laplacian of wavefunction for DMC walkers.
    real(dp), dimension(:,:), allocatable :: d2o

    !> Electron probabilities for DMC walkers.
    real(dp), dimension(:,:), allocatable :: peo_dmc

    !> Determinant wavefunction for DMC walkers.
    real(dp), dimension(:,:), allocatable :: psido_dmc

    !> Jastrow wavefunction for DMC walkers.
    real(dp), dimension(:,:), allocatable :: psijo_dmc

    !> Old electron velocities for DMC walkers.
    real(dp), dimension(:,:,:,:), allocatable :: vold_dmc

    !> Old electron positions for DMC walkers.
    real(dp), dimension(:,:,:,:), allocatable :: xold_dmc

    private
    public   :: delttn, enew, eold, nearestn, nearesto, pen, peo, psi2n
    public   :: psi2o, psido, psijo, rminn, rminno, rmino, rminon, rvminn
    public   :: rvminno, rvmino, rvminon, tjfn, tjfo, tjfoo, vnew, vold, xnew, xold
    public   :: allocate_config, deallocate_config
    public   :: d2o, peo_dmc, psido_dmc, psijo_dmc, vold_dmc, xold_dmc
    public   :: allocate_config_dmc, deallocate_config_dmc
    save
contains

    !> Allocates memory for VMC configuration arrays.
    subroutine allocate_config()
      use mstates_mod, only: MSTATES
      use system, only: nelec
      use multiple_geo, only: MFORCE
      use vmc_mod, only: nwftypejas
        implicit none
        if (.not. allocated(delttn)) allocate (delttn(nelec))
        if (.not. allocated(enew)) allocate (enew(MFORCE))
        if (.not. allocated(eold)) allocate (eold(MSTATES, MFORCE))
        if (.not. allocated(nearestn)) allocate (nearestn(nelec), source=0)
        if (.not. allocated(nearesto)) allocate (nearesto(nelec), source=0)
        if (.not. allocated(peo)) allocate (peo(MSTATES))
        peo = 0. ! Necessary since it is used in metrop_mov1_slat in vmc, but never set
        if (.not. allocated(psi2n)) allocate (psi2n(MFORCE))
        if (.not. allocated(psi2o)) allocate (psi2o(MSTATES, MFORCE))
        if (.not. allocated(psido)) allocate (psido(MSTATES))
        if (.not. allocated(psijo)) allocate (psijo(nwftypejas))
        if (.not. allocated(rminn)) allocate (rminn(nelec))
        if (.not. allocated(rminno)) allocate (rminno(nelec))
        if (.not. allocated(rmino)) allocate (rmino(nelec))
        if (.not. allocated(rminon)) allocate (rminon(nelec))
        if (.not. allocated(rvminn)) allocate (rvminn(3, nelec))
        if (.not. allocated(rvminno)) allocate (rvminno(3, nelec))
        if (.not. allocated(rvmino)) allocate (rvmino(3, nelec))
        if (.not. allocated(rvminon)) allocate (rvminon(3, nelec))
        if (.not. allocated(tjfo)) allocate (tjfo(MSTATES))
        if (.not. allocated(vnew)) allocate (vnew(3, nelec))
        if (.not. allocated(vold)) allocate (vold(3, nelec))
        if (.not. allocated(xnew)) allocate (xnew(3, nelec))
        if (.not. allocated(xold)) allocate (xold(3, nelec))
    end subroutine allocate_config

    !> Deallocates memory for VMC configuration arrays.
    subroutine deallocate_config()
        if (allocated(xold)) deallocate (xold)
        if (allocated(xnew)) deallocate (xnew)
        if (allocated(vold)) deallocate (vold)
        if (allocated(vnew)) deallocate (vnew)
        if (allocated(tjfo)) deallocate (tjfo)
        if (allocated(rvminon)) deallocate (rvminon)
        if (allocated(rvmino)) deallocate (rvmino)
        if (allocated(rvminno)) deallocate (rvminno)
        if (allocated(rvminn)) deallocate (rvminn)
        if (allocated(rminon)) deallocate (rminon)
        if (allocated(rmino)) deallocate (rmino)
        if (allocated(rminno)) deallocate (rminno)
        if (allocated(rminn)) deallocate (rminn)
        if (allocated(psido)) deallocate (psido)
        if (allocated(psi2o)) deallocate (psi2o)
        if (allocated(psi2n)) deallocate (psi2n)
        if (allocated(peo)) deallocate (peo)
        if (allocated(nearesto)) deallocate (nearesto)
        if (allocated(nearestn)) deallocate (nearestn)
        if (allocated(eold)) deallocate (eold)
        if (allocated(enew)) deallocate (enew)
        if (allocated(delttn)) deallocate (delttn)
    end subroutine deallocate_config

    !> Allocates memory for DMC configuration arrays.
    subroutine allocate_config_dmc()
      use dmc_mod, only: mwalk
      use multiple_geo, only: MFORCE
      use system, only: nelec

      implicit none

      if (.not. allocated(d2o)) allocate(d2o(mwalk,MFORCE))
      if (.not. allocated(peo_dmc)) allocate(peo_dmc(mwalk,MFORCE))
      if (.not. allocated(psido_dmc)) allocate(psido_dmc(mwalk,MFORCE))
      if (.not. allocated(psijo_dmc)) allocate(psijo_dmc(mwalk,MFORCE))
      if (.not. allocated(vold_dmc)) allocate(vold_dmc(3,nelec,mwalk,MFORCE))
      if (.not. allocated(xold_dmc)) allocate(xold_dmc(3,nelec,mwalk,MFORCE))
    end subroutine allocate_config_dmc

    !> Deallocates memory for DMC configuration arrays.
    subroutine deallocate_config_dmc()
      if (allocated(d2o)) deallocate(d2o)
      if (allocated(peo_dmc)) deallocate(peo_dmc)
      if (allocated(psido_dmc)) deallocate(psido_dmc)
      if (allocated(psijo_dmc)) deallocate(psijo_dmc)
      if (allocated(vold_dmc)) deallocate(vold_dmc)
      if (allocated(xold_dmc)) deallocate(xold_dmc)
    end subroutine deallocate_config_dmc
end module config

!> @brief Module for random number generator state.
!> @author CHAMP developers
!> @date 2022
!>
!> @details This module maintains the state variables for the random number generator
!> used in Monte Carlo sampling. It uses a linear congruential generator with seed values.
module random

    implicit none

    !> Random number generator seed array.
    integer  :: ll(4)

    !> Random number generator multiplier array.
    integer  :: mm(4)

    !> Switch to select random number generator (1=default, other values for alternatives).
    integer  :: switch_rng = 1

    data mm/502, 1521, 4071, 2107/
    data ll/0, 0, 0, 1/

    private
    public :: ll, mm, switch_rng
    save

end module random

!> @brief Module for Monte Carlo sampling statistics.
!> @author CHAMP developers
!> @date 2022
!>
!> @details This module tracks Monte Carlo and DMC sampling statistics including
!> acceptance rates, rejection rates, diffusion coefficients, and branching information.
module stats
      use precision_kinds, only: dp

    implicit none

    !> Maximum rejection ratio.
    real(dp) :: rejmax

    !> Acceptance ratio for DMC.
    real(dp) :: acc

    !> Accepted diffusion coefficient squared.
    real(dp) :: dfus2ac

    !> Unaccepted diffusion coefficient squared.
    real(dp) :: dfus2un

    !> Number of accepted moves.
    integer  :: nacc

    !> Number of branching events.
    integer  :: nbrnch

    !> Number of walker decreases.
    integer  :: nodecr

    !> Number of attempted moves.
    real(dp) :: trymove

    private
    public :: rejmax
    public :: acc, dfus2ac, dfus2un, nacc, nbrnch, nodecr, trymove
    save
end module stats

!> @brief Module for nodal distance tracking.
!> @author CHAMP developers
!> @date 2022
!>
!> @details This module tracks the cumulative distance to nodal surfaces
!> during Monte Carlo sampling, useful for analyzing nodal properties.
module tmpnode
      use precision_kinds, only: dp

    implicit none

    !> Cumulative sum of distances to nodal surfaces.
    real(dp) :: distance_node_sum

    private
    public :: distance_node_sum
    save
end module tmpnode

!> @brief Main module for Monte Carlo sampling.
!> @author CHAMP developers
!> @date 2020
!>
!> @details This module provides a unified interface for allocating and deallocating
!> all sampling-related data structures including configurations, random number state,
!> and statistics.
module m_sampling
contains

!> Allocates all sampling arrays.
subroutine allocate_m_sampling()
    use config, only: allocate_config

    call allocate_config()
end subroutine allocate_m_sampling

!> Deallocates all sampling arrays.
subroutine deallocate_m_sampling()
    use config, only: deallocate_config

    call deallocate_config()
end subroutine deallocate_m_sampling
end module 
