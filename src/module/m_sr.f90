!> @brief Module for Stochastic Reconfiguration (SR) method parameters.
!> @author CHAMP developers
!> @date 2025
!>
!> @details This module defines dimension parameters for the SR optimization method,
!> which is used for wavefunction optimization in quantum Monte Carlo calculations.
module sr_mod

    implicit none

    !> Maximum number of parameters to optimize.
    integer :: mparm

    !> mobs
    integer :: mobs

    !> Maximum number of configurations.
    integer :: mconf

    !> Flag for SR rescaling method.
    integer :: i_sr_rescale

    !> Zero-variance zero-bias parameter flag?.
    integer :: izvzb

    private
    public :: mparm, mobs, mconf
    public :: izvzb, i_sr_rescale
    save
end module sr_mod

!> @brief Module for SR indexing counters.
!> @author CHAMP developers
!> @date 2025
!>
!> @details This module contains index counters for SR calculations.
module sr_index

    implicit none

    !> jelo
    integer :: jelo

    !> jelo2
    integer :: jelo2

    !> jelohfj
    integer :: jelohfj

    private
    public :: jelo, jelo2, jelohfj
    save
end module sr_index

!> @brief Module for Stochastic Reconfiguration matrices and arrays.
!> @author CHAMP developers
!> @date 2025
!>
!> @details This module contains the main data structures for SR optimization including
!> local energies, SR matrices, overlap matrices, and observables for multiple states.
!> It handles the numerical computation of parameter updates in wavefunction optimization.
module sr_mat_n
      use mstates_mod, only: MSTATES
      use precision_kinds, only: dp
      use sr_mod, only: mparm, mobs, mconf

    implicit none

    !> Local energies for all configurations and states.
    real(dp), dimension(:, :), allocatable :: elocal

    !> SR Hamiltonian matrix.
    real(dp), dimension(:, :), allocatable :: h_sr

    !> SR Hamiltonian penalty matrix for constraints.
    real(dp), dimension(:, :), allocatable :: h_sr_penalty

    !> Lambda parameters for state orthogonalization (flattened).
    real(dp), dimension(:), allocatable :: isr_lambda

    !> Lambda matrix for state coupling.
    real(dp), dimension(:, :), allocatable :: sr_lambda

    !> Current state index for SR optimization.
    integer :: sr_state

    !> Orthogonality constraint flag.
    integer :: ortho=0

    !> jefj
    integer :: jefj

    !> jfj
    integer :: jfj

    !> jhfj
    integer :: jhfj

    !> nconf_n
    integer :: nconf_n

    !> Diagonal elements of overlap matrix S.
    real(dp), dimension(:, :), allocatable :: s_diag

    !> Inverse diagonal elements of overlap matrix S.
    real(dp), dimension(:, :), allocatable :: s_ii_inv

    !> SR overlap matrix times Hamiltonian.
    real(dp), dimension(:, :), allocatable :: sr_ho

    !> SR overlap matrix for parameters and configurations.
    real(dp), dimension(:, :, :), allocatable :: sr_o

    !> Configuration weights for all states.
    real(dp), dimension(:, :), allocatable :: wtg

    !> obs_tot
    real(dp), dimension(:, :), allocatable :: obs_tot

    private
    public :: elocal, h_sr, jefj, jfj, jhfj, nconf_n, s_diag, s_ii_inv, sr_ho, sr_o, wtg, obs_tot, isr_lambda, sr_lambda, sr_state, h_sr_penalty, ortho
    ! public :: obs
    public :: allocate_sr_mat_n, deallocate_sr_mat_n
    save
contains

    !> Allocates memory for SR matrices and arrays.
    subroutine allocate_sr_mat_n()
      use mstates_mod, only: MSTATES
      use optwf_func, only: ifunc_omega
      use sr_mod, only: mparm, mobs, mconf, izvzb, i_sr_rescale
        if (.not. allocated(elocal)) allocate (elocal(mconf, MSTATES))
        if (.not. allocated(h_sr)) allocate (h_sr(mparm, MSTATES))
        if (.not. allocated(h_sr_penalty)) allocate (h_sr_penalty(mparm,MSTATES))
        if (.not. allocated(isr_lambda)) allocate (isr_lambda(MSTATES*(MSTATES-1)/2))
        if (.not. allocated(sr_lambda)) allocate (sr_lambda(MSTATES,MSTATES))
        if (.not. allocated(s_diag)) allocate (s_diag(mparm, MSTATES))
        if (.not. allocated(s_ii_inv)) allocate (s_ii_inv(mparm, MSTATES))
        if (.not. allocated(sr_ho)) allocate (sr_ho(mparm, mconf))
        if (.not. allocated(sr_o)) allocate (sr_o(mparm, mconf, MSTATES))
        if (.not. allocated(wtg)) allocate (wtg(mconf, MSTATES))
        if (.not. allocated(obs_tot)) allocate (obs_tot(mobs, MSTATES))
    end subroutine allocate_sr_mat_n

    !> Deallocates memory for SR matrices and arrays.
    subroutine deallocate_sr_mat_n()
        if (allocated(obs_tot)) deallocate (obs_tot)
        if (allocated(wtg)) deallocate (wtg)
        if (allocated(sr_o)) deallocate (sr_o)
        if (allocated(sr_ho)) deallocate (sr_ho)
        if (allocated(s_ii_inv)) deallocate (s_ii_inv)
        if (allocated(s_diag)) deallocate (s_diag)
        if (allocated(h_sr)) deallocate (h_sr)
        if (allocated(h_sr_penalty)) deallocate (h_sr_penalty)
        if (allocated(isr_lambda)) deallocate (isr_lambda)
        if (allocated(sr_lambda)) deallocate (sr_lambda)
        if (allocated(elocal)) deallocate (elocal)
    end subroutine deallocate_sr_mat_n

end module sr_mat_n

!> @brief Main module for Stochastic Reconfiguration method.
!> @author CHAMP developers
!> @date 2025
!>
!> @details This module provides a unified interface for allocating and deallocating
!> all SR-related data structures. The Stochastic Reconfiguration method is an efficient
!> optimization algorithm for variational parameters in quantum Monte Carlo calculations.
module m_sr
contains

!> Allocates all SR method arrays.
subroutine allocate_m_sr()
      use sr_mat_n, only: allocate_sr_mat_n

    call allocate_sr_mat_n()
end subroutine allocate_m_sr

!> Deallocates all SR method arrays.
subroutine deallocate_m_sr()
      use sr_mat_n, only: deallocate_sr_mat_n

    call deallocate_sr_mat_n()
end subroutine deallocate_m_sr
end module 
