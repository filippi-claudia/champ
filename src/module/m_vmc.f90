!> @brief Module for Variational Monte Carlo (VMC) parameters and array dimensions.
!> @author CHAMP developers
!> @date 2025
!>
!> @details This module manages dimensions and indexing arrays for VMC calculations.
!> It handles:
!> - Orbital and Jastrow wavefunction type mappings
!> - Matrix dimensions for Slater determinants
!> - State-to-orbital and state-to-Jastrow index conversions
!> - Basis set and center type parameters
!>
!> @note Slater matrices are dimensioned assuming (nelec/2)^2 for equal numbers
!> of up and down spins. For spin-polarized calculations, dimensions must be
!> adjusted accordingly.
!>
!> @note For Jastrow4, neqsx=2*(nordj-1) is sufficient.
!> For Jastrow3, neqsx=2*nordj should be sufficient.
module vmc_mod
    use precision_kinds, only: dp

    implicit none

    !> Number of states per orbital wavefunction type.
    integer, dimension(:), allocatable :: nstoo

    !> Number of states per Jastrow wavefunction type.
    integer, dimension(:), allocatable :: nstoj

    !> Jastrow type to state mapping array.
    integer, dimension(:, :), allocatable :: jtos

    !> Orbital type to state mapping array.
    integer, dimension(:, :), allocatable :: otos

    !> State to Jastrow type mapping array.
    integer, dimension(:), allocatable :: stoj

    !> State to orbital type mapping array.
    integer, dimension(:), allocatable :: stoo

    !> State to backflow-Jastrow index mapping array.
    integer, dimension(:), allocatable :: stobjx

    !> Backflow-Jastrow index to orbital type mapping array.
    integer, dimension(:), allocatable :: bjxtoo

    !> Backflow-Jastrow index to Jastrow type mapping array.
    integer, dimension(:), allocatable :: bjxtoj

    !> norb_tot :: total number of orbitals
    integer :: norb_tot

    !> Maximum of 3 and the total number of center types.
    integer :: nctyp3x

    !> Dimension of Slater matrix for up-spin electrons (nup*nup).
    integer :: nmat_dim

    !> Dimension for electron pair indices (nelec*(nelec-1)/2).
    integer :: nmat_dim2

    !> Number of orbital wavefunction types.
    integer :: nwftypeorb

    !> Number of Jastrow wavefunction types.
    integer :: nwftypejas

    !> Maximum number of states per Jastrow type.
    integer :: nstojmax

    !> Maximum number of states per orbital type.
    integer :: nstoomax

    !> Number of backflow-Jastrow indices.
    integer :: nbjx

    !> Total number of states across all orbital types.
    integer :: nstoo_tot

    !> Total number of states across all Jastrow types.
    integer :: nstoj_tot

    !> Extra Jastrow parameter count.
    integer :: extraj

    !> Extra orbital parameter count.
    integer :: extrao

    !> Number of terms in expansion.
    integer :: mterms

    !> Three times the total number of centers (3*ncent_tot).
    integer :: ncent3

    !> Max number of coefficients parameter (fixed at 5).
    integer, parameter :: NCOEF = 5

    !> Maximum number of excited states (fixed at 10).
    integer, parameter :: MEXCIT = 10

    private
    public :: norb_tot, nctyp3x
    public :: nmat_dim, nmat_dim2

    public :: mterms, nwftypejas, nwftypeorb, nstojmax, nstoomax, nbjx, nstoj_tot, nstoo_tot
    public :: nstoo, nstoj, jtos, otos, stoj, stoo, stobjx, extrao, extraj, bjxtoo, bjxtoj

    public :: ncent3, NCOEF, MEXCIT
    public :: set_vmc_size

    save
contains

    !> Sets the dimensions for VMC arrays based on system properties.
    subroutine set_vmc_size
      use system,  only: ncent_tot,nctype_tot,ndn,nelec,nup

        nmat_dim = nup*nup
        nmat_dim2 = nelec*(nelec - 1)/2
        nctyp3x = max(3, nctype_tot)
        ncent3 = 3*ncent_tot

    end subroutine set_vmc_size
end module vmc_mod
