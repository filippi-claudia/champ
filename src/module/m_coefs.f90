!> @brief Module for wavefunction expansion and basis set dimension parameters.
!> @author CHAMP developers
!> @date 2020
!>
!> @details This module stores fundamental dimension parameters for the wavefunction
!> expansion in terms of atomic orbital (AO) basis functions. It defines the size
!> of the basis set and the maximum number of additional orbitals that can be used
!> for orbital optimization or excited state calculations.
!>
!> Key parameters:
!> - nbasis: Total number of atomic orbital basis functions
!> - next_max: Maximum number of additional orbitals for optimization
!>
!> @note These parameters are set during input file parsing and remain constant
!> throughout the calculation. They determine array dimensions in many modules.
module coefs
    implicit none

    !> Total number of atomic orbital (AO) basis functions in the expansion
    integer :: nbasis

    !> Maximum number of additional orbitals for orbital optimization or excited states
    integer :: next_max

    save
end module coefs
