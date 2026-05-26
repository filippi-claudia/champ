!> @brief Module for radial probability distributions of spin-up and spin-down electrons.
!> @author CHAMP developers
!>
!> @details This module stores per-radial-bin probability distributions
!> for spin-up and spin-down electrons used in DMC HDF5 restart files.
module denupdn
    use precision_kinds, only: dp

    implicit none

    !> Radial probability for spin-up electrons, dimension (nrad)
    real(dp), dimension(:), allocatable :: rprobup

    !> Radial probability for spin-down electrons, dimension (nrad)
    real(dp), dimension(:), allocatable :: rprobdn

    private
    public :: rprobup, rprobdn
    public :: allocate_denupdn, deallocate_denupdn
    save

contains

    !> Allocates spin-resolved radial probability arrays using nrad from vmc_mod.
    subroutine allocate_denupdn()
        use vmc_mod, only: nrad
        if (.not. allocated(rprobup)) allocate(rprobup(nrad))
        if (.not. allocated(rprobdn)) allocate(rprobdn(nrad))
        rprobup = 0.0_dp
        rprobdn = 0.0_dp
    end subroutine allocate_denupdn

    !> Deallocates spin-resolved radial probability arrays.
    subroutine deallocate_denupdn()
        if (allocated(rprobdn)) deallocate(rprobdn)
        if (allocated(rprobup)) deallocate(rprobup)
    end subroutine deallocate_denupdn

end module denupdn
