!> @brief Module for radial step-size acceptance statistics in VMC/DMC.
!> @author CHAMP developers
!>
!> @details This module stores per-radial-bin statistics for step-size
!> acceptance/rejection and kinetic energy accumulation used in
!> HDF5 restart files.
module step
    use precision_kinds, only: dp

    implicit none

    !> Radial probability distribution, dimension (nrad)
    real(dp), dimension(:), allocatable :: rprob

    !> Cumulative kinetic energy per radial bin, dimension (nrad)
    real(dp), dimension(:), allocatable :: ekin

    !> Cumulative kinetic energy squared per radial bin, dimension (nrad)
    real(dp), dimension(:), allocatable :: ekin2

    !> Number of accepted moves per radial bin, dimension (nrad)
    real(dp), dimension(:), allocatable :: suc

    !> Truncated forward-bias per radial bin, dimension (nrad)
    real(dp), dimension(:), allocatable :: trunfb

    !> Number of attempted moves per radial bin, dimension (nrad)
    real(dp), dimension(:), allocatable :: try

    private
    public :: rprob, ekin, ekin2, suc, trunfb, try
    public :: allocate_step, deallocate_step
    save

contains

    !> Allocates radial step-size arrays using nrad from vmc_mod.
    subroutine allocate_step()
        use vmc_mod, only: nrad
        if (.not. allocated(rprob))  allocate(rprob(nrad))
        if (.not. allocated(ekin))   allocate(ekin(nrad))
        if (.not. allocated(ekin2))  allocate(ekin2(nrad))
        if (.not. allocated(suc))    allocate(suc(nrad))
        if (.not. allocated(trunfb)) allocate(trunfb(nrad))
        if (.not. allocated(try))    allocate(try(nrad))
        rprob  = 0.0_dp
        ekin   = 0.0_dp
        ekin2  = 0.0_dp
        suc    = 0.0_dp
        trunfb = 0.0_dp
        try    = 0.0_dp
    end subroutine allocate_step

    !> Deallocates radial step-size arrays.
    subroutine deallocate_step()
        if (allocated(try))    deallocate(try)
        if (allocated(trunfb)) deallocate(trunfb)
        if (allocated(suc))    deallocate(suc)
        if (allocated(ekin2))  deallocate(ekin2)
        if (allocated(ekin))   deallocate(ekin)
        if (allocated(rprob))  deallocate(rprob)
    end subroutine deallocate_step

end module step
