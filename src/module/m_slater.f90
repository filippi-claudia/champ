module slater
    !> Arguments: d2dx2, ddx, fp, fpp, slmi

      use precision_kinds, only: dp
      use vmc_mod, only: nmat_dim

    implicit none

    real(dp), dimension(:, :), allocatable :: d2dx2 !(MELEC, nwftypeorb)
    real(dp), dimension(:, :, :), allocatable :: ddx !(3,MELEC, nwftypeorb)
    real(dp), dimension(:, :, :, :), allocatable :: fp !(3,nmat_dim,2, nwftypeorb)
    real(dp), dimension(:, :, :), allocatable :: fpp !(nmat_dim,2, nwftypeorb)
    real(dp), dimension(:, :, :), allocatable :: slmi !(nmat_dim,2, nwftypeorb)
    !> DMC extra variables:
    real(dp), dimension(:,:), allocatable :: fpd !(3,nmat_dim)
    real(dp), dimension(:), allocatable :: fppd !(nmat_dim)
    real(dp), dimension(:), allocatable :: fppu !(nmat_dim)
    real(dp), dimension(:,:), allocatable :: fpu !(3,nmat_dim)
    real(dp), dimension(:, :, :), allocatable :: cdet !(MDET,MSTATES,MWF)
    integer, dimension(:, :), allocatable :: iwundet !(MDET,2)
    real(dp), dimension(:), allocatable :: cdet_equiv !(MDET)
    real(dp), dimension(:), allocatable :: dcdet_equiv !(MDET)
    real(dp), dimension(:, :, :), allocatable :: coef !(MBASIS,norb_tot,MWF)

    integer :: ndet
    integer :: norb
    integer :: kref

    save

contains
    subroutine allocate_slater()
      use system,  only: nelec
      use vmc_mod, only: nmat_dim, nwftypeorb
        !use force_mod, only: MWF
        !use mstates_mod, only: MSTATES

        if (.not. allocated(d2dx2)) allocate(d2dx2(nelec, nwftypeorb))
        if (.not. allocated(ddx)) allocate(ddx(3, nelec, nwftypeorb))
        if (.not. allocated(fp)) allocate(fp(3, nmat_dim, 2, nwftypeorb))
        if (.not. allocated(fpp)) allocate(fpp(nmat_dim, 2, nwftypeorb))
        if (.not. allocated(slmi)) allocate(slmi(nmat_dim, 2, nwftypeorb))
        if (.not. allocated(fpd))  allocate(fpd(3,nmat_dim))
        if (.not. allocated(fppd)) allocate(fppd(nmat_dim))
        if (.not. allocated(fppu)) allocate(fppu(nmat_dim))
        if (.not. allocated(fpu))  allocate(fpu(3,nmat_dim))
        !if (.not. allocated(cdet)) allocate (cdet(MDET, MSTATES, MWF))
        if (.not. allocated(iwundet)) allocate (iwundet(ndet, 2))
        if (.not. allocated(cdet_equiv)) allocate (cdet_equiv(ndet))
        if (.not. allocated(dcdet_equiv)) allocate (dcdet_equiv(ndet))
        !if (.not. allocated(coef)) allocate (coef(MBASIS, norb_tot, MWF))

    end subroutine allocate_slater

    subroutine deallocate_slater()
        if (allocated(slmi)) deallocate(slmi)
        if (allocated(fpp)) deallocate(fpp)
        if (allocated(fp)) deallocate(fp)
        if (allocated(ddx)) deallocate(ddx)
        if (allocated(d2dx2)) deallocate(d2dx2)
        if (allocated(fpd))  deallocate(fpd)
        if (allocated(fppd)) deallocate(fppd)
        if (allocated(fppu)) deallocate(fppu)
        if (allocated(fpu))  deallocate(fpu)
        if (allocated(cdet)) deallocate (cdet)
        if (allocated(iwundet)) deallocate(iwundet)
        if (allocated(dcdet_equiv)) deallocate (dcdet_equiv)
        if (allocated(cdet_equiv)) deallocate (cdet_equiv)
        if (allocated(coef)) deallocate (coef)

    end subroutine deallocate_slater

end module slater
