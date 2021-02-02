module mstates_mod
    !> Arguments: MSTATES, MDETCSFX
    integer, parameter :: MSTATES = 2
    integer, parameter :: MDETCSFX = 20

    private
    public :: MSTATES, MDETCSFX
    save
end module mstates_mod

module mstates_ctrl
    !> Arguments: iefficiency, nstates_psig, iguiding

    integer :: iefficiency
    integer :: iguiding
    integer :: nstates_psig

    private
    public :: iefficiency, nstates_psig, iguiding
    save
end module mstates_ctrl

module wfsec
    !> Arguments: iwf, iwftype, nwftype

    integer :: iwf
    integer, dimension(:), allocatable :: iwftype !(MFORCE)
    integer :: nwftype

    private
    public :: iwf, iwftype, nwftype
    public :: allocate_wfsec, deallocate_wfsec
    save
contains
    ! subroutine allocate_wfsec()
    !     use force_mod, only: MFORCE
    !     if (.not. allocated(iwftype)) allocate (iwftype(MFORCE))
    ! end subroutine allocate_wfsec

    subroutine deallocate_wfsec()
        if (allocated(iwftype)) deallocate (iwftype)
    end subroutine deallocate_wfsec

end module wfsec

module csfs
    !> Arguments: ccsf, cxdet, iadet, ibdet, icxdet, ncsf, nstates
    use precision_kinds, only: dp

    real(dp), dimension(:, :, :), allocatable :: ccsf !(MDET,MSTATES,MWF)
    real(dp), dimension(:), allocatable :: cxdet !(MDET*MDETCSFX)
    integer, dimension(:), allocatable :: iadet !(MDET)
    integer, dimension(:), allocatable :: ibdet !(MDET)
    integer, dimension(:), allocatable :: icxdet !(MDET*MDETCSFX)
    integer :: ncsf
    integer :: nstates

    private
    public   :: ccsf, cxdet, iadet, ibdet, icxdet, ncsf, nstates
    public :: allocate_csfs, deallocate_csfs
    save
contains
    subroutine allocate_csfs()
        use wfsec, only: nwftype
        use dets, only: ndet
        use precision_kinds, only: dp
        use mstates_mod, only: MDETCSFX

        if (.not. allocated(ccsf)) allocate (ccsf(ndet, nstates, nwftype))
        if (.not. allocated(cxdet)) allocate (cxdet(ndet*MDETCSFX))
        if (.not. allocated(iadet)) allocate (iadet(ndet))
        if (.not. allocated(ibdet)) allocate (ibdet(ndet))
        if (.not. allocated(icxdet)) allocate (icxdet(ndet*MDETCSFX))
    end subroutine allocate_csfs

    subroutine deallocate_csfs()
        if (allocated(icxdet)) deallocate (icxdet)
        if (allocated(ibdet)) deallocate (ibdet)
        if (allocated(iadet)) deallocate (iadet)
        if (allocated(cxdet)) deallocate (cxdet)
        if (allocated(ccsf)) deallocate (ccsf)
    end subroutine deallocate_csfs

end module csfs

module mstates2
    !> Arguments: effcum, effcm2
    use precision_kinds, only: dp
    use mstates_mod, only: MSTATES

    real(dp), dimension(:), allocatable :: effcm2 !(MSTATES)
    real(dp), dimension(:), allocatable :: effcum !(MSTATES)
    private

    public :: effcum, effcm2
    public :: allocate_mstates2, deallocate_mstates2
    save
contains
    subroutine allocate_mstates2()
        use csfs, only: nstates
        use precision_kinds, only: dp
        use mstates_mod, only: MSTATES
        if (.not. allocated(effcm2)) allocate (effcm2(nstates))
        if (.not. allocated(effcum)) allocate (effcum(nstates))
    end subroutine allocate_mstates2

    subroutine deallocate_mstates2()
        if (allocated(effcum)) deallocate (effcum)
        if (allocated(effcm2)) deallocate (effcm2)
    end subroutine deallocate_mstates2

end module mstates2

module mstates3
    !> Arguments: weights_g, iweight_g
    use precision_kinds, only: dp
    use mstates_mod, only: MSTATES

    integer, dimension(:), allocatable :: iweight_g !(MSTATES)
    real(dp), dimension(:), allocatable :: weights_g !(MSTATES)
    private

    public :: weights_g, iweight_g
    public :: allocate_mstates3, deallocate_mstates3
    save
contains
    subroutine allocate_mstates3()
        use csfs, only: nstates
        use precision_kinds, only: dp
        use mstates_mod, only: MSTATES
        if (.not. allocated(iweight_g)) allocate (iweight_g(nstates))
        if (.not. allocated(weights_g)) allocate (weights_g(nstates))
    end subroutine allocate_mstates3

    subroutine deallocate_mstates3()
        if (allocated(weights_g)) deallocate (weights_g)
        if (allocated(iweight_g)) deallocate (iweight_g)
    end subroutine deallocate_mstates3

end module mstates3

subroutine allocate_m_mstates()
    use mstates2, only: allocate_mstates2
    use mstates3, only: allocate_mstates3

    call allocate_mstates2()
    call allocate_mstates3()
end subroutine allocate_m_mstates
