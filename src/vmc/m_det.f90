module dets
    use precision_kinds, only: dp
    real(dp), dimension(:, :, :), allocatable :: cdet !(MDET,MSTATES,MWF)
    integer :: ndet

    private
    public   :: cdet, ndet
    public :: allocate_dets, deallocate_dets
    save
contains
    ! subroutine allocate_dets()
    !     use force_mod, only: MWF
    !     use precision_kinds, only: dp
    !     use vmc_mod, only: MDET
    !     use mstates_mod, only: MSTATES
    !     if (.not. allocated(cdet)) allocate (cdet(MDET, MSTATES, MWF))
    ! end subroutine allocate_dets

    subroutine deallocate_dets()
        if (allocated(cdet)) deallocate (cdet)
    end subroutine deallocate_dets

end module dets

module dets_equiv
    !> Arguments: cdet_equiv, dcdet_equiv
    use precision_kinds, only: dp

    real(dp), dimension(:), allocatable :: cdet_equiv !(MDET)
    real(dp), dimension(:), allocatable :: dcdet_equiv !(MDET)

    private
    public   ::  cdet_equiv, dcdet_equiv
    public :: allocate_dets_equiv, deallocate_dets_equiv
    save
contains
    subroutine allocate_dets_equiv()
        use dets, only: ndet
        use vmc_mod, only: MDET
        use precision_kinds, only: dp

        if (.not. allocated(cdet_equiv)) allocate (cdet_equiv(MDET))
        if (.not. allocated(dcdet_equiv)) allocate (dcdet_equiv(MDET))

    end subroutine allocate_dets_equiv

    subroutine deallocate_dets_equiv()
        if (allocated(dcdet_equiv)) deallocate (dcdet_equiv)
        if (allocated(cdet_equiv)) deallocate (cdet_equiv)
    end subroutine deallocate_dets_equiv

end module dets_equiv
