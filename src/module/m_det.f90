module dets
    use precision_kinds, only: dp

    implicit none

    real(dp), dimension(:, :, :), allocatable :: cdet !(MDET,MSTATES,MWF)
    integer :: ndet
    integer :: nmap

    private
    public   :: cdet, ndet, nmap
!    public :: allocate_dets
    public :: deallocate_dets
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

module csfs
    !> Arguments: ccsf, cxdet, iadet, ibdet, icxdet, ncsf, nstates
    use precision_kinds, only: dp

    real(dp), dimension(:, :, :), allocatable :: ccsf !(MDET,MSTATES,MWF)
    real(dp), dimension(:), allocatable :: cxdet !(nmap)
    integer, dimension(:), allocatable :: iadet !(MDET)
    integer, dimension(:), allocatable :: ibdet !(MDET)
    integer, dimension(:), allocatable :: icxdet !(nmap)
    integer :: ncsf
    integer :: nstates

    private
    public   :: ccsf, cxdet, iadet, ibdet, icxdet, ncsf, nstates
    public :: allocate_csfs, deallocate_csfs
    save
contains
    subroutine allocate_csfs()
        use wfsec, only: nwftype
        use mstates_mod, only: MSTATES
        use dets, only: ndet, nmap
        if (.not. allocated(ccsf)) allocate (ccsf(ndet, MSTATES, nwftype))
        if (.not. allocated(cxdet)) allocate (cxdet(nmap))
        if (.not. allocated(iadet)) allocate (iadet(ndet))
        if (.not. allocated(ibdet)) allocate (ibdet(ndet))
        if (.not. allocated(icxdet)) allocate (icxdet(nmap))
    end subroutine allocate_csfs

    subroutine deallocate_csfs()
        if (allocated(icxdet)) deallocate (icxdet)
        if (allocated(ibdet)) deallocate (ibdet)
        if (allocated(iadet)) deallocate (iadet)
        if (allocated(cxdet)) deallocate (cxdet)
        if (allocated(ccsf)) deallocate (ccsf)
    end subroutine deallocate_csfs

end module csfs

module dets_equiv
    !> Arguments: cdet_equiv, dcdet_equiv
    use precision_kinds, only: dp

    implicit none

    real(dp), dimension(:), allocatable :: cdet_equiv !(MDET)
    real(dp), dimension(:), allocatable :: dcdet_equiv !(MDET)

    private
    public   ::  cdet_equiv, dcdet_equiv
    public :: allocate_dets_equiv, deallocate_dets_equiv
    save
contains
    subroutine allocate_dets_equiv()
        use dets, only: ndet

        if (.not. allocated(cdet_equiv)) allocate (cdet_equiv(ndet))
        if (.not. allocated(dcdet_equiv)) allocate (dcdet_equiv(ndet))

    end subroutine allocate_dets_equiv

    subroutine deallocate_dets_equiv()
        if (allocated(dcdet_equiv)) deallocate (dcdet_equiv)
        if (allocated(cdet_equiv)) deallocate (cdet_equiv)
    end subroutine deallocate_dets_equiv

end module dets_equiv
