module dets
    implicit none

    integer :: nmap

    save
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
      use dets,    only: nmap
      use mstates_mod, only: MSTATES
      use multiple_geo, only: nwftype
      use slater,  only: ndet
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

