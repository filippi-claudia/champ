module forcepar
    !> Arguments: deltot, istrech, nforce
    use precision_kinds, only: dp

    real(dp), dimension(:), allocatable :: deltot
    integer :: istrech
    integer :: nforce

    private
    public   ::  deltot, istrech, nforce
    public :: allocate_forcepar, deallocate_forcepar
    save
contains
    subroutine allocate_forcepar()
        use precision_kinds, only: dp
        if (.not. allocated(deltot)) allocate (deltot(nforce))
    end subroutine allocate_forcepar

    subroutine deallocate_forcepar()
        if (allocated(deltot)) deallocate (deltot)
    end subroutine deallocate_forcepar

end module forcepar
