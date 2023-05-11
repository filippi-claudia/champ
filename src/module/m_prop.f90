module properties
    !> Arguments: MAXPROP
    integer, parameter :: MAXPROP = 6
    private
    public :: MAXPROP
    save
end module properties

module prp000
    !> Arguments: ipropprt, iprop, nprop

    integer :: iprop
    integer :: ipropprt
    integer :: nprop

    private
    public :: ipropprt, iprop, nprop
    save
end module prp000

module prp001
    !> Arguments: vprop
      use precision_kinds, only: dp
      use properties, only: MAXPROP

    real(dp), dimension(:), allocatable :: vprop !(MAXPROP)

    private
    public :: vprop
    public :: allocate_prp001, deallocate_prp001
    save
contains
    subroutine allocate_prp001()
      use properties, only: MAXPROP
        if (.not. allocated(vprop)) allocate (vprop(MAXPROP))
    end subroutine allocate_prp001

    subroutine deallocate_prp001()
        if (allocated(vprop)) deallocate (vprop)
    end subroutine deallocate_prp001

end module prp001

module prp002
    !> Arguments: vprop_old
      use dmc_mod, only: mwalk
      use precision_kinds, only: dp
      use properties, only: MAXPROP

    real(dp), dimension(:, :), allocatable :: vprop_old !(MAXPROP,mwalk)
    real(dp), dimension(:), allocatable :: vprop_old2 !(MAXPROP)

    private
    public :: vprop_old, vprop_old2
    public :: allocate_prp002, deallocate_prp002
    save
contains
    subroutine allocate_prp002()
      use properties, only: MAXPROP
        if (.not. allocated(vprop_old)) allocate (vprop_old(MAXPROP, mwalk))
        if (.not. allocated(vprop_old2)) allocate (vprop_old2(MAXPROP))
    end subroutine allocate_prp002

    subroutine deallocate_prp002()
        if (allocated(vprop_old2)) deallocate (vprop_old2)
        if (allocated(vprop_old)) deallocate (vprop_old)
    end subroutine deallocate_prp002

end module prp002

module prp003
    !> Arguments: vprop_cm2, vprop_cum, cc_nuc, vprop_sum
      use precision_kinds, only: dp
      use properties, only: MAXPROP

    real(dp), dimension(:), allocatable :: cc_nuc !(3)
    real(dp), dimension(:), allocatable :: vprop_cm2 !(MAXPROP)
    real(dp), dimension(:), allocatable :: vprop_cum !(MAXPROP)
    real(dp), dimension(:), allocatable :: vprop_sum !(MAXPROP)

    private
    public :: vprop_cm2, vprop_cum, cc_nuc, vprop_sum
    public :: allocate_prp003, deallocate_prp003
    save
contains
    subroutine allocate_prp003()
      use properties, only: MAXPROP
        if (.not. allocated(cc_nuc)) allocate (cc_nuc(3))
        if (.not. allocated(vprop_cm2)) allocate (vprop_cm2(MAXPROP))
        if (.not. allocated(vprop_cum)) allocate (vprop_cum(MAXPROP))
        if (.not. allocated(vprop_sum)) allocate (vprop_sum(MAXPROP))
    end subroutine allocate_prp003

    subroutine deallocate_prp003()
        if (allocated(vprop_sum)) deallocate (vprop_sum)
        if (allocated(vprop_cum)) deallocate (vprop_cum)
        if (allocated(vprop_cm2)) deallocate (vprop_cm2)
        if (allocated(cc_nuc)) deallocate (cc_nuc)
    end subroutine deallocate_prp003

end module prp003

module m_prop
contains
subroutine allocate_m_prop()
      use prp001, only: allocate_prp001
      use prp002, only: allocate_prp002
      use prp003, only: allocate_prp003

    call allocate_prp001()
    call allocate_prp002()
    call allocate_prp003()
end subroutine allocate_m_prop

subroutine deallocate_m_prop()
      use prp001, only: deallocate_prp001
      use prp002, only: deallocate_prp002
      use prp003, only: deallocate_prp003

    call deallocate_prp001()
    call deallocate_prp002()
    call deallocate_prp003()
end subroutine deallocate_m_prop
end module
