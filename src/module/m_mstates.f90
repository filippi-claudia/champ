module mstates_mod
    !> Arguments: MSTATES
    integer, parameter :: MSTATES = 3

    private
    public :: MSTATES
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
        use mstates_mod, only: MSTATES
        if (.not. allocated(effcm2)) allocate (effcm2(MSTATES), source=0.0_dp)
        if (.not. allocated(effcum)) allocate (effcum(MSTATES), source=0.0_dp)
    end subroutine allocate_mstates2

    subroutine deallocate_mstates2()
        if (allocated(effcum)) deallocate (effcum)
        if (allocated(effcm2)) deallocate (effcm2)
    end subroutine deallocate_mstates2

end module mstates2

module mstates3
    !> Arguments: weights_g, iweight_g
    use precision_kinds, only: dp

    integer, dimension(:), allocatable :: iweight_g !(MSTATES)
    real(dp), dimension(:), allocatable :: weights_g !(MSTATES)
    private

    public :: weights_g, iweight_g
    public :: allocate_mstates3, deallocate_mstates3
    save
contains
    subroutine allocate_mstates3()
        use mstates_mod, only: MSTATES
        if (.not. allocated(iweight_g)) allocate (iweight_g(MSTATES), source=0)
        if (.not. allocated(weights_g)) allocate (weights_g(MSTATES), source=0.0_dp)
    end subroutine allocate_mstates3

    subroutine deallocate_mstates3()
        if (allocated(weights_g)) deallocate (weights_g)
        if (allocated(iweight_g)) deallocate (iweight_g)
    end subroutine deallocate_mstates3

end module mstates3

module m_mstates
contains
subroutine allocate_m_mstates()
    use mstates2, only: allocate_mstates2
    use mstates3, only: allocate_mstates3

    call allocate_mstates2()
    call allocate_mstates3()
end subroutine allocate_m_mstates

subroutine deallocate_m_mstates()
    use mstates2, only: deallocate_mstates2
    use mstates3, only: deallocate_mstates3

    call deallocate_mstates2()
    call deallocate_mstates3()
end subroutine deallocate_m_mstates
end module 
