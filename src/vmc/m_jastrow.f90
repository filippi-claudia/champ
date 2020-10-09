module jasn
    !> Arguments: d2ijn, d2n, fijn, fjn, fsn, fsumn
    use precision_kinds, only: dp
    use vmc_mod, only: MELEC

    real(dp), dimension(:, :), allocatable :: d2ijn !(MELEC,MELEC)
    real(dp) :: d2n
    real(dp), dimension(:, :, :), allocatable :: fijn !(3,MELEC,MELEC)
    real(dp), dimension(:, :), allocatable :: fjn !(3,MELEC)
    real(dp), dimension(:, :), allocatable :: fsn !(MELEC,MELEC)
    real(dp) :: fsumn

    private
    public :: d2ijn, d2n, fijn, fjn, fsn, fsumn
    public :: allocate_jasn, deallocate_jasn
    save
contains
    subroutine allocate_jasn()
        use precision_kinds, only: dp
        use vmc_mod, only: MELEC
        if (.not. allocated(d2ijn)) allocate (d2ijn(MELEC, MELEC))
        if (.not. allocated(fijn)) allocate (fijn(3, MELEC, MELEC))
        if (.not. allocated(fjn)) allocate (fjn(3, MELEC))
        if (.not. allocated(fsn)) allocate (fsn(MELEC, MELEC))
    end subroutine allocate_jasn

    subroutine deallocate_jasn()
        if (allocated(fsn)) deallocate (fsn)
        if (allocated(fjn)) deallocate (fjn)
        if (allocated(fijn)) deallocate (fijn)
        if (allocated(d2ijn)) deallocate (d2ijn)
    end subroutine deallocate_jasn

end module jasn

module jaso
    !> Arguments: d2ijo, d2o, fijo, fjo, fso, fsumo
    use precision_kinds, only: dp
    use vmc_mod, only: MELEC

    real(dp), dimension(:, :), allocatable :: d2ijo !(MELEC,MELEC)
    real(dp) :: d2o
    real(dp), dimension(:, :, :), allocatable :: fijo !(3,MELEC,MELEC)
    real(dp), dimension(:, :), allocatable :: fjo !(3,MELEC)
    real(dp), dimension(:, :), allocatable :: fso !(MELEC,MELEC)
    real(dp) :: fsumo

    private
    public :: d2ijo, d2o, fijo, fjo, fso, fsumo
    public :: allocate_jaso, deallocate_jaso
    save
contains
    subroutine allocate_jaso()
        use precision_kinds, only: dp
        use vmc_mod, only: MELEC
        if (.not. allocated(d2ijo)) allocate (d2ijo(MELEC, MELEC))
        if (.not. allocated(fijo)) allocate (fijo(3, MELEC, MELEC))
        if (.not. allocated(fjo)) allocate (fjo(3, MELEC))
        if (.not. allocated(fso)) allocate (fso(MELEC, MELEC))
    end subroutine allocate_jaso

    subroutine deallocate_jaso()
        if (allocated(fso)) deallocate (fso)
        if (allocated(fjo)) deallocate (fjo)
        if (allocated(fijo)) deallocate (fijo)
        if (allocated(d2ijo)) deallocate (d2ijo)
    end subroutine deallocate_jaso

end module jaso

module jaspar
    !> Arguments: nspin1, nspin2, sspin, sspinn, is
    use precision_kinds, only: dp

    integer :: is
    integer :: nspin1
    integer :: nspin2
    real(dp) :: sspin
    real(dp) :: sspinn

    private
    public   :: nspin1, nspin2, sspin, sspinn, is
    save
end module jaspar

module jaspar1
    !> Arguments: cjas1, cjas2
    use force_mod, only: MWF
    use precision_kinds, only: dp

    real(dp), dimension(:), allocatable :: cjas1 !(MWF)
    real(dp), dimension(:), allocatable :: cjas2 !(MWF)

    private
    public   ::  cjas1, cjas2
    public :: allocate_jaspar1, deallocate_jaspar1
    save
contains
    subroutine allocate_jaspar1()
        use force_mod, only: MWF
        use precision_kinds, only: dp
        if (.not. allocated(cjas1)) allocate (cjas1(MWF))
        if (.not. allocated(cjas2)) allocate (cjas2(MWF))
    end subroutine allocate_jaspar1

    subroutine deallocate_jaspar1()
        if (allocated(cjas2)) deallocate (cjas2)
        if (allocated(cjas1)) deallocate (cjas1)
    end subroutine deallocate_jaspar1

end module jaspar1

module jaspar2
    !> Arguments: a1, a2
    use force_mod, only: MWF
    use precision_kinds, only: dp

    real(dp), dimension(:, :, :), allocatable :: a1 !(83,3,MWF)
    real(dp), dimension(:, :, :), allocatable :: a2 !(83,3,MWF)

    private
    public :: a1, a2
    public :: allocate_jaspar2, deallocate_jaspar2
    save
contains
    subroutine allocate_jaspar2()
        use force_mod, only: MWF
        use precision_kinds, only: dp
        if (.not. allocated(a1)) allocate (a1(83, 3, MWF))
        if (.not. allocated(a2)) allocate (a2(83, 3, MWF))
    end subroutine allocate_jaspar2

    subroutine deallocate_jaspar2()
        if (allocated(a2)) deallocate (a2)
        if (allocated(a1)) deallocate (a1)
    end subroutine deallocate_jaspar2

end module jaspar2

module jaspar3
    !> Arguments: a, b, c, fck, nord, scalek
    use force_mod, only: MWF
    use precision_kinds, only: dp
    use vmc_mod, only: MCTYPE
    use vmc_mod, only: MORDJ1

    real(dp), dimension(:, :), allocatable :: a !(MORDJ1,MWF)
    real(dp), dimension(:, :, :), allocatable :: b !(MORDJ1,2,MWF)
    real(dp), dimension(:, :, :), allocatable :: c !(83,MCTYPE,MWF)
    real(dp), dimension(:, :, :), allocatable :: fck !(15,MCTYPE,MWF)
    integer :: nord
    real(dp), dimension(:), allocatable :: scalek !(MWF)

    private
    public :: a, b, c, fck, nord, scalek
    public :: allocate_jaspar3, deallocate_jaspar3
    save
contains
    subroutine allocate_jaspar3()
        use force_mod, only: MWF
        use precision_kinds, only: dp
        use vmc_mod, only: MCTYPE
        use vmc_mod, only: MORDJ1
        if (.not. allocated(a)) allocate (a(MORDJ1, MWF))
        ! if (.not. allocated(b)) allocate (b(MORDJ1, 2, MWF))
        ! if (.not. allocated(c)) allocate (c(83, MCTYPE, MWF))
        if (.not. allocated(fck)) allocate (fck(15, MCTYPE, MWF))
        ! if (.not. allocated(scalek)) allocate (scalek(MWF))
    end subroutine allocate_jaspar3

    subroutine deallocate_jaspar3()
        if (allocated(scalek)) deallocate (scalek)
        if (allocated(fck)) deallocate (fck)
        if (allocated(c)) deallocate (c)
        if (allocated(b)) deallocate (b)
        if (allocated(a)) deallocate (a)
    end subroutine deallocate_jaspar3

end module jaspar3

module jaspar4
    !> Arguments: a4, norda, nordb, nordc
    use force_mod, only: MWF
    use precision_kinds, only: dp
    use vmc_mod, only: MCTYPE
    use vmc_mod, only: MORDJ1

    real(dp), dimension(:, :, :), allocatable :: a4 !(MORDJ1,MCTYPE,MWF)
    integer :: norda
    integer :: nordb
    integer :: nordc

    private
    public :: a4, norda, nordb, nordc
    public :: allocate_jaspar4, deallocate_jaspar4
    save
contains
    ! subroutine allocate_jaspar4()
    !     use force_mod, only: MWF
    !     use precision_kinds, only: dp
    !     use vmc_mod, only: MCTYPE
    !     use vmc_mod, only: MORDJ1
    !     if (.not. allocated(a4)) allocate (a4(MORDJ1, MCTYPE, MWF))
    ! end subroutine allocate_jaspar4

    subroutine deallocate_jaspar4()
        if (allocated(a4)) deallocate (a4)
    end subroutine deallocate_jaspar4

end module jaspar4

module jaspar6
    !> Arguments: asymp_jasa, asymp_jasb, asymp_r, c1_jas6, c1_jas6i, c2_jas6, cutjas, cutjasi
    use precision_kinds, only: dp
    use vmc_mod, only: MCTYPE

    real(dp), dimension(:), allocatable :: asymp_jasa !(MCTYPE)
    real(dp), dimension(:), allocatable :: asymp_jasb !(2)
    real(dp) :: asymp_r
    real(dp) :: c1_jas6
    real(dp) :: c1_jas6i
    real(dp) :: c2_jas6
    real(dp) :: cutjas
    real(dp) :: cutjasi

    private
    public :: asymp_jasa, asymp_jasb, asymp_r, c1_jas6, c1_jas6i, c2_jas6, cutjas, cutjasi
    public :: allocate_jaspar6, deallocate_jaspar6
    save
contains
    subroutine allocate_jaspar6()
        use precision_kinds, only: dp
        use vmc_mod, only: MCTYPE
        if (.not. allocated(asymp_jasa)) allocate (asymp_jasa(MCTYPE))
        if (.not. allocated(asymp_jasb)) allocate (asymp_jasb(2))
    end subroutine allocate_jaspar6

    subroutine deallocate_jaspar6()
        if (allocated(asymp_jasb)) deallocate (asymp_jasb)
        if (allocated(asymp_jasa)) deallocate (asymp_jasa)
    end subroutine deallocate_jaspar6

end module jaspar6

module jaspointer
    !> Arguments: npoint, npointa
    use vmc_mod, only: MCTYP3X

    integer, dimension(:), allocatable :: npoint !(MCTYP3X)
    integer, dimension(:), allocatable :: npointa !(3*MCTYP3X)

    private
    public :: npoint, npointa
    public :: allocate_jaspointer, deallocate_jaspointer
    save
contains
    subroutine allocate_jaspointer()
        use vmc_mod, only: MCTYP3X
        if (.not. allocated(npoint)) allocate (npoint(MCTYP3X))
        if (.not. allocated(npointa)) allocate (npointa(3*MCTYP3X))
    end subroutine allocate_jaspointer

    subroutine deallocate_jaspointer()
        if (allocated(npointa)) deallocate (npointa)
        if (allocated(npoint)) deallocate (npoint)
    end subroutine deallocate_jaspointer

end module jaspointer

subroutine allocate_m_jastrow()
    use jasn, only: allocate_jasn
    use jaso, only: allocate_jaso
    use jaspar1, only: allocate_jaspar1
    use jaspar2, only: allocate_jaspar2
    use jaspar3, only: allocate_jaspar3
    use jaspar4, only: allocate_jaspar4
    use jaspar6, only: allocate_jaspar6
    use jaspointer, only: allocate_jaspointer

    call allocate_jasn()
    call allocate_jaso()
    call allocate_jaspar1()
    call allocate_jaspar2()
    call allocate_jaspar3()
    ! call allocate_jaspar4()
    call allocate_jaspar6()
    call allocate_jaspointer()
end subroutine allocate_m_jastrow
