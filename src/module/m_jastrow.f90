module jastrow_update
    !> Arguments: d2ijn, d2n, fijn, fjn, fsn, fsumn
    use precision_kinds, only: dp

    implicit none

    real(dp), dimension(:, :), allocatable :: d2ijn !(MELEC,MELEC)
    real(dp) :: d2n
    real(dp), dimension(:, :, :), allocatable :: fijn !(3,MELEC,MELEC)
    real(dp), dimension(:, :), allocatable :: fjn !(3,MELEC)
    real(dp), dimension(:, :), allocatable :: fsn !(MELEC,MELEC)
    real(dp) :: fsumn

    real(dp), dimension(:, :), allocatable :: d2ijo !(MELEC,MELEC)
    real(dp) :: d2o
    real(dp), dimension(:, :, :), allocatable :: fijo !(3,MELEC,MELEC)
    real(dp), dimension(:, :), allocatable :: fjo !(3,MELEC)
    real(dp), dimension(:, :), allocatable :: fso !(MELEC,MELEC)
    real(dp) :: fsumo
    !> DMC
    real(dp) :: d2jo

    private
    public :: d2ijo, d2o, fijo, fjo, fso, fsumo, d2jo
    public :: allocate_jaso, deallocate_jaso
    public :: d2ijn, d2n, fijn, fjn, fsn, fsumn
    public :: allocate_jasn, deallocate_jasn
    save
contains
    subroutine allocate_jasn()
      use system, only: nelec
        if (.not. allocated(d2ijn)) allocate (d2ijn(nelec, nelec))
        if (.not. allocated(fijn)) allocate (fijn(3, nelec, nelec))
        if (.not. allocated(fjn)) allocate (fjn(3, nelec))
        if (.not. allocated(fsn)) allocate (fsn(nelec, nelec))
    end subroutine allocate_jasn

    subroutine deallocate_jasn()
        if (allocated(fsn)) deallocate (fsn)
        if (allocated(fjn)) deallocate (fjn)
        if (allocated(fijn)) deallocate (fijn)
        if (allocated(d2ijn)) deallocate (d2ijn)
    end subroutine deallocate_jasn

    subroutine allocate_jaso()
      use system, only: nelec
        if (.not. allocated(d2ijo)) allocate (d2ijo(nelec, nelec))
        if (.not. allocated(fijo)) allocate (fijo(3, nelec, nelec))
        if (.not. allocated(fjo)) allocate (fjo(3, nelec))
        if (.not. allocated(fso)) allocate (fso(nelec, nelec))
    end subroutine allocate_jaso

    subroutine deallocate_jaso()
        if (allocated(fso)) deallocate (fso)
        if (allocated(fjo)) deallocate (fjo)
        if (allocated(fijo)) deallocate (fijo)
        if (allocated(d2ijo)) deallocate (d2ijo)
    end subroutine deallocate_jaso

end module jastrow_update

module jaspar4
    !> Arguments: a4, norda, nordb, nordc
    use precision_kinds, only: dp

    implicit none

    real(dp), dimension(:, :, :), allocatable :: a4 !(nordj1,nctype_tot,MWF)
    integer :: norda
    integer :: nordb
    integer :: nordc

    private
    public :: norda, nordb, nordc
    save
contains
    ! subroutine allocate_jaspar4()
    !     use multiple_geo, only: MWF
    !     use precision_kinds, only: dp
    !     use vmc_mod, only: nordj1
    !     if (.not. allocated(a4)) allocate (a4(nordj1, nctype_tot, MWF))
    ! end subroutine allocate_jaspar4

    subroutine deallocate_jaspar4()
        if (allocated(a4)) deallocate (a4)
    end subroutine deallocate_jaspar4

end module jaspar4

module jaspar6
    !> Arguments: asymp_jasa, asymp_jasb, asymp_r, c1_jas6, c1_jas6i, c2_jas6, cutjas, cutjasi
    use precision_kinds, only: dp

    implicit none

    real(dp) :: asymp_r
    real(dp) :: c1_jas6
    real(dp) :: c1_jas6i
    real(dp) :: c2_jas6
    real(dp) :: cutjas
    real(dp) :: cutjasi

    private
    public :: asymp_r, c1_jas6, c1_jas6i, c2_jas6, cutjas, cutjasi
    save
contains
    subroutine allocate_jaspar6()
        use system, only: nctype_tot
        if (.not. allocated(asymp_jasa)) allocate (asymp_jasa(nctype_tot))
        if (.not. allocated(asymp_jasb)) allocate (asymp_jasb(2))
    end subroutine allocate_jaspar6

    subroutine deallocate_jaspar6()
        if (allocated(asymp_jasb)) deallocate (asymp_jasb)
        if (allocated(asymp_jasa)) deallocate (asymp_jasa)
    end subroutine deallocate_jaspar6

end module jaspar6

module jaspointer
    !> Arguments: npoint, npointa

    implicit none

    integer, dimension(:), allocatable :: npoint
    integer, dimension(:), allocatable :: npointa

    private
    public :: npoint, npointa
    public :: allocate_jaspointer, deallocate_jaspointer
    save
contains
    subroutine allocate_jaspointer()
        use vmc_mod, only: nctyp3x
        if (.not. allocated(npoint)) allocate (npoint(nctyp3x), source=0)
        if (.not. allocated(npointa)) allocate (npointa(3*nctyp3x), source=0)
    end subroutine allocate_jaspointer

    subroutine deallocate_jaspointer()
        if (allocated(npointa)) deallocate (npointa)
        if (allocated(npoint)) deallocate (npoint)
    end subroutine deallocate_jaspointer

end module jaspointer



module jastrow
    use precision_kinds, only: dp
    implicit none

    ! From contr2
    integer :: ijas
    integer :: isc
    integer :: ianalyt_lap

    ! From jaspar
    integer :: is
    integer :: nspin1
    integer :: nspin2
    real(dp) :: sspin
    real(dp) :: sspinn

    ! From jaspar3
    real(dp), dimension(:, :, :), allocatable :: b !(nordj1,2,MWF)
    real(dp), dimension(:, :, :), allocatable :: c !(83,MCTYPE,MWF)
    real(dp), dimension(:)      , allocatable :: scalek !(MWF)

    ! From jaspar4
    real(dp), dimension(:, :, :), allocatable :: a4 !(nordj1,nctype_tot,MWF)
    integer :: norda
    integer :: nordb
    integer :: nordc

    ! From jaspar6
    real(dp), dimension(:)      , allocatable :: asymp_jasa !(MCTYPE)
    real(dp), dimension(:)      , allocatable :: asymp_jasb !(2)

    ! From vmc_mod
    integer :: nordj
    integer :: nordj1   ! nordj+1
    integer :: neqsx    ! 6*nordj
    
    save
contains

subroutine allocate_m_jastrow()
    use jastrow_update, only: allocate_jasn
    use jastrow_update, only: allocate_jaso
    use jaspointer, only: allocate_jaspointer

    use multiple_geo, only: MWF
    use system,      only: nctype_tot

    implicit none

    if (.not. allocated(a4))         allocate (a4(nordj1, nctype_tot, MWF))
    if (.not. allocated(b))          allocate (b(nordj1, 2, MWF))
    if (.not. allocated(c))          allocate (c(83, nctype_tot, MWF))
    if (.not. allocated(scalek))     allocate (scalek(MWF))

    call allocate_jasn()
    call allocate_jaso()
    call allocate_jaspointer()
end subroutine allocate_m_jastrow

subroutine allocate_jaspar6()
    use system,      only: nctype_tot
    implicit none

    if (.not. allocated(asymp_jasa)) allocate (asymp_jasa(nctype_tot))
    if (.not. allocated(asymp_jasb)) allocate (asymp_jasb(2))
end subroutine allocate_jaspar6

subroutine deallocate_m_jastrow()
    use jasn, only: deallocate_jasn
    use jaso, only: deallocate_jaso
    use jaspointer, only: deallocate_jaspointer

    implicit none

    if (allocated(a4))         deallocate (a4)
    if (allocated(asymp_jasa)) deallocate (asymp_jasa)
    if (allocated(asymp_jasb)) deallocate (asymp_jasb)
    if (allocated(b))          deallocate (b)
    if (allocated(c))          deallocate (c)
    if (allocated(scalek))     deallocate (scalek)

    call deallocate_jasn()
    call deallocate_jaso()
    call deallocate_jaspointer()
end subroutine deallocate_m_jastrow

end module jastrow

module cuspmat4
    !> Arguments: d, nterms
    use precision_kinds, only: dp

    implicit none

    real(dp), dimension(:, :), allocatable :: d !(neqsx,mterms)
    integer, dimension(:), allocatable :: iwc4 !(neqsx)
    integer :: nterms
    private

    public :: d, iwc4, nterms
    public :: allocate_cuspmat4, deallocate_cuspmat4
    save
contains
    subroutine allocate_cuspmat4()
        use vmc_mod, only: mterms
      use jastrow, only: neqsx
      use jastrow, only: neqsx
        if (.not. allocated(d)) allocate (d(neqsx, mterms))
        if (.not. allocated(iwc4)) allocate (iwc4(neqsx))
    end subroutine allocate_cuspmat4

    subroutine deallocate_cuspmat4()
        if (allocated(iwc4)) deallocate (iwc4)
        if (allocated(d)) deallocate (d)
    end subroutine deallocate_cuspmat4

end module cuspmat4

! module m_jastrow
! contains

! subroutine allocate_m_jastrow()
!     use jastrow_update, only: allocate_jasn
!     use jastrow_update, only: allocate_jaso
!     use jaspar3, only: allocate_jaspar3
! !    use jaspar4, only: allocate_jaspar4
!     use jaspar6, only: allocate_jaspar6
!     use jaspointer, only: allocate_jaspointer

!     implicit none

!     call allocate_jasn()
!     call allocate_jaso()
!     call allocate_jaspar3()
!     ! call allocate_jaspar4()
!     call allocate_jaspar6()
!     call allocate_jaspointer()
! end subroutine allocate_m_jastrow


! subroutine deallocate_m_jastrow()
!     use jastrow_update, only: deallocate_jasn
!     use jastrow_update, only: deallocate_jaso
!     use jaspar3, only: deallocate_jaspar3
!     use jaspar4, only: deallocate_jaspar4
!     use jaspar6, only: deallocate_jaspar6
!     use jaspointer, only: deallocate_jaspointer

!     implicit none

!     call deallocate_jasn()
!     call deallocate_jaso()
!     call deallocate_jaspar3()
!     call deallocate_jaspar4()
!     call deallocate_jaspar6()
!     call deallocate_jaspointer()
! end subroutine deallocate_m_jastrow
! end module
