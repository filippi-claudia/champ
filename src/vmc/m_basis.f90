module basis
    !> Arguments: zex, betaq, n1s, n2s, n2p, n3s, n3p, n3dzr, n3dx2, n3dxy, n3dxz, n3dyz, n4s, n4p,
    !> n4fxxx, n4fyyy, n4fzzz, n4fxxy, n4fxxz, n4fyyx, n4fyyz, n4fzzx, n4fzzy, n4fxyz, nsa, npa,
    !> ndzra, ndz2a, ndxya, ndxza, ndyza
    use precision_kinds, only: dp

    implicit none

    !  ncent  = number of centers
    !  zex    = screening constants for each basis function
    !  cent   = positions of each center
    !  pecent = potential energy of the centers
    !  znuc   = charge of the nuclei (centers)
    !  n1s    = number of 1s functions at each center
    !  n2s    = number of 2s functions at each center
    !  n2p    = number of 2p functions of each type at each center
    !  n3s    = number of 3s functions at each center
    !  n3p    = number of 3p functions of each type at each center
    !  n3dzr  = number of z**2-r**2 d functions at each center
    !  n3dx2  = number of x**2-y**2 d functions at each center
    !  n3dxy  = number of xy d functions at each center
    !  n3dxz  = number of xz d functions at each center
    !  n3dyz  = number of yz d functions at each center
    !  n4s    = number of 4s functions at each center
    !  n4p    = number of 4p functions of each type at each center

    real(dp), dimension(:, :), allocatable :: zex !(MBASIS,MWF)
    real(dp) :: betaq
    integer, dimension(:), allocatable :: n1s !(MCTYPE)
    integer, dimension(:), allocatable :: n2s !(MCTYPE)
    integer, dimension(:, :), allocatable :: n2p !(3,MCTYPE)
    integer, dimension(:), allocatable :: n3s !(MCTYPE)
    integer, dimension(:, :), allocatable :: n3p !(3,MCTYPE)
    integer, dimension(:), allocatable :: n3dzr !(MCTYPE)
    integer, dimension(:), allocatable :: n3dx2 !(MCTYPE)
    integer, dimension(:), allocatable :: n3dxy !(MCTYPE)
    integer, dimension(:), allocatable :: n3dxz !(MCTYPE)
    integer, dimension(:), allocatable :: n3dyz !(MCTYPE)
    integer, dimension(:), allocatable :: n4s !(MCTYPE)
    integer, dimension(:, :), allocatable :: n4p !(3,MCTYPE)
    integer, dimension(:), allocatable :: n4fxxx !(MCTYPE)
    integer, dimension(:), allocatable :: n4fyyy !(MCTYPE)
    integer, dimension(:), allocatable :: n4fzzz !(MCTYPE)
    integer, dimension(:), allocatable :: n4fxxy !(MCTYPE)
    integer, dimension(:), allocatable :: n4fxxz !(MCTYPE)
    integer, dimension(:), allocatable :: n4fyyx !(MCTYPE)
    integer, dimension(:), allocatable :: n4fyyz !(MCTYPE)
    integer, dimension(:), allocatable :: n4fzzx !(MCTYPE)
    integer, dimension(:), allocatable :: n4fzzy !(MCTYPE)
    integer, dimension(:), allocatable :: n4fxyz !(MCTYPE)
    integer, dimension(:), allocatable :: nsa !(MCTYPE)
    integer, dimension(:, :), allocatable :: npa !(3,MCTYPE)
    integer, dimension(:), allocatable :: ndzra !(MCTYPE)
    integer, dimension(:), allocatable :: ndz2a !(MCTYPE)
    integer, dimension(:), allocatable :: ndxya !(MCTYPE)
    integer, dimension(:), allocatable :: ndxza !(MCTYPE)
    integer, dimension(:), allocatable :: ndx2a !(MCTYPE)
    integer, dimension(:), allocatable :: ndyza !(MCTYPE)

    private
    public :: zex, betaq, n1s, n2s, n2p, n3s, n3p, n3dzr, n3dx2, n3dxy, n3dxz, n3dyz
    public :: n4s, n4p, n4fxxx, n4fyyy, n4fzzz, n4fxxy, n4fxxz, n4fyyx, n4fyyz
    public :: n4fzzx, n4fzzy, n4fxyz, nsa, npa, ndzra, ndz2a, ndxya, ndxza, ndyza
    public :: ndx2a
    public :: allocate_basis, deallocate_basis
    save
contains
    subroutine allocate_basis()

        ! if (.not. allocated(zex)) allocate (zex(MBASIS, MWF), source=0.0_dp)
        ! if (.not. allocated(n1s)) allocate (n1s(MCTYPE), source=0.0_dp)
        ! if (.not. allocated(n2s)) allocate (n2s(MCTYPE), source=0.0_dp)
        ! if (.not. allocated(n2p)) allocate (n2p(3, MCTYPE), source=0.0_dp)
        ! if (.not. allocated(n3s)) allocate (n3s(MCTYPE), source=0.0_dp)
        ! if (.not. allocated(n3p)) allocate (n3p(3, MCTYPE), source=0.0_dp)
        ! if (.not. allocated(n3dzr)) allocate (n3dzr(MCTYPE), source=0.0_dp)
        ! if (.not. allocated(n3dx2)) allocate (n3dx2(MCTYPE), source=0.0_dp)
        ! if (.not. allocated(n3dxy)) allocate (n3dxy(MCTYPE), source=0.0_dp)
        ! if (.not. allocated(n3dxz)) allocate (n3dxz(MCTYPE), source=0.0_dp)
        ! if (.not. allocated(n3dyz)) allocate (n3dyz(MCTYPE), source=0.0_dp)
        ! if (.not. allocated(n4s)) allocate (n4s(MCTYPE), source=0.0_dp)
        ! if (.not. allocated(n4p)) allocate (n4p(3, MCTYPE), source=0.0_dp)
        ! if (.not. allocated(n4fxxx)) allocate (n4fxxx(MCTYPE), source=0.0_dp)
        ! if (.not. allocated(n4fyyy)) allocate (n4fyyy(MCTYPE), source=0.0_dp)
        ! if (.not. allocated(n4fzzz)) allocate (n4fzzz(MCTYPE), source=0.0_dp)
        ! if (.not. allocated(n4fxxy)) allocate (n4fxxy(MCTYPE), source=0.0_dp)
        ! if (.not. allocated(n4fxxz)) allocate (n4fxxz(MCTYPE), source=0.0_dp)
        ! if (.not. allocated(n4fyyx)) allocate (n4fyyx(MCTYPE), source=0.0_dp)
        ! if (.not. allocated(n4fyyz)) allocate (n4fyyz(MCTYPE), source=0.0_dp)
        ! if (.not. allocated(n4fzzx)) allocate (n4fzzx(MCTYPE), source=0.0_dp)
        ! if (.not. allocated(n4fzzy)) allocate (n4fzzy(MCTYPE), source=0.0_dp)
        ! if (.not. allocated(n4fxyz)) allocate (n4fxyz(MCTYPE), source=0.0_dp)
        ! if (.not. allocated(nsa)) allocate (nsa(MCTYPE), source=0.0_dp)
        ! if (.not. allocated(npa)) allocate (npa(3, MCTYPE), source=0.0_dp)
        ! if (.not. allocated(ndzra)) allocate (ndzra(MCTYPE), source=0.0_dp)
        ! if (.not. allocated(ndz2a)) allocate (ndz2a(MCTYPE), source=0.0_dp)
        ! if (.not. allocated(ndxya)) allocate (ndxya(MCTYPE), source=0.0_dp)
        ! if (.not. allocated(ndxza)) allocate (ndxza(MCTYPE), source=0.0_dp)
        ! if (.not. allocated(ndx2a)) allocate (ndx2a(MCTYPE), source=0.0_dp)
        ! if (.not. allocated(ndyza)) allocate (ndyza(MCTYPE), source=0.0_dp)
    end subroutine allocate_basis

    subroutine deallocate_basis()
        if (allocated(ndyza)) deallocate (ndyza)
        if (allocated(ndx2a)) deallocate (ndx2a)
        if (allocated(ndxza)) deallocate (ndxza)
        if (allocated(ndxya)) deallocate (ndxya)
        if (allocated(ndz2a)) deallocate (ndz2a)
        if (allocated(ndzra)) deallocate (ndzra)
        if (allocated(npa)) deallocate (npa)
        if (allocated(nsa)) deallocate (nsa)
        if (allocated(n4fxyz)) deallocate (n4fxyz)
        if (allocated(n4fzzy)) deallocate (n4fzzy)
        if (allocated(n4fzzx)) deallocate (n4fzzx)
        if (allocated(n4fyyz)) deallocate (n4fyyz)
        if (allocated(n4fyyx)) deallocate (n4fyyx)
        if (allocated(n4fxxz)) deallocate (n4fxxz)
        if (allocated(n4fxxy)) deallocate (n4fxxy)
        if (allocated(n4fzzz)) deallocate (n4fzzz)
        if (allocated(n4fyyy)) deallocate (n4fyyy)
        if (allocated(n4fxxx)) deallocate (n4fxxx)
        if (allocated(n4p)) deallocate (n4p)
        if (allocated(n4s)) deallocate (n4s)
        if (allocated(n3dyz)) deallocate (n3dyz)
        if (allocated(n3dxz)) deallocate (n3dxz)
        if (allocated(n3dxy)) deallocate (n3dxy)
        if (allocated(n3dx2)) deallocate (n3dx2)
        if (allocated(n3dzr)) deallocate (n3dzr)
        if (allocated(n3p)) deallocate (n3p)
        if (allocated(n3s)) deallocate (n3s)
        if (allocated(n2p)) deallocate (n2p)
        if (allocated(n2s)) deallocate (n2s)
        if (allocated(n1s)) deallocate (n1s)
        if (allocated(zex)) deallocate (zex)
    end subroutine deallocate_basis

end module basis

module numbas_mod
    !> Arguments: MRWF_PTS, MRWF

    implicit none

    integer, parameter :: MRWF_PTS = 4000
    integer, parameter :: MRWF = 200
    private
    public :: MRWF, MRWF_PTS
    save

end module numbas_mod

module numexp
    !> Arguments: ae, ce

    use numbas_mod, only: MRWF
    use force_mod, only: MFORCE
    use precision_kinds, only: dp
    use vmc_mod, only: NCOEF

    implicit none

    real(dp), dimension(:, :, :, :), allocatable :: ae !(2,MRWF,MCTYPE,MFORCE)
    real(dp), dimension(:, :, :, :), allocatable :: ce !(NCOEF,MRWF,MCTYPE,MFORCE)

    private
    public :: ae, ce
    public :: allocate_numexp, deallocate_numexp
    save
contains
    subroutine allocate_numexp()
        use atom, only: nctype_tot
        use numbas_mod, only: MRWF
        use force_mod, only: MFORCE
        use vmc_mod, only: NCOEF
        if (.not. allocated(ae)) allocate (ae(2, MRWF, nctype_tot, MFORCE), source=0.0_dp)
        if (.not. allocated(ce)) allocate (ce(NCOEF, MRWF, nctype_tot, MFORCE), source=0.0_dp)
    end subroutine allocate_numexp

    subroutine deallocate_numexp()
        if (allocated(ce)) deallocate (ce)
        if (allocated(ae)) deallocate (ae)
    end subroutine deallocate_numexp

end module numexp

module numbas
    !> Arguments: arg, d2rwf, igrid, iwrwf, nr, nrbas, numr, r0, rwf

    use numbas_mod, only: MRWF, MRWF_PTS
    use precision_kinds, only: dp

    implicit none

    real(dp), dimension(:), allocatable :: arg !(MCTYPE)
    real(dp), dimension(:, :, :, :), allocatable :: d2rwf !(MRWF_PTS,MRWF,MCTYPE,MWF)
    integer, dimension(:), allocatable :: igrid !(MCTYPE)
    integer, dimension(:, :), allocatable :: iwrwf !(MBASIS,MCTYPE)
    integer, dimension(:), allocatable :: nr !(MCTYPE)
    integer, dimension(:), allocatable :: nrbas !(MCTYPE)
    integer :: numr
    real(dp), dimension(:), allocatable :: r0 !(MCTYPE)
    real(dp), dimension(:, :, :, :), allocatable :: rwf !(MRWF_PTS,MRWF,MCTYPE,MWF)

    private
    public :: arg, d2rwf, igrid, iwrwf, nr, nrbas, numr, r0, rwf
    public :: allocate_numbas, deallocate_numbas
    save
contains
    subroutine allocate_numbas()
        use wfsec, only: nwftype
        use coefs, only: nbasis
        use atom, only: nctype_tot
        use numbas_mod, only: MRWF, MRWF_PTS
        if (.not. allocated(arg)) allocate (arg(nctype_tot), source=0.0_dp)
        if (.not. allocated(d2rwf)) allocate (d2rwf(MRWF_PTS, MRWF, nctype_tot, nwftype), source=0.0_dp)
        if (.not. allocated(igrid)) allocate (igrid(nctype_tot), source=0)
        if (.not. allocated(iwrwf)) allocate (iwrwf(nbasis, nctype_tot), source=0)
        if (.not. allocated(nr)) allocate (nr(nctype_tot), source=0)
        if (.not. allocated(nrbas)) allocate (nrbas(nctype_tot), source=0)
        if (.not. allocated(r0)) allocate (r0(nctype_tot), source=0.0_dp)
        if (.not. allocated(rwf)) allocate (rwf(MRWF_PTS, MRWF, nctype_tot, nwftype), source=0.0_dp)
    end subroutine allocate_numbas

    subroutine deallocate_numbas()
        if (allocated(rwf)) deallocate (rwf)
        if (allocated(r0)) deallocate (r0)
        if (allocated(nrbas)) deallocate (nrbas)
        if (allocated(nr)) deallocate (nr)
        if (allocated(iwrwf)) deallocate (iwrwf)
        if (allocated(igrid)) deallocate (igrid)
        if (allocated(d2rwf)) deallocate (d2rwf)
        if (allocated(arg)) deallocate (arg)
    end subroutine deallocate_numbas

end module numbas

module numbas1
    !> Arguments: iwlbas, nbastyp

    implicit none

    integer, dimension(:, :), allocatable :: iwlbas !(MBASIS,MCTYPE)
    integer, dimension(:), allocatable :: nbastyp !(MCTYPE)

    private
    public :: iwlbas, nbastyp
    public :: allocate_numbas1, deallocate_numbas1
    save
contains
    subroutine allocate_numbas1()
        use coefs, only: nbasis
        use atom, only: nctype_tot
        if (.not. allocated(iwlbas)) allocate (iwlbas(nbasis, nctype_tot), source=0)
        ! if (.not. allocated(nbastyp)) allocate (nbastyp(nctype_tot), source=0.0_dp)
    end subroutine allocate_numbas1

    subroutine deallocate_numbas1()
        if (allocated(nbastyp)) deallocate (nbastyp)
        if (allocated(iwlbas)) deallocate (iwlbas)
    end subroutine deallocate_numbas1

end module numbas1

module numbas2
    !> Arguments: ibas0, ibas1

    implicit none

    integer, dimension(:), allocatable :: ibas0 !(MCENT)
    integer, dimension(:), allocatable :: ibas1 !(MCENT)

    private
    public :: ibas0, ibas1
    public :: allocate_numbas2, deallocate_numbas2
    save
contains
    subroutine allocate_numbas2()
        use atom, only: ncent_tot
        if (.not. allocated(ibas0)) allocate (ibas0(ncent_tot), source=0)
        if (.not. allocated(ibas1)) allocate (ibas1(ncent_tot), source=0)
    end subroutine allocate_numbas2

    subroutine deallocate_numbas2()
        if (allocated(ibas1)) deallocate (ibas1)
        if (allocated(ibas0)) deallocate (ibas0)
    end subroutine deallocate_numbas2

end module numbas2

subroutine allocate_m_basis()
    ! use basis, only: allocate_basis
    use numexp, only: allocate_numexp
    use numbas, only: allocate_numbas
    use numbas1, only: allocate_numbas1
    use numbas2, only: allocate_numbas2

    implicit none

    ! call allocate_basis()
    call allocate_numexp()
    call allocate_numbas()
    call allocate_numbas1()
    call allocate_numbas2()
end subroutine allocate_m_basis

subroutine deallocate_m_basis()
    use basis, only: deallocate_basis
    use numexp, only: deallocate_numexp
    use numbas, only: deallocate_numbas
    use numbas1, only: deallocate_numbas1
    use numbas2, only: deallocate_numbas2

    implicit none

    call deallocate_basis()
    call deallocate_numexp()
    call deallocate_numbas()
    call deallocate_numbas1()
    call deallocate_numbas2()
end subroutine deallocate_m_basis
