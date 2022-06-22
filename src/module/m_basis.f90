module basis
    !>  ncent  = number of centers
    !>  zex    = screening constants for each basis function
    !>  betaq  =
    !>  cent   = positions of each center
    !>
    !>  ns     = number of 1s functions at each center
    !>
    !>  npx    = number of 2px functions of each type at each center
    !>  npy    = number of 2py functions of each type at each center
    !>  npz    = number of 2pz functions of each type at each center
    !>
    !>  ndxx   = number of XY d functions at each center
    !>  ndxy   = number of XZ d functions at each center
    !>  ndxz   = number of YZ d functions at each center
    !>  ndyy   = number of YY d functions at each center
    !>  ndyz   = number of YZ d functions at each center
    !>  ndzz   = number of ZZ d functions at each center
    !>
    !>  nfxxx  = number of XYX f functions at each center
    !>  nfxxy  = number of XZX f functions at each center
    !>  nfxxz  = number of YZX f functions at each center
    !>  nfxyy  = number of YYY f functions at each center
    !>  nfxyz  = number of YZY f functions at each center
    !>  nfxzz  = number of ZZX f functions at each center
    !>  nfyyy  = number of YYY f functions at each center
    !>  nfyyz  = number of YZY f functions at each center
    !>  nfyzz  = number of ZZY f functions at each center
    !>  nfzzz  = number of ZZZ f functions at each center



    use precision_kinds, only: dp

    implicit none

    real(dp), dimension(:, :), allocatable :: zex !(MBASIS,MWF)
    real(dp) :: betaq
    integer, dimension(:), allocatable :: ns !(MCTYPE)
    integer, dimension(:), allocatable :: npx !(MCTYPE)
    integer, dimension(:), allocatable :: npy !(MCTYPE)
    integer, dimension(:), allocatable :: npz !(MCTYPE)
    integer, dimension(:), allocatable :: ndxx !(MCTYPE)
    integer, dimension(:), allocatable :: ndxy !(MCTYPE)
    integer, dimension(:), allocatable :: ndxz !(MCTYPE)
    integer, dimension(:), allocatable :: ndyy !(MCTYPE)
    integer, dimension(:), allocatable :: ndyz !(MCTYPE)
    integer, dimension(:), allocatable :: ndzz !(MCTYPE)
    integer, dimension(:), allocatable :: nfxxx !(MCTYPE)
    integer, dimension(:), allocatable :: nfxxy !(MCTYPE)
    integer, dimension(:), allocatable :: nfxxz !(MCTYPE)
    integer, dimension(:), allocatable :: nfxyy !(MCTYPE)
    integer, dimension(:), allocatable :: nfxyz !(MCTYPE)
    integer, dimension(:), allocatable :: nfxzz !(MCTYPE)
    integer, dimension(:), allocatable :: nfyyy !(MCTYPE)
    integer, dimension(:), allocatable :: nfyyz !(MCTYPE)
    integer, dimension(:), allocatable :: nfyzz !(MCTYPE)
    integer, dimension(:), allocatable :: nfzzz !(MCTYPE)
    integer, dimension(:), allocatable :: ngxxxx !(MCTYPE)
    integer, dimension(:), allocatable :: ngxxxy !(MCTYPE)
    integer, dimension(:), allocatable :: ngxxxz !(MCTYPE)
    integer, dimension(:), allocatable :: ngxxyy !(MCTYPE)
    integer, dimension(:), allocatable :: ngxxyz !(MCTYPE)
    integer, dimension(:), allocatable :: ngxxzz !(MCTYPE)
    integer, dimension(:), allocatable :: ngxyyy !(MCTYPE)
    integer, dimension(:), allocatable :: ngxyyz !(MCTYPE)
    integer, dimension(:), allocatable :: ngxyzz !(MCTYPE)
    integer, dimension(:), allocatable :: ngxzzz !(MCTYPE)
    integer, dimension(:), allocatable :: ngyyyy !(MCTYPE)
    integer, dimension(:), allocatable :: ngyyyz !(MCTYPE)
    integer, dimension(:), allocatable :: ngyyzz !(MCTYPE)
    integer, dimension(:), allocatable :: ngyzzz !(MCTYPE)
    integer, dimension(:), allocatable :: ngzzzz !(MCTYPE)

    private
    public :: zex, betaq
    public :: ns, npx, npy, npz, ndxx, ndxy, ndxz, ndyy, ndyz, ndzz
    public :: nfxxx, nfxxy, nfxxz, nfxyy, nfxyz, nfxzz, nfyyy, nfyyz, nfyzz, nfzzz
    public :: ngxxxx, ngxxxy, ngxxxz, ngxxyy, ngxxyz, ngxxzz, ngxyyy, ngxyyz, ngxyzz, ngxzzz, ngyyyy, ngyyyz, ngyyzz, ngyzzz, ngzzzz
    public :: allocate_basis, deallocate_basis
    save
contains
    subroutine allocate_basis()

        ! if (.not. allocated(zex)) allocate(zex(nbasis,MWF))
        ! if (.not. allocated(ns)) allocate(ns(nctype))
        ! if (.not. allocated(np)) allocate(np(3,nctype))
        ! if (.not. allocated(ndxx)) allocate(ndxx(nctype))
        ! if (.not. allocated(ndxy)) allocate(ndxy(nctype))
        ! if (.not. allocated(ndxz)) allocate(ndxz(nctype))
        ! if (.not. allocated(ndyy)) allocate(ndyy(nctype))
        ! if (.not. allocated(ndyz)) allocate(ndyz(nctype))
        ! if (.not. allocated(ndzz)) allocate(ndzz(nctype))
        ! if (.not. allocated(nfxxx)) allocate(nfxxx(nctype))
        ! if (.not. allocated(nfxxy)) allocate(nfxxy(nctype))
        ! if (.not. allocated(nfxxz)) allocate(nfxxz(nctype))
        ! if (.not. allocated(nfxyy)) allocate(nfxyy(nctype))
        ! if (.not. allocated(nfxyz)) allocate(nfxyz(nctype))
        ! if (.not. allocated(nfxzz)) allocate(nfxzz(nctype))
        ! if (.not. allocated(nfyyy)) allocate(nfyyy(nctype))
        ! if (.not. allocated(nfyyz)) allocate(nfyyz(nctype))
        ! if (.not. allocated(nfyzz)) allocate(nfyzz(nctype))
        ! if (.not. allocated(nfzzz)) allocate(nfzzz(nctype))

    end subroutine allocate_basis

    subroutine deallocate_basis()

        if (allocated(ngzzzz)) deallocate (ngzzzz)
        if (allocated(ngyzzz)) deallocate (ngyzzz)
        if (allocated(ngyyzz)) deallocate (ngyyzz)
        if (allocated(ngyyyz)) deallocate (ngyyyz)
        if (allocated(ngyyyy)) deallocate (ngyyyy)
        if (allocated(ngxzzz)) deallocate (ngxzzz)
        if (allocated(ngxyzz)) deallocate (ngxyzz)
        if (allocated(ngxyyz)) deallocate (ngxyyz)
        if (allocated(ngxyyy)) deallocate (ngxyyy)
        if (allocated(ngxxzz)) deallocate (ngxxzz)
        if (allocated(ngxxyz)) deallocate (ngxxyz)
        if (allocated(ngxxyy)) deallocate (ngxxyy)
        if (allocated(ngxxxz)) deallocate (ngxxxz)
        if (allocated(ngxxxy)) deallocate (ngxxxy)
        if (allocated(ngxxxx)) deallocate (ngxxxx)

        if (allocated(nfzzz)) deallocate (nfzzz)
        if (allocated(nfyzz)) deallocate (nfyzz)
        if (allocated(nfyyz)) deallocate (nfyyz)
        if (allocated(nfyyy)) deallocate (nfyyy)
        if (allocated(nfxzz)) deallocate (nfxzz)
        if (allocated(nfxyz)) deallocate (nfxyz)
        if (allocated(nfxyy)) deallocate (nfxyy)
        if (allocated(nfxxz)) deallocate (nfxxz)
        if (allocated(nfxxy)) deallocate (nfxxy)
        if (allocated(nfxxx)) deallocate (nfxxx)

        if (allocated(ndzz)) deallocate (ndzz)
        if (allocated(ndyz)) deallocate (ndyz)
        if (allocated(ndyy)) deallocate (ndyy)
        if (allocated(ndxz)) deallocate (ndxz)
        if (allocated(ndxy)) deallocate (ndxy)
        if (allocated(ndxx)) deallocate (ndxx)

        if (allocated(npz)) deallocate (npz)
        if (allocated(npy)) deallocate (npy)
        if (allocated(npx)) deallocate (npx)
        if (allocated(ns)) deallocate (ns)

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
    real(dp), dimension(:, :, :, :), allocatable :: ab !(2,MRWF,MCTYPE,MFORCE)
    real(dp), dimension(:, :, :, :), allocatable :: ce !(NCOEF,MRWF,MCTYPE,MFORCE)

    private
    public :: ae, ab, ce
    public :: allocate_numexp, deallocate_numexp
    save
contains
    subroutine allocate_numexp()
        use atom, only: nctype_tot
        use numbas_mod, only: MRWF
        use force_mod, only: MFORCE
        use vmc_mod, only: NCOEF
        if (.not. allocated(ae)) allocate (ae(2, MRWF, nctype_tot, MFORCE))
        if (.not. allocated(ab)) allocate (ab(2, MRWF, nctype_tot, MFORCE))
        if (.not. allocated(ce)) allocate (ce(NCOEF, MRWF, nctype_tot, MFORCE))
    end subroutine allocate_numexp

    subroutine deallocate_numexp()
        if (allocated(ce)) deallocate (ce)
        if (allocated(ab)) deallocate (ab)
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
    real(dp), dimension(:,:), allocatable :: rmax !(nrbas, MCTYPE)
    real(dp), dimension(:, :, :, :), allocatable :: rwf !(MRWF_PTS,MRWF,MCTYPE,MWF)

    private
    public :: arg, d2rwf, igrid, iwrwf, nr, nrbas, numr, r0, rwf, rmax
    public :: allocate_numbas, deallocate_numbas
    save
contains
    subroutine allocate_numbas()
        use wfsec, only: nwftype
        use coefs, only: nbasis
        use atom, only: nctype_tot
        use numbas_mod, only: MRWF, MRWF_PTS
        if (.not. allocated(arg)) allocate (arg(nctype_tot))
        if (.not. allocated(d2rwf)) allocate (d2rwf(MRWF_PTS, MRWF, nctype_tot, nwftype))
        if (.not. allocated(igrid)) allocate (igrid(nctype_tot), source=0)
        if (.not. allocated(iwrwf)) allocate (iwrwf(nbasis, nctype_tot), source=0)
        if (.not. allocated(nr)) allocate (nr(nctype_tot), source=0)
        if (.not. allocated(nrbas)) allocate (nrbas(nctype_tot), source=0)
        if (.not. allocated(r0)) allocate (r0(nctype_tot))
        if (.not. allocated(rmax)) allocate (rmax(MRWF,nctype_tot), source=0.0d0) ! This source is needed.
        if (.not. allocated(rwf)) allocate (rwf(MRWF_PTS, MRWF, nctype_tot, nwftype))
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
        ! if (.not. allocated(nbastyp)) allocate (nbastyp(nctype_tot))
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

module m_basis
contains
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
end module

! subroutines required by the trexio modules
module m_trexio_basis
    contains
    double precision function gnorm(exponent, l)
        use precision_kinds,    only: dp
        implicit none
        real(dp), intent (in)       :: exponent
        integer, intent (in)    :: l
        real(dp), parameter     :: pi = 4.0d0*atan(1.0d0)
        !real(dp), parameter     :: pi = 3.1415926535897932

        gnorm = 1.0d0

        if (l .eq. 0) then
            gnorm = (2.d0*exponent)**(3.d0/4.d0)*2.d0*(1.d0/(pi**(1.d0/4.d0)))
        elseif (l .eq. 1) then
            gnorm = (2.d0*exponent)**(5.d0/4.d0)*dsqrt(8.d0/3.d0)*(1.d0/(pi**(1.d0/4.d0)))
        elseif (l .eq. 2) then
            gnorm = (2.d0*exponent)**(7.d0/4.d0)*dsqrt(16.d0/15.d0)*(1.d0/(pi**(1.d0/4.d0)))
        elseif (l .eq. 3) then
            gnorm = (2.d0*exponent)**(9.d0/4.d0)*dsqrt(32.d0/105.d0)*(1.d0/(pi**(1.d0/4.d0)))
        elseif (l .eq. 4) then
            gnorm = (2.d0*exponent)**(11.d0/4.d0)*dsqrt(64.d0/945.d0)*(1.d0/(pi**(1.d0/4.d0)))
        endif

    end function gnorm

end module
