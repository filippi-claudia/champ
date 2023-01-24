module basis
    !>  ncent  = number of centers
    !>  zex    = screening constants for each basis function
    !>  betaq  =
    !>  cent   = positions of each center
    !>  ns      = number of s functions at each center
    !>  np      = number of px functions of each type at each center
    !>  nd      = number of d functions at each center
    !>  nf      = number of f functions at each center

      use precision_kinds, only: dp

    implicit none

    real(dp), dimension(:, :), allocatable :: zex !(MBASIS,MWF)
    real(dp) :: betaq
    integer, dimension(:), allocatable :: ns !(MCTYPE)
    integer, dimension(:), allocatable :: np !(MCTYPE)
    integer, dimension(:), allocatable :: nd !(MCTYPE)
    integer, dimension(:), allocatable :: nf !(MCTYPE)
    integer, dimension(:), allocatable :: ng !(MCTYPE)

    private
    public :: zex, betaq
    public :: ns, np, nd, nf, ng
    public :: allocate_basis, deallocate_basis
    save
contains
    subroutine allocate_basis()

        ! if (.not. allocated(zex)) allocate(zex(nbasis,MWF))
        ! if (.not. allocated(ns)) allocate(ns(nctype))
        ! if (.not. allocated(np)) allocate(np(3,nctype))
        ! if (.not. allocated(nd)) allocate(nd(nctype))
        ! if (.not. allocated(nf)) allocate(nf(nctype))
        ! if (.not. allocated(ng)) allocate(ng(nctype))

    end subroutine allocate_basis

    subroutine deallocate_basis()

        if (allocated(ng)) deallocate (ng)
        if (allocated(nf)) deallocate (nf)
        if (allocated(nd)) deallocate (nd)
        if (allocated(np)) deallocate (np)
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

      use multiple_geo, only: MFORCE
      use numbas_mod, only: MRWF
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
      use multiple_geo, only: MFORCE
      use numbas_mod, only: MRWF
      use system,  only: nctype_tot
      use vmc_mod, only: NCOEF
        if (.not. allocated(ae)) allocate (ae(2, MRWF, nctype_tot, MFORCE))
        if (.not. allocated(ce)) allocate (ce(NCOEF, MRWF, nctype_tot, MFORCE))
    end subroutine allocate_numexp

    subroutine deallocate_numexp()
        if (allocated(ce)) deallocate (ce)
        if (allocated(ae)) deallocate (ae)
    end subroutine deallocate_numexp

end module numexp

module numbas
    !> Arguments: arg, d2rwf, igrid, iwrwf, nr, nrbas, numr, r0, rwf

      use numbas_mod, only: MRWF,MRWF_PTS
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
    real(dp), dimension(:,:), allocatable :: rmaxwf !(nrbas, MCTYPE)
    real(dp), dimension(:, :, :, :), allocatable :: rwf !(MRWF_PTS,MRWF,MCTYPE,MWF)

    private
    public :: arg, d2rwf, igrid, iwrwf, nr, nrbas, numr, r0, rwf, rmaxwf
    public :: allocate_numbas, deallocate_numbas
    save
contains
    subroutine allocate_numbas()
      use coefs,   only: nbasis
      use multiple_geo, only: nwftype
      use numbas_mod, only: MRWF,MRWF_PTS
      use system,  only: nctype_tot
        if (.not. allocated(arg)) allocate (arg(nctype_tot))
        if (.not. allocated(d2rwf)) allocate (d2rwf(MRWF_PTS, MRWF, nctype_tot, nwftype))
        if (.not. allocated(igrid)) allocate (igrid(nctype_tot), source=0)
        if (.not. allocated(iwrwf)) allocate (iwrwf(nbasis, nctype_tot), source=0)
        if (.not. allocated(nr)) allocate (nr(nctype_tot), source=0)
        if (.not. allocated(nrbas)) allocate (nrbas(nctype_tot), source=0)
        if (.not. allocated(r0)) allocate (r0(nctype_tot))
        if (.not. allocated(rmaxwf)) allocate (rmaxwf(MRWF,nctype_tot), source=0.0d0) ! This source is needed.
        if (.not. allocated(rwf)) allocate (rwf(MRWF_PTS, MRWF, nctype_tot, nwftype))
    end subroutine allocate_numbas

    subroutine deallocate_numbas()
        if (allocated(rmaxwf)) deallocate (rmaxwf)
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
      use coefs,   only: nbasis
      use system,  only: nctype_tot
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
      use system,  only: ncent_tot
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
      use numbas,  only: allocate_numbas
      use numbas1, only: allocate_numbas1
      use numbas2, only: allocate_numbas2
      use numexp,  only: allocate_numexp

    implicit none

    ! call allocate_basis()
    call allocate_numexp()
    call allocate_numbas()
    call allocate_numbas1()
    call allocate_numbas2()
end subroutine allocate_m_basis

subroutine deallocate_m_basis()
      use basis,   only: deallocate_basis
      use numbas,  only: deallocate_numbas
      use numbas1, only: deallocate_numbas1
      use numbas2, only: deallocate_numbas2
      use numexp,  only: deallocate_numexp

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
      use coefs,   only: nbasis
    implicit none

    integer, dimension(5)       :: slm_per_l = (/1, 3, 6, 10, 15/) !s,p,d,f,g
    integer                     :: basis_num_shell
    integer, allocatable        :: index_slm(:)             !(nbasis)
    integer, allocatable        :: ao_radial_index(:)       !(nbasis)
    integer, allocatable        :: num_rad_per_cent(:)      !(ncent_tot)
    integer, allocatable        :: num_ao_per_cent(:)       !(ncent_tot)
    integer, allocatable        :: basis_shell_ang_mom(:)   !(nshell)

    private
    public :: slm_per_l, index_slm, num_rad_per_cent, num_ao_per_cent
    public :: basis_num_shell, basis_shell_ang_mom
    public :: ao_radial_index
    public :: gnorm

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
