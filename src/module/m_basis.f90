!> Module that contains the basis set information.
!> @var zex Screening constants for each basis function.
!> @var betaq
!> @var ns Number of s functions at each center.
!> @var np Number of px functions of each type at each center.
!> @var nd Number of d functions at each center.
!> @var nf Number of f functions at each center.
!> @var ng Number of g functions at each center.
!> @subroutine allocate_basis Allocates memory for the basis set.
!> @subroutine deallocate_basis Deallocates memory for the basis set.
!> @author Ravindra Shinde
module basis
    use precision_kinds, only: dp

    implicit none

    !> Screening constants for each basis function.
    real(dp), dimension(:, :), allocatable :: zex

    !> Beta coefficient parameter.
    real(dp) :: betaq

    !> Number of s functions at each center.
    integer, dimension(:), allocatable :: ns

    !> Number of px functions at each center.
    integer, dimension(:), allocatable :: np

    !> Number of d functions at each center.
    integer, dimension(:), allocatable :: nd

    !> Number of f functions at each center.
    integer, dimension(:), allocatable :: nf

    !> Number of g functions at each center.
    integer, dimension(:), allocatable :: ng

    private
    public :: zex, betaq
    public :: ns, np, nd, nf, ng
    public :: allocate_basis, deallocate_basis
    save
contains
    !> Allocates memory for the basis set.
    !> @note deprecated
    subroutine allocate_basis()
        ! Allocate basis is called while reading the basis set from the input file.
    end subroutine allocate_basis

    !> Deallocates memory for the basis set.
    !> This subroutine deallocates memory for the basis set parameters
    !> to avoid memory leaks after calculations are finished.
    subroutine deallocate_basis()

        if (allocated(ng)) deallocate (ng)
        if (allocated(nf)) deallocate (nf)
        if (allocated(nd)) deallocate (nd)
        if (allocated(np)) deallocate (np)
        if (allocated(ns)) deallocate (ns)
        if (allocated(zex)) deallocate (zex)
    end subroutine deallocate_basis

end module basis

!> Module defining constants for numerical basis representation.
!> Provides parameters used in various numerical modules for managing basis sets.
module numbas_mod
    !> Arguments: MRWF_PTS, MRWF

    implicit none

    !> Maximum number of wavefunction points.
    integer, parameter :: MRWF_PTS = 4000

    !> Maximum number of wavefunctions.
    integer, parameter :: MRWF = 200

    private
    public :: MRWF, MRWF_PTS
    save

end module numbas_mod

!> Module numexp
module numexp

      use multiple_geo, only: MFORCE
      use numbas_mod, only: MRWF
      use precision_kinds, only: dp
      use vmc_mod, only: NCOEF

    implicit none

    !> Array for storing the ae coefficients.
    real(dp), dimension(:, :, :, :), allocatable :: ae !(2,MRWF,MCTYPE,MFORCE)

    !> Array for storing the ce coefficients.
    real(dp), dimension(:, :, :, :), allocatable :: ce !(NCOEF,MRWF,MCTYPE,MFORCE)

    private
    public :: ae, ce
    public :: allocate_numexp, deallocate_numexp
    save
contains
    !> Allocates memory for the ae and ce arrays.
    subroutine allocate_numexp()
      use multiple_geo, only: MFORCE
      use numbas_mod, only: MRWF
      use system,  only: nctype_tot
      use vmc_mod, only: NCOEF
        if (.not. allocated(ae)) allocate (ae(2, MRWF, nctype_tot, MFORCE))
        if (.not. allocated(ce)) allocate (ce(NCOEF, MRWF, nctype_tot, MFORCE))
    end subroutine allocate_numexp

    !> Deallocates memory for the ae and ce arrays.
    subroutine deallocate_numexp()
        if (allocated(ce)) deallocate (ce)
        if (allocated(ae)) deallocate (ae)
    end subroutine deallocate_numexp

end module numexp

!> Module numbas for numerical basis information
module numbas
    !> Arguments: arg, d2rwf, igrid, iwrwf, nr, nrbas, numr, r0, rwf

      use numbas_mod, only: MRWF,MRWF_PTS
      use precision_kinds, only: dp

    implicit none

    !> Array of basis function arguments.
    real(dp), dimension(:), allocatable :: arg !(MCTYPE)

    !> Array of second derivatives of the wavefunctions.
    real(dp), dimension(:, :, :, :), allocatable :: d2rwf !(MRWF_PTS,MRWF,MCTYPE,MWF)

    !> Array of grid indices.
    integer, dimension(:), allocatable :: igrid !(MCTYPE)

    !> Array of wavefunction indices.
    integer, dimension(:, :), allocatable :: iwrwf !(MBASIS,MCTYPE)

    !> Array of number of radial points.
    integer, dimension(:), allocatable :: nr !(MCTYPE)

    !> Array of number of shells in basis functions.
    integer, dimension(:), allocatable :: nrbas !(MCTYPE)

    !> Number of radial points.
    integer :: numr

    !> Array r0
    real(dp), dimension(:), allocatable :: r0 !(MCTYPE)

    !> Array of wavefunction maximum radii for each basis function.
    real(dp), dimension(:,:), allocatable :: rmaxwf !(nrbas, MCTYPE)

    !> Array of wavefunctions.
    real(dp), dimension(:, :, :, :), allocatable :: rwf !(MRWF_PTS,MRWF,MCTYPE,MWF)

    private
    public :: arg, d2rwf, igrid, iwrwf, nr, nrbas, numr, r0, rwf, rmaxwf
    public :: allocate_numbas, deallocate_numbas
    save
contains
    !> Allocates memory for the numerical basis arrays.
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

    !> Deallocates memory for the numerical basis arrays.
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

!> Module numbas1 for numerical basis information
module numbas1
    !> Arguments: iwlbas, nbastyp

    implicit none

    !> Array of basis function indices.
    integer, dimension(:, :), allocatable :: iwlbas !(MBASIS,MCTYPE)

    !> Array of basis function types.
    integer, dimension(:), allocatable :: nbastyp !(MCTYPE)

    private
    public :: iwlbas, nbastyp
    public :: allocate_numbas1, deallocate_numbas1
    save
contains
    !> Allocates memory for the numerical basis arrays.
    subroutine allocate_numbas1()
      use coefs,   only: nbasis
      use system,  only: nctype_tot
        if (.not. allocated(iwlbas)) allocate (iwlbas(nbasis, nctype_tot), source=0)
        ! if (.not. allocated(nbastyp)) allocate (nbastyp(nctype_tot))
    end subroutine allocate_numbas1

    !> Deallocates memory for the numerical basis arrays.
    subroutine deallocate_numbas1()
        if (allocated(nbastyp)) deallocate (nbastyp)
        if (allocated(iwlbas)) deallocate (iwlbas)
    end subroutine deallocate_numbas1

end module numbas1

!> Module numbas2 for numerical basis information
module numbas2
    !> Arguments: ibas0, ibas1

    implicit none

    !> Array of basis function indices.
    integer, dimension(:), allocatable :: ibas0 !(MCENT)

    !> Array of basis function indices.
    integer, dimension(:), allocatable :: ibas1 !(MCENT)

    private
    public :: ibas0, ibas1
    public :: allocate_numbas2, deallocate_numbas2
    save
contains
    !> Allocates memory for the numerical basis arrays.
    subroutine allocate_numbas2()
      use system,  only: ncent_tot
        if (.not. allocated(ibas0)) allocate (ibas0(ncent_tot), source=0)
        if (.not. allocated(ibas1)) allocate (ibas1(ncent_tot), source=0)
    end subroutine allocate_numbas2

    !> Deallocates memory for the numerical basis arrays.
    subroutine deallocate_numbas2()
        if (allocated(ibas1)) deallocate (ibas1)
        if (allocated(ibas0)) deallocate (ibas0)
    end subroutine deallocate_numbas2

end module numbas2

!> Module m_basis for basis set information.
module m_basis
contains
!> Subroutine to allocate memory for the basis set.
subroutine allocate_m_basis()
      use numbas,  only: allocate_numbas
      use numbas1, only: allocate_numbas1
      use numbas2, only: allocate_numbas2
      use numexp,  only: allocate_numexp

    implicit none

    call allocate_numexp()
    call allocate_numbas()
    call allocate_numbas1()
    call allocate_numbas2()
end subroutine allocate_m_basis

!> Subroutine to deallocate memory for the basis set.
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

!> Module trexio basis for recreating basis grid from trexio file.
module m_trexio_basis
      use coefs,   only: nbasis
    implicit none

    !> Fixed array of slm per angular momentum.
    integer, dimension(5)       :: slm_per_l = (/1, 3, 6, 10, 15/) !s,p,d,f,g

    !> Total number of shells in the basis set.
    integer                     :: basis_num_shell

    !> Array of indices for the spherical harmonics.
    integer, allocatable        :: index_slm(:)             !(nbasis)

    !> Array of indices for the radial ao basis functions.
    integer, allocatable        :: ao_radial_index(:)       !(nbasis)

    !> Array of number of radial basis functions per center.
    integer, allocatable        :: num_rad_per_cent(:)      !(ncent_tot)

    !> Array of number of ao basis functions per center.
    integer, allocatable        :: num_ao_per_cent(:)       !(ncent_tot)

    !> Array of angular momentum for each shell.
    integer, allocatable        :: basis_shell_ang_mom(:)   !(nshell)

    private
    public :: slm_per_l, index_slm, num_rad_per_cent, num_ao_per_cent
    public :: basis_num_shell, basis_shell_ang_mom
    public :: ao_radial_index
    public :: gnorm

    contains

    !> Function to calculate the normalization constant for the basis functions.
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
