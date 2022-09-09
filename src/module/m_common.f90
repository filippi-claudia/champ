
module Bloc
    !> Arguments: b, tildem, xmat
    use precision_kinds, only: dp
    use vmc_mod, only: norb_tot
    use optwf_parms, only: nparmj

    implicit none

    real(dp), dimension(:, :), allocatable :: b !(norb_tot,MELEC)
    real(dp), dimension(:, :, :), allocatable :: tildem !(MELEC,norb_tot,2)
    real(dp), dimension(:, :), allocatable :: xmat !(MELEC**2,2)

    !> Former Bloc_da
    real(dp), dimension(:, :, :, :), allocatable :: b_da !(3,MELEC,norb_tot,MCENT)

    !> former Bloc_dj
    real(dp), dimension(:, :, :), allocatable :: b_dj !(norb_tot,MELEC,nparmj)

    private
    public :: b, tildem, xmat
    public :: b_da
    public :: b_dj
    public :: allocate_Bloc, deallocate_Bloc
    save
contains
    subroutine allocate_Bloc()
        use const, only: nelec
        use coefs, only: norb
        use atom, only: ncent_tot, ncent
        use vmc_mod, only: norb_tot
        use coefs, only: norb
        use optwf_parms, only: nparmj

        if (.not. allocated(b)) allocate (b(norb_tot, nelec))
        if (.not. allocated(tildem)) allocate (tildem(nelec, norb_tot, 2))
        if (.not. allocated(xmat)) allocate (xmat(nelec**2, 2))
        if (.not. allocated(b_da)) allocate (b_da(3, nelec, norb_tot, ncent_tot))
        if (.not. allocated(b_dj)) allocate (b_dj(norb_tot, nelec, nparmj))

    end subroutine allocate_Bloc

    subroutine deallocate_Bloc()
        if (allocated(b_dj)) deallocate (b_dj)
        if (allocated(b_da)) deallocate (b_da)
        if (allocated(xmat)) deallocate (xmat)
        if (allocated(tildem)) deallocate (tildem)
        if (allocated(b)) deallocate (b)
    end subroutine deallocate_Bloc

end module Bloc

module bparm
    !> Arguments: nocuspb, nspin2b

    implicit none

    integer :: nocuspb
    integer :: nspin2b

    private
    public :: nocuspb, nspin2b
    save
end module bparm

module casula
    !> Arguments: i_vpsp, icasula, t_vpsp
    use pseudo_mod, only: MPS_QUAD
    use precision_kinds, only: dp

    implicit none

    integer :: i_vpsp
    integer :: icasula
    real(dp), dimension(:, :, :), allocatable :: t_vpsp !(MCENT,MPS_QUAD,MELEC)

    private
    public :: i_vpsp, icasula, t_vpsp
    public :: allocate_casula, deallocate_casula
    save
contains
    subroutine allocate_casula()
        use const, only: nelec
        use atom, only: ncent_tot
        use pseudo_mod, only: MPS_QUAD
        if (.not. allocated(t_vpsp)) allocate (t_vpsp(ncent_tot, MPS_QUAD, nelec))
    end subroutine allocate_casula

    subroutine deallocate_casula()
        if (allocated(t_vpsp)) deallocate (t_vpsp)
    end subroutine deallocate_casula

end module casula

module distance_mod
    use precision_kinds, only: dp
    use vmc_mod, only: nmat_dim2

    implicit none

    real(dp), dimension(:, :, :), allocatable :: rshift !(3, MELEC, MCENT)
    real(dp), dimension(:, :), allocatable :: r_en !(MELEC, MCENT)
    real(dp), dimension(:, :, :), allocatable :: rvec_en !(3, MELEC, MCENT)
    real(dp), dimension(:), allocatable :: r_ee !(nmat_dim2)
    real(dp), dimension(:, :), allocatable :: rvec_ee !(3, nmat_dim2)

    private
    public :: rshift
    public :: r_en
    public :: rvec_en
    public :: r_ee
    public :: rvec_ee
    public :: allocate_distance_mod, deallocate_distance_mod
    save
contains
    subroutine allocate_distance_mod()
        use const, only: nelec
        use atom, only: ncent_tot
        use vmc_mod, only: nmat_dim2
        if (.not. allocated(r_en)) allocate (r_en(nelec, ncent_tot))
        if (.not. allocated(rvec_en)) allocate (rvec_en(3, nelec, ncent_tot))
        if (.not. allocated(r_ee)) allocate (r_ee(nmat_dim2))
        if (.not. allocated(rvec_ee)) allocate (rvec_ee(3, nmat_dim2))
        if (.not. allocated(rshift)) allocate (rshift(3, nelec, ncent_tot))
        rshift=0.d0
    end subroutine allocate_distance_mod

    subroutine deallocate_distance_mod()
        if (allocated(rvec_en)) deallocate (rvec_en)
        if (allocated(r_en)) deallocate (r_en)
        if (allocated(rvec_en)) deallocate (rvec_ee)
        if (allocated(r_en)) deallocate (r_ee)
        if (allocated(rshift)) deallocate (rshift)
    end subroutine deallocate_distance_mod

end module distance_mod

module distances_sav
    !> Arguments: r_ee_sav, r_en_sav, rshift_sav, rvec_ee_sav, rvec_en_sav
    use precision_kinds, only: dp

    implicit none

    real(dp), dimension(:), allocatable :: r_ee_sav !(MELEC)
    real(dp), dimension(:), allocatable :: r_en_sav !(MCENT)
    real(dp), dimension(:, :), allocatable :: rshift_sav !(3,MCENT)
    real(dp), dimension(:, :), allocatable :: rvec_ee_sav !(3,MELEC)
    real(dp), dimension(:, :), allocatable :: rvec_en_sav !(3,MCENT)

    private
    public   :: r_ee_sav, r_en_sav, rshift_sav, rvec_ee_sav, rvec_en_sav
    public :: allocate_distances_sav, deallocate_distances_sav
    save
contains
    subroutine allocate_distances_sav()
        use const, only: nelec
        use atom, only: ncent_tot
        if (.not. allocated(r_ee_sav)) allocate (r_ee_sav(nelec))
        if (.not. allocated(r_en_sav)) allocate (r_en_sav(ncent_tot))
        if (.not. allocated(rshift_sav)) allocate (rshift_sav(3, ncent_tot))
        if (.not. allocated(rvec_ee_sav)) allocate (rvec_ee_sav(3, nelec))
        if (.not. allocated(rvec_en_sav)) allocate (rvec_en_sav(3, ncent_tot))
    end subroutine allocate_distances_sav

    subroutine deallocate_distances_sav()
        if (allocated(rvec_en_sav)) deallocate (rvec_en_sav)
        if (allocated(rvec_ee_sav)) deallocate (rvec_ee_sav)
        if (allocated(rshift_sav)) deallocate (rshift_sav)
        if (allocated(r_en_sav)) deallocate (r_en_sav)
        if (allocated(r_ee_sav)) deallocate (r_ee_sav)
    end subroutine deallocate_distances_sav

end module distances_sav

module embed
    !> Never called

    implicit none

    ! tables for hartree and exchange potentials
    ! max interpolation order
    integer, parameter :: MXTABI = 6
    ! max number of  hartree, exchangetables
    integer, parameter :: MXHTAB = 1
    integer, parameter :: MXXTAB = 1
    integer, parameter :: MXTAB = MXHTAB + MXXTAB
    ! max number of gridpoints in hartree table
    integer, parameter :: MXHGRID = 1
    ! max number of gridpoints in exchange tables
    integer, parameter :: MXXGRID = 1

    private
    public :: MXTABI, MXTAB, MXHTAB, MXXTAB
    public :: MXHGRID, MXXGRID
    save
end module embed

module gauss_ecp
    !> only used in read_gauss ....
    !> Arguments: ecp_coef, ecp_exponent, necp_power, necp_term
    use pseudo_mod, only: MPS_L, MGAUSS
    use precision_kinds, only: dp

    implicit none

    real(dp), dimension(:, :, :), allocatable :: ecp_coef !(MGAUSS,MPS_L,MCTYPE)
    real(dp), dimension(:, :, :), allocatable :: ecp_exponent !(MGAUSS,MPS_L,MCTYPE)
    integer, dimension(:, :, :), allocatable :: necp_power !(MGAUSS,MPS_L,MCTYPE)
    integer, dimension(:, :), allocatable :: necp_term !(MPS_L,MCTYPE)

    private
    public   ::  ecp_coef, ecp_exponent, necp_power, necp_term
    public :: allocate_gauss_ecp, deallocate_gauss_ecp
    save
contains
    subroutine allocate_gauss_ecp()
        use atom, only: nctype_tot
        use pseudo_mod, only: MPS_L, MGAUSS
        if (.not. allocated(ecp_coef)) allocate (ecp_coef(MGAUSS, MPS_L, nctype_tot))
        if (.not. allocated(ecp_exponent)) allocate (ecp_exponent(MGAUSS, MPS_L, nctype_tot))
        if (.not. allocated(necp_power)) allocate (necp_power(MGAUSS, MPS_L, nctype_tot), source=0)
        if (.not. allocated(necp_term)) allocate (necp_term(MPS_L, nctype_tot), source=0)
    end subroutine allocate_gauss_ecp

    subroutine deallocate_gauss_ecp()
        if (allocated(necp_term)) deallocate (necp_term)
        if (allocated(necp_power)) deallocate (necp_power)
        if (allocated(ecp_exponent)) deallocate (ecp_exponent)
        if (allocated(ecp_coef)) deallocate (ecp_coef)
    end subroutine deallocate_gauss_ecp

end module gauss_ecp

module gradjerrb
    !> Arguments: nbj_current, ngrad_jas_bcum, ngrad_jas_blocks, njb_current

    implicit none

    integer :: nbj_current
    integer :: ngrad_jas_bcum
    integer :: ngrad_jas_blocks


    private
    public :: nbj_current, ngrad_jas_bcum, ngrad_jas_blocks

    save
end module gradjerrb

module insout
    !> something to do with grid ...
    !> Arguments: inout, inside

    implicit none

    integer :: inout
    integer :: inside

    private
    public :: inout, inside
    save
end module insout

module jd_scratch
    !> only for (jacobi) davidson in linear method
    !> Arguments: qr, rr
    use sr_mod, only: mparm
    use precision_kinds, only: dp

    implicit none

    real(dp), dimension(:), allocatable :: qr !(mparm)
    real(dp), dimension(:), allocatable :: rr !(mparm)

    private
    public :: qr, rr
    public :: allocate_jd_scratch, deallocate_jd_scratch
    save
contains
    subroutine allocate_jd_scratch()
        use sr_mod, only: mparm
        if (.not. allocated(qr)) allocate (qr(mparm))
        if (.not. allocated(rr)) allocate (rr(mparm))
    end subroutine allocate_jd_scratch

    subroutine deallocate_jd_scratch()
        if (allocated(rr)) deallocate (rr)
        if (allocated(qr)) deallocate (qr)
    end subroutine deallocate_jd_scratch

end module jd_scratch

module linear_norm
    !> Arguments: oav, ci_oav
    use optci, only: mxciterm
    use precision_kinds, only: dp

    implicit none

    real(dp), dimension(:), allocatable :: ci_oav !(mxciterm)

    private
    public :: ci_oav
    ! public :: oav
    public :: allocate_linear_norm, deallocate_linear_norm
    save
contains
    subroutine allocate_linear_norm()
        use optci, only: mxciterm
        if (.not. allocated(ci_oav)) allocate (ci_oav(mxciterm))
    end subroutine allocate_linear_norm

    subroutine deallocate_linear_norm()
        if (allocated(ci_oav)) deallocate (ci_oav)
    end subroutine deallocate_linear_norm

end module linear_norm

module multidet
    !> Arguments: iactv, irepcol_det, ireporb_det, ivirt, iwundet, kref, numrep_det

    implicit none

    integer, dimension(:), allocatable :: iactv !(2)
    integer, dimension(:, :, :), allocatable :: irepcol_det !(MELEC,MDET,2)
    integer, dimension(:, :, :), allocatable :: ireporb_det !(MELEC,MDET,2)
    integer, dimension(:), allocatable :: ivirt !(2)
    integer, dimension(:, :), allocatable :: iwundet !(MDET,2)
    integer :: kref
    integer :: kref_old
    integer :: ndet_req
    integer :: kchange
    integer :: kref_fixed
    integer, dimension(:, :), allocatable :: numrep_det !(MDET,2)
    integer, dimension(:, :), allocatable :: k_det !(MDET,2)
    integer, dimension(:, :), allocatable :: k_det2 !(MDET,2)
    integer, dimension(:, :), allocatable :: k_aux !(MDET,2)
    integer, dimension(:), allocatable :: ndetiab !(2)
    integer, dimension(:), allocatable :: ndetiab2 !(2)
    integer, dimension(:), allocatable :: ndetsingle !(2)

    private
    public :: iactv, irepcol_det, ireporb_det, ivirt, iwundet, kref, kref_fixed, numrep_det, k_det, ndetiab, ndet_req, k_det2, k_aux, ndetiab2, ndetsingle, kref_old, kchange
    public :: allocate_multidet, deallocate_multidet
    save
contains
    subroutine allocate_multidet()
        use const, only: nelec
        use dets, only: ndet
        if (.not. allocated(iactv)) allocate (iactv(2), source=0)
        if (.not. allocated(irepcol_det)) allocate (irepcol_det(nelec, ndet, 2), source=0)
        if (.not. allocated(ireporb_det)) allocate (ireporb_det(nelec, ndet, 2), source=0)
        if (.not. allocated(ivirt)) allocate (ivirt(2), source=0)
        if (.not. allocated(iwundet)) allocate (iwundet(ndet, 2), source=0)
        if (.not. allocated(numrep_det)) allocate (numrep_det(ndet, 2), source=0)
        if (.not. allocated(k_det)) allocate (k_det(ndet, 2), source=0)
        if (.not. allocated(ndetiab)) allocate (ndetiab(2), source=0)
        if (.not. allocated(k_det2)) allocate (k_det2(ndet, 2), source=0)
        if (.not. allocated(k_aux)) allocate (k_aux(ndet, 2), source=0)
        if (.not. allocated(ndetiab2)) allocate (ndetiab2(2), source=0)
        if (.not. allocated(ndetsingle)) allocate (ndetsingle(2), source=0)
    end subroutine allocate_multidet

    subroutine deallocate_multidet()
        if (allocated(numrep_det)) deallocate (numrep_det)
        if (allocated(iwundet)) deallocate (iwundet)
        if (allocated(ivirt)) deallocate (ivirt)
        if (allocated(ireporb_det)) deallocate (ireporb_det)
        if (allocated(irepcol_det)) deallocate (irepcol_det)
        if (allocated(iactv)) deallocate (iactv)
        if (allocated(k_det)) deallocate (k_det)
        if (allocated(ndetiab)) deallocate (ndetiab)
        if (allocated(k_det)) deallocate (k_det)
        if (allocated(k_det2)) deallocate (k_det2)
        if (allocated(ndetiab2)) deallocate (ndetiab2)
        if (allocated(ndetsingle)) deallocate (ndetsingle)
    end subroutine deallocate_multidet

end module multidet

module multimat
    !> Arguments: aa, wfmat
    use precision_kinds, only: dp

    implicit none

    real(dp), dimension(:, :, :), allocatable :: aa !(MELEC,norb_tot,2)
    real(dp), dimension(:, :, :), allocatable :: wfmat !(MDET,MEXCIT**2,2)

    private
    public :: aa, wfmat
    public :: allocate_multimat, deallocate_multimat
    save
contains
    subroutine allocate_multimat()
        use const, only: nelec
        use dets, only: ndet
        use coefs, only: norb
        use vmc_mod, only: norb_tot
        use vmc_mod, only: MEXCIT
        if (.not. allocated(aa)) allocate (aa(nelec, norb_tot, 2))
        if (.not. allocated(wfmat)) allocate (wfmat(ndet, MEXCIT**2, 2))
    end subroutine allocate_multimat

    subroutine deallocate_multimat()
        if (allocated(wfmat)) deallocate (wfmat)
        if (allocated(aa)) deallocate (aa)
    end subroutine deallocate_multimat

end module multimat

module multimatn
    !> Arguments: aan, wfmatn
    use precision_kinds, only: dp

    implicit none

    real(dp), dimension(:), allocatable :: aan !(MELEC*norb_tot)
    real(dp), dimension(:, :), allocatable :: wfmatn !(MDET, MEXCIT**2)

    private
    public :: aan, wfmatn
    public :: allocate_multimatn, deallocate_multimatn
    save
contains
    subroutine allocate_multimatn()
        use const, only: nelec
        use dets, only: ndet
        use vmc_mod, only: norb_tot
        use vmc_mod, only: MEXCIT
        if (.not. allocated(aan)) allocate (aan(nelec*norb_tot))
        if (.not. allocated(wfmatn)) allocate (wfmatn(ndet,MEXCIT**2))
    end subroutine allocate_multimatn

    subroutine deallocate_multimatn()
        if (allocated(wfmatn)) deallocate (wfmatn)
        if (allocated(aan)) deallocate (aan)
    end subroutine deallocate_multimatn

end module multimatn

module multislater
    !> Arguments: detiab
    use precision_kinds, only: dp

    implicit none

    real(dp), dimension(:, :), allocatable :: detiab !(MDET,2)
    !> DMC variables:


    private
    public :: detiab
    ! public :: detu, detd
    public :: allocate_multislater, deallocate_multislater
    save
contains
    subroutine allocate_multislater()
        use dets, only: ndet
        if (.not. allocated(detiab)) allocate(detiab(ndet, 2))
    end subroutine allocate_multislater

    subroutine deallocate_multislater()
        if (allocated(detiab)) deallocate (detiab)
    end subroutine deallocate_multislater

end module multislater

module multislatern
    !> Arguments: ddorbn, detn, dorbn, orbn
    use precision_kinds, only: dp
    use vmc_mod, only: norb_tot

    implicit none

    real(dp), dimension(:), allocatable :: ddorbn !(norb_tot)
    real(dp), dimension(:), allocatable :: detn !(MDET)
    real(dp), dimension(:, :), allocatable :: dorbn !(norb_tot, 3)
    real(dp), dimension(:), allocatable :: orbn !(norb_tot)
    private

    public ::  ddorbn, detn, dorbn, orbn
    public :: allocate_multislatern, deallocate_multislatern
    save
contains
    subroutine allocate_multislatern()
        use dets, only: ndet
        use coefs, only: norb
        use vmc_mod, only: norb_tot
        if (.not. allocated(ddorbn)) allocate (ddorbn(norb_tot))
        if (.not. allocated(detn)) allocate (detn(ndet))
        if (.not. allocated(dorbn)) allocate (dorbn(norb_tot, 3))
        if (.not. allocated(orbn)) allocate (orbn(norb_tot))
    end subroutine allocate_multislatern

    subroutine deallocate_multislatern()
        if (allocated(orbn)) deallocate (orbn)
        if (allocated(dorbn)) deallocate (dorbn)
        if (allocated(detn)) deallocate (detn)
        if (allocated(ddorbn)) deallocate (ddorbn)
    end subroutine deallocate_multislatern

end module multislatern

module m_icount
    !> Arguments: icount_ci, icount_orb, icount_prop

    implicit none

    integer :: icount_ci = 1
    integer :: icount_orb = 1
    integer :: icount_prop = 1

    private
    public :: icount_ci, icount_orb, icount_prop
    save
end module m_icount

module ncusp
    !> Never called !
    !> Arguments: ncnstr, ncuspc, nfock, nfockc, norbc

    implicit none

    integer :: ncnstr
    integer :: ncuspc
    integer :: nfock
    integer :: nfockc
    integer :: norbc

    private
    public :: ncnstr, ncuspc, nfock, nfockc, norbc
    save
end module ncusp

module orbval
    !> Arguments: ddorb, dorb, nadorb, ndetorb, orb
    use precision_kinds, only: dp
    use vmc_mod, only: norb_tot

    implicit none

    real(dp), dimension(:, :), allocatable :: ddorb !(norb_tot,MELEC)
    real(dp), dimension(:, :, :), allocatable :: dorb !(norb_tot,MELEC,3)
    integer :: nadorb
    integer :: ndetorb
    real(dp), dimension(:, :), allocatable :: orb !(MELEC,norb_tot)

    private
    public :: ddorb, dorb, nadorb, ndetorb, orb
    public :: allocate_orbval, deallocate_orbval
    save
contains
    subroutine allocate_orbval()
        use const, only: nelec
        use coefs, only: norb
        use precision_kinds, only: dp
        if (.not. allocated(ddorb)) allocate (ddorb(norb_tot, nelec))
        if (.not. allocated(dorb)) allocate (dorb(norb_tot, nelec, 3))
        if (.not. allocated(orb)) allocate (orb(nelec, norb_tot))
    end subroutine allocate_orbval

    subroutine deallocate_orbval()
        if (allocated(orb)) deallocate (orb)
        if (allocated(dorb)) deallocate (dorb)
        if (allocated(ddorb)) deallocate (ddorb)
    end subroutine deallocate_orbval

end module orbval

module phifun
    !> Arguments: d2phin, d2phin_all, d3phin, dphin, n0_ibasis, n0_ic, n0_nbasis, phin
    use precision_kinds, only: dp

    implicit none

    real(dp), dimension(:, :), allocatable :: d2phin !(MBASIS,MELEC)
    real(dp), dimension(:, :, :, :), allocatable :: d2phin_all !(3,3,MBASIS,MELEC)
    real(dp), dimension(:, :, :), allocatable :: d3phin !(3,MBASIS,MELEC)
    real(dp), dimension(:, :, :), allocatable :: dphin !(MBASIS,MELEC,3)
    integer, dimension(:, :), allocatable :: n0_ibasis !(MBASIS,MELEC)
    integer, dimension(:, :), allocatable :: n0_ic !(MBASIS,MELEC)
    integer, dimension(:), allocatable :: n0_nbasis !(MELEC)
    real(dp), dimension(:, :), allocatable :: phin !(MBASIS,MELEC)

    private
    public :: d2phin, d2phin_all, d3phin, dphin, n0_ibasis, n0_ic, n0_nbasis, phin
    public :: allocate_phifun, deallocate_phifun
    save
contains
    subroutine allocate_phifun()
        use const, only: nelec
        use coefs, only: nbasis

        if (.not. allocated(d2phin)) allocate (d2phin(nbasis, nelec))
        if (.not. allocated(d2phin_all)) allocate (d2phin_all(3, 3, nbasis, nelec))
        if (.not. allocated(d3phin)) allocate (d3phin(3, nbasis, nelec))
        if (.not. allocated(dphin)) allocate (dphin(nbasis, nelec, 3))
        if (.not. allocated(n0_ibasis)) allocate (n0_ibasis(nbasis, nelec), source=0)
        if (.not. allocated(n0_ic)) allocate (n0_ic(nbasis, nelec), source=0)
        if (.not. allocated(n0_nbasis)) allocate (n0_nbasis(nelec), source=0)
        if (.not. allocated(phin)) allocate (phin(nbasis, nelec))
    end subroutine allocate_phifun

    subroutine deallocate_phifun()
        if (allocated(phin)) deallocate (phin)
        if (allocated(n0_nbasis)) deallocate (n0_nbasis)
        if (allocated(n0_ic)) deallocate (n0_ic)
        if (allocated(n0_ibasis)) deallocate (n0_ibasis)
        if (allocated(dphin)) deallocate (dphin)
        if (allocated(d3phin)) deallocate (d3phin)
        if (allocated(d2phin_all)) deallocate (d2phin_all)
        if (allocated(d2phin)) deallocate (d2phin)
    end subroutine deallocate_phifun

end module phifun

module qua
    !> Arguments: nquad, wq, xq, xq0, yq, yq0, zq, zq0
    use pseudo_mod, only: MPS_QUAD
    use precision_kinds, only: dp

    implicit none

    integer :: nquad
    real(dp), dimension(:), allocatable :: wq !(MPS_QUAD)
    real(dp), dimension(:), allocatable :: xq !(MPS_QUAD)
    real(dp), dimension(:), allocatable :: xq0 !(MPS_QUAD)
    real(dp), dimension(:), allocatable :: yq !(MPS_QUAD)
    real(dp), dimension(:), allocatable :: yq0 !(MPS_QUAD)
    real(dp), dimension(:), allocatable :: zq !(MPS_QUAD)
    real(dp), dimension(:), allocatable :: zq0 !(MPS_QUAD)

    private
    public :: nquad, wq, xq, xq0, yq, yq0, zq, zq0
    public :: allocate_qua, deallocate_qua
    save
contains
    subroutine allocate_qua()
        use pseudo_mod, only: MPS_QUAD
        if (.not. allocated(wq)) allocate (wq(MPS_QUAD))
        if (.not. allocated(xq)) allocate (xq(MPS_QUAD))
        if (.not. allocated(xq0)) allocate (xq0(MPS_QUAD))
        if (.not. allocated(yq)) allocate (yq(MPS_QUAD))
        if (.not. allocated(yq0)) allocate (yq0(MPS_QUAD))
        if (.not. allocated(zq)) allocate (zq(MPS_QUAD))
        if (.not. allocated(zq0)) allocate (zq0(MPS_QUAD))
    end subroutine allocate_qua

    subroutine deallocate_qua()
        if (allocated(zq0)) deallocate (zq0)
        if (allocated(zq)) deallocate (zq)
        if (allocated(yq0)) deallocate (yq0)
        if (allocated(yq)) deallocate (yq)
        if (allocated(xq0)) deallocate (xq0)
        if (allocated(xq)) deallocate (xq)
        if (allocated(wq)) deallocate (wq)
    end subroutine deallocate_qua

end module qua

module b_tmove
    !> Arguments: b_t, iskip
    use pseudo_mod, only: MPS_QUAD
    use precision_kinds, only: dp
    use vmc_mod, only: norb_tot

    implicit none

    real(dp), dimension(:, :, :, :), allocatable :: b_t !(norb_tot,MPS_QUAD,MCENT,MELEC)
    integer, dimension(:, :), allocatable :: iskip !(MELEC,MCENT)

    private
    public :: b_t, iskip
    public :: allocate_b_tmove, deallocate_b_tmove
    save
contains
    subroutine allocate_b_tmove()
        use const, only: nelec
        use atom, only: ncent_tot
        use vmc_mod, only: norb_tot
        use qua, only: nquad
        use contr3, only: mode
        if(index(mode,'dmc').ne.0) then
          if (.not. allocated(b_t)) allocate (b_t(norb_tot, nquad, ncent_tot, nelec))
        endif
        if (.not. allocated(iskip)) allocate (iskip(nelec, ncent_tot), source=0)
    end subroutine allocate_b_tmove

    subroutine deallocate_b_tmove()
        if (allocated(iskip)) deallocate (iskip)
        if (allocated(b_t)) deallocate (b_t)
    end subroutine deallocate_b_tmove

end module b_tmove

module scale_more
    !> Arguments: dd3
    use precision_kinds, only: dp

    implicit none

    real(dp) :: dd3

    private
    public :: dd3
    save
end module scale_more

module scratch
    !> Arguments: denergy_det, dtildem
    use precision_kinds, only: dp
    use vmc_mod, only: norb_tot

    implicit none

    real(dp), dimension(:, :), allocatable :: denergy_det !(MDET,2)
    real(dp), dimension(:, :, :), allocatable :: dtildem !(MELEC,norb_tot,2)

    private
    public :: denergy_det, dtildem
    public :: allocate_scratch, deallocate_scratch
    save
contains
    subroutine allocate_scratch()
        use const, only: nelec
        use dets, only: ndet
        use vmc_mod, only: norb_tot
        if (.not. allocated(denergy_det)) allocate (denergy_det(ndet, 2))
        if (.not. allocated(dtildem)) allocate (dtildem(nelec, norb_tot, 2))
    end subroutine allocate_scratch

    subroutine deallocate_scratch()
        if (allocated(dtildem)) deallocate (dtildem)
        if (allocated(denergy_det)) deallocate (denergy_det)
    end subroutine deallocate_scratch

end module scratch

module slater
    !> Arguments: d2dx2, ddx, fp, fpp, slmi

    use precision_kinds, only: dp
    use vmc_mod, only: nmat_dim

    implicit none

    real(dp), dimension(:), allocatable :: d2dx2 !(MELEC)
    real(dp), dimension(:, :), allocatable :: ddx !(3,MELEC)
    real(dp), dimension(:, :, :), allocatable :: fp !(3,nmat_dim,2)
    real(dp), dimension(:, :), allocatable :: fpp !(nmat_dim,2)
    real(dp), dimension(:, :), allocatable :: slmi !(nmat_dim,2)
    !> DMC extra variables:
    real(dp), dimension(:,:), allocatable :: fpd !(3,nmat_dim)
    real(dp), dimension(:), allocatable :: fppd !(nmat_dim)
    real(dp), dimension(:), allocatable :: fppu !(nmat_dim)
    real(dp), dimension(:,:), allocatable :: fpu !(3,nmat_dim)


    private
    public :: d2dx2, ddx, fp, fpp, slmi
    public :: fpd, fppd, fppu, fpu
    ! public :: slmui, slmdi
    public :: allocate_slater, deallocate_slater
    save

contains
    subroutine allocate_slater()
        use const, only: nelec
        use vmc_mod, only: nmat_dim
        if (.not. allocated(d2dx2)) allocate(d2dx2(nelec))
        if (.not. allocated(ddx)) allocate(ddx(3, nelec))
        if (.not. allocated(fp)) allocate(fp(3, nmat_dim, 2))
        if (.not. allocated(fpp)) allocate(fpp(nmat_dim, 2))
        if (.not. allocated(slmi)) allocate(slmi(nmat_dim, 2))
        if (.not. allocated(fpd))  allocate(fpd(3,nmat_dim))
        if (.not. allocated(fppd)) allocate(fppd(nmat_dim))
        if (.not. allocated(fppu)) allocate(fppu(nmat_dim))
        if (.not. allocated(fpu))  allocate(fpu(3,nmat_dim))

    end subroutine allocate_slater

    subroutine deallocate_slater()
        if (allocated(slmi)) deallocate(slmi)
        if (allocated(fpp)) deallocate(fpp)
        if (allocated(fp)) deallocate(fp)
        if (allocated(ddx)) deallocate(ddx)
        if (allocated(d2dx2)) deallocate(d2dx2)
        if (allocated(fpd))  deallocate(fpd)
        if (allocated(fppd)) deallocate(fppd)
        if (allocated(fppu)) deallocate(fppu)
        if (allocated(fpu))  deallocate(fpu)

    end subroutine deallocate_slater

end module slater

module slatn
    !> Arguments: slmin
    use precision_kinds, only: dp
    use vmc_mod, only: nmat_dim

    implicit none

    real(dp), dimension(:), allocatable :: slmin !(nmat_dim)

    private
    public :: slmin
    public :: allocate_slatn, deallocate_slatn
    save
contains
    subroutine allocate_slatn()
        use vmc_mod, only: nmat_dim
        if (.not. allocated(slmin)) allocate (slmin(nmat_dim))
    end subroutine allocate_slatn

    subroutine deallocate_slatn()
        if (allocated(slmin)) deallocate (slmin)
    end subroutine deallocate_slatn

end module slatn

module svd_mod
    ! Not used anywhere !
    !> Arguments:

    implicit none

    integer, parameter :: MBUF = 10000
    integer, parameter :: MXDIM = 3000
    private
    public :: MBUF, MXDIM
    save
end module svd_mod

module vardep
    !> Arguments: cdep, iwdepend, nvdepend
    use precision_kinds, only: dp

    real(dp), dimension(:, :, :), allocatable :: cdep !(neqsx,83,MCTYPE)
    integer, dimension(:, :, :), allocatable :: iwdepend !(neqsx,83,MCTYPE)
    integer, dimension(:, :), allocatable :: nvdepend !(neqsx,MCTYPE)

    private
    public :: cdep, iwdepend, nvdepend
    public :: allocate_vardep, deallocate_vardep
    save
contains
    subroutine allocate_vardep()
        use atom, only: nctype_tot
        use vmc_mod, only: neqsx
        if (.not. allocated(cdep)) allocate (cdep(neqsx, 83, nctype_tot))
        if (.not. allocated(iwdepend)) allocate (iwdepend(neqsx, 83, nctype_tot), source=0)
        if (.not. allocated(nvdepend)) allocate (nvdepend(neqsx, nctype_tot), source=0)
    end subroutine allocate_vardep

    subroutine deallocate_vardep()
        if (allocated(nvdepend)) deallocate (nvdepend)
        if (allocated(iwdepend)) deallocate (iwdepend)
        if (allocated(cdep)) deallocate (cdep)
    end subroutine deallocate_vardep

end module vardep

module velocity_jastrow
    !> Arguments: vj, vjn
    use precision_kinds, only: dp

    implicit none

    real(dp), dimension(:, :), allocatable :: vj !(3,MELEC)
    real(dp), dimension(:, :), allocatable :: vjn !(3,MELEC)

    private
    public :: vj, vjn
    public :: allocate_velocity_jastrow, deallocate_velocity_jastrow
    save
contains
    subroutine allocate_velocity_jastrow()
        use const, only: nelec
        if (.not. allocated(vj)) allocate (vj(3, nelec))
        if (.not. allocated(vjn)) allocate (vjn(3, nelec))
    end subroutine allocate_velocity_jastrow

    subroutine deallocate_velocity_jastrow()
        if (allocated(vjn)) deallocate (vjn)
        if (allocated(vj)) deallocate (vj)
    end subroutine deallocate_velocity_jastrow

end module velocity_jastrow

module ycompact
    !> Arguments: dymat, ymat
    use precision_kinds, only: dp
    use vmc_mod, only: norb_tot
    use mstates_mod, only: MSTATES

    implicit none

    real(dp), dimension(:, :, :, :), allocatable :: dymat !(norb_tot,MELEC,2,MSTATES)
    real(dp), dimension(:, :, :, :), allocatable :: ymat !(norb_tot,MELEC,2,MSTATES)

    private
    public :: dymat, ymat
    public :: allocate_ycompact, deallocate_ycompact
    save
contains
    subroutine allocate_ycompact()
        use const, only: nelec
        use vmc_mod, only: norb_tot
        use mstates_mod, only: MSTATES
        if (.not. allocated(dymat)) allocate (dymat(norb_tot, nelec, 2, MSTATES))
        if (.not. allocated(ymat)) allocate (ymat(norb_tot, nelec, 2, MSTATES))
    end subroutine allocate_ycompact

    subroutine deallocate_ycompact()
        if (allocated(ymat)) deallocate (ymat)
        if (allocated(dymat)) deallocate (dymat)
    end subroutine deallocate_ycompact

end module ycompact

module ycompactn
    !> Arguments: ymatn
    use precision_kinds, only: dp
    use vmc_mod, only: norb_tot
    use mstates_mod, only: MSTATES

    implicit none

    real(dp), dimension(:, :, :), allocatable :: ymatn !(norb_tot,MELEC,MSTATES)

    private
    public :: ymatn
    public :: allocate_ycompactn, deallocate_ycompactn
    save
contains
    subroutine allocate_ycompactn()
        use const, only: nelec
        use vmc_mod, only: norb_tot
        use mstates_mod, only: MSTATES
        if (.not. allocated(ymatn)) allocate (ymatn(norb_tot, nelec, MSTATES))
    end subroutine allocate_ycompactn

    subroutine deallocate_ycompactn()
        if (allocated(ymatn)) deallocate (ymatn)
    end subroutine deallocate_ycompactn

end module ycompactn

module zcompact
    !> Arguments: aaz, dzmat, emz, zmat
    use precision_kinds, only: dp
    use vmc_mod, only: norb_tot
    use mstates_mod, only: MSTATES

    implicit none

    real(dp), dimension(:, :, :, :), allocatable :: aaz !(MELEC,MELEC,2,MSTATES)
    real(dp), dimension(:, :, :, :), allocatable :: dzmat !(norb_tot,MELEC,2,MSTATES)
    real(dp), dimension(:, :, :, :), allocatable :: emz !(MELEC,MELEC,2,MSTATES)
    real(dp), dimension(:, :, :, :), allocatable :: zmat !(norb_tot,MELEC,2,MSTATES)

    private
    public :: aaz, dzmat, emz, zmat
    public :: allocate_zcompact, deallocate_zcompact
    save
contains
    subroutine allocate_zcompact()
        use const, only: nelec
        use vmc_mod, only: norb_tot
        use mstates_mod, only: MSTATES
        if (.not. allocated(aaz)) allocate (aaz(nelec, nelec, 2, MSTATES))
        if (.not. allocated(dzmat)) allocate (dzmat(norb_tot, nelec, 2, MSTATES))
        if (.not. allocated(emz)) allocate (emz(nelec, nelec, 2, MSTATES))
        if (.not. allocated(zmat)) allocate (zmat(norb_tot, nelec, 2, MSTATES))
    end subroutine allocate_zcompact

    subroutine deallocate_zcompact()
        if (allocated(zmat)) deallocate (zmat)
        if (allocated(emz)) deallocate (emz)
        if (allocated(dzmat)) deallocate (dzmat)
        if (allocated(aaz)) deallocate (aaz)
    end subroutine deallocate_zcompact

end module zcompact

module zmatrix
    !> Arguments: czcart, czint, czcart_ref, izcmat, izmatrix
    use precision_kinds, only: dp

    implicit none

    real(dp), dimension(:, :), allocatable :: czcart !(3,MCENT)
    real(dp), dimension(:, :), allocatable :: czint !(3,MCENT)
    real(dp), dimension(:, :), allocatable :: czcart_ref !(3,3)
    integer, dimension(:, :), allocatable :: izcmat !(3,MCENT)
    integer :: izmatrix

    private
    public :: czcart, czint, czcart_ref, izcmat, izmatrix
    public :: allocate_zmatrix, deallocate_zmatrix
    save
contains
    subroutine allocate_zmatrix()
        use atom, only: ncent_tot
        if (.not. allocated(czcart)) allocate (czcart(3, ncent_tot))
        if (.not. allocated(czint)) allocate (czint(3, ncent_tot))
        if (.not. allocated(czcart_ref)) allocate (czcart_ref(3, 3))
        if (.not. allocated(izcmat)) allocate (izcmat(3, ncent_tot), source=0)
    end subroutine allocate_zmatrix

    subroutine deallocate_zmatrix()
        if (allocated(izcmat)) deallocate (izcmat)
        if (allocated(czcart_ref)) deallocate (czcart_ref)
        if (allocated(czint)) deallocate (czint)
        if (allocated(czcart)) deallocate (czcart)
    end subroutine deallocate_zmatrix

end module zmatrix

module zmatrix_grad
    !> never called
    !> Arguments: transform_grd
    use precision_kinds, only: dp

    implicit none

    real(dp), dimension(:, :), allocatable :: transform_grd !(MCENT3,MCENT3)

    private
    public :: transform_grd
    public :: allocate_zmatrix_grad, deallocate_zmatrix_grad
    save
contains
    subroutine allocate_zmatrix_grad()
        use vmc_mod, only: ncent3

        if (.not. allocated(transform_grd)) allocate (transform_grd(ncent3, ncent3))
    end subroutine allocate_zmatrix_grad

    subroutine deallocate_zmatrix_grad()
        if (allocated(transform_grd)) deallocate (transform_grd)
    end subroutine deallocate_zmatrix_grad

end module zmatrix_grad

module m_common
contains
subroutine allocate_m_common()

    use atom, only: allocate_atom
    use b_tmove, only: allocate_b_tmove
    use Bloc, only: allocate_Bloc
    use casula, only: allocate_casula
    ! use coefs, only: allocate_coefs
    use csfs, only: allocate_csfs
    use cuspmat4, only: allocate_cuspmat4
    ! use dets, only: allocate_dets
    use dets_equiv, only: allocate_dets_equiv
    use distance_mod, only: allocate_distance_mod
    use distances_sav, only: allocate_distances_sav
    use gauss_ecp, only: allocate_gauss_ecp
    use jd_scratch, only: allocate_jd_scratch
    use linear_norm, only: allocate_linear_norm
    use multidet, only: allocate_multidet
    use multimat, only: allocate_multimat
    use multimatn, only: allocate_multimatn
    use multislater, only: allocate_multislater
    use multislatern, only: allocate_multislatern
    use orbval, only: allocate_orbval
    use phifun, only: allocate_phifun
    use qua, only: allocate_qua
    use scratch, only: allocate_scratch
    use slater, only: allocate_slater
    use slatn, only: allocate_slatn
    use vardep, only: allocate_vardep
    use velocity_jastrow, only: allocate_velocity_jastrow
    ! use wfsec, only: allocate_wfsec
    use ycompact, only: allocate_ycompact
    use ycompactn, only: allocate_ycompactn
    use zcompact, only: allocate_zcompact
    use zmatrix, only: allocate_zmatrix
    use zmatrix_grad, only: allocate_zmatrix_grad

    implicit none

    call allocate_atom()
    call allocate_b_tmove()
    call allocate_Bloc()
    call allocate_casula()
    ! call allocate_coefs()
    call allocate_csfs()
    call allocate_cuspmat4()
    ! call allocate_dets()
    call allocate_dets_equiv()
    call allocate_distance_mod()
    call allocate_distances_sav()
    call allocate_gauss_ecp()
    call allocate_jd_scratch()
    call allocate_linear_norm()
    call allocate_multidet()
    call allocate_multimat()
    call allocate_multimatn()
    call allocate_multislater()
    call allocate_multislatern()
    call allocate_orbval()
    call allocate_phifun()
    call allocate_qua()
    call allocate_scratch()
    call allocate_slater()
    call allocate_slatn()
    call allocate_vardep()
    call allocate_velocity_jastrow()
    ! call allocate_wfsec()
    call allocate_ycompact()
    call allocate_ycompactn()
    call allocate_zcompact()
    call allocate_zmatrix()
    call allocate_zmatrix_grad()
end subroutine allocate_m_common

subroutine deallocate_m_common()
    ! use atom, only: deallocate_atom
    use b_tmove, only: deallocate_b_tmove
    use Bloc, only: deallocate_Bloc
    use casula, only: deallocate_casula
    use coefs, only: deallocate_coefs
    use csfs, only: deallocate_csfs
    use cuspmat4, only: deallocate_cuspmat4
    use dets, only: deallocate_dets
    use dets_equiv, only: deallocate_dets_equiv
    use distance_mod, only: deallocate_distance_mod
    use distances_sav, only: deallocate_distances_sav
    use gauss_ecp, only: deallocate_gauss_ecp
    use jd_scratch, only: deallocate_jd_scratch
    use linear_norm, only: deallocate_linear_norm
    use multidet, only: deallocate_multidet
    use multimat, only: deallocate_multimat
    use multimatn, only: deallocate_multimatn
    use multislater, only: deallocate_multislater
    use multislatern, only: deallocate_multislatern
    use orbval, only: deallocate_orbval
    use phifun, only: deallocate_phifun
    use qua, only: deallocate_qua
    use scratch, only: deallocate_scratch
    use slater, only: deallocate_slater
    use slatn, only: deallocate_slatn
    use vardep, only: deallocate_vardep
    use velocity_jastrow, only: deallocate_velocity_jastrow
    use wfsec, only: deallocate_wfsec
    use ycompact, only: deallocate_ycompact
    use ycompactn, only: deallocate_ycompactn
    use zcompact, only: deallocate_zcompact
    use zmatrix, only: deallocate_zmatrix
    use zmatrix_grad, only: deallocate_zmatrix_grad

    implicit none

    ! call deallocate_atom()
    call deallocate_b_tmove()
    call deallocate_Bloc()
    call deallocate_casula()
    call deallocate_coefs()
    call deallocate_csfs()
    call deallocate_cuspmat4()
    call deallocate_dets()
    call deallocate_dets_equiv()
    call deallocate_distance_mod()
    call deallocate_distances_sav()
    call deallocate_gauss_ecp()
    call deallocate_jd_scratch()
    call deallocate_linear_norm()
    call deallocate_multidet()
    call deallocate_multimat()
    call deallocate_multimatn()
    call deallocate_multislater()
    call deallocate_multislatern()
    call deallocate_orbval()
    call deallocate_phifun()
    call deallocate_qua()
    call deallocate_scratch()
    call deallocate_slater()
    call deallocate_slatn()
    call deallocate_vardep()
    call deallocate_velocity_jastrow()
    call deallocate_wfsec()
    call deallocate_ycompact()
    call deallocate_ycompactn()
    call deallocate_zcompact()
    call deallocate_zmatrix()
    call deallocate_zmatrix_grad()
end subroutine deallocate_m_common
end module
