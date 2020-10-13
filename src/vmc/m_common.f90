module atom
    !> Arguments: znuc, cent, pecent, iwctype, nctype, ncent
    use precision_kinds, only: dp
    use vmc_mod, only: MCENT, MCTYPE

    real(dp), dimension(:, :), allocatable :: cent
    real(dp), dimension(:), allocatable :: znuc
    real(dp) :: pecent
    integer, dimension(:), allocatable :: iwctype
    integer :: nctype, ncent
    integer :: nctype_tot, ncent_tot

    private
    public   :: znuc, cent, pecent, iwctype, nctype, ncent, ncent_tot, nctype_tot
    public :: allocate_atom, deallocate_atom
    save
contains
    subroutine allocate_atom()
        use precision_kinds, only: dp
        use vmc_mod, only: MCENT, MCTYPE
        if (.not. allocated(cent)) allocate (cent(3, MCENT))
        if (.not. allocated(znuc)) allocate (znuc(MCTYPE))
        if (.not. allocated(iwctype)) allocate (iwctype(MCTYPE))
    end subroutine allocate_atom

    subroutine deallocate_atom()
        if (allocated(iwctype)) deallocate (iwctype)
        if (allocated(znuc)) deallocate (znuc)
        if (allocated(cent)) deallocate (cent)
    end subroutine deallocate_atom

end module atom

module ghostatom
    !> Arguments: newghostype, nghostcent

    integer :: newghostype
    integer :: nghostcent

    private
    public   :: newghostype, nghostcent
    save
end module ghostatom

module b_tmove
    !> Arguments: b_t, iskip
    use pseudo_mod, only: MPS_QUAD
    use precision_kinds, only: dp
    use vmc_mod, only: MELEC, MORB, MCENT

    real(dp), dimension(:, :, :, :), allocatable :: b_t !(MORB,MPS_QUAD,MCENT,MELEC)
    integer, dimension(:, :), allocatable :: iskip !(MELEC,MCENT)

    private
    public :: b_t, iskip
    public :: allocate_b_tmove, deallocate_b_tmove
    save
contains
    subroutine allocate_b_tmove()
        use pseudo_mod, only: MPS_QUAD
        use precision_kinds, only: dp
        use vmc_mod, only: MELEC, MORB, MCENT
        if (.not. allocated(b_t)) allocate (b_t(MORB, MPS_QUAD, MCENT, MELEC))
        if (.not. allocated(iskip)) allocate (iskip(MELEC, MCENT))
    end subroutine allocate_b_tmove

    subroutine deallocate_b_tmove()
        if (allocated(iskip)) deallocate (iskip)
        if (allocated(b_t)) deallocate (b_t)
    end subroutine deallocate_b_tmove

end module b_tmove

module Bloc
    !> Arguments: b, tildem, xmat
    use precision_kinds, only: dp
    use vmc_mod, only: MELEC, MORB, MCENT
    use optjas, only: MPARMJ

    real(dp), dimension(:, :), allocatable :: b !(MORB,MELEC)
    real(dp), dimension(:, :, :), allocatable :: tildem !(MELEC,MORB,2)
    real(dp), dimension(:, :), allocatable :: xmat !(MELEC**2,2)

    !> Former Bloc_da
    real(dp), dimension(:, :, :, :), allocatable :: b_da !(3,MELEC,MORB,MCENT)
    real(dp), dimension(:, :, :, :), allocatable :: db !(3,MELEC,MORB,MCENT)

    !> former Bloc_dj
    real(dp), dimension(:, :, :), allocatable :: b_dj !(MORB,MELEC,MPARMJ)

    private
    public :: b, tildem, xmat
    public :: b_da, db
    public :: b_dj
    public :: allocate_Bloc, deallocate_Bloc
    save
contains
    subroutine allocate_Bloc()
        use precision_kinds, only: dp
        use vmc_mod, only: MELEC, MORB, MCENT
        use optjas, only: MPARMJ
        if (.not. allocated(b)) allocate (b(MORB, MELEC))
        if (.not. allocated(tildem)) allocate (tildem(MELEC, MORB, 2))
        if (.not. allocated(xmat)) allocate (xmat(MELEC**2, 2))
        if (.not. allocated(b_da)) allocate (b_da(3, MELEC, MORB, MCENT))
        if (.not. allocated(db)) allocate (db(3, MELEC, MORB, MCENT))
        if (.not. allocated(b_dj)) allocate (b_dj(MORB, MELEC, MPARMJ))
    end subroutine allocate_Bloc

    subroutine deallocate_Bloc()
        if (allocated(b_dj)) deallocate (b_dj)
        if (allocated(db)) deallocate (db)
        if (allocated(b_da)) deallocate (b_da)
        if (allocated(xmat)) deallocate (xmat)
        if (allocated(tildem)) deallocate (tildem)
        if (allocated(b)) deallocate (b)
    end subroutine deallocate_Bloc

end module Bloc

module bparm
    !> Arguments: nocuspb, nspin2b
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
    use vmc_mod, only: MELEC, MCENT

    integer :: i_vpsp
    integer :: icasula
    real(dp), dimension(:, :, :), allocatable :: t_vpsp !(MCENT,MPS_QUAD,MELEC)

    private
    public :: i_vpsp, icasula, t_vpsp
    public :: allocate_casula, deallocate_casula
    save
contains
    subroutine allocate_casula()
        use pseudo_mod, only: MPS_QUAD
        use precision_kinds, only: dp
        use vmc_mod, only: MELEC, MCENT
        if (.not. allocated(t_vpsp)) allocate (t_vpsp(MCENT, MPS_QUAD, MELEC))
    end subroutine allocate_casula

    subroutine deallocate_casula()
        if (allocated(t_vpsp)) deallocate (t_vpsp)
    end subroutine deallocate_casula

end module casula

module chck
    !> Never called
    !> Arguments: bot
    use precision_kinds, only: dp

    real(dp) :: bot

    private
    public :: bot
    save
end module chck

module coefs
    !> need a better name is that the MO ?
    !> if yes can we put it in basis ?
    !> Arguments: coef, nbasis, norb
    use force_mod, only: MWF
    use precision_kinds, only: dp
    use vmc_mod, only: MORB, MBASIS

    real(dp), dimension(:, :, :), allocatable :: coef !(MBASIS,MORB,MWF)
    integer :: nbasis
    integer :: norb

    private
    public :: coef, nbasis, norb
    public :: allocate_coefs, deallocate_coefs
    save
contains
    ! subroutine allocate_coefs()
    !     use force_mod, only: MWF
    !     use precision_kinds, only: dp
    !     use vmc_mod, only: MORB, MBASIS
    !     if (.not. allocated(coef)) allocate (coef(MBASIS, MORB, MWF))
    ! end subroutine allocate_coefs

    subroutine deallocate_coefs()
        if (allocated(coef)) deallocate (coef)
    end subroutine deallocate_coefs

end module coefs

module csfs
    !> Arguments: ccsf, cxdet, iadet, ibdet, icxdet, ncsf, nstates
    use force_mod, only: MWF
    use precision_kinds, only: dp
    use vmc_mod, only: MDET
    use mstates_mod, only: MSTATES, MDETCSFX

    real(dp), dimension(:, :, :), allocatable :: ccsf !(MDET,MSTATES,MWF)
    real(dp), dimension(:), allocatable :: cxdet !(MDET*MDETCSFX)
    integer, dimension(:), allocatable :: iadet !(MDET)
    integer, dimension(:), allocatable :: ibdet !(MDET)
    integer, dimension(:), allocatable :: icxdet !(MDET*MDETCSFX)
    integer :: ncsf
    integer :: nstates

    private
    public   :: ccsf, cxdet, iadet, ibdet, icxdet, ncsf, nstates
    public :: allocate_csfs, deallocate_csfs
    save
contains
    subroutine allocate_csfs()
        use force_mod, only: MWF
        use precision_kinds, only: dp
        use vmc_mod, only: MDET
        use mstates_mod, only: MSTATES, MDETCSFX
        if (.not. allocated(ccsf)) allocate (ccsf(MDET, MSTATES, MWF))
        if (.not. allocated(cxdet)) allocate (cxdet(MDET*MDETCSFX))
        if (.not. allocated(iadet)) allocate (iadet(MDET))
        if (.not. allocated(ibdet)) allocate (ibdet(MDET))
        if (.not. allocated(icxdet)) allocate (icxdet(MDET*MDETCSFX))
    end subroutine allocate_csfs

    subroutine deallocate_csfs()
        if (allocated(icxdet)) deallocate (icxdet)
        if (allocated(ibdet)) deallocate (ibdet)
        if (allocated(iadet)) deallocate (iadet)
        if (allocated(cxdet)) deallocate (cxdet)
        if (allocated(ccsf)) deallocate (ccsf)
    end subroutine deallocate_csfs

end module csfs

module cuspmat
    !> Never called !
    !> Arguments: cm, ishe, iwc3, neqs
    use precision_kinds, only: dp
    use vmc_mod, only: NEQSX

    real(dp), dimension(:, :), allocatable :: cm !(NEQSX,NEQSX)
    integer :: ishe
    integer, dimension(:), allocatable :: iwc3 !(NEQSX)
    integer :: neqs

    private
    public :: cm, ishe, iwc3, neqs
    public :: allocate_cuspmat, deallocate_cuspmat
    save
contains
    subroutine allocate_cuspmat()
        use precision_kinds, only: dp
        use vmc_mod, only: NEQSX
        if (.not. allocated(cm)) allocate (cm(NEQSX, NEQSX))
        if (.not. allocated(iwc3)) allocate (iwc3(NEQSX))
    end subroutine allocate_cuspmat

    subroutine deallocate_cuspmat()
        if (allocated(iwc3)) deallocate (iwc3)
        if (allocated(cm)) deallocate (cm)
    end subroutine deallocate_cuspmat

end module cuspmat

module cuspmat4
    !> Arguments: d, icusp, nterms
    use vmc_mod, only: NEQSX, MTERMS
    use precision_kinds, only: dp

    real(dp), dimension(:, :), allocatable :: d !(NEQSX,MTERMS)
    integer, dimension(:), allocatable :: iwc4 !(NEQSX)
    integer :: nterms
    private

    public :: d, iwc4, nterms
    public :: allocate_cuspmat4, deallocate_cuspmat4
    save
contains
    subroutine allocate_cuspmat4()
        use vmc_mod, only: NEQSX, MTERMS
        use precision_kinds, only: dp
        if (.not. allocated(d)) allocate (d(NEQSX, MTERMS))
        if (.not. allocated(iwc4)) allocate (iwc4(NEQSX))
    end subroutine allocate_cuspmat4

    subroutine deallocate_cuspmat4()
        if (allocated(iwc4)) deallocate (iwc4)
        if (allocated(d)) deallocate (d)
    end subroutine deallocate_cuspmat4

end module cuspmat4

module dets
    !> Arguments: cdet, ndet
    use force_mod, only: MWF
    use precision_kinds, only: dp
    use vmc_mod, only: MDET
    use mstates_mod, only: MSTATES

    real(dp), dimension(:, :, :), allocatable :: cdet !(MDET,MSTATES,MWF)
    integer :: ndet

    private
    public   :: cdet, ndet
    public :: allocate_dets, deallocate_dets
    save
contains
    ! subroutine allocate_dets()
    !     use force_mod, only: MWF
    !     use precision_kinds, only: dp
    !     use vmc_mod, only: MDET
    !     use mstates_mod, only: MSTATES
    !     if (.not. allocated(cdet)) allocate (cdet(MDET, MSTATES, MWF))
    ! end subroutine allocate_dets

    subroutine deallocate_dets()
        if (allocated(cdet)) deallocate (cdet)
    end subroutine deallocate_dets

end module dets

module dets_equiv
    !> Arguments: cdet_equiv, dcdet_equiv
    use precision_kinds, only: dp
    use vmc_mod, only: MDET

    real(dp), dimension(:), allocatable :: cdet_equiv !(MDET)
    real(dp), dimension(:), allocatable :: dcdet_equiv !(MDET)

    private
    public   ::  cdet_equiv, dcdet_equiv
    public :: allocate_dets_equiv, deallocate_dets_equiv
    save
contains
    subroutine allocate_dets_equiv()
        use precision_kinds, only: dp
        use vmc_mod, only: MDET
        if (.not. allocated(cdet_equiv)) allocate (cdet_equiv(MDET))
        if (.not. allocated(dcdet_equiv)) allocate (dcdet_equiv(MDET))
    end subroutine allocate_dets_equiv

    subroutine deallocate_dets_equiv()
        if (allocated(dcdet_equiv)) deallocate (dcdet_equiv)
        if (allocated(cdet_equiv)) deallocate (cdet_equiv)
    end subroutine deallocate_dets_equiv

end module dets_equiv

module distance_mod
    use precision_kinds, only: dp
    use vmc_mod, only: MELEC, MCENT, MMAT_DIM2

    real(dp), dimension(:, :, :), allocatable :: rshift !(3, MELEC, MCENT)
    real(dp), dimension(:, :), allocatable :: r_en !(MELEC, MCENT)
    real(dp), dimension(:, :, :), allocatable :: rvec_en !(3, MELEC, MCENT)
    real(dp), dimension(:), allocatable :: r_ee !(MMAT_DIM2)
    real(dp), dimension(:, :), allocatable :: rvec_ee !(3, MMAT_DIM2)

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
        use precision_kinds, only: dp
        use vmc_mod, only: MELEC, MCENT, MMAT_DIM2
        if (.not. allocated(r_en)) allocate (r_en(MELEC, MCENT))
        if (.not. allocated(rvec_en)) allocate (rvec_en(3, MELEC, MCENT))
        if (.not. allocated(r_ee)) allocate (r_ee(MMAT_DIM2))
        if (.not. allocated(rvec_ee)) allocate (rvec_ee(3, MMAT_DIM2))
        if (.not. allocated(rshift)) allocate (rshift(3, MELEC, MCENT))
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
    use vmc_mod, only: MELEC, MCENT

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
        use precision_kinds, only: dp
        use vmc_mod, only: MELEC, MCENT
        if (.not. allocated(r_ee_sav)) allocate (r_ee_sav(MELEC))
        if (.not. allocated(r_en_sav)) allocate (r_en_sav(MCENT))
        if (.not. allocated(rshift_sav)) allocate (rshift_sav(3, MCENT))
        if (.not. allocated(rvec_ee_sav)) allocate (rvec_ee_sav(3, MELEC))
        if (.not. allocated(rvec_en_sav)) allocate (rvec_en_sav(3, MCENT))
    end subroutine allocate_distances_sav

    subroutine deallocate_distances_sav()
        if (allocated(rvec_en_sav)) deallocate (rvec_en_sav)
        if (allocated(rvec_ee_sav)) deallocate (rvec_ee_sav)
        if (allocated(rshift_sav)) deallocate (rshift_sav)
        if (allocated(r_en_sav)) deallocate (r_en_sav)
        if (allocated(r_ee_sav)) deallocate (r_ee_sav)
    end subroutine deallocate_distances_sav

end module distances_sav

module elec
    !> Arguments: ndn, nup

    integer :: ndn
    integer :: nup

    private
    public   :: ndn, nup
    save
end module elec

module embed
    !> Never called
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
    use vmc_mod, only: MCTYPE

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
        use pseudo_mod, only: MPS_L, MGAUSS
        use precision_kinds, only: dp
        use vmc_mod, only: MCTYPE
        if (.not. allocated(ecp_coef)) allocate (ecp_coef(MGAUSS, MPS_L, MCTYPE))
        if (.not. allocated(ecp_exponent)) allocate (ecp_exponent(MGAUSS, MPS_L, MCTYPE))
        if (.not. allocated(necp_power)) allocate (necp_power(MGAUSS, MPS_L, MCTYPE))
        if (.not. allocated(necp_term)) allocate (necp_term(MPS_L, MCTYPE))
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

    integer :: nbj_current
    integer :: ngrad_jas_bcum
    integer :: ngrad_jas_blocks
    integer :: njb_current

    private
    public :: nbj_current, ngrad_jas_bcum, ngrad_jas_blocks, njb_current
    save
end module gradjerrb

module insout
    !> something to do with grid ...
    !> Arguments: inout, inside

    integer :: inout
    integer :: inside

    private
    public :: inout, inside
    save
end module insout

module jd_scratch
    !> only for (jacobi) davidson in linear method
    !> Arguments: qr, rr
    use sr_mod, only: MPARM
    use precision_kinds, only: dp

    real(dp), dimension(:), allocatable :: qr !(MPARM)
    real(dp), dimension(:), allocatable :: rr !(MPARM)

    private
    public :: qr, rr
    public :: allocate_jd_scratch, deallocate_jd_scratch
    save
contains
    subroutine allocate_jd_scratch()
        use sr_mod, only: MPARM
        use precision_kinds, only: dp
        if (.not. allocated(qr)) allocate (qr(MPARM))
        if (.not. allocated(rr)) allocate (rr(MPARM))
    end subroutine allocate_jd_scratch

    subroutine deallocate_jd_scratch()
        if (allocated(rr)) deallocate (rr)
        if (allocated(qr)) deallocate (qr)
    end subroutine deallocate_jd_scratch

end module jd_scratch

module linear_norm
    !> Arguments: oav, ci_oav
    use optci, only: MXCITERM
    use precision_kinds, only: dp

    real(dp), dimension(:), allocatable :: oav !(MXCITERM)
    real(dp), dimension(:), allocatable :: ci_oav !(MXCITERM)

    private
    public :: oav, ci_oav
    public :: allocate_linear_norm, deallocate_linear_norm
    save
contains
    subroutine allocate_linear_norm()
        use optci, only: MXCITERM
        use precision_kinds, only: dp
        if (.not. allocated(oav)) allocate (oav(MXCITERM))
        if (.not. allocated(ci_oav)) allocate (ci_oav(MXCITERM))
    end subroutine allocate_linear_norm

    subroutine deallocate_linear_norm()
        if (allocated(ci_oav)) deallocate (ci_oav)
        if (allocated(oav)) deallocate (oav)
    end subroutine deallocate_linear_norm

end module linear_norm

module multidet
    !> Arguments: iactv, irepcol_det, ireporb_det, ivirt, iwundet, kref, numrep_det
    use vmc_mod, only: MELEC, MDET

    integer, dimension(:), allocatable :: iactv !(2)
    integer, dimension(:, :, :), allocatable :: irepcol_det !(MELEC,MDET,2)
    integer, dimension(:, :, :), allocatable :: ireporb_det !(MELEC,MDET,2)
    integer, dimension(:), allocatable :: ivirt !(2)
    integer, dimension(:, :), allocatable :: iwundet !(MDET,2)
    integer :: kref
    integer, dimension(:, :), allocatable :: numrep_det !(MDET,2)

    private
    public :: iactv, irepcol_det, ireporb_det, ivirt, iwundet, kref, numrep_det
    public :: allocate_multidet, deallocate_multidet
    save
contains
    subroutine allocate_multidet()
        use vmc_mod, only: MELEC, MDET
        if (.not. allocated(iactv)) allocate (iactv(2))
        ! if (.not. allocated(irepcol_det)) allocate (irepcol_det(MELEC, MDET, 2))
        ! if (.not. allocated(ireporb_det)) allocate (ireporb_det(MELEC, MDET, 2))
        if (.not. allocated(ivirt)) allocate (ivirt(2))
        ! if (.not. allocated(iwundet)) allocate (iwundet(MDET, 2))
        ! if (.not. allocated(numrep_det)) allocate (numrep_det(MDET, 2))
    end subroutine allocate_multidet

    subroutine deallocate_multidet()
        if (allocated(numrep_det)) deallocate (numrep_det)
        if (allocated(iwundet)) deallocate (iwundet)
        if (allocated(ivirt)) deallocate (ivirt)
        if (allocated(ireporb_det)) deallocate (ireporb_det)
        if (allocated(irepcol_det)) deallocate (irepcol_det)
        if (allocated(iactv)) deallocate (iactv)
    end subroutine deallocate_multidet

end module multidet

module multimat
    !> Arguments: aa, wfmat
    use precision_kinds, only: dp
    use vmc_mod, only: MELEC, MORB, MDET
    use vmc_mod, only: MEXCIT

    real(dp), dimension(:, :, :), allocatable :: aa !(MELEC,MORB,2)
    real(dp), dimension(:, :, :), allocatable :: wfmat !(MEXCIT**2,MDET,2)

    private
    public :: aa, wfmat
    public :: allocate_multimat, deallocate_multimat
    save
contains
    subroutine allocate_multimat()
        use precision_kinds, only: dp
        use vmc_mod, only: MELEC, MORB, MDET
        use vmc_mod, only: MEXCIT
        if (.not. allocated(aa)) allocate (aa(MELEC, MORB, 2))
        if (.not. allocated(wfmat)) allocate (wfmat(MEXCIT**2, MDET, 2))
    end subroutine allocate_multimat

    subroutine deallocate_multimat()
        if (allocated(wfmat)) deallocate (wfmat)
        if (allocated(aa)) deallocate (aa)
    end subroutine deallocate_multimat

end module multimat

module multimatn
    !> Arguments: aan, wfmatn
    use precision_kinds, only: dp
    use vmc_mod, only: MELEC, MORB, MDET
    use vmc_mod, only: MEXCIT

    real(dp), dimension(:, :), allocatable :: aan !(MELEC,MORB)
    real(dp), dimension(:, :), allocatable :: wfmatn !(MEXCIT**2,MDET)

    private
    public :: aan, wfmatn
    public :: allocate_multimatn, deallocate_multimatn
    save
contains
    subroutine allocate_multimatn()
        use precision_kinds, only: dp
        use vmc_mod, only: MELEC, MORB, MDET
        use vmc_mod, only: MEXCIT
        if (.not. allocated(aan)) allocate (aan(MELEC, MORB))
        if (.not. allocated(wfmatn)) allocate (wfmatn(MEXCIT**2, MDET))
    end subroutine allocate_multimatn

    subroutine deallocate_multimatn()
        if (allocated(wfmatn)) deallocate (wfmatn)
        if (allocated(aan)) deallocate (aan)
    end subroutine deallocate_multimatn

end module multimatn

module multislater
    !> Arguments: detiab
    use precision_kinds, only: dp
    use vmc_mod, only: MDET

    real(dp), dimension(:, :), allocatable :: detiab !(MDET,2)

    private
    public :: detiab
    public :: allocate_multislater, deallocate_multislater
    save
contains
    subroutine allocate_multislater()
        use precision_kinds, only: dp
        use vmc_mod, only: MDET
        if (.not. allocated(detiab)) allocate (detiab(MDET, 2))
    end subroutine allocate_multislater

    subroutine deallocate_multislater()
        if (allocated(detiab)) deallocate (detiab)
    end subroutine deallocate_multislater

end module multislater

module multislatern
    !> Arguments: ddorbn, detn, dorbn, orbn

    use precision_kinds, only: dp
    use vmc_mod, only: MORB, MDET

    real(dp), dimension(:), allocatable :: ddorbn !(MORB)
    real(dp), dimension(:), allocatable :: detn !(MDET)
    real(dp), dimension(:, :), allocatable :: dorbn !(3,MORB)
    real(dp), dimension(:), allocatable :: orbn !(MORB)
    private

    public ::  ddorbn, detn, dorbn, orbn
    public :: allocate_multislatern, deallocate_multislatern
    save
contains
    subroutine allocate_multislatern()
        use precision_kinds, only: dp
        use vmc_mod, only: MORB, MDET
        if (.not. allocated(ddorbn)) allocate (ddorbn(MORB))
        if (.not. allocated(detn)) allocate (detn(MDET))
        if (.not. allocated(dorbn)) allocate (dorbn(3, MORB))
        if (.not. allocated(orbn)) allocate (orbn(MORB))
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
    use vmc_mod, only: MELEC, MORB

    real(dp), dimension(:, :), allocatable :: ddorb !(MELEC,MORB)
    real(dp), dimension(:, :, :), allocatable :: dorb !(3,MELEC,MORB)
    integer :: nadorb
    integer :: ndetorb
    real(dp), dimension(:, :), allocatable :: orb !(MELEC,MORB)

    private
    public :: ddorb, dorb, nadorb, ndetorb, orb
    public :: allocate_orbval, deallocate_orbval
    save
contains
    subroutine allocate_orbval()
        use precision_kinds, only: dp
        use vmc_mod, only: MELEC, MORB
        if (.not. allocated(ddorb)) allocate (ddorb(MELEC, MORB))
        if (.not. allocated(dorb)) allocate (dorb(3, MELEC, MORB))
        if (.not. allocated(orb)) allocate (orb(MELEC, MORB))
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
    use vmc_mod, only: MELEC, MBASIS

    real(dp), dimension(:, :), allocatable :: d2phin !(MBASIS,MELEC)
    real(dp), dimension(:, :, :, :), allocatable :: d2phin_all !(3,3,MBASIS,MELEC)
    real(dp), dimension(:, :, :), allocatable :: d3phin !(3,MBASIS,MELEC)
    real(dp), dimension(:, :, :), allocatable :: dphin !(3,MBASIS,MELEC)
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
        use precision_kinds, only: dp
        use vmc_mod, only: MELEC, MBASIS
        if (.not. allocated(d2phin)) allocate (d2phin(MBASIS, MELEC))
        if (.not. allocated(d2phin_all)) allocate (d2phin_all(3, 3, MBASIS, MELEC))
        if (.not. allocated(d3phin)) allocate (d3phin(3, MBASIS, MELEC))
        if (.not. allocated(dphin)) allocate (dphin(3, MBASIS, MELEC))
        if (.not. allocated(n0_ibasis)) allocate (n0_ibasis(MBASIS, MELEC))
        if (.not. allocated(n0_ic)) allocate (n0_ic(MBASIS, MELEC))
        if (.not. allocated(n0_nbasis)) allocate (n0_nbasis(MELEC))
        if (.not. allocated(phin)) allocate (phin(MBASIS, MELEC))
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
        use precision_kinds, only: dp
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

module rlobxy
    !> only rlobx used in read input
    !> Arguments: rlobx, rloby, rloby2
    use precision_kinds, only: dp
    use vmc_mod, only: NSPLIN

    real(dp), dimension(:), allocatable :: rlobx !(NSPLIN)
    real(dp), dimension(:), allocatable :: rloby !(NSPLIN)
    real(dp), dimension(:), allocatable :: rloby2 !(NSPLIN)

    private
    public :: rlobx, rloby, rloby2
    public :: allocate_rlobxy, deallocate_rlobxy
    save
contains
    subroutine allocate_rlobxy()
        use precision_kinds, only: dp
        use vmc_mod, only: NSPLIN
        if (.not. allocated(rlobx)) allocate (rlobx(NSPLIN))
        if (.not. allocated(rloby)) allocate (rloby(NSPLIN))
        if (.not. allocated(rloby2)) allocate (rloby2(NSPLIN))
    end subroutine allocate_rlobxy

    subroutine deallocate_rlobxy()
        if (allocated(rloby2)) deallocate (rloby2)
        if (allocated(rloby)) deallocate (rloby)
        if (allocated(rlobx)) deallocate (rlobx)
    end subroutine deallocate_rlobxy

end module rlobxy

module scale_more
    !> Arguments: dd3
    use precision_kinds, only: dp

    real(dp) :: dd3

    private
    public :: dd3
    save
end module scale_more

module scratch
    !> Arguments: denergy_det, dtildem
    use precision_kinds, only: dp
    use vmc_mod, only: MELEC, MORB, MDET

    real(dp), dimension(:, :), allocatable :: denergy_det !(MDET,2)
    real(dp), dimension(:, :, :), allocatable :: dtildem !(MELEC,MORB,2)

    private
    public :: denergy_det, dtildem
    public :: allocate_scratch, deallocate_scratch
    save
contains
    subroutine allocate_scratch()
        use precision_kinds, only: dp
        use vmc_mod, only: MELEC, MORB, MDET
        if (.not. allocated(denergy_det)) allocate (denergy_det(MDET, 2))
        if (.not. allocated(dtildem)) allocate (dtildem(MELEC, MORB, 2))
    end subroutine allocate_scratch

    subroutine deallocate_scratch()
        if (allocated(dtildem)) deallocate (dtildem)
        if (allocated(denergy_det)) deallocate (denergy_det)
    end subroutine deallocate_scratch

end module scratch

module slater
    !> Arguments: d2dx2, ddx, fp, fpp, slmi

    use precision_kinds, only: dp
    use vmc_mod, only: MELEC, MMAT_DIM

    real(dp), dimension(:), allocatable :: d2dx2 !(MELEC)
    real(dp), dimension(:, :), allocatable :: ddx !(3,MELEC)
    real(dp), dimension(:, :, :), allocatable :: fp !(3,MMAT_DIM,2)
    real(dp), dimension(:, :), allocatable :: fpp !(MMAT_DIM,2)
    real(dp), dimension(:, :), allocatable :: slmi !(MMAT_DIM,2)

    private
    public :: d2dx2, ddx, fp, fpp, slmi
    public :: allocate_slater, deallocate_slater
    save
contains
    subroutine allocate_slater()
        use precision_kinds, only: dp
        use vmc_mod, only: MELEC, MMAT_DIM
        if (.not. allocated(d2dx2)) allocate (d2dx2(MELEC))
        if (.not. allocated(ddx)) allocate (ddx(3, MELEC))
        if (.not. allocated(fp)) allocate (fp(3, MMAT_DIM, 2))
        if (.not. allocated(fpp)) allocate (fpp(MMAT_DIM, 2))
        if (.not. allocated(slmi)) allocate (slmi(MMAT_DIM, 2))
    end subroutine allocate_slater

    subroutine deallocate_slater()
        if (allocated(slmi)) deallocate (slmi)
        if (allocated(fpp)) deallocate (fpp)
        if (allocated(fp)) deallocate (fp)
        if (allocated(ddx)) deallocate (ddx)
        if (allocated(d2dx2)) deallocate (d2dx2)
    end subroutine deallocate_slater

end module slater

module slatn
    !> Arguments: slmin
    use precision_kinds, only: dp
    use vmc_mod, only: MMAT_DIM

    real(dp), dimension(:), allocatable :: slmin !(MMAT_DIM)

    private
    public :: slmin
    public :: allocate_slatn, deallocate_slatn
    save
contains
    subroutine allocate_slatn()
        use precision_kinds, only: dp
        use vmc_mod, only: MMAT_DIM
        if (.not. allocated(slmin)) allocate (slmin(MMAT_DIM))
    end subroutine allocate_slatn

    subroutine deallocate_slatn()
        if (allocated(slmin)) deallocate (slmin)
    end subroutine deallocate_slatn

end module slatn

module svd_mod
    ! Not used anywhere !
    !> Arguments:
    integer, parameter :: MBUF = 10000
    integer, parameter :: MXDIM = 3000
    private
    public :: MBUF, MXDIM
    save
end module svd_mod

module vardep
    !> Arguments: cdep, iwdepend, nvdepend
    use precision_kinds, only: dp
    use vmc_mod, only: MCTYPE
    use vmc_mod, only: NEQSX

    real(dp), dimension(:, :, :), allocatable :: cdep !(NEQSX,83,MCTYPE)
    integer, dimension(:, :, :), allocatable :: iwdepend !(NEQSX,83,MCTYPE)
    integer, dimension(:, :), allocatable :: nvdepend !(NEQSX,MCTYPE)

    private
    public :: cdep, iwdepend, nvdepend
    public :: allocate_vardep, deallocate_vardep
    save
contains
    subroutine allocate_vardep()
        use precision_kinds, only: dp
        use vmc_mod, only: MCTYPE
        use vmc_mod, only: NEQSX
        if (.not. allocated(cdep)) allocate (cdep(NEQSX, 83, MCTYPE))
        if (.not. allocated(iwdepend)) allocate (iwdepend(NEQSX, 83, MCTYPE))
        if (.not. allocated(nvdepend)) allocate (nvdepend(NEQSX, MCTYPE))
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
    use vmc_mod, only: MELEC

    real(dp), dimension(:, :), allocatable :: vj !(3,MELEC)
    real(dp), dimension(:, :), allocatable :: vjn !(3,MELEC)

    private
    public :: vj, vjn
    public :: allocate_velocity_jastrow, deallocate_velocity_jastrow
    save
contains
    subroutine allocate_velocity_jastrow()
        use precision_kinds, only: dp
        use vmc_mod, only: MELEC
        if (.not. allocated(vj)) allocate (vj(3, MELEC))
        if (.not. allocated(vjn)) allocate (vjn(3, MELEC))
    end subroutine allocate_velocity_jastrow

    subroutine deallocate_velocity_jastrow()
        if (allocated(vjn)) deallocate (vjn)
        if (allocated(vj)) deallocate (vj)
    end subroutine deallocate_velocity_jastrow

end module velocity_jastrow

module wfsec
    use force_mod, only: MFORCE
    !> Arguments: iwf, iwftype, nwftype

    integer :: iwf
    integer, dimension(:), allocatable :: iwftype !(MFORCE)
    integer :: nwftype

    private
    public :: iwf, iwftype, nwftype
    public :: allocate_wfsec, deallocate_wfsec
    save
contains
    ! subroutine allocate_wfsec()
    !     use force_mod, only: MFORCE
    !     if (.not. allocated(iwftype)) allocate (iwftype(MFORCE))
    ! end subroutine allocate_wfsec

    subroutine deallocate_wfsec()
        if (allocated(iwftype)) deallocate (iwftype)
    end subroutine deallocate_wfsec

end module wfsec

module ycompact
    !> Arguments: dymat, ymat
    use precision_kinds, only: dp
    use vmc_mod, only: MELEC, MORB
    use mstates_mod, only: MSTATES

    real(dp), dimension(:, :, :, :), allocatable :: dymat !(MORB,MELEC,2,MSTATES)
    real(dp), dimension(:, :, :, :), allocatable :: ymat !(MORB,MELEC,2,MSTATES)

    private
    public :: dymat, ymat
    public :: allocate_ycompact, deallocate_ycompact
    save
contains
    subroutine allocate_ycompact()
        use precision_kinds, only: dp
        use vmc_mod, only: MELEC, MORB
        use mstates_mod, only: MSTATES
        if (.not. allocated(dymat)) allocate (dymat(MORB, MELEC, 2, MSTATES))
        if (.not. allocated(ymat)) allocate (ymat(MORB, MELEC, 2, MSTATES))
    end subroutine allocate_ycompact

    subroutine deallocate_ycompact()
        if (allocated(ymat)) deallocate (ymat)
        if (allocated(dymat)) deallocate (dymat)
    end subroutine deallocate_ycompact

end module ycompact

module ycompactn
    !> Arguments: ymatn
    use precision_kinds, only: dp
    use vmc_mod, only: MELEC, MORB
    use mstates_mod, only: MSTATES

    real(dp), dimension(:, :, :), allocatable :: ymatn !(MORB,MELEC,MSTATES)

    private
    public :: ymatn
    public :: allocate_ycompactn, deallocate_ycompactn
    save
contains
    subroutine allocate_ycompactn()
        use precision_kinds, only: dp
        use vmc_mod, only: MELEC, MORB
        use mstates_mod, only: MSTATES
        if (.not. allocated(ymatn)) allocate (ymatn(MORB, MELEC, MSTATES))
    end subroutine allocate_ycompactn

    subroutine deallocate_ycompactn()
        if (allocated(ymatn)) deallocate (ymatn)
    end subroutine deallocate_ycompactn

end module ycompactn

module zcompact
    !> Arguments: aaz, dzmat, emz, zmat
    use precision_kinds, only: dp
    use vmc_mod, only: MELEC, MORB
    use mstates_mod, only: MSTATES

    real(dp), dimension(:, :, :, :), allocatable :: aaz !(MELEC,MELEC,2,MSTATES)
    real(dp), dimension(:, :, :, :), allocatable :: dzmat !(MORB,MELEC,2,MSTATES)
    real(dp), dimension(:, :, :, :), allocatable :: emz !(MELEC,MELEC,2,MSTATES)
    real(dp), dimension(:, :, :, :), allocatable :: zmat !(MORB,MELEC,2,MSTATES)

    private
    public :: aaz, dzmat, emz, zmat
    public :: allocate_zcompact, deallocate_zcompact
    save
contains
    subroutine allocate_zcompact()
        use precision_kinds, only: dp
        use vmc_mod, only: MELEC, MORB
        use mstates_mod, only: MSTATES
        if (.not. allocated(aaz)) allocate (aaz(MELEC, MELEC, 2, MSTATES))
        if (.not. allocated(dzmat)) allocate (dzmat(MORB, MELEC, 2, MSTATES))
        if (.not. allocated(emz)) allocate (emz(MELEC, MELEC, 2, MSTATES))
        if (.not. allocated(zmat)) allocate (zmat(MORB, MELEC, 2, MSTATES))
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
    use vmc_mod, only: MCENT

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
        use precision_kinds, only: dp
        use vmc_mod, only: MCENT
        if (.not. allocated(czcart)) allocate (czcart(3, MCENT))
        if (.not. allocated(czint)) allocate (czint(3, MCENT))
        if (.not. allocated(czcart_ref)) allocate (czcart_ref(3, 3))
        if (.not. allocated(izcmat)) allocate (izcmat(3, MCENT))
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
    use vmc_mod, only: MCENT3

    real(dp), dimension(:, :), allocatable :: transform_grd !(MCENT3,MCENT3)

    private
    public :: transform_grd
    public :: allocate_zmatrix_grad, deallocate_zmatrix_grad
    save
contains
    subroutine allocate_zmatrix_grad()
        use precision_kinds, only: dp
        use vmc_mod, only: MCENT3
        if (.not. allocated(transform_grd)) allocate (transform_grd(MCENT3, MCENT3))
    end subroutine allocate_zmatrix_grad

    subroutine deallocate_zmatrix_grad()
        if (allocated(transform_grd)) deallocate (transform_grd)
    end subroutine deallocate_zmatrix_grad

end module zmatrix_grad

subroutine allocate_m_common()
    use atom, only: allocate_atom
    use b_tmove, only: allocate_b_tmove
    use Bloc, only: allocate_Bloc
    use casula, only: allocate_casula
    ! use coefs, only: allocate_coefs
    use csfs, only: allocate_csfs
    use cuspmat, only: allocate_cuspmat
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
    use rlobxy, only: allocate_rlobxy
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

    call allocate_atom()
    call allocate_b_tmove()
    call allocate_Bloc()
    call allocate_casula()
    ! call allocate_coefs()
    call allocate_csfs()
    call allocate_cuspmat()
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
    call allocate_rlobxy()
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
