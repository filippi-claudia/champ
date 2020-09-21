!> \brief File collecting all modules that replace common blocks.
!>
!> \author P. Lopez-Tarifa & F. Zapata NLeSC(2019)

 module atom
   !> Arguments: znuc, cent, pecent, iwctype, nctype, ncent
   use precision_kinds, only: dp
   use vmc, only: MCENT, MCTYPE

   real(dp) :: cent( 3, MCENT)
   real(dp) :: znuc( MCTYPE)
   real(dp) :: pecent
   integer  :: iwctype( MCENT), nctype, ncent

   private
   public   :: znuc, cent, pecent, iwctype, nctype, ncent
   save
 end module atom

 module ghostatom
   !> Arguments: newghostype, nghostcent

   integer  :: newghostype
   integer  :: nghostcent

   private
   public   :: newghostype, nghostcent
   save
 end module ghostatom

 module b_tmove
   !> Arguments: b_t, iskip
   use pseudo_mod, only: MPS_QUAD
   use precision_kinds, only: dp
   use vmc, only: MELEC, MORB, MCENT

   real(dp) :: b_t(MORB,MPS_QUAD,MCENT,MELEC)
   integer  :: iskip(MELEC,MCENT)

   private
   public :: b_t, iskip
   save
 end module b_tmove

 module Bloc
   !> Arguments: b, tildem, xmat
   use precision_kinds, only: dp
   use vmc, only: MELEC, MORB, MCENT
   use optjas, only: MPARMJ

   real(dp) :: b(MORB,MELEC)
   real(dp) :: tildem(MELEC,MORB,2)
   real(dp) :: xmat(MELEC**2,2)

   !> Former Bloc_da
   real(dp) :: b_da(3,MELEC,MORB,MCENT)
   real(dp) :: db(3,MELEC,MORB,MCENT)

   !> former Bloc_dj
   real(dp) :: b_dj(MORB,MELEC,MPARMJ)

   private
   public :: b, tildem, xmat
   public :: b_da, db
   public :: b_dj
   save
 end module Bloc

 module bparm
   !> Arguments: nocuspb, nspin2b
   integer  :: nocuspb
   integer  :: nspin2b

   private
   public :: nocuspb, nspin2b
   save
 end module bparm

 module casula
   !> Arguments: i_vpsp, icasula, t_vpsp
   use pseudo_mod, only: MPS_QUAD
   use precision_kinds, only: dp
   use vmc, only: MELEC, MCENT

   integer  :: i_vpsp
   integer  :: icasula
   real(dp) :: t_vpsp(MCENT,MPS_QUAD,MELEC)

   private
   public :: i_vpsp, icasula, t_vpsp
   save
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
   use vmc, only: MORB, MBASIS

   real(dp) :: coef(MBASIS,MORB,MWF)
   integer  :: nbasis
   integer  :: norb

   private
   public :: coef, nbasis, norb
   save
 end module coefs

 module csfs
   !> Arguments: ccsf, cxdet, iadet, ibdet, icxdet, ncsf, nstates
   use force_mod, only: MWF
   use precision_kinds, only: dp
   use vmc, only: MDET
   use mstates_mod, only: MSTATES, MDETCSFX

   real(dp) :: ccsf(MDET,MSTATES,MWF)
   real(dp) :: cxdet(MDET*MDETCSFX)
   integer  :: iadet(MDET)
   integer  :: ibdet(MDET)
   integer  :: icxdet(MDET*MDETCSFX)
   integer  :: ncsf
   integer  :: nstates

   private
   public   :: ccsf, cxdet, iadet, ibdet, icxdet, ncsf, nstates
   save
 end module csfs

 module cuspmat
   !> Never called !
   !> Arguments: cm, ishe, iwc3, neqs
   use precision_kinds, only: dp
   use vmc, only: NEQSX

   real(dp) :: cm(NEQSX,NEQSX)
   integer  :: ishe
   integer  :: iwc3(NEQSX)
   integer  :: neqs

   private
   public :: cm, ishe, iwc3, neqs
   save
 end module cuspmat

 module cuspmat4
   !> Arguments: d, icusp, nterms
   use vmc, only: NEQSX, MTERMS
   use precision_kinds, only: dp
   

    real(dp) :: d(NEQSX,MTERMS)
    integer  :: iwc4(NEQSX)
    integer  :: nterms
    private

    public :: d, iwc4, nterms
    save
 end module cuspmat4

 module dets
   !> Arguments: cdet, ndet
   use force_mod, only: MWF
   use precision_kinds, only: dp
   use vmc, only: MDET
   use mstates_mod, only: MSTATES

   real(dp) :: cdet(MDET,MSTATES,MWF)
   integer  :: ndet

   private
   public   :: cdet, ndet
   save
 end module dets

 module dets_equiv
  !> Arguments: cdet_equiv, dcdet_equiv
   use precision_kinds, only: dp
   use vmc, only: MDET

   real(dp) :: cdet_equiv(MDET)
   real(dp) :: dcdet_equiv(MDET)

   private
   public   ::  cdet_equiv, dcdet_equiv
   save
 end module dets_equiv

 module distances_sav
  !> Arguments: r_ee_sav, r_en_sav, rshift_sav, rvec_ee_sav, rvec_en_sav
  use precision_kinds, only: dp
  use vmc, only: MELEC, MCENT

  real(dp) :: r_ee_sav(MELEC)
  real(dp) :: r_en_sav(MCENT)
  real(dp) :: rshift_sav(3,MCENT)
  real(dp) :: rvec_ee_sav(3,MELEC)
  real(dp) :: rvec_en_sav(3,MCENT)

  private
  public   :: r_ee_sav, r_en_sav, rshift_sav, rvec_ee_sav, rvec_en_sav
  save
 end module distances_sav

 module elec
  !> Arguments: ndn, nup

  integer  :: ndn
  integer  :: nup

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
  integer, parameter :: MXHTAB=1
  integer, parameter :: MXXTAB=1
  integer, parameter :: MXTAB=MXHTAB+MXXTAB
  ! max number of gridpoints in hartree table
  integer, parameter :: MXHGRID=1
  ! max number of gridpoints in exchange tables
  integer, parameter :: MXXGRID=1

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
   use vmc, only: MCTYPE
 
   real(dp) :: ecp_coef(MGAUSS,MPS_L,MCTYPE)
   real(dp) :: ecp_exponent(MGAUSS,MPS_L,MCTYPE)
   integer  :: necp_power(MGAUSS,MPS_L,MCTYPE)
   integer  :: necp_term(MPS_L,MCTYPE)
 
   private
   public   ::  ecp_coef, ecp_exponent, necp_power, necp_term
   save
 end module gauss_ecp

 module gradjerrb
  !> Arguments: nbj_current, ngrad_jas_bcum, ngrad_jas_blocks, njb_current

   integer  :: nbj_current
   integer  :: ngrad_jas_bcum
   integer  :: ngrad_jas_blocks
   integer  :: njb_current

   private
   public :: nbj_current, ngrad_jas_bcum, ngrad_jas_blocks, njb_current
   save
 end module gradjerrb

 module insout
   !> something to do with grid ...
   !> Arguments: inout, inside

   integer  :: inout
   integer  :: inside

   private 
   public :: inout, inside 
   save
 end module insout

 module jd_scratch
   !> only for (jacobi) davidson in linear method
   !> Arguments: qr, rr
   use sr_mod, only: MPARM
   use precision_kinds, only: dp

   real(dp) :: qr(MPARM)
   real(dp) :: rr(MPARM)

   private 
   public :: qr, rr 
   save
 end module jd_scratch

 module linear_norm
   !> Arguments: oav, ci_oav
   use optci, only: MXCITERM
   use precision_kinds, only: dp

   real(dp) :: oav(MXCITERM)
   real(dp) :: ci_oav(MXCITERM) 

   private 
   public :: oav, ci_oav
   save
 end module linear_norm

 module multidet
   !> Arguments: iactv, irepcol_det, ireporb_det, ivirt, iwundet, kref, numrep_det
   use vmc, only: MELEC, MDET

   integer  :: iactv(2)
   integer  :: irepcol_det(MELEC,MDET,2)
   integer  :: ireporb_det(MELEC,MDET,2)
   integer  :: ivirt(2)
   integer  :: iwundet(MDET,2)
   integer  :: kref
   integer  :: numrep_det(MDET,2)

   private
   public :: iactv, irepcol_det, ireporb_det, ivirt, iwundet, kref, numrep_det
   save
 end module multidet

 module multimat
   !> Arguments: aa, wfmat
   use precision_kinds, only: dp
   use vmc, only: MELEC, MORB, MDET
   use vmc, only: MEXCIT

   real(dp) :: aa(MELEC,MORB,2)
   real(dp) :: wfmat(MEXCIT**2,MDET,2)

   private
   public :: aa, wfmat
   save
 end module multimat

 module multimatn
   !> Arguments: aan, wfmatn
   use precision_kinds, only: dp
   use vmc, only: MELEC, MORB, MDET
   use vmc, only: MEXCIT

   real(dp) :: aan(MELEC,MORB)
   real(dp) :: wfmatn(MEXCIT**2,MDET)

   private 
   public :: aan, wfmatn 
   save
 end module multimatn

 module multislater
   !> Arguments: detiab
   use precision_kinds, only: dp
   use vmc, only: MDET 

   real(dp) :: detiab(MDET,2)

   private
   public :: detiab
   save
 end module multislater

  module multislatern
    !> Arguments: ddorbn, detn, dorbn, orbn
 
    use precision_kinds, only: dp
    use vmc, only: MORB, MDET 
 
     real(dp) :: ddorbn(MORB)
     real(dp) :: detn(MDET)
     real(dp) :: dorbn(3,MORB)
     real(dp) :: orbn(MORB)
     private
 
     public ::  ddorbn, detn, dorbn, orbn
     save
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

   integer  :: ncnstr
   integer  :: ncuspc
   integer  :: nfock
   integer  :: nfockc
   integer  :: norbc

   private
   public :: ncnstr, ncuspc, nfock, nfockc, norbc
   save
 end module ncusp

 module orbval
   !> Arguments: ddorb, dorb, nadorb, ndetorb, orb
   use precision_kinds, only: dp
   use vmc, only: MELEC, MORB
 
   real(dp) :: ddorb(MELEC,MORB)
   real(dp) :: dorb(3,MELEC,MORB)
   integer  :: nadorb
   integer  :: ndetorb
   real(dp) :: orb(MELEC,MORB)

   private
   public :: ddorb, dorb, nadorb, ndetorb, orb
   save
 end module orbval


 module phifun
   !> Arguments: d2phin, d2phin_all, d3phin, dphin, n0_ibasis, n0_ic, n0_nbasis, phin
   use precision_kinds, only: dp
   use vmc, only: MELEC, MBASIS

   real(dp) :: d2phin(MBASIS,MELEC)
   real(dp) :: d2phin_all(3,3,MBASIS,MELEC)
   real(dp) :: d3phin(3,MBASIS,MELEC)
   real(dp) :: dphin(3,MBASIS,MELEC)
   integer  :: n0_ibasis(MBASIS,MELEC)
   integer  :: n0_ic(MBASIS,MELEC)
   integer  :: n0_nbasis(MELEC)
   real(dp) :: phin(MBASIS,MELEC)

   private
   public :: d2phin, d2phin_all, d3phin, dphin, n0_ibasis, n0_ic, n0_nbasis, phin
   save
 end module phifun

 module qua
   !> Arguments: nquad, wq, xq, xq0, yq, yq0, zq, zq0
   use pseudo_mod, only: MPS_QUAD
   use precision_kinds, only: dp

   integer  :: nquad
   real(dp) :: wq(MPS_QUAD)
   real(dp) :: xq(MPS_QUAD)
   real(dp) :: xq0(MPS_QUAD)
   real(dp) :: yq(MPS_QUAD)
   real(dp) :: yq0(MPS_QUAD)
   real(dp) :: zq(MPS_QUAD)
   real(dp) :: zq0(MPS_QUAD)

   private
   public :: nquad, wq, xq, xq0, yq, yq0, zq, zq0
   save
 end module qua

 module rlobxy
   !> only rlobx used in read input
   !> Arguments: rlobx, rloby, rloby2
   use precision_kinds, only: dp
   use vmc, only: NSPLIN

    real(dp) :: rlobx(NSPLIN)
    real(dp) :: rloby(NSPLIN)
    real(dp) :: rloby2(NSPLIN)

    private
    public :: rlobx, rloby, rloby2
    save
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
   use vmc, only: MELEC, MORB, MDET

   real(dp) :: denergy_det(MDET,2)
   real(dp) :: dtildem(MELEC,MORB,2)

   private
   public :: denergy_det, dtildem
   save
 end module scratch

 module slater
   !> Arguments: d2dx2, ddx, fp, fpp, slmi
 
   use precision_kinds, only: dp
   use vmc, only: MELEC, MMAT_DIM
 
   real(dp) :: d2dx2(MELEC)
   real(dp) :: ddx(3,MELEC)
   real(dp) :: fp(3,MMAT_DIM,2)
   real(dp) :: fpp(MMAT_DIM,2)
   real(dp) :: slmi(MMAT_DIM,2)

   private
   public :: d2dx2, ddx, fp, fpp, slmi
   save
 end module slater

 module slatn
   !> Arguments: slmin
   use precision_kinds, only: dp
   use vmc, only: MMAT_DIM

   real(dp) :: slmin(MMAT_DIM)

   private
   public :: slmin
   save
 end module slatn

 module svd_mod
  ! Not used anywhere !
  !> Arguments:
  integer, parameter :: MBUF=10000
  integer, parameter :: MXDIM=3000
  private 
  public :: MBUF, MXDIM 
  save 
 end module svd_mod 

 
 module vardep
   !> Arguments: cdep, iwdepend, nvdepend
   use precision_kinds, only: dp
   use vmc, only: MCTYPE
   use vmc, only: NEQSX

   real(dp) :: cdep(NEQSX,83,MCTYPE)
   integer  :: iwdepend(NEQSX,83,MCTYPE)
   integer  :: nvdepend(NEQSX,MCTYPE)

   private 
   public :: cdep, iwdepend, nvdepend 
   save
 end module vardep

 module velocity_jastrow
   !> Arguments: vj, vjn
   use precision_kinds, only: dp
   use vmc, only: MELEC

   real(dp) :: vj(3,MELEC)
   real(dp) :: vjn(3,MELEC)

   private 
   public :: vj, vjn 
   save
 end module velocity_jastrow
 
 module wfsec
   use force_mod, only: MFORCE
   !> Arguments: iwf, iwftype, nwftype

   integer  :: iwf
   integer  :: iwftype(MFORCE)
   integer  :: nwftype

   private
   public :: iwf, iwftype, nwftype
   save
 end module wfsec

 module ycompact
   !> Arguments: dymat, ymat
   use precision_kinds, only: dp
   use vmc, only: MELEC, MORB
   use mstates_mod, only: MSTATES

   real(dp) :: dymat(MORB,MELEC,2,MSTATES)
   real(dp) :: ymat(MORB,MELEC,2,MSTATES)

   private
   public :: dymat, ymat
   save
 end module ycompact

 module ycompactn
   !> Arguments: ymatn
   use precision_kinds, only: dp
   use vmc, only: MELEC, MORB
   use mstates_mod, only: MSTATES

   real(dp) :: ymatn(MORB,MELEC,MSTATES)

   private
   public :: ymatn
   save
 end module ycompactn

 module zcompact
   !> Arguments: aaz, dzmat, emz, zmat
   use precision_kinds, only: dp
   use vmc, only: MELEC, MORB
   use mstates_mod, only: MSTATES

   real(dp) :: aaz(MELEC,MELEC,2,MSTATES)
   real(dp) :: dzmat(MORB,MELEC,2,MSTATES)
   real(dp) :: emz(MELEC,MELEC,2,MSTATES)
   real(dp) :: zmat(MORB,MELEC,2,MSTATES)

   private
   public :: aaz, dzmat, emz, zmat
   save
 end module zcompact

 module zmatrix
   !> Arguments: czcart, czint, czcart_ref, izcmat, izmatrix  
   use precision_kinds, only: dp
   use vmc, only: MCENT

   real(dp) :: czcart(3,MCENT)
   real(dp) :: czint(3,MCENT)
   real(dp) :: czcart_ref(3,3)
   integer  :: izcmat(3,MCENT)
   integer  :: izmatrix

   private
   public :: czcart, czint, czcart_ref, izcmat, izmatrix 
   save
 end module zmatrix
 
 module zmatrix_grad
   !> never called
   !> Arguments: transform_grd
   use precision_kinds, only: dp
   use vmc, only: MCENT3

   real(dp) :: transform_grd(MCENT3,MCENT3)

   private 
   public :: transform_grd 
   save
 end module zmatrix_grad
