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

 module basis
   !> Arguments: zex, betaq, n1s, n2s, n2p, n3s, n3p, n3dzr, n3dx2, n3dxy, n3dxz, n3dyz, n4s, n4p, n4fxxx, n4fyyy, n4fzzz, n4fxxy, n4fxxz, n4fyyx, n4fyyz, n4fzzx, n4fzzy, n4fxyz, nsa, npa, ndzra, ndz2a, ndxya, ndxza, ndyza
   use force_mod, only: MWF
   use precision_kinds, only: dp
   use vmc, only: MBASIS, MCTYPE

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

   real(dp) :: zex(MBASIS,MWF)
   real(dp) :: betaq
   integer  :: n1s(MCTYPE)
   integer  :: n2s(MCTYPE)
   integer  :: n2p(3,MCTYPE)
   integer  :: n3s(MCTYPE)
   integer  :: n3p(3,MCTYPE)
   integer  :: n3dzr(MCTYPE)
   integer  :: n3dx2(MCTYPE)
   integer  :: n3dxy(MCTYPE)
   integer  :: n3dxz(MCTYPE)
   integer  :: n3dyz(MCTYPE)
   integer  :: n4s(MCTYPE)
   integer  :: n4p(3,MCTYPE)
   integer  :: n4fxxx(MCTYPE)
   integer  :: n4fyyy(MCTYPE)
   integer  :: n4fzzz(MCTYPE)
   integer  :: n4fxxy(MCTYPE)
   integer  :: n4fxxz(MCTYPE)
   integer  :: n4fyyx(MCTYPE)
   integer  :: n4fyyz(MCTYPE)
   integer  :: n4fzzx(MCTYPE)
   integer  :: n4fzzy(MCTYPE)
   integer  :: n4fxyz(MCTYPE)
   integer  :: nsa(MCTYPE)
   integer  :: npa(3,MCTYPE)
   integer  :: ndzra(MCTYPE)
   integer  :: ndz2a(MCTYPE)
   integer  :: ndxya(MCTYPE)
   integer  :: ndxza(MCTYPE)
   integer  :: ndx2a(MCTYPE)
   integer  :: ndyza(MCTYPE)

   private
   public :: zex, betaq, n1s, n2s, n2p, n3s, n3p, n3dzr, n3dx2, n3dxy, n3dxz, n3dyz
   public :: n4s, n4p, n4fxxx, n4fyyy, n4fzzz, n4fxxy, n4fxxz, n4fyyx, n4fyyz
   public :: n4fzzx, n4fzzy, n4fxyz, nsa, npa, ndzra, ndz2a, ndxya, ndxza, ndyza
   public :: ndx2a
   save
 end module basis

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
   use vmc, only: MELEC, MORB

   real(dp) :: b(MORB,MELEC)
   real(dp) :: tildem(MELEC,MORB,2)
   real(dp) :: xmat(MELEC**2,2) 

   private 
   public :: b, tildem, xmat
   save
 end module Bloc

 module Bloc_da
   !> Arguments: b_da, db
   use precision_kinds, only: dp
   use vmc, only: MELEC, MORB, MCENT

   real(dp) :: b_da(3,MELEC,MORB,MCENT)
   real(dp) :: db(3,MELEC,MORB,MCENT) 

   private
   public :: b_da, db
   save
 end module Bloc_da
  
 module Bloc_dj
   !> Arguments: b_dj
   use optjas, only: MPARMJ
   use precision_kinds, only: dp
   use vmc, only: MELEC, MORB

   real(dp) :: b_dj(MORB,MELEC,MPARMJ)

   private
   public :: b_dj
   save
 end module Bloc_dj

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
  !> Arguments: bot
  use precision_kinds, only: dp

   real(dp) :: bot

   private
   public :: bot
   save
 end module chck

 module coefs
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

 module config
   !> Arguments: delttn, enew, eold, nearestn, nearesto, pen, peo, psi2n, psi2o, psido, psijo, rminn, rminno, rmino, rminon, rvminn, rvminno, rvmino, rvminon, tjfn, tjfo, tjfoo, vnew, vold, xnew, xold
   use force_mod, only: MFORCE
   use precision_kinds, only: dp
   use vmc, only: MELEC
   use mstates_mod, only: MSTATES

   real(dp) :: delttn(MELEC)
   real(dp) :: enew(MFORCE)
   real(dp) :: eold(MSTATES,MFORCE)
   integer  :: nearestn(MELEC)
   integer  :: nearesto(MELEC)
   real(dp) :: pen
   real(dp) :: peo(MSTATES)
   real(dp) :: psi2n(MFORCE)
   real(dp) :: psi2o(MSTATES,MFORCE)
   real(dp) :: psido(MSTATES)
   real(dp) :: psijo
   real(dp) :: rminn(MELEC)
   real(dp) :: rminno(MELEC)
   real(dp) :: rmino(MELEC)
   real(dp) :: rminon(MELEC)
   real(dp) :: rvminn(3,MELEC)
   real(dp) :: rvminno(3,MELEC)
   real(dp) :: rvmino(3,MELEC)
   real(dp) :: rvminon(3,MELEC)
   real(dp) :: tjfn
   real(dp) :: tjfo(MSTATES)
   real(dp) :: tjfoo
   real(dp) :: vnew(3,MELEC)
   real(dp) :: vold(3,MELEC)
   real(dp) :: xnew(3,MELEC)
   real(dp) :: xold(3,MELEC)

   private
   public   :: delttn, enew, eold, nearestn, nearesto, pen, peo, psi2n
   public   :: psi2o, psido, psijo, rminn, rminno, rmino, rminon, rvminn
   public   :: rvminno, rvmino, rvminon, tjfn, tjfo, tjfoo, vnew, vold, xnew, xold
   save
 end module config

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

 module da_energy_sumcum
   !> Arguments: da_energy_cm2, da_energy_cum, da_energy_sum, da_psi_cum, da_psi_sum
   use precision_kinds, only: dp
   use vmc, only: MCENT

   real(dp) :: da_energy_cm2(3,MCENT)
   real(dp) :: da_energy_cum(3,MCENT)
   real(dp) :: da_energy_sum(3,MCENT)
   real(dp) :: da_psi_cum(3,MCENT)
   real(dp) :: da_psi_sum(3,MCENT)

   private 
   public :: da_energy_cm2, da_energy_cum, da_energy_sum, da_psi_cum, da_psi_sum 
   save
 end module da_energy_sumcum

 module da_jastrow4val
   !> Arguments: da_d2j, da_j, da_vj
   use precision_kinds, only: dp
   use vmc, only: MELEC, MCENT

   real(dp) :: da_d2j(3,MELEC,MCENT)
   real(dp) :: da_j(3,MELEC,MCENT)
   real(dp) :: da_vj(3,3,MELEC,MCENT)

   private
   public   ::  da_d2j, da_j, da_vj
   save
 end module da_jastrow4val

 module da_orbval
   !> Arguments: da_d2orb, da_dorb, da_orb
   use precision_kinds, only: dp
   use vmc, only: MELEC, MORB, MCENT

   real(dp) :: da_d2orb(3,MELEC,MORB,MCENT)
   real(dp) :: da_dorb(3,3,MELEC,MORB,MCENT)
   real(dp) :: da_orb(3,MELEC,MORB,MCENT)

   private
   public   ::  da_d2orb, da_dorb, da_orb
   save
 end module da_orbval

 module da_pseudo
   !> Arguments: da_pecent, da_vps, da_nonloc  

   use pseudo_mod, only: MPS_L
   use precision_kinds, only: dp
   use vmc, only: MELEC, MCENT


   real(dp) :: da_pecent( 3, MCENT), da_vps( 3, MELEC, MCENT, MPS_L)
   real(dp) :: da_nonloc( 3, MCENT)= 0.0D0 

   private
   public   :: da_pecent, da_vps, da_nonloc 
   save
 end module da_pseudo 
 
 module da_energy_now
   !> Arguments: da_energy, da_psi
   use precision_kinds, only: dp
   use vmc, only: MCENT
 
   real(dp) :: da_energy(3,MCENT)
   real(dp) :: da_psi(3,MCENT)
 
   private
   public   ::  da_energy, da_psi
   save
 end module da_energy_now

 module deloc_dj_m
   !> Arguments: denergy
   use optjas, only: MPARMJ
   use precision_kinds, only: dp
   use mstates_mod, only: MSTATES

   real(dp) :: denergy(MPARMJ,MSTATES)

   private 
   public :: denergy 
   save
 end module deloc_dj_m

 module denergy_det_m
   !> Arguments: denergy_det
   use precision_kinds, only: dp
   use vmc, only: MDET

    real(dp) :: denergy_det(MDET,2)

    private 
    public :: denergy_det 
    save
 end module denergy_det_m

 module denupdn
   !> Arguments: rprobdn, rprobup
   use precision_kinds, only: dp
   use vmc, only: nrad

   real(dp) :: rprobdn(nrad)
   real(dp) :: rprobup(nrad)

   private
   public   ::  rprobdn, rprobup 
   save
 end module denupdn

 module derivjas
   !> Arguments: d2g, g, go, gvalue
   use optjas, only: MPARMJ
   use precision_kinds, only: dp
   use vmc, only: MELEC

   real(dp) :: d2g(MPARMJ)
   real(dp) :: g(3,MELEC,MPARMJ)
   real(dp) :: go(MELEC,MELEC,MPARMJ)
   real(dp) :: gvalue(MPARMJ)

   private
   public   :: d2g, g, go, gvalue
   save
 end module derivjas

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

 module dorb_m
   !> Arguments: iworbd
   use vmc, only: MELEC, MDET

   integer  :: iworbd(MELEC,MDET)

   private 
   public :: iworbd 
   save
 end module dorb_m

 module elec
  !> Arguments: ndn, nup

  integer  :: ndn
  integer  :: nup

  private
  public   :: ndn, nup
  save
 end module elec

 module embed
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


 module estcum
   !> Arguments: ecum, ecum1, iblk, pecum, r2cum, tjfcum, tpbcum, avcum
   use force_mod, only: MFORCE
   use precision_kinds, only: dp
   use mstates_mod, only: MSTATES

   real(dp) :: ecum(MSTATES,MFORCE)
   real(dp) :: ecum1(MSTATES)
   integer  :: iblk
   real(dp) :: pecum(MSTATES)
   real(dp) :: r2cum
   real(dp) :: tjfcum(MSTATES)
   real(dp) :: tpbcum(MSTATES)
   real(dp) :: avcum(MSTATES*3) 

   private
   public   ::  ecum, ecum1, iblk, pecum, r2cum, tjfcum, tpbcum, avcum
   save
 end module estcum
 
 module estsig
   !> Arguments: ecm21s, ecum1s
   use precision_kinds, only: dp
   use mstates_mod, only: MSTATES

   real(dp) :: ecm21s(MSTATES)
   real(dp) :: ecum1s(MSTATES)

   private
   public   ::  ecm21s, ecum1s
   save
 end module estsig

 module estsum
  !> Arguments: acc, esum, esum1, pesum, r2sum, tjfsum, tpbsum
  use force_mod, only: MFORCE
  use precision_kinds, only: dp
  use mstates_mod, only: MSTATES

  real(dp) :: acc
  real(dp) :: esum(MSTATES,MFORCE)
  real(dp) :: esum1(MSTATES)
  real(dp) :: pesum(MSTATES)
  real(dp) :: r2sum
  real(dp) :: tjfsum(MSTATES)
  real(dp) :: tpbsum(MSTATES)

  private
  public   :: acc, esum, esum1, pesum, r2sum, tjfsum, tpbsum
  save
 end module estsum

 module estpsi
   !> Arguments: apsi, aref, detref
   use precision_kinds, only: dp
   use mstates_mod, only: MSTATES

   real(dp) :: apsi(MSTATES)
   real(dp) :: aref
   real(dp) :: detref(2)

   private
   public   ::  apsi, aref, detref 
   save
 end module estpsi

 module est2cm
   !> Arguments: ecm2, ecm21, pecm2, r2cm2, tjfcm2, tpbcm2, avcm2
   use force_mod, only: MFORCE
   use precision_kinds, only: dp
   use mstates_mod, only: MSTATES

   real(dp) :: ecm2(MSTATES,MFORCE)
   real(dp) :: ecm21(MSTATES)
   real(dp) :: pecm2(MSTATES)
   real(dp) :: r2cm2
   real(dp) :: tjfcm2(MSTATES)
   real(dp) :: tpbcm2(MSTATES)
   real(dp) :: avcm2(MSTATES*3) 

   private
   public   :: ecm2, ecm21, pecm2, r2cm2, tjfcm2, tpbcm2, avcm2 
   save
 end module est2cm


 module gauss_ecp
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

 module grdnthes
   !> Arguments: hessian_zmat
   use precision_kinds, only: dp
   use vmc, only: MCENT

   real(dp) :: hessian_zmat(3,MCENT)

   private
   public   ::  hessian_zmat 
   save
 end module grdnthes

 module grdntsmv
   !> Arguments: igrdaidx, igrdcidx, igrdmv
   use force_mod, only: MFORCE
   use vmc, only: MCENT

    integer  :: igrdaidx(MFORCE)
    integer  :: igrdcidx(MFORCE)
    integer  :: igrdmv(3,MCENT)

    private 
    public :: igrdaidx, igrdcidx, igrdmv 
    save
 end module grdntsmv

 module grdntspar
   !> Arguments: delgrdba, delgrdbl, delgrdda, delgrdxyz, igrdtype, ngradnts
   use precision_kinds, only: dp

   real(dp) :: delgrdba
   real(dp) :: delgrdbl
   real(dp) :: delgrdda
   real(dp) :: delgrdxyz
   integer  :: igrdtype
   integer  :: ngradnts

   private 
   public :: delgrdba, delgrdbl, delgrdda, delgrdxyz, igrdtype, ngradnts 
   save
 end module grdntspar
 
 

 module header
   !> Arguments: date, title

   character*20 title
   character*24 date

   private 
   public :: date, title 
   save
 end module header

 module ijasnonlin
   !> Arguments: d1d2a, d1d2b, d2d2a, d2d2b
   use precision_kinds, only: dp
   use vmc, only: MCTYPE

   real(dp) :: d1d2a(MCTYPE)
   real(dp) :: d1d2b(2)
   real(dp) :: d2d2a(MCTYPE)
   real(dp) :: d2d2b(2)

   private 
   public :: d1d2a, d1d2b, d2d2a, d2d2b 
   save
 end module ijasnonlin

 module inputflags
   !> Arguments: iznuc,igeometry,ibasis_num,ilcao,iexponents,
   !             ideterminants,ijastrow_parameter, ioptorb_def,ilattice,
   !             ici_def,iforces,icsfs,imstates,igradients,icharge_efield,
   !             imultideterminants,ioptorb_mixvirt,imodify_zmat,izmatrix_check,
   !             ihessian_zmat 

   integer :: iznuc
   integer :: igeometry 
   integer :: ibasis_num
   integer :: ilcao
   integer :: iexponents                              
   integer :: ideterminants
   integer :: ijastrow_parameter
   integer :: ioptorb_def
   integer :: ilattice
   integer :: ici_def
   integer :: iforces
   integer :: icsfs
   integer :: imstates
   integer :: igradients
   integer :: icharge_efield
   integer :: imultideterminants
   integer :: ioptorb_mixvirt
   integer :: imodify_zmat
   integer :: izmatrix_check
   integer :: ihessian_zmat

   private
   public :: iznuc,igeometry,ibasis_num,ilcao, iexponents
   public :: ideterminants,ijastrow_parameter, ioptorb_def,ilattice
   public :: ici_def,iforces,icsfs,imstates,igradients,icharge_efield
   public :: imultideterminants,ioptorb_mixvirt,imodify_zmat,izmatrix_check
   public :: ihessian_zmat
   save
 end module inputflags

 module insout
   !> Arguments: inout, inside

   integer  :: inout
   integer  :: inside

   private 
   public :: inout, inside 
   save
 end module insout

 

 module jd_scratch
   !> Arguments: qr, rr
   use sr_mod, only: MPARM
   use precision_kinds, only: dp

   real(dp) :: qr(MPARM)
   real(dp) :: rr(MPARM)

   private 
   public :: qr, rr 
   save
 end module jd_scratch

 module kinet
   !> Arguments: dtdx2n, dtdx2o
   use precision_kinds, only: dp
   use vmc, only: MELEC

   real(dp) :: dtdx2n(MELEC)
   real(dp) :: dtdx2o(MELEC)

   private 
   public :: dtdx2n, dtdx2o 
   save
 end module kinet

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

 module method_opt
   !> Arguments: method
 
   character*20 :: method

   private
   public :: method
   save
 end module method_opt

 module mix_jas_ci
   !> Arguments: de_o_ci, dj_de_ci, dj_o_ci, dj_oe_ci
   use optjas, only: MPARMJ
   use precision_kinds, only: dp
   use vmc, only: MDET

   real(dp) :: de_o_ci(MPARMJ,MDET)
   real(dp) :: dj_de_ci(MPARMJ,MDET)
   real(dp) :: dj_o_ci(MPARMJ,MDET)
   real(dp) :: dj_oe_ci(MPARMJ,MDET)

   private 
   public :: de_o_ci, dj_de_ci, dj_o_ci, dj_oe_ci 
   save
 end module mix_jas_ci

 module mix_jas_orb
   !> Arguments: de_o, dj_ho, dj_o, dj_oe
   use optorb_mod, only: MXREDUCED
   use optjas, only: MPARMJ
   use precision_kinds, only: dp
   use mstates_mod, only: MSTATES
 
   real(dp) :: de_o(MPARMJ,MXREDUCED,MSTATES)
   real(dp) :: dj_ho(MPARMJ,MXREDUCED,MSTATES)
   real(dp) :: dj_o(MPARMJ,MXREDUCED,MSTATES)
   real(dp) :: dj_oe(MPARMJ,MXREDUCED,MSTATES)
 
   private 
   public :: de_o, dj_ho, dj_o, dj_oe 
   save
 end module mix_jas_orb

 module mix_orb_ci
   !> Arguments: ci_de_o, ci_o_ho, ci_o_o, ci_o_oe
   use optorb_mod, only: MXREDUCED
   use optci, only: MXCITERM
   use precision_kinds, only: dp

   real(dp) :: ci_de_o(MXCITERM,MXREDUCED)
   real(dp) :: ci_o_ho(MXCITERM,MXREDUCED)
   real(dp) :: ci_o_o(MXCITERM,MXREDUCED)
   real(dp) :: ci_o_oe(MXCITERM,MXREDUCED)

   private 
   public :: ci_de_o, ci_o_ho, ci_o_o, ci_o_oe 
   save
 end module mix_orb_ci

 


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

 module numbas_mod
  !> Arguments: MRWF_PTS, MRWF
  integer, parameter :: MRWF_PTS=4000
  integer, parameter :: MRWF=200
  private 
  public :: MRWF, MRWF_PTS
  save 
 end module numbas_mod 

 module numexp
   !> Arguments: ae, ce
   use numbas_mod, only: MRWF
   use force_mod, only: MFORCE
   use precision_kinds, only: dp
   use vmc, only: MCTYPE
   use vmc, only: NCOEF
 
   real(dp) :: ae(2,MRWF,MCTYPE,MFORCE)
   real(dp) :: ce(NCOEF,MRWF,MCTYPE,MFORCE)
 
   private 
   public :: ae, ce 
   save
 end module numexp

 module ncusp
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

 module numbas
   !> Arguments: arg, d2rwf, igrid, iwrwf, nr, nrbas, numr, r0, rwf
   use numbas_mod, only: MRWF, MRWF_PTS
   use force_mod, only: MWF
   use precision_kinds, only: dp
   use vmc, only: MBASIS, MCTYPE

   real(dp) :: arg(MCTYPE)
   real(dp) :: d2rwf(MRWF_PTS,MRWF,MCTYPE,MWF)
   integer  :: igrid(MCTYPE)
   integer  :: iwrwf(MBASIS,MCTYPE)
   integer  :: nr(MCTYPE)
   integer  :: nrbas(MCTYPE)
   integer  :: numr
   real(dp) :: r0(MCTYPE)
   real(dp) :: rwf(MRWF_PTS,MRWF,MCTYPE,MWF)

   private
   public :: arg, d2rwf, igrid, iwrwf, nr, nrbas, numr, r0, rwf
   save
 end module numbas

 module numbas1
   !> Arguments: iwlbas, nbastyp
   use vmc, only: MBASIS, MCTYPE

   integer  :: iwlbas(MBASIS,MCTYPE)
   integer  :: nbastyp(MCTYPE)

   private
   public :: iwlbas, nbastyp
   save
 end module numbas1

 module numbas2
   !> Arguments: ibas0, ibas1
   use vmc, only: MCENT

   integer  :: ibas0(MCENT)
   integer  :: ibas1(MCENT)

   private
   public :: ibas0, ibas1
   save
 end module numbas2

 module orbital_num_lag
   !> Arguments: denom
   use precision_kinds, only: dp
   use grid_lagrange_mod, only: LAGSTART, LAGEND

   real(dp) :: denom(LAGSTART:LAGEND,3)
   real(dp) :: step_inv(3,3)

   private
   public :: denom, step_inv
   save
end module orbital_num_lag


 module pars
   !> Arguments: Z, a00, a20, a21, c0000, c1110, c2000, eps_fock, xm1, xm12, xm2, xma, xms
   use precision_kinds, only: dp

   real(dp) :: Z
   real(dp) :: a00
   real(dp) :: a20
   real(dp) :: a21
   real(dp) :: c0000
   real(dp) :: c1110
   real(dp) :: c2000
   real(dp) :: eps_fock
   real(dp) :: xm1
   real(dp) :: xm12
   real(dp) :: xm2
   real(dp) :: xma
   real(dp) :: xms

   private
   public :: Z, a00, a20, a21, c0000, c1110, c2000, eps_fock, xm1, xm12, xm2, xma, xms
   save
 end module pars

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

 module rnyucm
   !> Arguments: ll, lm

   integer  :: ll(4)
   integer  :: mm(4)
   data mm / 502,1521,4071,2107/
   data ll /   0,   0,   0,   1/

   private
   public :: ll, mm
   save
end module rnyucm

 module sa_check
   !> Arguments: energy_all, energy_err_all
   use precision_kinds, only: dp
   use mstates_mod, only: MSTATES

   real(dp) :: energy_all(MSTATES)
   real(dp) :: energy_err_all(MSTATES)

   private
   public :: energy_all, energy_err_all
   save
 end module sa_check

 module sa_weights
   !> Arguments: iweight, nweight, weights
   use precision_kinds, only: dp
   use mstates_mod, only: MSTATES

   integer  :: iweight(MSTATES)
   integer  :: nweight
   real(dp) :: weights(MSTATES)

   private
   public :: iweight, nweight, weights
   save
 end module sa_weights

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

 module slatn
   !> Arguments: slmin
   use precision_kinds, only: dp
   use vmc, only: MMAT_DIM

   real(dp) :: slmin(MMAT_DIM)

   private
   public :: slmin
   save
 end module slatn

 module spc
   !> Arguments: nsf, num

   integer  :: nsf
   integer  :: num(50)

   private
   public :: nsf, num
   save
end module spc

module spc1
   !> Arguments: csf, qsf, rsf
   use precision_kinds, only: dp

   real(dp) :: csf(750,4,50)
   real(dp) :: qsf(50,3)
   real(dp) :: rsf(50)

   private
   public :: csf, qsf, rsf
   save
end module spc1

module spc2
   !> Arguments: nxyz, sfxyz, usf
   use precision_kinds, only: dp

   integer  :: nxyz
   real(dp) :: sfxyz(5000,4)
   real(dp) :: usf(5000,3)

   private
   public :: nxyz, sfxyz, usf
   save
 end module spc2

 module stats
   !> Arguments: rejmax
   use precision_kinds, only: dp

   real(dp) :: rejmax

   private
   public :: rejmax
   save
 end module stats

 module step
   !> Arguments: ekin, ekin2, rprob, suc, trunfb, try
   use precision_kinds, only: dp
   use vmc, only: nrad

   real(dp) :: ekin(nrad)
   real(dp) :: ekin2(nrad)
   real(dp) :: rprob(nrad)
   real(dp) :: suc(nrad)
   real(dp) :: trunfb(nrad)
   real(dp) :: try(nrad)

   private
   public :: ekin, ekin2, rprob, suc, trunfb, try
   save
 end module step

 module svd_mod
  ! Not used anywhere !
  !> Arguments:
  integer, parameter :: MBUF=10000
  integer, parameter :: MXDIM=3000
  private 
  public :: MBUF, MXDIM 
  save 
 end module svd_mod 

 module tempor
   !> Arguments: dist_nn
   use precision_kinds, only: dp

   real(dp) :: dist_nn

   private
   public :: dist_nn
   save
 end module tempor

 module tempor_test
   !> Arguments: c_imag, c_real, igvec_dft, iwgvec, ngg, ngvec_dft, rkvec_tmp, rkvec_tmp2
   use ewald_mod, only: IVOL_RATIO
   use ewald_mod, only: NGVEC_BIGX
   use precision_kinds, only: dp

   real(dp) :: c_imag(NGVEC_BIGX)
   real(dp) :: c_real(NGVEC_BIGX)
   integer  :: igvec_dft(3,NGVEC_BIGX)
   integer  :: iwgvec(NGVEC_BIGX)
   integer  :: ngg(IVOL_RATIO)
   integer  :: ngvec_dft
   real(dp) :: rkvec_tmp(3)
   real(dp) :: rkvec_tmp2(3)

   private
   public :: c_imag, c_real, igvec_dft, iwgvec, ngg, ngvec_dft, rkvec_tmp, rkvec_tmp2
   save
 end module tempor_test

 module test
   !> Arguments: f, vbare_coul, vbare_jas, vbare_psp
   use ewald_mod, only: NGNORM_BIGX
   use ewald_mod, only: NGNORM_SIM_BIGX
   use precision_kinds, only: dp

   real(dp) :: f
   real(dp) :: vbare_coul(NGNORM_SIM_BIGX)
   real(dp) :: vbare_jas(NGNORM_SIM_BIGX)
   real(dp) :: vbare_psp(NGNORM_BIGX)

   private
   public :: f, vbare_coul, vbare_jas, vbare_psp
   save
 end module test

 module tmpnode
   !> Arguments: distance_node_sum
   use precision_kinds, only: dp

   real(dp) :: distance_node_sum

   private
   public :: distance_node_sum
   save
 end module tmpnode

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
   !> Arguments: transform_grd
   use precision_kinds, only: dp
   use vmc, only: MCENT3

   real(dp) :: transform_grd(MCENT3,MCENT3)

   private 
   public :: transform_grd 
   save
 end module zmatrix_grad
