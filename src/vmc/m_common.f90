!> \brief File collecting all modules that replace common blocks.  
!>
!> \author P. Lopez-Tarifa & F. Zapata NLeSC(2019)

 module precision_kinds
   ! named constants for 4, 2, and 1 byte integers:
   integer, parameter :: &
        i4b = selected_int_kind(9), &
        i2b = selected_int_kind(4), &
        i1b = selected_int_kind(2)
   ! single, double and quadruple precision reals:
   integer, parameter :: &
        sp = kind(1.0), &
        dp = selected_real_kind(2 * precision(1.0_sp)), &
        qp = selected_real_kind(2 * precision(1.0_dp))
 end module precision_kinds

 module atom
   !> Arguments: znuc, cent, pecent, iwctype, nctype, ncent
   use precision_kinds, only: dp  
   include 'vmc.h'

   real(dp) :: cent( 3, MCENT)
   real(dp) :: znuc( MCTYPE)
   real(dp) :: pecent
   integer  :: iwctype( MCENT), nctype, ncent

   private
   public   :: znuc, cent, pecent, iwctype, nctype, ncent 
   save
 end module atom
   
 module config
   !> Arguments: delttn, enew, eold, nearestn, nearesto, pen, peo, psi2n, psi2o, psido, psijo, rminn, rminno, rmino, rminon, rvminn, rvminno, rvmino, rvminon, tjfn, tjfo, tjfoo, vnew, vold, xnew, xold
   use precision_kinds, only: dp
   include 'vmc.h'
   include 'force.h'
   include 'mstates.h'

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

 module const
   !> Arguments: pi, hb, etrial, delta, deltai, fbias, nelec, imetro, ipr
   use precision_kinds, only: dp
   include 'vmc.h'

   real(dp) :: delta
   real(dp) :: deltai
   real(dp) :: etrial
   real(dp) :: fbias
   real(dp) :: hb
   integer  :: imetro
   integer  :: ipr
   integer  :: nelec
   real(dp) :: pi

   private
   public   :: pi, hb, etrial, delta, deltai, fbias, nelec, imetro, ipr
   save
 end module const

 module contrl_per
   !> Arguments: iperiodic, ibasis 

   integer  :: iperiodic, ibasis

   private
   public   :: iperiodic, ibasis
   save
 end module contrl_per
 
 module csfs
   !> Arguments: ccsf, cxdet, iadet, ibdet, icxdet, ncsf, nstates
   use precision_kinds, only: dp
   include 'vmc.h'
   include 'mstates.h'
   include 'force.h'

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
 
 module da_jastrow4val
   !> Arguments: da_d2j, da_j, da_vj
   use precision_kinds, only: dp
   include 'vmc.h'

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
   include 'vmc.h'

   real(dp) :: da_d2orb(3,MELEC,MORB,MCENT)
   real(dp) :: da_dorb(3,3,MELEC,MORB,MCENT)
   real(dp) :: da_orb(3,MELEC,MORB,MCENT)

   private
   public   ::  da_d2orb, da_dorb, da_orb
   save
 end module da_orbval

 module da_pseudo
   !> Arguments: da_pecent, da_vps, da_nonloc  

   use precision_kinds, only: dp  

   include 'vmc.h'
   include 'pseudo.h'

   real(dp) :: da_pecent( 3, MCENT), da_vps( 3, MELEC, MCENT, MPS_L)
   real(dp) :: da_nonloc( 3, MCENT)= 0.0D0 

   private
   public   :: da_pecent, da_vps, da_nonloc 
   save
 end module da_pseudo 
 
 module da_energy_now
   !> Arguments: da_energy, da_psi
   use precision_kinds, only: dp
   include 'vmc.h'
 
   real(dp) :: da_energy(3,MCENT)
   real(dp) :: da_psi(3,MCENT)
 
   private
   public   ::  da_energy, da_psi
   save
 end module da_energy_now

 module denupdn
   !> Arguments: rprobdn, rprobup
   use precision_kinds, only: dp
   include 'vmc.h'

   real(dp) :: rprobdn(nrad)
   real(dp) :: rprobup(nrad)

   private
   public   ::  rprobdn, rprobup 
   save
 end module denupdn

 module derivjas
   !> Arguments: d2g, g, go, gvalue
   use precision_kinds, only: dp
   include 'vmc.h'
   include 'optjas.h'

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
   use precision_kinds, only: dp
   include 'vmc.h'
   include 'force.h'
   include 'mstates.h'

   real(dp) :: cdet(MDET,MSTATES,MWF)
   integer  :: ndet

   private
   public   :: cdet, ndet
   save
 end module dets

 module dets_equiv
  !> Arguments: cdet_equiv, dcdet_equiv
   use precision_kinds, only: dp
   include 'vmc.h'

   real(dp) :: cdet_equiv(MDET)
   real(dp) :: dcdet_equiv(MDET)

   private
   public   ::  cdet_equiv, dcdet_equiv
   save
 end module dets_equiv

 module distances_sav
  !> Arguments: r_ee_sav, r_en_sav, rshift_sav, rvec_ee_sav, rvec_en_sav
  use precision_kinds, only: dp
  include 'vmc.h'

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
  use precision_kinds, only: dp
  include 'vmc.h'

  integer  :: ndn
  integer  :: nup

  private
  public   :: ndn, nup
  save
 end module elec

 module estcum
   !> Arguments: ecum, ecum1, iblk, pecum, r2cum, tjfcum, tpbcum, avcum
   use precision_kinds, only: dp
   include 'mstates.h'
   include 'force.h'

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
   include 'mstates.h'

   real(dp) :: ecm21s(MSTATES)
   real(dp) :: ecum1s(MSTATES)

   private
   public   ::  ecm21s, ecum1s
   save
 end module estsig

 module estsum
  !> Arguments: acc, esum, esum1, pesum, r2sum, tjfsum, tpbsum
  use precision_kinds, only: dp
  include 'force.h'
  include 'mstates.h'

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
   include 'mstates.h'

   real(dp) :: apsi(MSTATES)
   real(dp) :: aref
   real(dp) :: detref(2)

   private
   public   ::  apsi, aref, detref 
   save
 end module estpsi

 module est2cm
   !> Arguments: ecm2, ecm21, pecm2, r2cm2, tjfcm2, tpbcm2, avcm2
   use precision_kinds, only: dp
   include 'mstates.h'
   include 'force.h'

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

 module ewald
   !> Arguments: b_coul, b_coul_sim, y_coul, y_coul_sim
   use precision_kinds, only: dp
   include 'ewald.h'

   real(dp) :: b_coul(NCOEFX)
   real(dp) :: b_coul_sim(NCOEFX)
   real(dp) :: y_coul(NGNORMX)
   real(dp) :: y_coul_sim(NGNORM_SIMX)

   private
   public   ::  b_coul, b_coul_sim, y_coul, y_coul_sim
   save
 end module ewald

 module ewald_basis
   !> Arguments: vps_basis_fourier
   use precision_kinds, only: dp
   include 'ewald.h'

   real(dp) :: vps_basis_fourier(NGNORM_BIGX)

   private
   public   :: vps_basis_fourier
   save
 end module ewald_basis

 module force_analy 
   !> Arguments: iforce_analy 

   integer  :: iforce_analy 

   private
   public   :: iforce_analy 
   save
 end module force_analy 

 module forcest
   !> Arguments: fcm2, fcum
   use precision_kinds, only: dp
   include 'force.h'
   include 'mstates.h'

   real(dp) :: fcm2(MSTATES,MFORCE)
   real(dp) :: fcum(MSTATES,MFORCE)

   private
   public   ::  fcm2, fcum
   save
 end module forcest

 module forcestr
   !> Arguments: delc
   use precision_kinds, only: dp
   include 'vmc.h'
   include 'force.h'
 
   real(dp) :: delc(3,MCENT,MFORCE)
 
   private
   public   ::  delc 
   save
 end module forcestr

 module forcewt
   !> Arguments: wcum, wsum
   use precision_kinds, only: dp
   include 'mstates.h'
   include 'force.h'

   real(dp) :: wcum(MSTATES,MFORCE)
   real(dp) :: wsum(MSTATES,MFORCE)

   private
   public   ::  wcum, wsum 
   save
end module forcewt

 module force_dmc
   !> Arguments: itausec, nwprod
   use precision_kinds, only: dp

   integer  :: itausec
   integer  :: nwprod

   private
   public   ::   itausec, nwprod
   save
 end module force_dmc

 module force_fin
   !> Arguments: da_energy_ave, da_energy_err
   use precision_kinds, only: dp
   include 'vmc.h'

   real(dp) :: da_energy_ave(3,MCENT)
   real(dp) :: da_energy_err(3)

   private
   public   :: da_energy_ave, da_energy_err
   save
 end module force_fin

 module force_mat_n
   !> Arguments: force_o
   use precision_kinds, only: dp
   include 'vmc.h'
   include 'sr.h'

   real(dp) :: force_o(6*MCENT,MCONF)

   private
   public   ::  force_o 
   save
 end module force_mat_n  

 module forcepar
   !> Arguments: deltot, istrech, nforce
   use precision_kinds, only: dp
   include 'force.h'

   real(dp) :: deltot(MFORCE)
   integer  :: istrech
   integer  :: nforce

   private
   public   ::  deltot, istrech, nforce 
   save
 end module forcepar

 module gauss_ecp
   !> Arguments: ecp_coef, ecp_exponent, necp_power, necp_term
   use precision_kinds, only: dp
   include 'vmc.h'
   include 'pseudo.h'
 
   real(dp) :: ecp_coef(MGAUSS,MPS_L,MCTYPE)
   real(dp) :: ecp_exponent(MGAUSS,MPS_L,MCTYPE)
   integer  :: necp_power(MGAUSS,MPS_L,MCTYPE)
   integer  :: necp_term(MPS_L,MCTYPE)
 
   private
   public   ::  ecp_coef, ecp_exponent, necp_power, necp_term
   save
 end module gauss_ecp

 module ghostatom
   !> Arguments: newghostype, nghostcent
   use precision_kinds, only: dp

   integer  :: newghostype
   integer  :: nghostcent

   private
   public   :: newghostype, nghostcent
   save
 end module ghostatom

 module gradhessjo
   !> Arguments: d1d2a_old, d1d2b_old, d2d2a_old, d2d2b_old, denergy_old, gvalue_old
   use precision_kinds, only: dp
   include 'vmc.h'
   include 'mstates.h'
   include 'optjas.h'

   real(dp) :: d1d2a_old(MCTYPE)
   real(dp) :: d1d2b_old(2)
   real(dp) :: d2d2a_old(MCTYPE)
   real(dp) :: d2d2b_old(2)
   real(dp) :: denergy_old(MPARMJ,MSTATES)
   real(dp) :: gvalue_old(MPARMJ)

   private
   public   ::  d1d2a_old, d1d2b_old, d2d2a_old, d2d2b_old, denergy_old, gvalue_old
   save
 end module gradhessjo

 module gradhess_ci
   !> Arguments: grad_ci, h_ci, s_ci
   use precision_kinds, only: dp
   include 'optci.h'

   real(dp) :: grad_ci(MXCITERM)
   real(dp) :: h_ci(MXCITERM,MXCIREDUCED)
   real(dp) :: s_ci(MXCITERM,MXCIREDUCED)

   private
   public   ::  grad_ci, h_ci, s_ci
   save
 end module gradhess_ci

 module gradhess_jas
   !> Arguments: grad_jas, h_jas, s_jas
   use precision_kinds, only: dp
   include 'optjas.h'
 
   real(dp) :: grad_jas(MPARMJ)
   real(dp) :: h_jas(MPARMJ,MPARMJ)
   real(dp) :: s_jas(MPARMJ,MPARMJ)
 
   private
   public   ::  grad_jas, h_jas, s_jas
   save
 end module gradhess_jas
 
 module gradhess_mix_jas_ci
   !> Arguments: h_mix_jas_ci, s_mix_jas_ci
   use precision_kinds, only: dp
   include 'optjas.h'
   include 'optci.h'

   real(dp) :: h_mix_jas_ci(2*MPARMJ,MXCITERM)
   real(dp) :: s_mix_jas_ci(MPARMJ,MXCITERM)

   private
   public   ::  h_mix_jas_ci, s_mix_jas_ci
   save
 end module gradhess_mix_jas_ci

 module gradhess_mix_jas_orb
   !> Arguments: h_mix_jas_orb, s_mix_jas_orb
   use precision_kinds, only: dp
   include 'optorb.h'
   include 'optjas.h'

   real(dp) :: h_mix_jas_orb(2*MPARMJ,MXREDUCED)
   real(dp) :: s_mix_jas_orb(MPARMJ,MXREDUCED)

   private
   public   ::  h_mix_jas_orb, s_mix_jas_orb
   save
 end module gradhess_mix_jas_orb

 module gradhess_mix_orb_ci
   !> Arguments: h_mix_ci_orb, s_mix_ci_orb
   use precision_kinds, only: dp
   include 'optorb.h'
   include 'optci.h'

   real(dp) :: h_mix_ci_orb(2*MXCITERM,MXREDUCED)
   real(dp) :: s_mix_ci_orb(MXCITERM,MXREDUCED)

   private
   public   ::  h_mix_ci_orb, s_mix_ci_orb
   save
 end module gradhess_mix_orb_ci

 module gradjerr
   !> Arguments: dj_bsum, dj_e_bsum, dj_e_save, dj_save, e_bsum, grad_jas_bcm2, grad_jas_bcum
   use precision_kinds, only: dp
   include 'optjas.h'
   include 'mstates.h'

   real(dp) :: dj_bsum(MPARMJ,MSTATES)
   real(dp) :: dj_e_bsum(MPARMJ,MSTATES)
   real(dp) :: dj_e_save(MPARMJ,MSTATES)
   real(dp) :: dj_save(MPARMJ,MSTATES)
   real(dp) :: e_bsum(MSTATES)
   real(dp) :: grad_jas_bcm2(MPARMJ,MSTATES)
   real(dp) :: grad_jas_bcum(MPARMJ,MSTATES)

   private
   public   ::  dj_bsum, dj_e_bsum, dj_e_save, dj_save, e_bsum, grad_jas_bcm2, grad_jas_bcum
   save
 end module gradjerr

 module grdnthes
   !> Arguments: hessian_zmat
   use precision_kinds, only: dp
   include 'vmc.h'

   real(dp) :: hessian_zmat(3,MCENT)

   private
   public   ::  hessian_zmat 
   save
 end module grdnthes

 module grdntsmv
   !> Arguments: igrdaidx, igrdcidx, igrdmv
   use precision_kinds, only: dp
   include 'vmc.h'
   include 'force.h'

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
   use precision_kinds, only: dp

   character*20 title
   character*24 date

   private 
   public :: date, title 
   save
 end module header

 module ijasnonlin
   !> Arguments: d1d2a, d1d2b, d2d2a, d2d2b
   use precision_kinds, only: dp
   include 'vmc.h'

   real(dp) :: d1d2a(MCTYPE)
   real(dp) :: d1d2b(2)
   real(dp) :: d2d2a(MCTYPE)
   real(dp) :: d2d2b(2)

   private 
   public :: d1d2a, d1d2b, d2d2a, d2d2b 
   save
 end module ijasnonlin

 module insout
   !> Arguments: inout, inside
   use precision_kinds, only: dp
   include 'vmc.h'

   integer  :: inout
   integer  :: inside

   private 
   public :: inout, inside 
   save
 end module insout

 module jaso
   !> Arguments: d2ijo, d2o, fijo, fjo, fso, fsumo
   use precision_kinds, only: dp
   include 'vmc.h'

    real(dp) :: d2ijo(MELEC,MELEC)
    real(dp) :: d2o
    real(dp) :: fijo(3,MELEC,MELEC)
    real(dp) :: fjo(3,MELEC)
    real(dp) :: fso(MELEC,MELEC)
    real(dp) :: fsumo

    private 
    public :: d2ijo, d2o, fijo, fjo, fso, fsumo 
    save
 end module jaso

 module jaspar
   !> Arguments: nspin1, nspin2, sspin, sspinn, is
   use precision_kinds, only: dp

   integer  :: is
   integer  :: nspin1
   integer  :: nspin2
   real(dp) :: sspin
   real(dp) :: sspinn

   private
   public   :: nspin1, nspin2, sspin, sspinn, is
   save
 end module jaspar

 module jaspar1
   !> Arguments: cjas1, cjas2
   use precision_kinds, only: dp
   include 'force.h'

   real(dp) :: cjas1(MWF)
   real(dp) :: cjas2(MWF)

   private
   public   ::  cjas1, cjas2
   save
 end module jaspar1

 module jaspar2
   !> Arguments: a1, a2
   use precision_kinds, only: dp
   include 'force.h'

   real(dp) :: a1(83,3,MWF)
   real(dp) :: a2(83,3,MWF)

   private 
   public :: a1, a2 
   save
 end module jaspar2

 module jaspar3
   !> Arguments: a, b, c, fck, nord, scalek
   use precision_kinds, only: dp
   include 'vmc.h'
   include 'force.h'

   real(dp) :: a(MORDJ1,MWF)
   real(dp) :: b(MORDJ1,2,MWF)
   real(dp) :: c(83,MCTYPE,MWF)
   real(dp) :: fck(15,MCTYPE,MWF)
   integer  :: nord
   real(dp) :: scalek(MWF)

   private 
   public :: a, b, c, fck, nord, scalek 
   save
 end module jaspar3

 module jaspar4
   !> Arguments: a4, norda, nordb, nordc
   use precision_kinds, only: dp
   include 'vmc.h'
   include 'force.h'

   real(dp) :: a4(MORDJ1,MCTYPE,MWF)
   integer  :: norda
   integer  :: nordb
   integer  :: nordc

   private 
   public :: a4, norda, nordb, nordc 
   save
 end module jaspar4

 module jaspar6
   !> Arguments: asymp_jasa, asymp_jasb, asymp_r, c1_jas6, c1_jas6i, c2_jas6, cutjas, cutjasi
   use precision_kinds, only: dp
   include 'vmc.h'

   real(dp) :: asymp_jasa(MCTYPE)
   real(dp) :: asymp_jasb(2)
   real(dp) :: asymp_r
   real(dp) :: c1_jas6
   real(dp) :: c1_jas6i
   real(dp) :: c2_jas6
   real(dp) :: cutjas
   real(dp) :: cutjasi

   private 
   public :: asymp_jasa, asymp_jasb, asymp_r, c1_jas6, c1_jas6i, c2_jas6, cutjas, cutjasi 
   save
 end module jaspar6

 module jaspointer
   !> Arguments: npoint, npointa
   use precision_kinds, only: dp
   include 'vmc.h'

   integer  :: npoint(MCTYP3X)
   integer  :: npointa(3*MCTYP3X)

   private 
   public :: npoint, npointa 
   save
 end module jaspointer

 module jd_scratch
   !> Arguments: qr, rr
   use precision_kinds, only: dp
   include 'sr.h'

   real(dp) :: qr(MPARM)
   real(dp) :: rr(MPARM)

   private 
   public :: qr, rr 
   save
 end module jd_scratch

 module kinet
   !> Arguments: dtdx2n, dtdx2o
   use precision_kinds, only: dp
   include 'vmc.h'

   real(dp) :: dtdx2n(MELEC)
   real(dp) :: dtdx2o(MELEC)

   private 
   public :: dtdx2n, dtdx2o 
   save
 end module kinet

 module linear_norm
   !> Arguments: oav, ci_oav
   use precision_kinds, only: dp
   include 'optci.h'

   real(dp) :: oav(MXCITERM)
   real(dp) :: ci_oav(MXCITERM) 

   private 
   public :: oav, ci_oav
   save
 end module linear_norm

 module mix_jas_ci
   !> Arguments: de_o_ci, dj_de_ci, dj_o_ci, dj_oe_ci
   use precision_kinds, only: dp
   include 'vmc.h'
   include 'optjas.h'

   real(dp) :: de_o_ci(MPARMJ,MDET)
   real(dp) :: dj_de_ci(MPARMJ,MDET)
   real(dp) :: dj_o_ci(MPARMJ,MDET)
   real(dp) :: dj_oe_ci(MPARMJ,MDET)

   private 
   public :: de_o_ci, dj_de_ci, dj_o_ci, dj_oe_ci 
   save
 end module mix_jas_ci

 module mmpol_hpsi
   !> Arguments: eek_pol, peQMdp, peQMq
   use precision_kinds, only: dp
   include 'mmpol.h'

   real(dp) :: eek_pol(3,MCHMM)
   real(dp) :: peQMdp
   real(dp) :: peQMq

   private 
   public :: eek_pol, peQMdp, peQMq
   save
 end module mmpol_hpsi

 module mmpolo
   !> Arguments: cmmpolo, dmmpolo, eeko
   use precision_kinds, only: dp
   include 'mmpol.h'

   real(dp) :: cmmpolo
   real(dp) :: dmmpolo
   real(dp) :: eeko(3,MCHMM)

   private 
   public :: cmmpolo, dmmpolo, eeko
   save
 end module mmpolo

 module mpiconf
   !> Arguments: idtask, nproc, wid
    integer  :: idtask
    integer  :: nproc
    logical  :: wid 

    private 
    public :: idtask, nproc, wid 
    save
 end module mpiconf

 module multidet
   !> Arguments: iactv, irepcol_det, ireporb_det, ivirt, iwundet, kref, numrep_det
   use precision_kinds, only: dp
   include 'vmc.h'

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

 module ncusp
   !> Arguments: ncnstr, ncuspc, nfock, nfockc, norbc
   use precision_kinds, only: dp

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
   use precision_kinds, only: dp
   include 'vmc.h'
   include 'force.h'
   include 'numbas.h'

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
   use precision_kinds, only: dp
   include 'vmc.h'

   integer  :: iwlbas(MBASIS,MCTYPE)
   integer  :: nbastyp(MCTYPE)

   private
   public :: iwlbas, nbastyp
   save
 end module numbas1

 module numbas2
   !> Arguments: ibas0, ibas1
   use precision_kinds, only: dp
   include 'vmc.h'

   integer  :: ibas0(MCENT)
   integer  :: ibas1(MCENT)

   private
   public :: ibas0, ibas1
   save
 end module numbas2

 module optorb_mix
   !> Arguments: iwmix_virt, norbopt, norbvirt
   use precision_kinds, only: dp
   include 'vmc.h'

   integer  :: iwmix_virt(MORB,MORB)
   integer  :: norbopt
   integer  :: norbvirt

   private
   public :: iwmix_virt, norbopt, norbvirt
   save
 end module optorb_mix
