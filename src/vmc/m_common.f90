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

 module b_tmove
   !> Arguments: b_t, iskip
   use precision_kinds, only: dp
   include 'vmc.h'
   include 'pseudo.h'

   real(dp) :: b_t(MORB,MPS_QUAD,MCENT,MELEC)
   integer  :: iskip(MELEC,MCENT)

   private
   public :: b_t, iskip
   save
 end module b_tmove

 module Bloc
   !> Arguments: b, tildem, xmat, xmatd, xmatu
   use precision_kinds, only: dp
   include 'vmc.h'

   real(dp) :: b(MORB,MELEC)
   real(dp) :: tildem(MELEC,MORB,2)
   real(dp) :: xmat(MELEC**2,2) 
   real(dp) :: xmatd(MELEC**2)
   real(dp) :: xmatu(MELEC**2)

   private 
   public :: b, tildem, xmat, xmatd, xmatu 
   save
 end module Bloc

 module Bloc_da
   !> Arguments: b_da, db
   use precision_kinds, only: dp
   include 'vmc.h'

   real(dp) :: b_da(3,MELEC,MORB,MCENT)
   real(dp) :: db(3,MELEC,MORB,MCENT) 

   private
   public :: b_da, db
   save
 end module Bloc_da
   
 module Bloc_dj
   !> Arguments: b_dj
   use precision_kinds, only: dp
   include 'vmc.h'
   include 'optjas.h'

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
   use precision_kinds, only: dp
   include 'vmc.h'
   include 'pseudo.h'

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
  include 'vmc.h'

   real(dp) :: bot

   private
   public :: bot
   save
 end module chck

 module coefs
   !> Arguments: coef, nbasis, norb
   use precision_kinds, only: dp
   include 'vmc.h'
   include 'force.h'

   real(dp) :: coef(MBASIS,MORB,MWF)
   integer  :: nbasis
   integer  :: norb

   private
   public :: coef, nbasis, norb
   save
 end module coefs

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

 module const2
  !> Arguments: deltar, deltat
  use precision_kinds, only: dp

   real(dp) :: deltar
   real(dp) :: deltat

   private
   public :: deltar, deltat
   save
 end module const2

 module constant
  !> Arguments: twopi
  use precision_kinds, only: dp

  real(dp) :: twopi

  private
  public :: twopi
  save
 end module constant

 module contrl
  !> Arguments: idump, irstar, isite, n_conf, nblk, nblkeq, nconf_new, nstep
  use precision_kinds, only: dp
  include 'vmc.h'

   integer  :: idump
   integer  :: irstar
   integer  :: isite
   integer  :: n_conf
   integer  :: nblk
   integer  :: nblkeq
   integer  :: nconf_new
   integer  :: nstep

   private 
   public :: idump, irstar, isite, n_conf, nblk, nblkeq, nconf_new, nstep 
   save
 end module contrl

 module contr2
   !> Arguments: i3body, ianalyt_lap, iaver, icusp, icusp2, ifock, ijas, irewgt, isc, istrch

    integer  :: i3body
    integer  :: ianalyt_lap
    integer  :: iaver
    integer  :: icusp
    integer  :: icusp2
    integer  :: ifock
    integer  :: ijas
    integer  :: irewgt
    integer  :: isc
    integer  :: istrch

    private
    public :: i3body, ianalyt_lap, iaver, icusp, icusp2, ifock, ijas, irewgt, isc, istrch
    save
 end module contr2

 module contr3
  !> Arguments: mode

  character*12 :: mode

  private
  public :: mode
  save
end module contr3

 module contrl_per
   !> Arguments: iperiodic, ibasis 

   integer  :: iperiodic, ibasis

   private
   public   :: iperiodic, ibasis
   save
 end module contrl_per
 
 module contrldmc
   !> Arguments: iacc_rej, icross, icuspg, icut_br, icut_e, idiv_v, idmc, ipq, itau_eff, nfprod, rttau, tau, taueff, tautot
   use precision_kinds, only: dp
   include 'force.h'

   integer  :: iacc_rej
   integer  :: icross
   integer  :: icuspg
   integer  :: icut_br
   integer  :: icut_e
   integer  :: idiv_v
   integer  :: idmc
   integer  :: ipq
   integer  :: itau_eff
   integer  :: nfprod
   real(dp) :: rttau
   real(dp) :: tau
   real(dp) :: taueff(MFORCE)
   real(dp) :: tautot

   private
   public :: iacc_rej, icross, icuspg, icut_br, icut_e, idiv_v, idmc, ipq, itau_eff, nfprod, rttau, tau, taueff, tautot
   save
 end module contrldmc

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

 module cuspmat
   !> Arguments: cm, ishe, iwc3, neqs
   use precision_kinds, only: dp
   include 'vmc.h'

   real(dp) :: cm(NEQSX,NEQSX)
   integer  :: ishe
   integer  :: iwc3(NEQSX)
   integer  :: neqs

   private 
   public :: cm, ishe, iwc3, neqs 
   save
 end module cuspmat

 module da_energy_ave_m
   !> Arguments: da_energy_cm2, da_energy_cum, da_energy_sum, da_psi_cum, da_psi_sum
   use precision_kinds, only: dp
   include 'vmc.h'

   real(dp) :: da_energy_cm2(3,MCENT)
   real(dp) :: da_energy_cum(3,MCENT)
   real(dp) :: da_energy_sum(3,MCENT)
   real(dp) :: da_psi_cum(3,MCENT)
   real(dp) :: da_psi_sum(3,MCENT)

   private 
   public :: da_energy_cm2, da_energy_cum, da_energy_sum, da_psi_cum, da_psi_sum 
   save
 end module da_energy_ave_m 

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

 module deloc_dj_m
   !> Arguments: denergy
   use precision_kinds, only: dp
   include 'optjas.h'
   include 'mstates.h'

   real(dp) :: denergy(MPARMJ,MSTATES)

   private 
   public :: denergy 
   save
 end module deloc_dj_m

 module denergy_det_m
   !> Arguments: denergy_det
   use precision_kinds, only: dp
   include 'vmc.h'

    real(dp) :: denergy_det(MDET,2)

    private 
    public :: denergy_det 
    save
 end module denergy_det_m

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

 module dorb_m
   !> Arguments: iworbd
   use precision_kinds, only: dp
   include 'vmc.h'

   integer  :: iworbd(MELEC,MDET)

   private 
   public :: iworbd 
   save
end module dorb_m

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

 module gradhessj
   !> Arguments: d2j, d2j_e, de, de_de, de_e, dj, dj_de, dj_dj, dj_dj_e, dj_e, dj_e2, e2
   use precision_kinds, only: dp
   include 'mstates.h'
   include 'optjas.h'

   real(dp) :: d2j(MPARMJ,MPARMJ,MSTATES)
   real(dp) :: d2j_e(MPARMJ,MPARMJ,MSTATES)
   real(dp) :: de(MPARMJ,MSTATES)
   real(dp) :: de_de(MPARMJ,MPARMJ,MSTATES)
   real(dp) :: de_e(MPARMJ,MSTATES)
   real(dp) :: dj(MPARMJ,MSTATES)
   real(dp) :: dj_de(MPARMJ,MPARMJ,MSTATES)
   real(dp) :: dj_dj(MPARMJ,MPARMJ,MSTATES)
   real(dp) :: dj_dj_e(MPARMJ,MPARMJ,MSTATES)
   real(dp) :: dj_e(MPARMJ,MSTATES)
   real(dp) :: dj_e2(MPARMJ,MSTATES)
   real(dp) :: e2(MPARMJ,MSTATES)

   private 
   public :: d2j, d2j_e, de, de_de, de_e, dj, dj_de, dj_dj, dj_dj_e, dj_e, dj_e2, e2 
   save
 end module gradhessj

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

 module gradjerrb
  !> Arguments: nbj_current, ngrad_jas_bcum, ngrad_jas_blocks, njb_current
  use precision_kinds, only: dp
  include 'vmc.h'

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

 module jasn
   !> Arguments: d2ijn, d2n, fijn, fjn, fsn, fsumn
   use precision_kinds, only: dp
   include 'vmc.h'

   real(dp) :: d2ijn(MELEC,MELEC)
   real(dp) :: d2n
   real(dp) :: fijn(3,MELEC,MELEC)
   real(dp) :: fjn(3,MELEC)
   real(dp) :: fsn(MELEC,MELEC)
   real(dp) :: fsumn

   private 
   public :: d2ijn, d2n, fijn, fjn, fsn, fsumn 
   save
 end module jasn

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

 module multimat
   !> Arguments: aa, wfmat
   use precision_kinds, only: dp
   include 'vmc.h'

   real(dp) :: aa(MELEC,MORB,2)
   real(dp) :: wfmat(MEXCIT**2,MDET,2)

   private 
   public :: aa, wfmat 
   save
 end module multimat

 module multimatn
   !> Arguments: aan, wfmatn
   use precision_kinds, only: dp
   include 'vmc.h'

   real(dp) :: aan(MELEC,MORB)
   real(dp) :: wfmatn(MEXCIT**2,MDET)

   private 
   public :: aan, wfmatn 
   save
 end module multimatn

 module numexp
   !> Arguments: ae, ce
   use precision_kinds, only: dp
   include 'vmc.h'
   include 'force.h'
   include 'numbas.h'
 
   real(dp) :: ae(2,MRWF,MCTYPE,MFORCE)
   real(dp) :: ce(NCOEF,MRWF,MCTYPE,MFORCE)
 
   private 
   public :: ae, ce 
   save
 end module numexp

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

 module optorb
   !> Arguments: dmat_diag, irrep, orb_energy
   use precision_kinds, only: dp
   include 'vmc.h'

   real(dp) :: dmat_diag(MORB)
   integer  :: irrep(MORB)
   real(dp) :: orb_energy(MORB)

   private 
   public :: dmat_diag, irrep, orb_energy 
   save
 end module optorb

 module optorb_cblock   ! from optorb.h 

   ! norbterm: number of terms (possibly after a transformation)
   ! norbprim: number of primitive terms (determinant ratios)
   integer :: norbterm
   integer :: norbprim

   ! PLT: From old common block /orb004/
   integer :: nefp_blocks
   integer :: nb_current
   integer :: norb_f_bcum

   ! reduced correlation matrix pointers
   ! threshold in terms of std dev. , limit for keeping operators
   ! if iuse_trafo: linearly transformed operators sampled instead of primitive
   !     replacement operators 
   ! PLT: From old common blocks /orb006/ and /orb008/.
   integer :: isample_cmat
   integer :: nreduced 
   integer :: iuse_trafiuse_trafoo 

   ! Dumping block averages for error analysis. 
   ! PLT: From old common blocks /orb009/ and /orb010/.
   integer :: idump_blockav 
   integer :: iorbsample 
   integer :: ns_current

   ! Printing flags:
   integer :: iorbprt
   integer :: iorbprt_sav

   private 
   public :: norbterm, norbprim
   public :: nefp_blocks, nb_current, norb_f_bcum
   public :: isample_cmat, nreduced, iuse_trafiuse_trafoo 
   public :: idump_blockav, iorbsample, ns_current 
   public :: iorbprt, iorbprt_sav
   save
 end module optorb_cblock 

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

 module optwf_contrl
   !> Arguments: ioptci, ioptjas, ioptorb, nparm

    integer  :: ioptci
    integer  :: ioptjas
    integer  :: ioptorb
    integer  :: nparm

    private
    public :: ioptci, ioptjas, ioptorb, nparm
    save
 end module optwf_contrl

 module optwf_corsam
   !> Arguments: add_diag_tmp, energy, energy_err, force, force_err
   use precision_kinds, only: dp
   include 'force.h'

   real(dp) :: add_diag(MFORCE)
   real(dp) :: add_diag_tmp(MFORCE)
   real(dp) :: energy(MFORCE)
   real(dp) :: energy_err(MFORCE)
   real(dp) :: force(MFORCE)
   real(dp) :: force_err(MFORCE)

   private
   public :: add_diag, add_diag_tmp, energy, energy_err, force, force_err
   save
 end module optwf_corsam

 module optwf_func
   !> Arguments: ifunc_omega, omega, omega_hes
   use precision_kinds, only: dp

   integer  :: ifunc_omega
   real(dp) :: omega
   real(dp) :: omega_hes

   private
   public :: ifunc_omega, omega, omega_hes
   save
 end module optwf_func

 module optwf_nparmj
   !> Arguments: nparma, nparmb, nparmc, nparmf
   include 'vmc.h'

   integer  :: nparma(MCTYP3X)
   integer  :: nparmb(3)
   integer  :: nparmc(MCTYPE)
   integer  :: nparmf(MCTYPE)

   private
   public :: nparma, nparmb, nparmc, nparmf
   save
 end module optwf_nparmj

 module optwf_parms
   !> Arguments: nparmd, nparme, nparmg, nparmj, nparml, nparms

   integer  :: nparmd
   integer  :: nparme
   integer  :: nparmg
   integer  :: nparmj
   integer  :: nparml
   integer  :: nparms

   private
   public :: nparmd, nparme, nparmg, nparmj, nparml, nparms
   save
 end module optwf_parms

 module optwf_wjas
   !> Arguments: iwjasa, iwjasb, iwjasc, iwjasf
   use precision_kinds, only: dp
   include 'vmc.h'

   integer  :: iwjasa(83,MCTYP3X)
   integer  :: iwjasb(83,3)
   integer  :: iwjasc(83,MCTYPE)
   integer  :: iwjasf(15,MCTYPE)

   private
   public :: iwjasa, iwjasb, iwjasc, iwjasf
   save
 end module optwf_wjas

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

 module pcm_force
   !> Arguments: sch_s
   use precision_kinds, only: dp
   include 'pcm.h'
   include 'force.h'

   real(dp) :: sch_s(MCHS,MFORCE)

   private
   public :: sch_s
   save
 end module pcm_force

 module pcm_hpsi
   !> Arguments: enfpcm, pepcms, pepcmv, qopcm
   use precision_kinds, only: dp
   include 'pcm.h'

   real(dp) :: enfpcm(MCHS)
   real(dp) :: pepcms
   real(dp) :: pepcmv
   real(dp) :: qopcm

   private
   public :: enfpcm, pepcms, pepcmv, qopcm
   save
 end module pcm_hpsi

 module pcm_num_spl2
   !> Arguments: bc, wk
   use precision_kinds, only: dp

   real(dp) :: bc
   real(dp) :: wk

   private
   public :: bc, wk
   save
 end module pcm_num_spl2

 module pcm_xv_new
   !> Arguments: xv_new
   use precision_kinds, only: dp
   include 'pcm.h'

   real(dp) :: xv_new(3,MCHV)

   private
   public :: xv_new
   save
 end module pcm_xv_new

 module pcmo
   !> Arguments: enfpcmo, qopcmo, spcmo, vpcmo
   use precision_kinds, only: dp
   include 'pcm.h'

   real(dp) :: enfpcmo(MCHS)
   real(dp) :: qopcmo
   real(dp) :: spcmo
   real(dp) :: vpcmo

   private
   public :: enfpcmo, qopcmo, spcmo, vpcmo
   save
 end module pcmo

 module periodic
   !> Arguments: cutg, cutg_big, cutg_sim, cutg_sim_big, cutr, cutr_sim, glatt, glatt_inv, glatt_sim, gnorm, gnorm_sim, gvec, gvec_sim, igmult, igmult_sim, igvec, igvec_sim, ireal_imag, isrange, k_inv, kvec, nband, ncoef_per, ng1d, ng1d_sim, ngnorm, ngnorm_big, ngnorm_orb, ngnorm_sim, ngnorm_sim_big, ngvec, ngvec_big, ngvec_orb, ngvec_sim, ngvec_sim_big, nkvec, np, npoly, rknorm, rkvec, rkvec_shift, rlatt, rlatt_inv, rlatt_sim, rlatt_sim_inv, vcell, vcell_sim, znuc2_sum, znuc_sum
   use precision_kinds, only: dp
   include 'vmc.h'
   include 'ewald.h'

   real(dp) :: cutg
   real(dp) :: cutg_big
   real(dp) :: cutg_sim
   real(dp) :: cutg_sim_big
   real(dp) :: cutr
   real(dp) :: cutr_sim
   real(dp) :: glatt(3,3)
   real(dp) :: glatt_inv(3,3)
   real(dp) :: glatt_sim(3,3)
   real(dp) :: gnorm(NGNORM_BIGX)
   real(dp) :: gnorm_sim(NGNORM_SIM_BIGX)
   real(dp) :: gvec(3,NGVEC_BIGX)
   real(dp) :: gvec_sim(3,NGVEC_SIM_BIGX)
   integer  :: igmult(NGNORM_BIGX)
   integer  :: igmult_sim(NGNORM_SIM_BIGX)
   integer  :: igvec(3,NGVEC_BIGX)
   integer  :: igvec_sim(3,NGVEC_SIM_BIGX)
   integer  :: ireal_imag(MORB)
   integer  :: isrange
   integer  :: k_inv(IVOL_RATIO)
   integer  :: kvec(3,IVOL_RATIO)
   integer  :: nband(IVOL_RATIO)
   integer  :: ncoef_per
   integer  :: ng1d(3)
   integer  :: ng1d_sim(3)
   integer  :: ngnorm
   integer  :: ngnorm_big
   integer  :: ngnorm_orb
   integer  :: ngnorm_sim
   integer  :: ngnorm_sim_big
   integer  :: ngvec
   integer  :: ngvec_big
   integer  :: ngvec_orb
   integer  :: ngvec_sim
   integer  :: ngvec_sim_big
   integer  :: nkvec
   integer  :: np
   integer  :: npoly
   real(dp) :: rknorm(IVOL_RATIO)
   real(dp) :: rkvec(3,IVOL_RATIO)
   real(dp) :: rkvec_shift(3)
   real(dp) :: rlatt(3,3)
   real(dp) :: rlatt_inv(3,3)
   real(dp) :: rlatt_sim(3,3)
   real(dp) :: rlatt_sim_inv(3,3)
   real(dp) :: vcell
   real(dp) :: vcell_sim
   real(dp) :: znuc2_sum
   real(dp) :: znuc_sum

   private
   public :: cutg, cutg_big, cutg_sim, cutg_sim_big, cutr, cutr_sim, glatt, glatt_inv
   public :: glatt_sim, gnorm, gnorm_sim, gvec, gvec_sim, igmult, igmult_sim, igvec
   public :: igvec_sim, ireal_imag, isrange, k_inv, kvec, nband, ncoef_per, ng1d, ng1d_sim
   public :: ngnorm, ngnorm_big, ngnorm_orb, ngnorm_sim, ngnorm_sim_big, ngvec, ngvec_big
   public :: ngvec_orb, ngvec_sim, ngvec_sim_big, nkvec, np, npoly, rknorm, rkvec, rkvec_shift
   public :: rlatt, rlatt_inv, rlatt_sim, rlatt_sim_inv, vcell, vcell_sim, znuc2_sum, znuc_sum
   save
 end module periodic

 module phifun
   !> Arguments: d2phin, d2phin_all, d3phin, dphin, n0_ibasis, n0_ic, n0_nbasis, phin
   use precision_kinds, only: dp
   include 'vmc.h'

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

 module pseudo
   !> Arguments: lpot, nloc, vps, vpso
   use precision_kinds, only: dp
   include 'vmc.h'
   include 'pseudo.h'
   include 'force.h'

   integer  :: lpot(MCTYPE)
   integer  :: nloc
   real(dp) :: vps(MELEC,MCENT,MPS_L)
   real(dp) :: vpso(MELEC,MCENT,MPS_L,MFORCE)

   private 
   public :: lpot, nloc, vps, vpso 
   save
end module pseudo

 module pseudo_champ
   !> Arguments: igrid_ps, rmax_coul, rmax_nloc
   use precision_kinds, only: dp
   include 'vmc.h'

   integer  :: igrid_ps(MCTYPE)
   real(dp) :: rmax_coul(MCTYPE)
   real(dp) :: rmax_nloc(MCTYPE)

   private
   public :: igrid_ps, rmax_coul, rmax_nloc
   save
 end module pseudo_champ

 module pseudo_fahy
   !> Arguments: drad, dradl, nlrad, npotl, potl, ptnlc, rcmax
   use precision_kinds, only: dp
   include 'vmc.h'
   include 'pseudo.h'

   real(dp) :: drad(MCTYPE)
   real(dp) :: dradl(MCTYPE)
   integer  :: nlrad(MCTYPE)
   integer  :: npotl(MCTYPE)
   real(dp) :: potl(MPS_GRID,MCTYPE)
   real(dp) :: ptnlc(MPS_GRID,MCTYPE,MPS_L)
   real(dp) :: rcmax(MCTYPE)

   private
   public :: drad, dradl, nlrad, npotl, potl, ptnlc, rcmax
   save
 end module pseudo_fahy

 module pseudo_tm
   !> Arguments: arg, arg_ps, d2pot, nr_ps, r0, r0_ps, rmax, rmax_ps, vpseudo
   use precision_kinds, only: dp
   include 'vmc.h'
   include 'pseudo.h'

    real(dp) :: arg(MCTYPE)
    real(dp) :: arg_ps(MCTYPE)
    real(dp) :: d2pot(MPS_GRID,MCTYPE,MPS_L)
    integer  :: nr_ps(MCTYPE)
    real(dp) :: r0(MCTYPE)
    real(dp) :: r0_ps(MCTYPE)
    real(dp) :: rmax(MCTYPE)
    real(dp) :: rmax_ps(MCTYPE)
    real(dp) :: vpseudo(MPS_GRID,MCTYPE,MPS_L)

    private
    public :: arg, arg_ps, d2pot, nr_ps, r0, r0_ps, rmax, rmax_ps, vpseudo
    save
 end module pseudo_tm

 module pworbital
   !> Arguments: c_im, c_ip, c_rm, c_rp, icmplx, isortg, isortk, ngorb
   use precision_kinds, only: dp
   include 'vmc.h'
   include 'ewald.h'

   real(dp) :: c_im(NGVECX,MORB)
   real(dp) :: c_ip(NGVECX,MORB)
   real(dp) :: c_rm(NGVECX,MORB)
   real(dp) :: c_rp(NGVECX,MORB)
   integer  :: icmplx
   integer  :: isortg(NGVECX,MORB)
   integer  :: isortk(IVOL_RATIO)
   integer  :: ngorb(MORB)

   private
   public :: c_im, c_ip, c_rm, c_rp, icmplx, isortg, isortk, ngorb
   save
 end module pworbital

 module rlobxy
   !> Arguments: rlobx, rloby, rloby2
   use precision_kinds, only: dp
   include 'vmc.h'

    real(dp) :: rlobx(NSPLIN)
    real(dp) :: rloby(NSPLIN)
    real(dp) :: rloby2(NSPLIN)

    private
    public :: rlobx, rloby, rloby2
    save
 end module rlobxy

 module sa_check
   !> Arguments: energy_all, energy_err_all
   use precision_kinds, only: dp
   include 'mstates.h'

   real(dp) :: energy_all(MSTATES)
   real(dp) :: energy_err_all(MSTATES)

   private
   public :: energy_all, energy_err_all
   save
 end module sa_check

 module sa_weights
   !> Arguments: iweight, nweight, weights
   use precision_kinds, only: dp
   include 'mstates.h'

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
   include 'vmc.h'

   real(dp) :: denergy_det(MDET,2)
   real(dp) :: dtildem(MELEC,MORB,2)

   private
   public :: denergy_det, dtildem
   save
 end module scratch

 module slatn
   !> Arguments: slmin
   use precision_kinds, only: dp
   include 'vmc.h'

   real(dp) :: slmin(MMAT_DIM)

   private
   public :: slmin
   save
 end module slatn

 module sr_index
   !> Arguments: jelo, jelo2, jelohfj

   integer  :: jelo
   integer  :: jelo2
   integer  :: jelohfj

   private
   public :: jelo, jelo2, jelohfj
   save
 end module sr_index

 module sr_mat_n
   !> Arguments: elocal, h_sr, jefj, jfj, jhfj, nconf, obs, s_diag, s_ii_inv, sr_ho, sr_o, wtg, obs_tot
   use precision_kinds, only: dp
   include 'sr.h'
   include 'mstates.h'

   real(dp) :: elocal(MCONF,MSTATES)
   real(dp) :: h_sr(MPARM)
   integer  :: jefj
   integer  :: jfj
   integer  :: jhfj
   integer  :: nconf
   real(dp) :: obs(MOBS,MSTATES)
   real(dp) :: s_diag(MPARM,MSTATES)
   real(dp) :: s_ii_inv(MPARM)
   real(dp) :: sr_ho(MPARM,MCONF)
   real(dp) :: sr_o(MPARM,MCONF)
   real(dp) :: wtg(MCONF,MSTATES)
   real(dp) :: obs_tot(MOBS,MSTATES)

   private
   public :: elocal, h_sr, jefj, jfj, jhfj, nconf, obs, s_diag, s_ii_inv, sr_ho, sr_o, wtg, obs_tot
   save
 end module sr_mat_n

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
   include 'vmc.h'

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
   use precision_kinds, only: dp
   include 'ewald.h'

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
   use precision_kinds, only: dp
   include 'ewald.h'

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
   include 'vmc.h'

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
   include 'vmc.h'

   real(dp) :: vj(3,MELEC)
   real(dp) :: vjn(3,MELEC)

   private 
   public :: vj, vjn 
   save
 end module velocity_jastrow

 module wfsec
   !> Arguments: iwf, iwftype, nwftype
   use precision_kinds, only: dp
   include 'force.h'

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
   include 'vmc.h'
   include 'mstates.h'

   real(dp) :: dymat(MORB,MELEC,2,MSTATES)
   real(dp) :: ymat(MORB,MELEC,2,MSTATES)

   private
   public :: dymat, ymat
   save
 end module ycompact

 module ycompactn
   !> Arguments: ymatn
   use precision_kinds, only: dp
   include 'vmc.h'
   include 'mstates.h'

   real(dp) :: ymatn(MORB,MELEC,MSTATES)

   private
   public :: ymatn
   save
 end module ycompactn

 module zcompact
   !> Arguments: aaz, dzmat, emz, zmat
   use precision_kinds, only: dp
   include 'vmc.h'
   include 'mstates.h'

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
   include 'vmc.h'

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
   include 'vmc.h'

   real(dp) :: transform_grd(MCENT3,MCENT3)

   private 
   public :: transform_grd 
   save
 end module zmatrix_grad
 
 
