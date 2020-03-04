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

 module ghostatom
   !> Arguments: newghostype, nghostcent
   use precision_kinds, only: dp
   include 'vmc.h'

   integer  :: newghostype
   integer  :: nghostcent

   private
   public   :: newghostype, nghostcent
   save
 end module ghostatom

 module jaspar
   !> Arguments: nspin1, nspin2, sspin, sspinn, is
   use precision_kinds, only: dp
   include 'vmc.h'

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


