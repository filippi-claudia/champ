module mmpol_mod
  !> Arguments MCHMM, mmpolfile_sites, mmpolfile_chmm
  
  integer, parameter :: MCHMM=1
  character*80 mmpolfile_sites
  character*80 mmpolfile_chmm
  private 
  public :: MCHMM, mmpolfile_sites, mmpolfile_chmm
  save
 end module mmpol_mod

 module mmpol_cntrl
   !> Arguments: isites_mmpol, immpolprt, icall_mm, ich_mmpol, immpol

   integer  :: icall_mm
   integer  :: ich_mmpol
   integer  :: immpol
   integer  :: immpolprt
   integer  :: isites_mmpol

   private
   public :: isites_mmpol, immpolprt, icall_mm, ich_mmpol, immpol
   save
 end module mmpol_cntrl

 module mmpol_dipol
   !> Arguments: dipo, alfa
   use mmpol_mod, only: MCHMM
   use precision_kinds, only: dp

   real(dp) :: alfa(MCHMM)
   real(dp) :: dipo(3,MCHMM)

   private
   public :: dipo, alfa
   save
 end module mmpol_dipol

 module mmpol_hpsi
   !> Arguments: eek_pol, peQMdp, peQMq
   use mmpol_mod, only: MCHMM
   use precision_kinds, only: dp

   real(dp) :: eek_pol(3,MCHMM)
   real(dp) :: peQMdp
   real(dp) :: peQMq

   private 
   public :: eek_pol, peQMdp, peQMq
   save
 end module mmpol_hpsi

 module mmpol_parms
   !> Arguments: x_mmpol, nchmm, chmm, rqq
   use mmpol_mod, only: MCHMM
   use precision_kinds, only: dp
 
   real(dp) :: chmm(MCHMM)
   integer  :: nchmm
   real(dp) :: rqq(MCHMM,MCHMM)
   real(dp) :: x_mmpol(3,MCHMM)

   private
   public :: x_mmpol, nchmm, chmm, rqq
   save
 end module mmpol_parms
 
 module mmpolo
   !> Arguments: cmmpolo, dmmpolo, eeko
   use mmpol_mod, only: MCHMM
   use precision_kinds, only: dp

   real(dp) :: cmmpolo
   real(dp) :: dmmpolo
   real(dp) :: eeko(3,MCHMM)

   private 
   public :: cmmpolo, dmmpolo, eeko
   save
 end module mmpolo

 module mmpol_ahpol
   !> Arguments: ah_pol, bh_pol
   use mmpol_mod, only: MCHMM
   use precision_kinds, only: dp

   real(dp) :: ah_pol(3*MCHMM,3*MCHMM)
   real(dp) :: bh_pol(3*MCHMM)

   private
   public :: ah_pol, bh_pol
   save
 end module mmpol_ahpol

 module mmpol_averages
   !> Arguments: cmmpol_cum, cmmpol_cm2, eek2_cum, dmmpol_sum, eek1_cm2, eek_sum, eek2_cm2, cmmpol_sum, dmmpol_cum, dmmpol_cm2, eek3_cum, eek1_cum, eek3_cm2
   use mmpol_mod, only: MCHMM
   use precision_kinds, only: dp
 
    real(dp) :: cmmpol_cm2
    real(dp) :: cmmpol_cum
    real(dp) :: cmmpol_sum
    real(dp) :: dmmpol_cm2
    real(dp) :: dmmpol_cum
    real(dp) :: dmmpol_sum
    real(dp) :: eek1_cm2(MCHMM)
    real(dp) :: eek1_cum(MCHMM)
    real(dp) :: eek2_cm2(MCHMM)
    real(dp) :: eek2_cum(MCHMM)
    real(dp) :: eek3_cm2(MCHMM)
    real(dp) :: eek3_cum(MCHMM)
    real(dp) :: eek_sum(3,MCHMM)
 
    private
    public :: cmmpol_cum, cmmpol_cm2, eek2_cum, dmmpol_sum, eek1_cm2, eek_sum
    public :: eek2_cm2, cmmpol_sum, dmmpol_cum, dmmpol_cm2, eek3_cum, eek1_cum, eek3_cm2
    save
 end module mmpol_averages

 module mmpol_fdc
   !> Arguments: a_cutoff, rcolm, screen1, screen2
   use mmpol_mod, only: MCHMM
   use precision_kinds, only: dp
 
   real(dp) :: a_cutoff
   real(dp) :: rcolm
   real(dp) :: screen1(MCHMM,MCHMM)
   real(dp) :: screen2(MCHMM,MCHMM)
 
   private
   public :: a_cutoff, rcolm, screen1, screen2
   save
 end module mmpol_fdc
 
 module mmpol_field
   !> Arguments: eqk_pol, enk_pol
   use mmpol_mod, only: MCHMM
   use precision_kinds, only: dp
 
   real(dp) :: enk_pol(3,MCHMM)
   real(dp) :: eqk_pol(3,MCHMM)
 
   private
   public :: eqk_pol, enk_pol
   save
 end module mmpol_field

 module mmpol_inds
   use mmpol_mod, only: MCHMM
   !> Arguments: inds_pol
 
   integer  :: inds_pol(MCHMM)
 
   private
   public :: inds_pol
   save
 end module mmpol_inds

 module mmpol_pot
   !> Arguments: peqq, pepol_dp, pepol_q, penu_q, peq_dp, penu_dp, u_dd, u_self
   use precision_kinds, only: dp

   real(dp) :: penu_dp
   real(dp) :: penu_q
   real(dp) :: pepol_dp
   real(dp) :: pepol_q
   real(dp) :: peq_dp
   real(dp) :: peqq
   real(dp) :: u_dd
   real(dp) :: u_self

   private
   public :: peqq, pepol_dp, pepol_q, penu_q, peq_dp, penu_dp, u_dd, u_self
   save
 end module mmpol_pot