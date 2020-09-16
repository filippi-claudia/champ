 module force_mod
  integer, parameter :: MFORCE=3
  integer, parameter :: MFORCE_WT_PRD=1000
  integer, parameter :: MWF=3

  private 
  public :: MFORCE, MFORCE_WT_PRD, MWF
  save
 end module force_mod

  module force_analy 
   !> Arguments: iforce_analy, iuse_zmat, alfgeo 
   use precision_kinds, only: dp

   integer  :: iforce_analy 
   integer  :: iuse_zmat
   real(dp) :: alfgeo

   private
   public   :: iforce_analy, iuse_zmat, alfgeo 
   save
 end module force_analy 

  module forcest
   !> Arguments: fcm2, fcum
   use force_mod, only: MFORCE
   use precision_kinds, only: dp
   use mstates_mod, only: MSTATES

   real(dp) :: fcm2(MSTATES,MFORCE)
   real(dp) :: fcum(MSTATES,MFORCE)

   private
   public   ::  fcm2, fcum
   save
 end module forcest

  module forcestr
   !> Arguments: delc
   use force_mod, only: MFORCE
   use precision_kinds, only: dp
   use vmc, only: MCENT
 
   real(dp) :: delc(3,MCENT,MFORCE)
 
   private
   public   ::  delc 
   save
 end module forcestr

  module forcewt
   !> Arguments: wcum, wsum
   use force_mod, only: MFORCE
   use precision_kinds, only: dp
   use mstates_mod, only: MSTATES

   real(dp) :: wcum(MSTATES,MFORCE)
   real(dp) :: wsum(MSTATES,MFORCE)

   private
   public   ::  wcum, wsum 
   save
end module forcewt

module force_dmc
   !> Arguments: itausec, nwprod

   integer  :: itausec
   integer  :: nwprod

   private
   public   ::   itausec, nwprod
   save
 end module force_dmc

 module force_fin
   !> Arguments: da_energy_ave, da_energy_err
   use precision_kinds, only: dp
   use vmc, only: MCENT

   real(dp) :: da_energy_ave(3,MCENT)
   real(dp) :: da_energy_err(3)

   private
   public   :: da_energy_ave, da_energy_err
   save
 end module force_fin

  module force_mat_n
   !> Arguments: force_o
   use sr_mod, only: MCONF
   use precision_kinds, only: dp
   use vmc, only: MCENT

   real(dp) :: force_o(6*MCENT,MCONF)

   private
   public   ::  force_o 
   save
 end module force_mat_n  

 module forcepar
   !> Arguments: deltot, istrech, nforce
   use force_mod, only: MFORCE
   use precision_kinds, only: dp

   real(dp) :: deltot(MFORCE)
   integer  :: istrech
   integer  :: nforce

   private
   public   ::  deltot, istrech, nforce 
   save
 end module forcepar