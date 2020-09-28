 module da_energy_sumcum
   !> Arguments: da_energy_cm2, da_energy_cum, da_energy_sum, da_psi_cum, da_psi_sum
   use precision_kinds, only: dp
   use vmc_mod, only: MCENT

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
   use vmc_mod, only: MELEC, MCENT

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
   use vmc_mod, only: MELEC, MORB, MCENT

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
   use vmc_mod, only: MELEC, MCENT


   real(dp) :: da_pecent( 3, MCENT), da_vps( 3, MELEC, MCENT, MPS_L)
   real(dp) :: da_nonloc( 3, MCENT)= 0.0D0 

   private
   public   :: da_pecent, da_vps, da_nonloc 
   save
 end module da_pseudo 
 
 module da_energy_now
   !> Arguments: da_energy, da_psi
   use precision_kinds, only: dp
   use vmc_mod, only: MCENT
 
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
   use vmc_mod, only: MDET

    real(dp) :: denergy_det(MDET,2)

    private 
    public :: denergy_det 
    save
 end module denergy_det_m

 module denupdn
   !> Arguments: rprobdn, rprobup
   use precision_kinds, only: dp
   use vmc_mod, only: nrad

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
   use vmc_mod, only: MELEC

   real(dp) :: d2g(MPARMJ)
   real(dp) :: g(3,MELEC,MPARMJ)
   real(dp) :: go(MELEC,MELEC,MPARMJ)
   real(dp) :: gvalue(MPARMJ)

   private
   public   :: d2g, g, go, gvalue
   save
 end module derivjas

  module dorb_m
   !> Arguments: iworbd
   use vmc_mod, only: MELEC, MDET

   integer  :: iworbd(MELEC,MDET)

   private 
   public :: iworbd 
   save
 end module dorb_m

  module ijasnonlin
   !> Arguments: d1d2a, d1d2b, d2d2a, d2d2b
   use precision_kinds, only: dp
   use vmc_mod, only: MCTYPE

   real(dp) :: d1d2a(MCTYPE)
   real(dp) :: d1d2b(2)
   real(dp) :: d2d2a(MCTYPE)
   real(dp) :: d2d2b(2)

   private 
   public :: d1d2a, d1d2b, d2d2a, d2d2b 
   save
 end module ijasnonlin