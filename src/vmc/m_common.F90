!> \brief File collecting all modules that replace common blocks.  
!>
!> \author P. Lopez-Tarifa NLeSC(2019)
! Victor's suggestion
!module vmc
!  include 'vmc.h'
!end module

module atom
!> Arguments: znuc, cent, pecent, iwctype, nctype, ncent
   include 'vmc.h'
   real*8 :: cent( 3, MCENT) 
   real*8 :: znuc( MCTYPE)
   integer:: pecent
   integer:: iwctype( MCENT), nctype, ncent
  SAVE
end module 

module contrl_per
!> Arguments: iperiodic, ibasis 
  integer:: iperiodic, ibasis
  SAVE
end module contrl_per

module da_pseudo
!> Arguments: da_pecent, da_vps, da_nonloc  
  include 'vmc.h'
  include 'pseudo.h'
  real*8 :: da_pecent( 3, MCENT), da_vps( 3, MELEC, MCENT, MPS_L)
  real*8 :: da_nonloc( 3, MCENT)= 0.0D0 
 SAVE
end module 

module force_analy 
!> Arguments: iforce_analy 
  integer :: iforce_analy 
  save
end module 
