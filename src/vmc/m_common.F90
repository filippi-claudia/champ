!> \brief File collecting all modules that replace common blocks.  
!>
!> \author P. Lopez-Tarifa NLeSC(2019)
MODULE contrl_per
  INTEGER:: iperiodic,ibasis
  SAVE
END MODULE contrl_per

MODULE force_analy 
  INTEGER:: iforce_analy 
  SAVE
END MODULE 

MODULE da_pseudo
   include 'vmc.h'
   include 'pseudo.h'
   real*8 :: da_pecent(3,MCENT), da_vps(3,MELEC,MCENT,MPS_L)
   real*8 :: da_nonloc(3,MCENT)=0.0D0 
  SAVE
END MODULE 

MODULE atom
   include 'vmc.h'
   integer znuc(MCTYPE),cent(3,MCENT), pecent
   integer iwctype(MCENT),nctype,ncent
  SAVE
END MODULE 