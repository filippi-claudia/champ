      subroutine determinant_psig(psid,psig)

      implicit real*8(a-h,o-z)
      include 'vmc.h'
      include 'force.h'
      include 'optci.h'
      include 'mstates.h'

      common /const/ pi,hb,etrial,delta,deltai,fbias,nelec,imetro,ipr
      common /contrl_per/ iperiodic,ibasis

      common /csfs/ ccsf(MDET,MSTATES,MWF),cxdet(MDET*MDETCSFX)
     &,icxdet(MDET*MDETCSFX),iadet(MDET),ibdet(MDET),ncsf,nstates

      dimension psid(*)

      psig=0
      do 200 i=1,nstates
        istate=iweight_g(i)
  200   psig=psig+weights_g(i)*psid(istate)*psid(istate)

      psig=dsqrt(psig)

      return
      end
