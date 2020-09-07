      subroutine determinant_psig(psid,psig)

      use csfs, only: nstates
      use mstates_mod, only: MSTATES, MDETCSFX

      use mstates3, only: iweight_g, weights_g

      implicit real*8(a-h,o-z)



      include 'vmc.h'
      include 'force.h'
      include 'optci.h'



      dimension psid(*)

      psig=0
      do 200 i=1,nstates
        istate=iweight_g(i)
  200   psig=psig+weights_g(i)*psid(istate)*psid(istate)

      psig=dsqrt(psig)

      return
      end
