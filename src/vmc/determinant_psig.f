      subroutine determinant_psig(psid,psig)

      use const, only: pi, hb, etrial, delta, deltai, fbias, nelec, imetro, ipr
      use csfs, only: ccsf, cxdet, iadet, ibdet, icxdet, ncsf, nstates
      use contrl_per, only: iperiodic, ibasis

      implicit real*8(a-h,o-z)


      include 'vmc.h'
      include 'force.h'
      include 'optci.h'
      include 'mstates.h'



      dimension psid(*)

      psig=0
      do 200 i=1,nstates
        istate=iweight_g(i)
  200   psig=psig+weights_g(i)*psid(istate)*psid(istate)

      psig=dsqrt(psig)

      return
      end
