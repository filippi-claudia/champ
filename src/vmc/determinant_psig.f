      subroutine determinant_psig(psid,psig)

      use csfs, only: nstates
      use config, only: anormo
      use mstates3, only: iweight_g, weights_g

      implicit real*8(a-h,o-z)

      dimension psid(*)

      psig=0.0d0
      do i=1,nstates
         istate=iweight_g(i)
         psig=psig+weights_g(i)*psid(istate)*psid(istate)/anormo(istate)
      enddo

      psig=dsqrt(psig)

      end subroutine
