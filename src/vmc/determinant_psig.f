      subroutine determinant_psig(psid,psig)

      use csfs, only: nstates

      use mstates3, only: iweight_g, weights_g

      use precision_kinds, only: dp
      implicit none

      integer :: i, istate
      real(dp) :: psig
      real(dp), dimension(*) :: psid







      psig=0
      do 200 i=1,nstates
        istate=iweight_g(i)
  200   psig=psig+weights_g(i)*psid(istate)*psid(istate)

      psig=dsqrt(psig)

      return
      end
