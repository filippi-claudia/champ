      module determinant_psig_mod
      contains
      subroutine determinant_psig(psid,psig)

      use csfs, only: nstates

      use mstates3, only: iweight_g, weights_g

      use precision_kinds, only: dp
      implicit none

      integer :: i, istate
      real(dp) :: psig
      real(dp), dimension(*) :: psid







      psig=0
      do i=1,nstates
        istate=iweight_g(i)
        psig=psig+weights_g(i)*psid(istate)*psid(istate)/anormo(istate)
c missing a exp(2*psij(istate)) term here compared to ortho
      enddo

      psig=dsqrt(psig)

      return
      end
      end module
