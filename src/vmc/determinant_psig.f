      module determinant_psig_mod
      contains
      subroutine determinant_psig(psid,psij,psig)

      use csfs, only: nstates, anormo

      use mstates3, only: iweight_g, weights_g

      use precision_kinds, only: dp
      implicit none

      integer :: i, istate
      real(dp) :: psig
      real(dp), dimension(*) :: psid
      real(dp), dimension(*) :: psij
    
      psig=0
      do i=1,nstates !STU state map psij
        istate=iweight_g(i)
        psig=psig+weights_g(i)*psid(istate)*psid(istate)*exp(2*psij(istate))/1.0d0!anormo(istate)
      enddo

      psig=dsqrt(psig)

      return
      end
      end module
