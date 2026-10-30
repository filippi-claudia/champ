module determinant_psig_mod
contains
      subroutine determinant_psig(psid,psij,psi2g)

      use csfs, only: nstates, anormo
      use mstates3, only: iweight_g, weights_g
      use vmc_mod, only: stoj
      use precision_kinds, only: dp
      use contrl_file, only: ounit
      implicit none

      integer :: i, ijas1, istate
      real(dp) :: psi2g, ratio
      real(dp), dimension(*) :: psid
      real(dp), dimension(*) :: psij

! psig computed with respect to state 1
      ijas1=stoj(iweight_g(1))

      psi2g=weights_g(1)/anormo(1)
      do i=2,nstates
        istate=iweight_g(i)
        ratio=psid(istate)/psid(1)*exp(psij(stoj(istate))-psij(ijas1))
        psi2g=psi2g+weights_g(i)/anormo(i)*ratio*ratio
      enddo

      return
      end
end module
