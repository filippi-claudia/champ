c-----------------------------------------------------------------------

c compute gauss-pseudopotential for electron iel
      subroutine getvps_gauss(rvec_en,r_en,iel)

      use vmc_mod, only: MELEC, MCENT
      use atom, only: znuc, iwctype, ncent, ncent_tot
      use const, only: nelec
      use pseudo, only: lpot, vps

      use da_pseudo, only: da_vps

      use precision_kinds, only: dp
      implicit none

      integer :: ic, ict, iel, k, l
      real(dp) :: dvpot, r, reni, ri, ri2
      real(dp) :: vpot
      real(dp), dimension(nelec,ncent_tot) :: r_en
      real(dp), dimension(3,nelec,ncent_tot) :: rvec_en






      do 10 ic=1,ncent
        ict=iwctype(ic)

        r=max(1.0d-10,r_en(iel,ic))
        ri=1.d0/r
        ri2=ri*ri
c local potential
        vpot=0.d0
        dvpot=0.d0
        call gauss_pot(r,lpot(ict),ict,vpot,dvpot)
        vpot=vpot-znuc(ict)*ri
        vps(iel,ic,lpot(ict))=vpot
        do 5 k=1,3
          reni=rvec_en(k,iel,ic)*ri
    5     da_vps(k,iel,ic,lpot(ict))=-(dvpot+znuc(ict)*ri2)*reni
c non-local pseudopotential
        do 10 l=1,lpot(ict)-1
         vpot=0.d0
         dvpot=0.d0
         call gauss_pot(r,l,ict,vpot,dvpot)
         vps(iel,ic,l)=vpot
         do 10 k=1,3
          reni=rvec_en(k,iel,ic)*ri
          da_vps(k,iel,ic,l)=-dvpot*reni
   10 continue

c     do ic=1,ncent
c       write(6,*) 'HELLO_GAUSS',da_vps(1,iel,ic,lpot(iwctype(ic)))
c     enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine gauss_pot(r,l,ict,vpot,dvpot)
      use gauss_ecp, only: ecp_coef, ecp_exponent, necp_power, necp_term

      use precision_kinds, only: dp
      implicit none

      integer :: i, ict, l
      real(dp) :: dp_val, dv, dvpot, e, p
      real(dp) :: r, rsq, v, vpot


      v = 0.d0
      dv = 0.d0
      rsq=r**2

      do i=1, necp_term(l,ict)
       if(necp_power(i,l,ict).ne.2)then
        p = r**(necp_power(i,l,ict)-2)
        dp_val = (necp_power(i,l,ict)-2)*p/r
       else
        p = 1.d0
        dp_val= 0.d0
       endif
       if (ecp_exponent(i,l,ict)*rsq.gt.7.0D2) then
          e = 0.0d0
       else
          e = ecp_coef(i,l,ict)*exp(-ecp_exponent(i,l,ict)*rsq)
       endif
       v = v + p*e
       dv = dv + (dp_val -2*p*ecp_exponent(i,l,ict)*r)*e
      enddo

      vpot=v
      dvpot=dv
      return
      end
