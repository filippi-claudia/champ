      subroutine bxmatrix(kref,xmatu,xmatd,b)

      use const, only: pi, hb, etrial, delta, deltai, fbias, nelec, imetro, ipr
      use dets, only: cdet, ndet
      use elec, only: ndn, nup
      use coefs, only: coef, nbasis, norb
      use dorb_m, only: iworbd
      implicit real*8(a-h,o-z)

      include 'vmc.h'
      include 'force.h'
      include 'mstates.h'
      include 'pseudo.h'

      common /slater/ slmi(MMAT_DIM,2)
     &,fp(3,MMAT_DIM,2)
     &,fpp(MMAT_DIM,2)
     &,ddx(3,MELEC),d2dx2(MELEC)

      dimension b(MORB,MELEC),btemp(MELEC**2,2),xmatu(MELEC**2),xmatd(MELEC**2),work(MELEC)

      do 110 iab=1,2
        if(iab.eq.1) then
          iel=0
          nel=nup
         else
          iel=nup
          nel=ndn
        endif
        ish=-nel
        do 110 i=1,nel
          ish=ish+nel
          do 110 j=1,nel
  110       btemp(j+ish,iab)=b(iworbd(j+iel,kref),i+iel)

      call multiply_slmi_mderiv_simple(nup,btemp(1,1),work,slmi(1,1),xmatu)
      call multiply_slmi_mderiv_simple(ndn,btemp(1,2),work,slmi(1,2),xmatd)

      return
      end
c-----------------------------------------------------------------------
