      subroutine bxmatrix(kref,xmatu,xmatd,b)

      use const, only: pi, hb, etrial, delta, deltai, fbias, nelec, imetro, ipr
      use dets, only: cdet, ndet
      use elec, only: ndn, nup
      implicit real*8(a-h,o-z)



      include 'vmc.h'
      include 'force.h'
      include 'mstates.h'
      include 'pseudo.h'

      common /coefs/ coef(MBASIS,MORB,MWF),nbasis,norb
      common /slater/ slmui(MMAT_DIM),slmdi(MMAT_DIM)
     &,fpu(3,MMAT_DIM),fpd(3,MMAT_DIM)
     &,fppu(MMAT_DIM),fppd(MMAT_DIM)
     &,ddx(3,MELEC),d2dx2(MELEC)
      common /multislater/ detu(MDET),detd(MDET)


      common /dorb/ iworbd(MELEC,MDET)

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

      call multiply_slmi_mderiv_simple(nup,btemp(1,1),work,slmui(1),xmatu)
      call multiply_slmi_mderiv_simple(ndn,btemp(1,2),work,slmdi(1),xmatd)

      return
      end
c-----------------------------------------------------------------------
