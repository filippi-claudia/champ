      subroutine bxmatrix(kref,xmatu,xmatd,b)

      use vmc_mod, only: MELEC, MORB
      use vmc_mod, only: MMAT_DIM
      use elec, only: ndn, nup
      use dorb_m, only: iworbd
      use coefs, only: norb
      use slater, only: d2dx2, ddx, fp, fpp, slmi
      use const, only: nelec

      implicit real*8(a-h,o-z)




      dimension b(MORB,nelec),btemp(nelec**2,2),xmatu(nelec**2),xmatd(nelec**2),work(nelec)

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
