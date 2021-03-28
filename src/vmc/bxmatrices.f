      subroutine bxmatrix(kref,xmatu,xmatd,b)

      use vmc_mod, only: MELEC, MORB
      use vmc_mod, only: MMAT_DIM
      use elec, only: ndn, nup
      use dorb_m, only: iworbd
      use slater, only: d2dx2, ddx, fp, fpp, slmi

      implicit real*8(a-h,o-z)

      dimension b(MORB,MELEC),btemp(MELEC**2,2),xmatu(MELEC**2),xmatd(MELEC**2),work(MELEC)

c     RLPB (here we need a loop)
      kstate=1

      do iab=1,2
         if(iab.eq.1) then
            iel=0
            nel=nup
         else
            iel=nup
            nel=ndn
         endif
         ish=-nel
         do i=1,nel
            ish=ish+nel
            do j=1,nel
               btemp(j+ish,iab)=b(iworbd(j+iel,kref),i+iel)
            enddo
         enddo
      enddo

      call multiply_slmi_mderiv_simple(nup,btemp(1,1),work,slmi(1,kstate,1),xmatu)
      call multiply_slmi_mderiv_simple(ndn,btemp(1,2),work,slmi(1,kstate,2),xmatd)

      end subroutine
