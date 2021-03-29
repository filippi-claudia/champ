      subroutine bxmatrix(kref,xmatu,xmatd,b,istate)

      use vmc_mod, only: MELEC, MORB
      use vmc_mod, only: MMAT_DIM
      use elec, only: ndn, nup
      use dorb_m, only: iworbd
      use slater, only: d2dx2, ddx, fp, fpp, slmi

      implicit real*8(a-h,o-z)

      dimension b(MORB,MELEC),btemp(MELEC**2,2),xmatu(MELEC**2),xmatd(MELEC**2),work(MELEC)

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

      call multiply_slmi_mderiv_simple(nup,btemp(1,1),work,slmi(1,1,istate),xmatu)
      call multiply_slmi_mderiv_simple(ndn,btemp(1,2),work,slmi(1,2,istate),xmatd)

      end subroutine
