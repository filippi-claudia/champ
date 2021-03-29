      subroutine detsav(iel,iflag)
c     Written by Claudia Filippi, modified by RLPB
      use force_mod, only: MFORCE, MFORCE_WT_PRD, MWF
      use vmc_mod, only: MELEC, MORB, MBASIS, MDET, MCENT, MCTYPE, MCTYP3X
      use vmc_mod, only: NSPLIN, nrad, MORDJ, MORDJ1, MMAT_DIM, MMAT_DIM2, MMAT_DIM20
      use vmc_mod, only: radmax, delri
      use vmc_mod, only: NEQSX, MTERMS
      use vmc_mod, only: MCENT3, NCOEF, MEXCIT
      use csfs, only: nstates
      use dets, only: ndet
      use elec, only: ndn, nup
      use multidet, only: ivirt, kref, numrep_det
      use slatn, only: slmin
      use ycompact, only: ymat
      use ycompactn, only: ymatn
      use coefs, only: norb
      use dorb_m, only: iworbd
      use multimat, only: aa, wfmat
      use multimatn, only: aan, wfmatn
      use multislatern, only: ddorbn, detn, dorbn, orbn
      use orbval, only: ddorb, dorb, nadorb, ndetorb, orb
      use slater, only: d2dx2, ddx, fp, fpp, slmi
      use multislater, only: detiab

      implicit real*8(a-h,o-z)

      if(iel.le.nup) then
         iab=1
         nel=nup
         ish=0
      else
         iab=2
         nel=ndn
         ish=nup
      endif

      ikel=nel*(iel-ish-1)
      do istate=1,nstates   
         do j=1,nel*nel
            slmi(j,iab,istate)=slmin(j,istate)
         enddo
         do j=ivirt(iab),norb
            do i=1,nel
               ymat(j,i,iab,istate)=ymatn(j,i,istate)
               aa(i,j,iab,istate)=aan(i,j,istate)
            enddo
         enddo

         do k=1,ndet
            if(k.eq.kref) cycle
            ndim=numrep_det(k,iab)
            do i=1,ndim*ndim
               wfmat(i,k,iab,istate)=wfmatn(i,k,istate)
            enddo
         enddo

         do j=1,nel
            fp(1,j+ikel,iab,istate)=dorbn(1,iworbd(j+ish,kref),istate)
            fp(2,j+ikel,iab,istate)=dorbn(2,iworbd(j+ish,kref),istate)
            fp(3,j+ikel,iab,istate)=dorbn(3,iworbd(j+ish,kref),istate)
         enddo

         do k=1,ndet
            detiab(k,iab,istate)=detn(k,istate)
         enddo

         do iorb=1,norb
            orb(iel,iorb,istate)=orbn(iorb,istate)
            do kk=1,3
               dorb(kk,iel,iorb,istate)=dorbn(kk,iorb,istate)
            enddo
         enddo
      enddo

      end subroutine
