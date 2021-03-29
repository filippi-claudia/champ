      subroutine multideterminante(iel)

      use force_mod, only: MFORCE, MFORCE_WT_PRD, MWF
      use vmc_mod, only: MELEC, MORB, MBASIS, MDET, MCENT, MCTYPE, MCTYP3X
      use vmc_mod, only: NSPLIN, nrad, MORDJ, MORDJ1, MMAT_DIM, MMAT_DIM2, MMAT_DIM20
      use vmc_mod, only: radmax, delri
      use vmc_mod, only: NEQSX, MTERMS
      use vmc_mod, only: MCENT3, NCOEF, MEXCIT
      use csfs, only: nstates
      use dets, only: ndet
      use elec, only: ndn, nup
      use multidet, only: irepcol_det, ireporb_det, ivirt, iwundet, kref, numrep_det
      use slatn, only: slmin
      use ycompactn, only: ymatn
      use coefs, only: norb
      use multimatn, only: aan, wfmatn
      use multislatern, only: ddorbn, detn, dorbn, orbn
      use orbval, only: ddorb, dorb, nadorb, ndetorb, orb
      use multislater, only: detiab
      use mstates_mod, only: MSTATES

      implicit real*8(a-h,o-z)

      parameter (one=1.d0,half=0.5d0)
      dimension gmat(MELEC,MORB,3),gmatn(MEXCIT**2,3)
      dimension b(MORB,3),ddx_mdet(3)
      dimension orb_sav(MORB,MSTATES)
      dimension dum1(MSTATES)

      if(ndet.eq.1) return

      do istate=1,nstates
         iab=1
         nel=nup
         ish=0
         if(iel.gt.nup) then
            iab=2
            nel=ndn
            ish=nup
         endif

c     temporarely copy orbn to orb
         do iorb=1,norb
            orb_sav(iorb,istate)=orb(iel,iorb,istate)
            orb(iel,iorb,istate)=orbn(iorb,istate)
         enddo

         do jrep=ivirt(iab),norb
            do irep=1,nel
               dum1(istate)=0.0d0
               do i=1,nel
                  dum1(istate)=dum1(istate)+slmin(irep+(i-1)*nel,istate)*orb(i+ish,jrep,istate)
               enddo
               aan(irep,jrep,istate)=dum1(istate)
            enddo
         enddo

c     compute wave function 
         do k=1,ndet
            if(k.ne.kref) then
               if(iwundet(k,iab).eq.k) then
                  ndim=numrep_det(k,iab)
                  jj=0
                  do jrep=1,ndim
                     jorb=ireporb_det(jrep,k,iab)
                     do irep=1,ndim
                        iorb=irepcol_det(irep,k,iab)
                        jj=jj+1
                        wfmatn(jj,k,istate)=aan(iorb,jorb,istate)
                     enddo
                  enddo
                  call matinv(wfmatn(1,k,istate),ndim,det)
                  detn(k,istate)=det
               else
                  index_det=iwundet(k,iab)
                  detn(k,istate)=detn(index_det,istate)
               endif
            endif
         enddo

         do k=1,ndet
            if(k.ne.kref.and.iwundet(k,iab).ne.kref) then
               detn(k,istate)=detn(k,istate)*detn(kref,istate)
            endif
         enddo

         if(iab.eq.1) call compute_ymat(iab,detn(1,istate),detiab(1,2,istate),
     &        wfmatn(1,1,istate),ymatn(1,1,istate),istate)
         if(iab.eq.2) call compute_ymat(iab,detiab(1,1,istate),detn(1,istate),
     &        wfmatn(1,1,istate),ymatn(1,1,istate),istate)

         do iorb=1,norb
            orb(iel,iorb,istate)=orb_sav(iorb,istate)
         enddo
      enddo

      end subroutine

c-----------------------------------------------------------------------

      subroutine multideterminante_grad(iel,dorb,detratio,slmi,aa,wfmat,ymat,velocity)

      use force_mod, only: MFORCE, MFORCE_WT_PRD, MWF
      use vmc_mod, only: MELEC, MORB, MBASIS, MDET, MCENT, MCTYPE, MCTYP3X
      use vmc_mod, only: NSPLIN, nrad, MORDJ, MORDJ1, MMAT_DIM, MMAT_DIM2, MMAT_DIM20
      use vmc_mod, only: radmax, delri
      use vmc_mod, only: NEQSX, MTERMS
      use vmc_mod, only: MCENT3, NCOEF, MEXCIT
      use dets, only: ndet
      use elec, only: ndn, nup
      use multidet, only: iactv, ivirt, kref
      use coefs, only: norb
      use dorb_m, only: iworbd

      implicit real*8(a-h,o-z)

      parameter (one=1.d0,half=0.5d0)
      dimension slmi(MMAT_DIM)
      dimension aa(MELEC,MORB),wfmat(MEXCIT**2,MDET),ymat(MORB,MELEC)
      dimension b(MORB,3),dorb(3,MORB)
      dimension gmat(MELEC,MORB,3)
      dimension velocity(3)

      do k=1,3
         velocity(k)=0.0d0
      enddo

      if(ndet.eq.1) return

      if(iel.le.nup) then
         iab=1
         nel=nup
         ish=0
      else
         iab=2
         nel=ndn
         ish=nup
      endif

      jel=iel-ish
c     TMP to fix
      do kk=1,3
         do iorb=1,norb
            b(iorb,kk)=dorb(kk,iorb)
         enddo
      enddo

      do kk=1,3
         do jrep=ivirt(iab),norb
            dum=0.0d0
            do j=1,nel
               dum=dum+b(iworbd(j+ish,kref),kk)*aa(j,jrep)
            enddo
            dum=b(jrep,kk)-dum
            do irep=iactv(iab),nel
               gmat(irep,jrep,kk)=dum*slmi(irep+(jel-1)*nel)
            enddo
         enddo
      enddo

      do kk=1,3
         dum=0.0d0
         do jrep=ivirt(iab),norb
            do irep=iactv(iab),nel
               dum=dum+ymat(jrep,irep)*gmat(irep,jrep,kk)
            enddo
         enddo
         velocity(kk)=dum*detratio
      enddo

      end subroutine
