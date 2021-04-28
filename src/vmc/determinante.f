      subroutine determinante(iel,x,rvec_en,r_en,iflag)

      use force_mod, only: MFORCE, MFORCE_WT_PRD, MWF
      use vmc_mod, only: MELEC, MORB, MBASIS, MDET, MCENT, MCTYPE, MCTYP3X
      use vmc_mod, only: NSPLIN, nrad, MORDJ, MORDJ1, MMAT_DIM, MMAT_DIM2, MMAT_DIM20
      use vmc_mod, only: radmax, delri
      use vmc_mod, only: NEQSX, MTERMS
      use vmc_mod, only: MCENT3, NCOEF, MEXCIT
      use elec, only: ndn, nup
      use multidet, only: kref
      use slatn, only: slmin
      use dorb_m, only: iworbd
      use multislatern, only: ddorbn, detn, dorbn, orbn 
      use slater, only: d2dx2, ddx, fp, fpp, slmi
      use multislater, only: detiab
      use csfs, only: nstates

      implicit real*8(a-h,o-z)

      dimension x(3,*),rvec_en(3,MELEC,MCENT),r_en(MELEC,MCENT)

      call orbitalse(iel,x,rvec_en,r_en,iflag)

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
         ratio_kref=0.0d0
         do j=1,nel
            ratio_kref=ratio_kref+slmi(j+ikel,iab,istate)*orbn(iworbd(j+ish,kref),istate)
         enddo

         detn(kref,istate)=detiab(kref,iab,istate)*ratio_kref

         if(ratio_kref.eq.0.d0) return

         do i=1,nel
            if(i+ish.ne.iel) then
               ik=nel*(i-1)
               sum=0.0d0
               do j=1,nel
                  sum=sum+slmi(j+ik,iab,istate)*orbn(iworbd(j+ish,kref),istate)
               enddo
               sum=sum/ratio_kref
               do j=1,nel
                  slmin(j+ik,istate)=slmi(j+ik,iab,istate)-slmi(j+ikel,iab,istate)*sum
               enddo
            endif
         enddo

         do j=1,nel
            slmin(j+ikel,istate)=slmi(j+ikel,iab,istate)/ratio_kref
         enddo
      enddo

      end subroutine

c-----------------------------------------------------------------------

      subroutine compute_determinante_grad(iel,psig,psid,vd,iflag_move)
      use force_mod, only: MFORCE, MFORCE_WT_PRD, MWF
      use vmc_mod, only: MELEC, MORB, MBASIS, MDET, MCENT, MCTYPE, MCTYP3X
      use vmc_mod, only: NSPLIN, nrad, MORDJ, MORDJ1, MMAT_DIM, MMAT_DIM2, MMAT_DIM20
      use vmc_mod, only: radmax, delri
      use vmc_mod, only: NEQSX, MTERMS
      use vmc_mod, only: MCENT3, NCOEF, MEXCIT
      use mstates_mod, only: MSTATES
      use csfs, only: nstates
      use elec, only: nup
      use multidet, only: kref
      use slatn, only: slmin
      use ycompact, only: ymat
      use ycompactn, only: ymatn
      use coefs, only: norb
      use multimat, only: aa, wfmat
      use multimatn, only: aan, wfmatn
      use velocity_jastrow, only: vj, vjn
      use mstates_ctrl, only: iguiding
      use mstates3, only: iweight_g, weights_g
      use multislatern, only: ddorbn, detn, dorbn, orbn 
      use orbval, only: ddorb, dorb, nadorb, ndetorb, orb
      use slater, only: d2dx2, ddx, fp, fpp, slmi
      use multislater, only: detiab
      use config, only: anormo

      implicit real*8(a-h,o-z)

      dimension psid(*),vd(3),vref(3,MSTATES),vd_s(3),dorb_tmp(3,MORB,MSTATES)
      dimension ymat_tmp(MORB,MELEC,MSTATES)
      save ymat_tmp

      if(iel.le.nup) then
         iab=1
      else
         iab=2
      endif

      psi2g=psig*psig
      psi2gi=1.0d0/psi2g

c     All quantities saved (old) avaliable
      if(iflag_move.eq.1) then
         do istate=1,nstates
            do kk=1,3
               do iorb=1,norb
                  dorb_tmp(kk,iorb,istate)=dorb(kk,iel,iorb,istate)
               enddo
            enddo
            call determinante_ref_grad(iel,slmi(1,iab,istate),dorb_tmp(1,1,istate),vref(1,istate))
         enddo

         if(iguiding.eq.0) then
            detratio=detiab(kref,1,1)*detiab(kref,2,1)/psid(1)
            call multideterminante_grad(iel,dorb_tmp(1,1,1),detratio,slmi(1,iab,1),
     &           aa(1,1,iab,1),wfmat(1,1,iab,1),ymat(1,1,iab,1),vd)
            do kk=1,3
               vd(kk)=vd(kk)+vref(kk,1)
            enddo
         else
            do kk=1,3
               vd(kk)=0.0d0
            enddo
            do i=1,nstates
               istate=iweight_g(i)
               detratio=detiab(kref,1,istate)*detiab(kref,2,istate)/psid(istate)
               call multideterminante_grad(iel,dorb_tmp(1,1,istate),detratio,slmi(1,iab,istate),
     &              aa(1,1,iab,istate),wfmat(1,1,iab,istate),ymat(1,1,iab,istate),vd_s)
               do kk=1,3
                  vd(kk)=vd(kk)+weights_g(i)*psid(istate)*psid(istate)
     &                   *(vd_s(kk)+vref(kk,istate))/anormo(istate)
               enddo
            enddo
            vd(1)=vd(1)*psi2gi
            vd(2)=vd(2)*psi2gi
            vd(3)=vd(3)*psi2gi
         endif
         
         vd(1)=vj(1,iel)+vd(1)
         vd(2)=vj(2,iel)+vd(2)
         vd(3)=vj(3,iel)+vd(3)

c     Within single-electron move - quantities of electron iel not saved 
      elseif(iflag_move.eq.0) then
         do istate=1,nstates
            call determinante_ref_grad(iel,slmin(1,istate),dorbn(1,1,istate),vref(1,istate))
	 enddo
         if(iguiding.eq.0) then
            if(iab.eq.1) then
               detratio=detn(kref,1)*detiab(kref,2,1)/psid(1)
            else
               detratio=detiab(kref,1,1)*detn(kref,1)/psid(1)
            endif
            call multideterminante_grad(iel,dorbn(1,1,1),
     &           detratio,slmin(1,1),aan(1,1,1),wfmatn(1,1,1),ymatn(1,1,1),vd)
            do kk=1,3
               vd(kk)=vd(kk)+vref(kk,1)
            enddo
         else
            do kk=1,3
               vd(kk)=0.0d0
            enddo
            do i=1,nstates
               istate=iweight_g(i)
               if(iab.eq.1) then
                  detratio=detn(kref,istate)*detiab(kref,2,istate)/psid(istate)
               else
                  detratio=detiab(kref,1,istate)*detn(kref,istate)/psid(istate)
               endif
               call multideterminante_grad(iel,dorbn(1,1,istate),detratio,slmin(1,istate),
     &              aan(1,1,istate),wfmatn(1,1,istate),ymatn(1,1,istate),vd_s)
               do kk=1,3
                  vd(kk)=vd(kk)+weights_g(i)*psid(istate)*psid(istate)
     &                   *(vd_s(kk)+vref(kk,istate))/anormo(istate)
               enddo
            enddo
            vd(1)=vd(1)*psi2gi
            vd(2)=vd(2)*psi2gi
            vd(3)=vd(3)*psi2gi
         endif

         vd(1)=vjn(1,iel)+vd(1)
         vd(2)=vjn(2,iel)+vd(2)
         vd(3)=vjn(3,iel)+vd(3)
      else
c     Within single-electron move - iel not equal to electron moved - quantities of electron iel not saved 
         do kk=1,3
            do istate=1,nstates
               do iorb=1,norb
                  dorb_tmp(kk,iorb,istate)=dorb(kk,iel,iorb,istate)
               enddo
            enddo
         enddo
c     iel has same spin as electron moved
         if(iflag_move.eq.2) then
            do istate=1,nstates
               if(iab.eq.1) then
                  detratio=detn(kref,istate)*detiab(kref,2,istate)/psid(istate)
               else
                  detratio=detiab(kref,1,istate)*detn(kref,istate)/psid(istate)
               endif
               call determinante_ref_grad(iel,slmin(1,istate),dorb_tmp(1,1,istate),vref(1,istate))
               call multideterminante_grad(iel,dorb_tmp(1,1,istate),
     &              detratio,slmin(1,istate),aan(1,1,istate),wfmatn(1,1,istate),ymatn(1,1,istate),vd)
            enddo
c     iel has different spin than the electron moved
         else
            do istate=1,nstates
               if(iab.eq.1) then
                  detratio=detiab(kref,1,istate)*detn(kref,istate)/psid(istate)
               else
                  detratio=detn(kref,istate)*detiab(kref,2,istate)/psid(istate)
               endif
               call determinante_ref_grad(iel,slmi(1,iab,istate),dorb_tmp(1,1,istate),vref(1,istate))
               if(iel.eq.1) then
                  call compute_ymat(1,detiab(1,1,istate),detn(1,istate),
     &                 wfmat(1,1,1,istate),ymat_tmp(1,1,istate),1)
               endif
               if(iel.eq.nup+1) then
                  call compute_ymat(2,detn(1,istate),detiab(1,2,istate),
     &                 wfmat(1,1,2,istate),ymat_tmp(1,1,istate),1)
               endif
               call multideterminante_grad(iel,dorb_tmp(1,1,istate),detratio,slmi(1,iab,istate),
     &              aa(1,1,iab,istate),wfmat(1,1,iab,istate),ymat_tmp(1,1,istate),vd)
            enddo
         endif

         vd(1)=vjn(1,iel)+vd(1)+vref(1,istate)
         vd(2)=vjn(2,iel)+vd(2)+vref(2,istate)
         vd(3)=vjn(3,iel)+vd(3)+vref(3,istate)
      endif

      end subroutine

c-----------------------------------------------------------------------

      subroutine determinante_ref_grad(iel,slmi,dorb,ddx_ref)

      use force_mod, only: MFORCE, MFORCE_WT_PRD, MWF
      use vmc_mod, only: MELEC, MORB, MBASIS, MDET, MCENT, MCTYPE, MCTYP3X
      use vmc_mod, only: NSPLIN, nrad, MORDJ, MORDJ1, MMAT_DIM, MMAT_DIM2, MMAT_DIM20
      use vmc_mod, only: radmax, delri
      use vmc_mod, only: NEQSX, MTERMS
      use vmc_mod, only: MCENT3, NCOEF, MEXCIT
      use elec, only: ndn, nup
      use multidet, only: kref
      use dorb_m, only: iworbd

      implicit real*8(a-h,o-z)

      dimension slmi(MMAT_DIM),dorb(3,MORB)
      dimension ddx_ref(3)

      ddx_ref(1)=0.0d0
      ddx_ref(2)=0.0d0
      ddx_ref(3)=0.0d0

      if(iel.le.nup) then
         ish=0
         jel=iel
         nel=nup
      else
         ish=nup
         jel=iel-nup
         nel=ndn
      endif

      ik=(jel-1)*nel
      do j=1,nel
         ddx_ref(1)=ddx_ref(1)+slmi(j+ik)*dorb(1,iworbd(j+ish,kref))
         ddx_ref(2)=ddx_ref(2)+slmi(j+ik)*dorb(2,iworbd(j+ish,kref))
         ddx_ref(3)=ddx_ref(3)+slmi(j+ik)*dorb(3,iworbd(j+ish,kref))
      enddo

      end subroutine
