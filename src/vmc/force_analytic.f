      subroutine compute_force(psid,denergy)

      use atom, only: znuc, cent, pecent, iwctype, nctype, ncent
      use const, only: pi, hb, etrial, delta, deltai, fbias, nelec, imetro, ipr
      use da_energy_now, only: da_energy, da_psi
      use da_jastrow4val, only: da_d2j, da_j, da_vj
      use da_orbval, only: da_d2orb, da_dorb, da_orb
      use elec, only: ndn, nup
      use multidet, only: iactv, irepcol_det, ireporb_det, ivirt, iwundet, kref, numrep_det
      use coefs, only: coef, nbasis, norb
      use dorb_m, only: iworbd

      implicit real*8(a-h,o-z)

      include 'vmc.h'
      include 'force.h'

      common /slater/ slmi(MMAT_DIM,2)
     &,fpu(3,MMAT_DIM),fpd(3,MMAT_DIM)
     &,fppu(MMAT_DIM),fppd(MMAT_DIM)
     &,ddx(3,MELEC),d2dx2(MELEC)
      common /multislater/ detu(MDET),detd(MDET)
      common /orbval/ orb(MELEC,MORB),dorb(3,MELEC,MORB),ddorb(MELEC,MORB),ndetorb,nadorb

      dimension da_psi_ref(3,MCENT)

      call compute_da_psi(psid,da_psi_ref)
      call compute_da_energy(psid,denergy)

      do 800 ic=1,ncent
        do 800 k=1,3
          da_psi(k,ic)=da_psi(k,ic)+da_psi_ref(k,ic)
          do 800 i=1,nelec
 800        da_psi(k,ic)=da_psi(k,ic)+da_j(k,i,ic)

c     write(6,*)'da_ref',((da_psi_ref(l,ic),l=1,3),ic=1,ncent)
c     write(6,*) 'da_psi',((da_psi(k,ic),k=1,3),ic=1,ncent)

      return
      end
c-----------------------------------------------------------------------
      subroutine compute_da_psi(psid,da_psi_ref)
      use atom, only: znuc, cent, pecent, iwctype, nctype, ncent

      use const, only: pi, hb, etrial, delta, deltai, fbias, nelec, imetro, ipr
      use da_energy_now, only: da_energy, da_psi
      use da_jastrow4val, only: da_d2j, da_j, da_vj
      use da_orbval, only: da_d2orb, da_dorb, da_orb

      use elec, only: ndn, nup
      use multidet, only: iactv, irepcol_det, ireporb_det, ivirt, iwundet, kref, numrep_det

      use ycompact, only: dymat, ymat
      use zcompact, only: aaz, dzmat, emz, zmat

      use coefs, only: coef, nbasis, norb
      use dorb_m, only: iworbd
      implicit real*8(a-h,o-z)










      include 'vmc.h'
      include 'force.h'
      include 'mstates.h'


      common /slater/ slmi(MMAT_DIM,2)
     &,fpu(3,MMAT_DIM),fpd(3,MMAT_DIM)
     &,fppu(MMAT_DIM),fppd(MMAT_DIM)
     &,ddx(3,MELEC),d2dx2(MELEC)
      common /multislater/ detu(MDET),detd(MDET)




      common /orbval/ orb(MELEC,MORB),dorb(3,MELEC,MORB),ddorb(MELEC,MORB),ndetorb,nadorb



      dimension b_a(MORB,MELEC),b_kref(MELEC*MELEC),tildem_a(MELEC,MORB),xmat(MELEC*MELEC,2),work(MELEC)
      dimension da_psi_ref(3,MCENT)

      do 400 ic=1,ncent
        do 400 k=1,3

          do i=1,nelec
            do iorb=1,norb
              b_a(iorb,i)=da_orb(k,i,iorb,ic)
            enddo
          enddo

          trace=0
          da_psi_ref(k,ic)=0
          do iab=1,2

            if(iab.eq.1) then
              ish=0
              nel=nup
             else
              ish=nup
              nel=ndn
            endif

            ii=-nel
            do 110 i=1,nel
              ii=ii+nel
              do 110 j=1,nel
  110           b_kref(j+ii)=b_a(iworbd(j+ish,kref),i+ish)
          
c compute force for reference determinant
            do 120 i=1,nel
              do 120 j=1,nel
  120           da_psi_ref(k,ic)=da_psi_ref(k,ic)+slmi(i+(j-1)*nel,iab)*b_kref(i+(j-1)*nel)

            do jrep=ivirt(iab),norb
              do irep=1,nel
                trace=trace+zmat(jrep,irep,iab,1)*b_a(jrep,irep+ish)
              enddo
            enddo

            do jrep=1,nel
              jorb=iworbd(jrep+ish,kref)
              do irep=1,nel
                trace=trace-aaz(jrep,irep,iab,1)*b_a(jorb,irep+ish)
              enddo
            enddo

c enddo iab
          enddo

          da_psi(k,ic)=trace*detu(kref)*detd(kref)/psid

 400  continue

c     do 800 ic=1,ncent
c       do 800 k=1,3
c         da_psi(k,ic)=da_psi(k,ic)+da_psi_ref(k,ic)
c         do 800 i=1,nelec
c800        da_psi(k,ic)=da_psi(k,ic)+da_j(k,i,ic)

c     if(ipr.gt.3) write(6,*)'da_ref',((da_psi_ref(l,ic),l=1,3),ic=1,ncent)
c     if(ipr.gt.3) write(6,*)'da_psi',((da_psi(l,ic),l=1,3),ic=1,ncent)

      return
      end
c-----------------------------------------------------------------------
      subroutine compute_da_energy(psid,denergy)
      use atom, only: znuc, cent, pecent, iwctype, nctype, ncent
      use const, only: pi, hb, etrial, delta, deltai, fbias, nelec, imetro, ipr
      use da_energy_now, only: da_energy, da_psi
      use da_jastrow4val, only: da_d2j, da_j, da_vj
      use da_orbval, only: da_d2orb, da_dorb, da_orb
      use elec, only: ndn, nup
      use multidet, only: iactv, irepcol_det, ireporb_det, ivirt, iwundet, kref, numrep_det
      use zcompact, only: aaz, dzmat, emz, zmat
      use Bloc_da, only: b_da, db
      use coefs, only: coef, nbasis, norb
      use Bloc, only: b, tildem, xmat, xmatd, xmatu
      use dorb_m, only: iworbd
      use multimat, only: aa, wfmat

      implicit real*8(a-h,o-z)

      include 'vmc.h'
      include 'pseudo.h'
      include 'mstates.h'
      include 'force.h'

      common /slater/ slmi(MMAT_DIM,2)
     &,fpu(3,MMAT_DIM),fpd(3,MMAT_DIM)
     &,fppu(MMAT_DIM),fppd(MMAT_DIM)
     &,ddx(3,MELEC),d2dx2(MELEC)
      common /multislater/ detu(MDET),detd(MDET)




      common /orbval/ orb(MELEC,MORB),dorb(3,MELEC,MORB),ddorb(MELEC,MORB),ndetorb,nadorb
      common /velocity_jastrow/vj(3,MELEC),vjn(3,MELEC)


      common /pseudo/ vps(MELEC,MCENT,MPS_L),vpso(MELEC,MCENT,MPS_L,MFORCE)
     &,lpot(MCTYPE),nloc

      common /da_pseudo/ da_pecent(3,MCENT),da_vps(3,MELEC,MCENT,MPS_L),
     & da_nonloc(3,MCENT)



      dimension da_energy_ref(3,MCENT)

      do 400 ic=1,ncent
        do 400 k=1,3

          trace=0
          da_energy_ref(k,ic)=0
          do iab=1,2

            if(iab.eq.1) then
              ish=0
              nel=nup
             else
              ish=nup
              nel=ndn
            endif

c compute force for reference determinant
            do 120 i=1,nel
              do 120 j=1,nel
                jorb=iworbd(j+ish,kref)
  120           da_energy_ref(k,ic)=da_energy_ref(k,ic)+slmi(j+(i-1)*nel,iab)*b_da(k,i+ish,jorb,ic)
     &                                                 -da_orb(k,i+ish,jorb,ic)*xmat(i+(j-1)*nel,iab)
            do jrep=ivirt(iab),norb
              do irep=1,nel
                trace=trace+zmat(jrep,irep,iab,1)*b_da(k,irep+ish,jrep,ic) 
     &                     +dzmat(jrep,irep,iab,1)*da_orb(k,irep+ish,jrep,ic)
              enddo
            enddo

            do jrep=1,nel
              jorb=iworbd(jrep+ish,kref)
              do irep=1,nel
                trace=trace-emz(jrep,irep,iab,1)*da_orb(k,irep+ish,jorb,ic)
     &                     -aaz(jrep,irep,iab,1)*b_da(k,irep+ish,jorb,ic)
              enddo
            enddo
c enddo iab
          enddo

          da_energy(k,ic)=trace*detu(kref)*detd(kref)/psid
  400 continue

      do 800 ic=1,ncent
        ict=iwctype(ic)

        do 800 k=1,3
       
          da_other_kin=0.d0
          da_other_pot=da_pecent(k,ic)
          do 410 i=1,nelec
            da_other_kin=da_other_kin+da_d2j(k,i,ic)
     &               +2*(vj(1,i)*da_vj(k,1,i,ic)+vj(2,i)*da_vj(k,2,i,ic)+vj(3,i)*da_vj(k,3,i,ic))
  410       da_other_pot=da_other_pot+da_vps(k,i,ic,lpot(ict))

          da_energy(k,ic)=da_energy(k,ic)+da_energy_ref(k,ic)-hb*da_other_kin+da_other_pot
     &                   -denergy*da_psi(k,ic)

c complete da_psi
  800 continue

c     write(6,*)'da_energy',((da_energy(l,ic),l=1,3),ic=1,ncent)

      return
      end
c-----------------------------------------------------------------------
      subroutine force_analy_init(iflag)
      use atom, only: znuc, cent, pecent, iwctype, nctype, ncent
      use da_energy_now, only: da_energy, da_psi
      use da_energy_ave_m, only: da_energy_cm2, da_energy_cum, da_energy_sum, da_psi_cum, da_psi_sum

      use force_analy, only: iforce_analy
      implicit real*8(a-h,o-z)




      include 'vmc.h'




      if(iforce_analy.eq.0) return

      do 10 ic=1,ncent
        do 10 k=1,3
          da_psi_sum(k,ic)=0.0d0
  10      da_energy_sum(k,ic)=0.0d0

      if(iflag.gt.0) return

      do 20 ic=1,ncent
        do 20 k=1,3
          da_psi_cum(k,ic)=0.0d0
          da_energy_cum(k,ic)=0.0d0
  20      da_energy_cm2(k,ic)=0.0d0

      return
      end

c-----------------------------------------------------------------------
      subroutine force_analy_sum(p,q,eloc,eloco)
      use atom, only: znuc, cent, pecent, iwctype, nctype, ncent
      use da_energy_now, only: da_energy, da_psi
      use da_energy_ave_m, only: da_energy_cm2, da_energy_cum, da_energy_sum, da_psi_cum, da_psi_sum

      use force_analy, only: iforce_analy
      implicit real*8(a-h,o-z)




      include 'vmc.h'




      if(iforce_analy.eq.0) return

      do 10 ic=1,ncent
        do 10 k=1,3
          da_energy(k,ic)=da_energy(k,ic)+2*eloc*da_psi(k,ic)
          da_psi_sum(k,ic)= da_psi_sum(k,ic)+p*da_psi(k,ic)
  10      da_energy_sum(k,ic)= da_energy_sum(k,ic)+p*da_energy(k,ic)

      return
      end
c-----------------------------------------------------------------------
      subroutine force_analy_cum(wsum,eave,wcum)
      use atom, only: znuc, cent, pecent, iwctype, nctype, ncent
      use da_energy_ave_m, only: da_energy_cm2, da_energy_cum, da_energy_sum, da_psi_cum, da_psi_sum

      use force_analy, only: iforce_analy
      implicit real*8(a-h,o-z)



      include 'vmc.h'



      if(iforce_analy.eq.0) return

      do 10 ic=1,ncent
        do 10 k=1,3
          da_energy_now=(da_energy_sum(k,ic)-2*eave*da_psi_sum(k,ic))/wsum
          da_energy_cm2(k,ic)=da_energy_cm2(k,ic)+wsum*da_energy_now**2
          da_psi_cum(k,ic)=da_psi_cum(k,ic)+da_psi_sum(k,ic)
  10      da_energy_cum(k,ic)=da_energy_cum(k,ic)+da_energy_sum(k,ic)

      return
      end
c-----------------------------------------------------------------------
      subroutine force_analy_fin(wcum,iblk,eave)
      use atom, only: znuc, cent, pecent, iwctype, nctype, ncent
      use force_fin, only: da_energy_ave, da_energy_err
      use contrl, only: idump, irstar, isite, n_conf, nblk, nblkeq, nconf_new, nstep
      use da_energy_ave_m, only: da_energy_cm2, da_energy_cum, da_energy_sum, da_psi_cum, da_psi_sum

      use force_analy, only: iforce_analy
      implicit real*8(a-h,o-z)




      include 'vmc.h'





      err(x,x2)=dsqrt(abs(x2/wcum-(x/wcum)**2)/iblk)

      if(iforce_analy.eq.0) return

      rtpass=dsqrt(wcum)

      open(80,file='force_analytic',form='formatted',status='unknown')
      do 20 ic=1,ncent
        do 10 k=1,3
          da_energy_ave(k,ic)=(da_energy_cum(k,ic)-2*eave*da_psi_cum(k,ic))/wcum
   10     da_energy_err(k)=err(da_energy_ave(k,ic),da_energy_cm2(k,ic))
   20   write(80,'(i5,1p6e14.5)') ic,(da_energy_ave(k,ic),k=1,3),(da_energy_err(k),k=1,3)

       ! TODO JF this is included in the treatment of internal
       ! coordinates, remove this when finished
       !call transform_grad_zmat(da_energy_ave)

      return
      end
c-----------------------------------------------------------------------
      subroutine force_analy_dump(iu)
      use atom, only: znuc, cent, pecent, iwctype, nctype, ncent
      use contrl, only: idump, irstar, isite, n_conf, nblk, nblkeq, nconf_new, nstep
      use da_energy_ave_m, only: da_energy_cm2, da_energy_cum, da_energy_sum, da_psi_cum, da_psi_sum

      use force_analy, only: iforce_analy
      implicit real*8(a-h,o-z)



      include 'vmc.h'




      if(iforce_analy.eq.0) return

      write(iu) ((da_energy_cum(k,ic),da_psi_cum(k,ic),da_energy_cm2(k,ic),k=1,3),ic=1,ncent)

      return
      end
c-----------------------------------------------------------------------
      subroutine force_analy_rstrt(iu)
      use atom, only: znuc, cent, pecent, iwctype, nctype, ncent
      use contrl, only: idump, irstar, isite, n_conf, nblk, nblkeq, nconf_new, nstep
      use da_energy_ave_m, only: da_energy_cm2, da_energy_cum, da_energy_sum, da_psi_cum, da_psi_sum

      use force_analy, only: iforce_analy
      implicit real*8(a-h,o-z)



      include 'vmc.h'




      if(iforce_analy.eq.0) return

      read(iu) ((da_energy_cum(k,ic),da_psi_cum(k,ic),da_energy_cm2(k,ic),k=1,3),ic=1,ncent)

      return
      end
c-----------------------------------------------------------------------
      subroutine force_analy_save
      implicit real*8(a-h,o-z)

      return
      end
