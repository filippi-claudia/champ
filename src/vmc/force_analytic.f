      subroutine compute_force(psid,denergy)
      use vmc_mod, only: MCENT
      use mstates_mod, only: MSTATES
      use atom, only: ncent
      use const, only: nelec
      use da_jastrow4val, only: da_j
      use da_energy_now, only: da_psi
      use csfs, only: nstates

      implicit real*8(a-h,o-z)

      dimension da_psi_ref(3,MCENT,MSTATES)

      do istate=1,nstates
         call compute_da_psi(psid,da_psi_ref,istate)
         call compute_da_energy(psid,denergy,istate)

         do ic=1,ncent
            do k=1,3
               da_psi(k,ic,istate)=da_psi(k,ic,istate)+da_psi_ref(k,ic,istate)
               do i=1,nelec
                  da_psi(k,ic,istate)=da_psi(k,ic,istate)+da_j(k,i,ic)
               enddo
            enddo
         enddo
      enddo

      end subroutine

c-----------------------------------------------------------------------

      subroutine compute_da_psi(psid,da_psi_ref,istate)
      use vmc_mod, only: MELEC, MORB, MDET, MCENT
      use vmc_mod, only: MMAT_DIM
      use atom, only: ncent
      use const, only: nelec, ipr
      use da_energy_now, only: da_psi
      use da_jastrow4val, only: da_j
      use da_orbval, only: da_orb
      use elec, only: ndn, nup
      use multidet, only: ivirt, kref
      use zcompact, only: aaz, zmat
      use coefs, only: norb
      use dorb_m, only: iworbd
      use orbval, only: ddorb, dorb, nadorb, ndetorb, orb
      use slater, only: d2dx2, ddx, fp, fpp, slmi
      use multislater, only: detiab

      implicit real*8(a-h,o-z)

      dimension b_a(MORB,MELEC),b_kref(MELEC*MELEC),tildem_a(MELEC,MORB)
      dimension da_psi_ref(3,MCENT)

      do ic=1,ncent
         do k=1,3
            do i=1,nelec
               do iorb=1,norb
                  b_a(iorb,i)=da_orb(k,i,iorb,ic,istate)
               enddo
            enddo
            trace=0.0d0
            da_psi_ref(k,ic)=0.0d0
            do iab=1,2
               if(iab.eq.1) then
                  ish=0
                  nel=nup
               else
                  ish=nup
                  nel=ndn
               endif
               ii=-nel
               do i=1,nel
                  ii=ii+nel
                  do j=1,nel
                     b_kref(j+ii)=b_a(iworbd(j+ish,kref),i+ish)
                  enddo
               enddo
c     compute force for reference determinant
               do i=1,nel
                  do j=1,nel
                     da_psi_ref(k,ic)=da_psi_ref(k,ic)+
     &                    slmi(i+(j-1)*nel,istate,iab)*b_kref(i+(j-1)*nel)
                  enddo
               enddo
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
c     enddo iab
            enddo
            da_psi(k,ic,istate)=trace*detiab(kref,istate,1)*detiab(kref,istate,2)/psid
         enddo
      enddo

      end subroutine

c-----------------------------------------------------------------------

      subroutine compute_da_energy(psid,denergy,istate)
      use vmc_mod, only: MELEC, MORB, MDET, MCENT
      use vmc_mod, only: MMAT_DIM
      use atom, only: iwctype, ncent
      use const, only: hb, nelec
      use da_energy_now, only: da_energy, da_psi
      use da_jastrow4val, only: da_d2j, da_vj
      use da_orbval, only: da_orb
      use elec, only: ndn, nup
      use multidet, only: ivirt, kref
      use zcompact, only: aaz, dzmat, emz, zmat
      use Bloc, only: b_da
      use coefs, only: norb
      use Bloc, only: xmat
      use dorb_m, only: iworbd
      use pseudo, only: lpot
      use da_pseudo, only: da_pecent, da_vps
      use velocity_jastrow, only: vj
      use orbval, only: ddorb, dorb, nadorb, ndetorb, orb
      use slater, only: d2dx2, ddx, fp, fpp, slmi
      use multislater, only: detiab

      implicit real*8(a-h,o-z)

      dimension da_energy_ref(3,MCENT)

      do ic=1,ncent
         do k=1,3
            trace=0.0d0
            da_energy_ref(k,ic)=0.0d0
            do iab=1,2
               if(iab.eq.1) then
                  ish=0
                  nel=nup
               else
                  ish=nup
                  nel=ndn
               endif
c     compute force for reference determinant
               do i=1,nel
                  do j=1,nel
                     jorb=iworbd(j+ish,kref)
                     da_energy_ref(k,ic)=da_energy_ref(k,ic)
     &                    +slmi(j+(i-1)*nel,istate,iab)*b_da(k,i+ish,jorb,ic,istate)
     &                    -da_orb(k,i+ish,jorb,ic,istate)*xmat(i+(j-1)*nel,istate,iab)
                  enddo
               enddo
               do jrep=ivirt(iab),norb
                  do irep=1,nel
                     trace=trace+zmat(jrep,irep,iab,1)*b_da(k,irep+ish,jrep,ic,istate) 
     &                    +dzmat(jrep,irep,iab,1)*da_orb(k,irep+ish,jrep,ic,istate)
                  enddo
               enddo
               do jrep=1,nel
                  jorb=iworbd(jrep+ish,kref)
                  do irep=1,nel
                     trace=trace-emz(jrep,irep,iab,1)*da_orb(k,irep+ish,jorb,ic,istate)
     &                    -aaz(jrep,irep,iab,1)*b_da(k,irep+ish,jorb,ic,istate)
                  enddo
               enddo
c     enddo iab
            enddo
            da_energy(k,ic,istate)=trace*detiab(kref,istate,1)*detiab(kref,istate,2)/psid
         enddo
      enddo

      do ic=1,ncent
         ict=iwctype(ic)
         do k=1,3
            da_other_kin=0.d0
            da_other_pot=da_pecent(k,ic)
            do i=1,nelec
               da_other_kin=da_other_kin+da_d2j(k,i,ic)
     &              +2*(vj(1,i)*da_vj(k,1,i,ic)+vj(2,i)*da_vj(k,2,i,ic)+vj(3,i)*da_vj(k,3,i,ic))
               da_other_pot=da_other_pot+da_vps(k,i,ic,lpot(ict))
            enddo
            da_energy(k,ic,istate)=da_energy(k,ic,istate)+da_energy_ref(k,ic)-hb*da_other_kin+da_other_pot
     &           -denergy*da_psi(k,ic,istate)
c     complete da_psi
         enddo
      enddo

      end subroutine

c-----------------------------------------------------------------------

      subroutine force_analy_init(iflag)
      use atom, only: ncent
      use da_energy_sumcum, only: da_energy_cm2, da_energy_cum, da_energy_sum, da_psi_cum, da_psi_sum
      use force_analy, only: iforce_analy
      use csfs, only: nstates

      implicit real*8(a-h,o-z)


      if(iforce_analy.eq.0) return

      do istate=1,nstates
         do ic=1,ncent
            do k=1,3
               da_psi_sum(k,ic,istate)=0.0d0
               da_energy_sum(k,ic,istate)=0.0d0
            enddo
         enddo
      enddo

      if(iflag.gt.0) return

      do istate=1,nstates
         do ic=1,ncent
            do k=1,3
               da_psi_cum(k,ic,istate)=0.0d0
               da_energy_cum(k,ic,istate)=0.0d0
               da_energy_cm2(k,ic,istate)=0.0d0
            enddo
         enddo
      enddo

      end subroutine

c-----------------------------------------------------------------------

      subroutine force_analy_sum(p,q,eloc,eloco,istate)
      use atom, only: ncent
      use da_energy_now, only: da_energy, da_psi
      use da_energy_sumcum, only: da_energy_sum, da_psi_sum
      use force_analy, only: iforce_analy

      implicit real*8(a-h,o-z)

      if(iforce_analy.eq.0) return

      do ic=1,ncent
         do k=1,3
            da_energy(k,ic,istate)=da_energy(k,ic,istate)+2*eloc*da_psi(k,ic,istate)
            da_psi_sum(k,ic,istate)= da_psi_sum(k,ic,istate)+p*da_psi(k,ic,istate)
            da_energy_sum(k,ic,istate)= da_energy_sum(k,ic,istate)+p*da_energy(k,ic,istate)
         enddo
      enddo

      end subroutine

c-----------------------------------------------------------------------

      subroutine force_analy_cum(wsum,eave,wcum,istate)
      use atom, only: ncent
      use da_energy_sumcum, only: da_energy_cm2, da_energy_cum, da_energy_sum, da_psi_cum, da_psi_sum
      use force_analy, only: iforce_analy

      implicit real*8(a-h,o-z)

      if(iforce_analy.eq.0) return

      do ic=1,ncent
         do k=1,3
            da_energy_now=(da_energy_sum(k,ic,istate)-2*eave*da_psi_sum(k,ic,istate))/wsum
            da_energy_cm2(k,ic,istate)=da_energy_cm2(k,ic,istate)+wsum*da_energy_now**2
            da_psi_cum(k,ic,istate)=da_psi_cum(k,ic,istate)+da_psi_sum(k,ic,istate)
            da_energy_cum(k,ic,istate)=da_energy_cum(k,ic,istate)+da_energy_sum(k,ic,istate)
         enddo
      enddo

      end subroutine

c-----------------------------------------------------------------------

      subroutine force_analy_fin(wcum,iblk,eave,istate)
      use atom, only: ncent
      use force_fin, only: da_energy_ave, da_energy_err
      use da_energy_sumcum, only: da_energy_cm2, da_energy_cum, da_psi_cum
      use force_analy, only: iforce_analy

      implicit real*8(a-h,o-z)

      err(x,x2)=dsqrt(abs(x2/wcum-(x/wcum)**2)/iblk)

      if(iforce_analy.eq.0) return

      rtpass=dsqrt(wcum)

      open(80,file='force_analytic',form='formatted',status='unknown')

      do ic=1,ncent
         do k=1,3
            da_energy_ave(k,ic,istate)=(da_energy_cum(k,ic,istate)-2*eave*da_psi_cum(k,ic,istate))/wcum
            da_energy_err(k,istate)=err(da_energy_ave(k,ic,istate),da_energy_cm2(k,ic,istate))
         enddo
c        RLPB fix this to print all states
         write(80,'(i5,1p6e14.5)') ic,(da_energy_ave(k,ic,1),k=1,3),(da_energy_err(k,1),k=1,3)
      enddo

      end subroutine

c-----------------------------------------------------------------------

      subroutine force_analy_dump(iu)
      use atom, only: ncent
      use da_energy_sumcum, only: da_energy_cm2, da_energy_cum, da_psi_cum
      use force_analy, only: iforce_analy

      implicit real*8(a-h,o-z)

      if(iforce_analy.eq.0) return

c     RLPB need to extend to all states
      istate=1

      write(iu) ((da_energy_cum(k,ic,istate),da_psi_cum(k,ic,istate),
     &     da_energy_cm2(k,ic,istate),k=1,3),ic=1,ncent)

      end subroutine

c-----------------------------------------------------------------------

      subroutine force_analy_rstrt(iu)
      use atom, only: ncent
      use da_energy_sumcum, only: da_energy_cm2, da_energy_cum, da_psi_cum
      use force_analy, only: iforce_analy

      implicit real*8(a-h,o-z)

      if(iforce_analy.eq.0) return
      
c     RLPB need to extend to all states
      istate=1

      read(iu) ((da_energy_cum(k,ic,istate),da_psi_cum(k,ic,istate),
     &     da_energy_cm2(k,ic,istate),k=1,3),ic=1,ncent)

      end subroutine

c-----------------------------------------------------------------------

      subroutine force_analy_save
      implicit real*8(a-h,o-z)

      return

      end subroutine
