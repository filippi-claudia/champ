      module force_analytic
      contains
      subroutine compute_force(psid,denergy)

      use contrl_file, only: ounit
      use da_energy_now, only: da_psi
      use da_jastrow4val, only: da_j
      use precision_kinds, only: dp
      use system,  only: ncent,ncent_tot,nelec
      implicit none

      integer :: i, ic, k
      real(dp) :: denergy, psid
      real(dp), dimension(3, ncent_tot) :: da_psi_ref

      call compute_da_psi(psid,da_psi_ref)
      call compute_da_energy(psid,denergy)

      do ic=1,ncent
        do k=1,3
          da_psi(k,ic)=da_psi(k,ic)+da_psi_ref(k,ic)
          do i=1,nelec
            da_psi(k,ic)=da_psi(k,ic)+da_j(k,i,ic)
          enddo
        enddo
      enddo

c     write(ounit,*)'da_ref',((da_psi_ref(l,ic),l=1,3),ic=1,ncent)
c     write(ounit,*) 'da_psi',((da_psi(k,ic),k=1,3),ic=1,ncent)

      return
      end
c-----------------------------------------------------------------------
      subroutine compute_da_psi(psid,da_psi_ref)

      use da_energy_now, only: da_psi
      use da_jastrow4val, only: da_j
      use da_orbval, only: da_orb
      use dorb_m,  only: iworbd
      use multidet, only: ivirt
      use multislater, only: detiab
      use precision_kinds, only: dp
      use slater,  only: kref,norb,slmi
      use system,  only: ncent,ncent_tot,ndn,nelec,nup
      use vmc_mod, only: norb_tot
      use zcompact, only: aaz,zmat

      implicit none

      integer :: i, iab, ic, ii, iorb
      integer :: irep, ish, j, jorb
      integer :: jrep, k, nel
      real(dp) :: c800, psid, trace
      real(dp), dimension(norb_tot, nelec) :: b_a
      real(dp), dimension(nelec*nelec) :: b_kref
      real(dp), dimension(nelec, norb_tot) :: tildem_a
      real(dp), dimension(3, ncent_tot) :: da_psi_ref

      do ic=1,ncent
        do k=1,3

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
            do i=1,nel
              ii=ii+nel
              do j=1,nel
                b_kref(j+ii)=b_a(iworbd(j+ish,kref),i+ish)
              enddo
            enddo

c compute force for reference determinant
            do i=1,nel
              do j=1,nel
                da_psi_ref(k,ic)=da_psi_ref(k,ic)+slmi(i+(j-1)*nel,iab)*b_kref(i+(j-1)*nel)
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

c enddo iab
          enddo

          da_psi(k,ic)=trace*detiab(kref,1)*detiab(kref,2)/psid

        enddo
      enddo

c     do 800 ic=1,ncent
c       do 800 k=1,3
c         da_psi(k,ic)=da_psi(k,ic)+da_psi_ref(k,ic)
c         do 800 i=1,nelec
c800        da_psi(k,ic)=da_psi(k,ic)+da_j(k,i,ic)

c     if(ipr.gt.3) write(ounit,*)'da_ref',((da_psi_ref(l,ic),l=1,3),ic=1,ncent)
c     if(ipr.gt.3) write(ounit,*)'da_psi',((da_psi(l,ic),l=1,3),ic=1,ncent)

      return
      end
c-----------------------------------------------------------------------
      subroutine compute_da_energy(psid,denergy)

      use Bloc,    only: b_da,xmat
      use constants, only: hb
      use da_energy_now, only: da_energy,da_psi
      use da_jastrow4val, only: da_d2j,da_vj
      use da_orbval, only: da_orb
      use da_pseudo, only: da_pecent,da_vps
      use dorb_m,  only: iworbd
      use multidet, only: ivirt
      use multislater, only: detiab
      use precision_kinds, only: dp
      use pseudo,  only: lpot
      use slater,  only: kref,norb,slmi
      use system,  only: iwctype,ncent,ncent_tot,ndn,nelec,nup
      use velocity_jastrow, only: vj
      use zcompact, only: aaz,dzmat,emz,zmat

      implicit none

      integer :: i, iab, ic, ict, irep
      integer :: ish, j, jorb, jrep
      integer :: k, nel
      real(dp) :: da_other_kin, da_other_pot, denergy, psid, trace
      real(dp), dimension(3, ncent_tot) :: da_energy_ref

      do ic=1,ncent
        do k=1,3

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
            do i=1,nel
              do j=1,nel
                jorb=iworbd(j+ish,kref)
                da_energy_ref(k,ic)=da_energy_ref(k,ic)+slmi(j+(i-1)*nel,iab)*b_da(k,i+ish,jorb,ic)
     &                                                 -da_orb(k,i+ish,jorb,ic)*xmat(i+(j-1)*nel,iab)
              enddo
            enddo
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

          da_energy(k,ic)=trace*detiab(kref,1)*detiab(kref,2)/psid
        enddo
      enddo

      do ic=1,ncent
        ict=iwctype(ic)

        do k=1,3

          da_other_kin=0.d0
          da_other_pot=da_pecent(k,ic)
          do i=1,nelec
            da_other_kin=da_other_kin+da_d2j(k,i,ic)
     &               +2*(vj(1,i)*da_vj(k,1,i,ic)+vj(2,i)*da_vj(k,2,i,ic)+vj(3,i)*da_vj(k,3,i,ic))
            da_other_pot=da_other_pot+da_vps(k,i,ic,lpot(ict))
          enddo

          da_energy(k,ic)=da_energy(k,ic)+da_energy_ref(k,ic)-hb*da_other_kin+da_other_pot
     &                   -denergy*da_psi(k,ic)

c complete da_psi
        enddo
      enddo

c     write(ounit,*)'da_energy',((da_energy(l,ic),l=1,3),ic=1,ncent)

      return
      end
c-----------------------------------------------------------------------
      subroutine force_analy_init(iflag)

      use da_energy_sumcum, only: da_energy_cm2,da_energy_cum
      use da_energy_sumcum, only: da_energy_sum,da_psi_cum,da_psi_sum
      use m_force_analytic, only: iforce_analy
      use system,  only: ncent

      implicit none

      integer :: ic, iflag, k

      if(iforce_analy.eq.0) return

      do ic=1,ncent
        do k=1,3
          da_psi_sum(k,ic)=0.0d0
          da_energy_sum(k,ic)=0.0d0
        enddo
      enddo

      if(iflag.gt.0) return

      do ic=1,ncent
        do k=1,3
          da_psi_cum(k,ic)=0.0d0
          da_energy_cum(k,ic)=0.0d0
          da_energy_cm2(k,ic)=0.0d0
        enddo
      enddo

      return
      end

c-----------------------------------------------------------------------
      subroutine force_analy_sum(p,q,eloc,eloco)

      use da_energy_now, only: da_energy,da_psi
      use da_energy_sumcum, only: da_energy_sum,da_psi_sum
      use m_force_analytic, only: iforce_analy
      use precision_kinds, only: dp
      use system,  only: ncent

      implicit none

      integer :: ic, k
      real(dp) :: eloc, eloco, p, q

      if(iforce_analy.eq.0) return

      do ic=1,ncent
        do k=1,3
          da_energy(k,ic)=da_energy(k,ic)+2*eloc*da_psi(k,ic)
          da_psi_sum(k,ic)= da_psi_sum(k,ic)+p*da_psi(k,ic)
          da_energy_sum(k,ic)= da_energy_sum(k,ic)+p*da_energy(k,ic)
        enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine force_analy_cum(wsum,eave,wcum)

      use da_energy_sumcum, only: da_energy_cm2,da_energy_cum
      use da_energy_sumcum, only: da_energy_sum,da_psi_cum,da_psi_sum
      use m_force_analytic, only: iforce_analy
      use precision_kinds, only: dp
      use system,  only: ncent

      implicit none

      integer :: ic, k
      real(dp) :: da_energy_now, eave, wcum, wsum

      if(iforce_analy.eq.0) return

      do ic=1,ncent
        do k=1,3
          da_energy_now=(da_energy_sum(k,ic)-2*eave*da_psi_sum(k,ic))/wsum
          da_energy_cm2(k,ic)=da_energy_cm2(k,ic)+wsum*da_energy_now**2
          da_psi_cum(k,ic)=da_psi_cum(k,ic)+da_psi_sum(k,ic)
          da_energy_cum(k,ic)=da_energy_cum(k,ic)+da_energy_sum(k,ic)
        enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine force_analy_fin(wcum,iblk,eave)

      use da_energy_sumcum, only: da_energy_cm2,da_energy_cum,da_psi_cum
      use m_force_analytic, only: da_energy_ave,iforce_analy
      use precision_kinds, only: dp
      use system,  only: ncent

      implicit none

      integer :: iblk, ic, k
      real(dp) :: eave, err, rtpass, wcum, x
      real(dp) :: x2
      real(dp), dimension(3) :: da_energy_err

      err(x,x2)=dsqrt(abs(x2/wcum-(x/wcum)**2)/iblk)

      if(iforce_analy.eq.0) return

      rtpass=dsqrt(wcum)

      open(80,file='force_analytic',form='formatted',status='unknown')
      do ic=1,ncent
        do k=1,3
          da_energy_ave(k,ic)=(da_energy_cum(k,ic)-2*eave*da_psi_cum(k,ic))/wcum
          da_energy_err(k)=err(da_energy_ave(k,ic),da_energy_cm2(k,ic))
        enddo
        write(80,'(i5,1p6e14.5)') ic,(da_energy_ave(k,ic),k=1,3),(da_energy_err(k),k=1,3)
      enddo

       ! TODO JF this is included in the treatment of internal
       ! coordinates, remove this when finished
       !call transform_grad_zmat(da_energy_ave)

      return
      end
c-----------------------------------------------------------------------
      subroutine force_analy_dump(iu)

      use da_energy_sumcum, only: da_energy_cm2,da_energy_cum,da_psi_cum
      use m_force_analytic, only: iforce_analy
      use system,  only: ncent

      implicit none

      integer :: ic, iu, k

      if(iforce_analy.eq.0) return

      write(iu) ((da_energy_cum(k,ic),da_psi_cum(k,ic),da_energy_cm2(k,ic),k=1,3),ic=1,ncent)

      return
      end
c-----------------------------------------------------------------------
      subroutine force_analy_rstrt(iu)

      use da_energy_sumcum, only: da_energy_cm2,da_energy_cum,da_psi_cum
      use m_force_analytic, only: iforce_analy
      use system,  only: ncent

      implicit none

      integer :: ic, iu, k

      if(iforce_analy.eq.0) return

      read(iu) ((da_energy_cum(k,ic),da_psi_cum(k,ic),da_energy_cm2(k,ic),k=1,3),ic=1,ncent)

      return
      end
c-----------------------------------------------------------------------
      subroutine force_analy_save
      implicit none

      return
      end
      end module
