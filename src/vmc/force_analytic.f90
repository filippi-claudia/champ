module force_analytic
contains
      subroutine compute_force(psid,denergy)

      use contrl_file, only: ounit
      use da_energy_now, only: da_psi
      use da_jastrow, only: da_j
      use precision_kinds, only: dp
      use system,  only: ncent,ncent_tot,nelec
      use contrl_file,    only: ounit
      implicit none

      integer :: i, ic, j, k
      real(dp) :: denergy, psid
      real(dp), dimension(3, ncent_tot) :: da_psi_ref

!     ! multistate indices were not added
      call compute_da_psi(psid,da_psi_ref)
      call compute_da_energy(psid,denergy)
      do ic=1,ncent
        do k=1,3
          da_psi(k,ic)=da_psi(k,ic)+da_psi_ref(k,ic)
          do i=1,nelec
            do j=1,i
              da_psi(k,ic)=da_psi(k,ic)+da_j(k,i,j,ic)
            enddo
          enddo
        enddo
      enddo
      return
      end
!-----------------------------------------------------------------------
      subroutine compute_da_psi(psid,da_psi_ref)

      use contrl_file, only: ounit
      use da_energy_now, only: da_psi
      use da_jastrow, only: da_j
      use da_orbval, only: da_orb
      use dorb_m,  only: iworbd
      use multidet, only: ivirt
      use multislater, only: detiab
      use precision_kinds, only: dp
      use slater,  only: kref,norb,slmi
      use system,  only: ncent,ncent_tot,ndn,nelec,nup
      use vmc_mod, only: norb_tot
      use zcompact, only: aaz,zmat
      use control, only: ipr

      implicit none

      integer :: i, iab, ic, ii, iorb
      integer :: irep, ish, j, jorb, l
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

! compute force for reference determinant
            do i=1,nel
              do j=1,nel
                da_psi_ref(k,ic)=da_psi_ref(k,ic)+slmi(i+(j-1)*nel,iab,1)*b_kref(i+(j-1)*nel)
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

! enddo iab
          enddo

          da_psi(k,ic)=trace*detiab(kref,1,1)*detiab(kref,2,1)/psid

        enddo
      enddo

!     do 800 ic=1,ncent
!       do 800 k=1,3
!         da_psi(k,ic)=da_psi(k,ic)+da_psi_ref(k,ic)
!         do 800 i=1,nelec
!800        da_psi(k,ic)=da_psi(k,ic)+da_j(k,i,ic)

!     if(ipr.gt.3) write(ounit,*)'da_ref',((da_psi_ref(l,ic),l=1,3),ic=1,ncent)
     if(ipr.gt.3) write(ounit,*)'da_psi',((da_psi(l,ic),l=1,3),ic=1,ncent)

      return
      end
!-----------------------------------------------------------------------
      subroutine compute_da_energy(psid,denergy)

      use Bloc,    only: b_da,xmat
      use constants, only: hb
      use da_energy_now, only: da_energy,da_psi
      use da_jastrow, only: da_d2j,da_vj
      use da_orbval, only: da_orb
      use da_pseudo, only: da_pecent,da_vps, da_pe_en
      use dorb_m,  only: iworbd
      use multidet, only: ivirt
      use multislater, only: detiab
      use precision_kinds, only: dp
      use pseudo,  only: lpot
      use slater,  only: kref,norb,slmi
      use system,  only: iwctype,ncent,ncent_tot,ndn,nelec,nup
      use velocity_jastrow, only: vj
      use zcompact, only: aaz,dzmat,emz,zmat
      use contrl_per, only: iperiodic
      use contrl_file,    only: ounit
      use control, only: ipr
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

! compute force for reference determinant
            do i=1,nel
              do j=1,nel
                jorb=iworbd(j+ish,kref)
                da_energy_ref(k,ic)=da_energy_ref(k,ic)+slmi(j+(i-1)*nel,iab,1)*b_da(k,i+ish,jorb,ic) &
                -da_orb(k,i+ish,jorb,ic)*xmat(i+(j-1)*nel,iab,1)
              enddo
            enddo
            do jrep=ivirt(iab),norb
              do irep=1,nel
                trace=trace+zmat(jrep,irep,iab,1)*b_da(k,irep+ish,jrep,ic) &
                           +dzmat(jrep,irep,iab,1)*da_orb(k,irep+ish,jrep,ic)
              enddo
            enddo

            do jrep=1,nel
              jorb=iworbd(jrep+ish,kref)
              do irep=1,nel
                trace=trace-emz(jrep,irep,iab,1)*da_orb(k,irep+ish,jorb,ic) &
                           -aaz(jrep,irep,iab,1)*b_da(k,irep+ish,jorb,ic)
              enddo
            enddo
! enddo iab
          enddo

          da_energy(k,ic)=trace*detiab(kref,1,1)*detiab(kref,2,1)/psid
        enddo
      enddo

      do ic=1,ncent
        ict=iwctype(ic)

        do k=1,3

          da_other_kin=da_d2j(k,ic)
          if (iperiodic.gt.0) then
            da_other_pot=da_pecent(k,ic) + da_pe_en(k,ic)
          else
            da_other_pot=da_pecent(k,ic)
          endif

          do i=1,nelec
            da_other_kin=da_other_kin + &
            +2*(vj(1,i,1)*da_vj(k,1,i,ic)+vj(2,i,1)*da_vj(k,2,i,ic)+vj(3,i,1)*da_vj(k,3,i,ic))
            da_other_pot=da_other_pot+da_vps(k,i,ic,lpot(ict))
          enddo
          da_energy(k,ic)=da_energy(k,ic)+da_energy_ref(k,ic)-hb*da_other_kin+da_other_pot &
                         -denergy*da_psi(k,ic)

! complete da_psi
        enddo
      enddo
      if(ipr.gt.3) then
    write(ounit,*)'da_energy',((da_energy(k,ic),k=1,3),ic=1,ncent)
    write(ounit,*)'da_energy',da_energy(2,1)
    write(ounit,*)'da_vj',da_vj(2,1,1,1)
    write(ounit,*)'da_d2j',da_d2j(2,1)
    write(ounit,*)'denergy',denergy
    write(ounit,*)'da_psi',da_psi(2,1)
      endif

      return
      end
!-----------------------------------------------------------------------
      subroutine force_analy_init(iflag)

      use da_energy_sumcum, only: da_energy_cm2,da_energy_cum
      use da_energy_sumcum, only: da_energy_sum,da_psi_cum,da_psi_sum
      use m_force_analytic, only: iforce_analy
      use system,  only: ncent
      use vd_mod, only: dmc_ivd, da_branch_sum, da_branch_cum
      use force_pth, only: PTH

      implicit none

      integer :: ic, iflag, k, iph

      if(iforce_analy.eq.0) return

      do iph=1,PTH
        do ic=1,ncent
           do k=1,3
              da_psi_sum(k,ic,iph)=0.0d0
              da_energy_sum(k,ic,iph)=0.0d0
              if (dmc_ivd.gt.0) da_branch_sum(k,ic,iph)=0.d0
           enddo
        enddo
      enddo

      if(iflag.gt.0) return

      do iph=1,PTH
        do ic=1,ncent
          do k=1,3
            da_psi_cum(k,ic,iph)=0.0d0
            da_energy_cum(k,ic,iph)=0.0d0
            da_energy_cm2(k,ic,iph)=0.0d0
            if (dmc_ivd.gt.0) da_branch_cum(k,ic,iph)=0.d0
          enddo
        enddo
      enddo

      return
      end

!-----------------------------------------------------------------------
      subroutine force_analy_sum(p,q,eloc,eloco)

      use da_energy_now, only: da_energy,da_psi
      use da_energy_sumcum, only: da_energy_sum,da_psi_sum
      use force_pth, only: PTH
      use m_force_analytic, only: iforce_analy
      use precision_kinds, only: dp
      use system,  only: ncent
      use vd_mod, only: dmc_ivd, da_branch_sum, da_branch
      use pathak_mod, only: ipathak, pnew

      implicit none

      integer :: ic, k, iph
      real(dp) :: eloc, eloco, p, q
      real(dp), dimension(3, ncent, PTH) :: da_energy_here

      if(iforce_analy.eq.0) return

      do iph=1,PTH
        do ic=1,ncent
          do k=1,3
            if (dmc_ivd.gt.0) then
              if (ipathak.gt.0) then
                da_energy_here(k,ic,iph)=da_energy(k,ic)*pnew(iph)+2*eloc*da_psi(k,ic)*pnew(iph)+eloc*da_branch(k,ic,iph)

                da_branch_sum(k,ic,iph) = da_branch_sum(k,ic,iph) + p*da_branch(k,ic,iph)
                da_psi_sum(k,ic,iph)=da_psi_sum(k,ic,iph)+p*da_psi(k,ic)*pnew(iph)
                da_energy_sum(k,ic,iph)= da_energy_sum(k,ic,iph)+p*da_energy_here(k,ic,iph)
              else
                da_energy_here(k,ic,iph)=da_energy(k,ic)+2*eloc*da_psi(k,ic)+eloc*da_branch(k,ic,iph)

                da_branch_sum(k,ic,iph) = da_branch_sum(k,ic,iph) + p*da_branch(k,ic,iph)
                da_psi_sum(k,ic,iph)= da_psi_sum(k,ic,iph)+p*da_psi(k,ic)
                da_energy_sum(k,ic,iph)= da_energy_sum(k,ic,iph)+p*da_energy_here(k,ic,iph)
              endif
            else
              if (ipathak.gt.0) then
                da_energy_here(k,ic,iph)=da_energy(k,ic)*pnew(iph)+2*eloc*da_psi(k,ic)*pnew(iph)

                da_psi_sum(k,ic,iph)=da_psi_sum(k,ic,iph)+p*da_psi(k,ic)*pnew(iph)
                da_energy_sum(k,ic,iph)=da_energy_sum(k,ic,iph)+p*da_energy_here(k,ic,iph)
              else
                da_energy_here(k,ic,iph)=da_energy(k,ic)+2*eloc*da_psi(k,ic)

                da_psi_sum(k,ic,iph)= da_psi_sum(k,ic,iph)+p*da_psi(k,ic)
                da_energy_sum(k,ic,iph)=da_energy_sum(k,ic,iph)+p*da_energy_here(k,ic,iph)
              endif
            endif
          enddo
        enddo
      enddo

      return
      end

!-----------------------------------------------------------------------

      subroutine force_analy_vd(ecutn, ecuto, e_cutoff, iw, iwmod)

      use branch, only: eold, eest
      use contrldmc, only: tau, ibranching_c, icut_e
      use da_energy_now, only: da_energy
      use force_pth, only: PTH
      use pathak_mod, only: ipathak, eps_pathak, pold, pnew, pathak
      use precision_kinds, only: dp
      use system,  only: ncent, nelec
      use vd_mod, only: da_branch, deriv_eold, esnake, ehist
      use velratio, only: fratio
      use da_energy_sumcum, only: da_energy_cum, da_energy_sum
      use estcum, only: wgcum
      use estsum, only: wgsum

      implicit none

      integer :: ic, k, iph
      integer :: iw, iwmod
      real(dp) :: ecuto, ecutn, e_cutoff
      real(dp) :: sqrt_pi_o2, sqrt_nelec, fratio_aux, fratio_aux2, da_av

      real(dp), dimension(3, ncent) :: deriv_energy_new

      real(dp), dimension(3, ncent) :: deriv_f_old
      real(dp), dimension(3, ncent) :: deriv_f_new

      real(dp), parameter :: zero = 0.d0
      real(dp), parameter :: one = 1.d0
      real(dp), parameter :: half = .5d0
      real(dp), parameter :: two = 2.d0
      real(dp), parameter :: small = 1.d-10

      sqrt_pi_o2 = 0.88622692545d0
      sqrt_nelec = dsqrt(dble(nelec))

      deriv_energy_new=da_energy
      
      if (icut_e.eq.3) then
        sqrt_nelec = 0.2d0*sqrt_nelec
        ! New first
        fratio_aux = max(ibranching_c * e_cutoff * dabs(eest - ecutn)/sqrt_nelec, 1e-9)
        do ic=1,ncent
          do k=1,3
            fratio_aux2 = (exp(-1* fratio_aux**2)/fratio_aux - ecuto/fratio_aux)*(ibranching_c * e_cutoff * ((eest - ecutn) * (- deriv_energy_new(k, ic))/dabs(eest - ecutn))/sqrt_nelec)
            deriv_f_new(k,ic) = deriv_energy_new(k,ic) * ecuto + ecutn * fratio_aux2 - eest * fratio_aux2 
          enddo
        enddo
        ! Now old
        fratio_aux = max(ibranching_c * e_cutoff * dabs(eest - eold(iw, 1))/sqrt_nelec, 1e-9)
        do ic=1,ncent
          do k=1,3
            fratio_aux2 = (exp(-1* fratio_aux**2)/fratio_aux - sqrt_pi_o2 * derf(fratio_aux)/fratio_aux/fratio_aux)*(ibranching_c * e_cutoff * ((eest - eold(iw, 1)) * (- deriv_eold(k,ic,iw))/dabs(eest - eold(iw, 1)))/sqrt_nelec)
            deriv_f_old(k,ic) = deriv_eold(k,ic,iw) * sqrt_pi_o2 * derf(fratio_aux)/fratio_aux + eold(iw,1) * fratio_aux2 - eest * fratio_aux2
          enddo
        enddo
      else
        if(dabs(ecutn-e_cutoff).lt.small) deriv_energy_new=zero
        if(dabs(ecuto-e_cutoff).lt.small) then
          do ic=1,ncent
            do k=1,3
              deriv_eold(k,ic,iw)=zero
            enddo
          enddo
        endif
        deriv_f_new(:,:)=deriv_energy_new(:,:)
        deriv_f_old(:,:)=deriv_eold(:,:,iw)
      endif

      do iph=1,PTH
        do ic=1,ncent
          do k=1,3
            if (ipathak.gt.0) then
              esnake(k,ic,iw,iph)=esnake(k,ic,iw,iph)+deriv_f_new(k,ic)*pnew(iph) &
      +deriv_f_old(k,ic)*pold(iw,iph)-ehist(k,ic,iw,iwmod,iph)
              ehist(k,ic,iw,iwmod,iph)=deriv_f_old(k,ic)*pold(iw,iph)+deriv_f_new(k,ic)*pnew(iph)

              da_branch(k,ic,iph)=-half*tau*esnake(k,ic,iw,iph)
            else
              esnake(k,ic,iw,iph)=esnake(k,ic,iw,iph)+deriv_f_new(k,ic) &
      +deriv_f_old(k,ic)-ehist(k,ic,iw,iwmod,iph)

              ehist(k,ic,iw,iwmod,iph)=deriv_f_old(k,ic)+deriv_f_new(k,ic)

              da_branch(k,ic,iph)=-half*tau*esnake(k,ic,iw,iph)
            endif
          enddo
        enddo
      enddo

      do ic=1,ncent
        do k=1,3
          deriv_eold(k,ic,iw)=deriv_energy_new(k,ic)
        enddo
      enddo
      if (ipathak.gt.0) then
        do iph=1,PTH
          pold(iw,iph)=pnew(iph)
        enddo
      endif

      return
      end
!-----------------------------------------------------------------------
      subroutine force_analy_cum(wsum,eave)

      use da_energy_sumcum, only: da_energy_cm2,da_energy_cum
      use da_energy_sumcum, only: da_energy_sum,da_psi_cum,da_psi_sum
      use force_pth, only: PTH
      use m_force_analytic, only: iforce_analy
      use precision_kinds, only: dp
      use system,  only: ncent
      use vd_mod, only: dmc_ivd, da_branch_sum, da_branch_cum
      use pathak_mod, only: ipathak, pnew

      implicit none

      integer :: ic, k, iph
      real(dp) :: eave, wcum, wsum
      real(dp), dimension(3, ncent, PTH) :: da_energy_now !(3,ncent,PTH)

      if(iforce_analy.eq.0) return

      do iph=1,PTH
        do ic=1,ncent
           do k=1,3
              if (dmc_ivd.gt.0) then
                da_energy_now(k,ic,iph)=(da_energy_sum(k,ic,iph)-2*eave*da_psi_sum(k,ic,iph)-eave*da_branch_sum(k,ic,iph))/wsum
                da_energy_cm2(k,ic,iph)=da_energy_cm2(k,ic,iph)+wsum*da_energy_now(k,ic,iph)**2

                da_branch_cum(k,ic,iph)=da_branch_cum(k,ic,iph)+da_branch_sum(k,ic,iph)
                da_psi_cum(k,ic,iph)=da_psi_cum(k,ic,iph)+da_psi_sum(k,ic,iph)
                da_energy_cum(k,ic,iph)=da_energy_cum(k,ic,iph)+da_energy_sum(k,ic,iph)
              else
                da_energy_now(k,ic,iph)=(da_energy_sum(k,ic,iph)-2*eave*da_psi_sum(k,ic,iph))/wsum
                da_energy_cm2(k,ic,iph)=da_energy_cm2(k,ic,iph)+wsum*da_energy_now(k,ic,iph)**2

                da_psi_cum(k,ic,iph)=da_psi_cum(k,ic,iph)+da_psi_sum(k,ic,iph)
                da_energy_cum(k,ic,iph)=da_energy_cum(k,ic,iph)+da_energy_sum(k,ic,iph)
              endif
            enddo
          enddo
        enddo

      return
      end
!-----------------------------------------------------------------------
      subroutine force_analy_fin(wcum,iblk_eff,eave)

      use da_energy_sumcum, only: da_energy_cm2,da_energy_cum,da_psi_cum
      use force_pth, only: PTH
      use m_force_analytic, only: da_energy_ave,iforce_analy
      use precision_kinds, only: dp
      use system,  only: ncent
      use vd_mod, only: dmc_ivd, da_branch_cum
      use pathak_mod, only: ipathak

      implicit none

      real(dp) :: iblk_eff
      integer :: iblk, ic, k, iph
      real(dp) :: eave, rtpass, wcum, x
      real(dp) :: x2
      real(dp), dimension(3) :: da_energy_err

      if(iforce_analy.eq.0) return

      rtpass=dsqrt(wcum)

      ! write(ounit, *) 'iblk_eff in force_analy_fin is ', iblk_eff

      open(80,file='force_analytic',form='formatted',status='unknown')
      do iph=1,PTH
        do ic=1,ncent
          do k=1,3
            if (dmc_ivd.gt.0) then
              da_energy_ave(k,ic,iph)=(da_energy_cum(k,ic,iph)-2*eave*da_psi_cum(k,ic,iph)-eave*da_branch_cum(k,ic,iph))/wcum
            else
              da_energy_ave(k,ic,iph)=(da_energy_cum(k,ic,iph)-2*eave*da_psi_cum(k,ic,iph))/wcum
            endif
            x = da_energy_ave(k,ic,iph)
            x2 = da_energy_cm2(k,ic,iph)
            da_energy_err(k)=dsqrt(abs(x2/wcum-x**2)/iblk_eff)
          enddo

          if (ipathak.gt.0) then
            write(80,'(i5,i5,1p6e14.5)')iph,ic,(da_energy_ave(k,ic,iph),k=1,3),(da_energy_err(k),k=1,3)
          else
            write(80,'(i5,1p6e14.5)') ic,(da_energy_ave(k,ic,iph),k=1,3),(da_energy_err(k),k=1,3)
          endif

        enddo
      enddo

       ! TODO JF this is included in the treatment of internal
       ! coordinates, remove this when finished
       !call transform_grad_zmat(da_energy_ave)

      return
      end
!-----------------------------------------------------------------------
      subroutine force_analy_dump(iu)

      use da_energy_sumcum, only: da_energy_cm2,da_energy_cum,da_psi_cum
      use force_pth, only: PTH
      use m_force_analytic, only: iforce_analy
      use system,  only: ncent

      implicit none

      integer :: ic, iu, k, iph

      if(iforce_analy.eq.0) return

      write(iu) (((da_energy_cum(k,ic,iph),da_psi_cum(k,ic,iph),da_energy_cm2(k,ic,iph),k=1,3),ic=1,ncent),iph=1,PTH)

      return
      end
!-----------------------------------------------------------------------
      subroutine force_analy_rstrt(iu)

      use da_energy_sumcum, only: da_energy_cm2,da_energy_cum,da_psi_cum
      use force_pth, only: PTH
      use m_force_analytic, only: iforce_analy
      use system,  only: ncent

      implicit none

      integer :: ic, iu, k, iph

      if(iforce_analy.eq.0) return

      read(iu) (((da_energy_cum(k,ic,iph),da_psi_cum(k,ic,iph),da_energy_cm2(k,ic,iph),k=1,3),ic=1,ncent),iph=1,PTH)

      return
      end
!-----------------------------------------------------------------------
      subroutine force_analy_save
      implicit none

      return
      end
end module
