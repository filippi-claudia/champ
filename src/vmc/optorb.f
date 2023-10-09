      module optorb_f_mod
      use error,   only: fatal_error

      contains
      subroutine optorb_deriv(psid,denergy,zmat,dzmat,emz,aaz,orbprim,eorbprim,istate)

      use vmc_mod, only: norb_tot, stoo, stobjx
      use system, only: ndn, nup
      use slater, only: kref
      use multidet, only: ivirt
      use optwf_control, only: ioptorb
      use orb_mat_022, only: ideriv
      use orb_mat_033, only: ideriv_ref,irepcol_ref
      use orbval,  only: orb      
      use Bloc, only: b, tildem
      use multimat, only: aa
      use multislater, only: detiab
      use optorb_cblock, only: norbterm
      use system, only: nelec
      use precision_kinds, only: dp
      use slater,  only: kref
      use system,  only: ndn,nelec,nup
      use vmc_mod, only: norb_tot

      implicit none

      integer :: i, iab, io, irep, ish, istate
      integer :: iterm, jo, nel, o, x
      real(dp) :: denergy, detratio, dorb_energy, dorb_energy_ref, dorb_psi
      real(dp) :: dorb_psi_ref, psid
      real(dp), dimension(norb_tot, nelec, 2) :: zmat
      real(dp), dimension(norb_tot, nelec, 2) :: dzmat
      real(dp), dimension(nelec, nelec, 2) :: emz
      real(dp), dimension(nelec, nelec, 2) :: aaz
      real(dp), dimension(*) :: orbprim
      real(dp), dimension(*) :: eorbprim


      if(ioptorb.eq.0) return


c     ns_current=ns_current+1
c     if(ns_current.ne.iorbsample) return
c ns_current reset in optorb_sum

      o=stoo(istate)
      x=stobjx(istate)
      detratio=detiab(kref,1,o)*detiab(kref,2,o)/psid
      do iterm=1,norbterm

        io=ideriv(1,iterm)
        jo=ideriv(2,iterm)

        dorb_psi_ref=0.0d0
        dorb_energy_ref=0.d0

        dorb_psi=0.d0
        dorb_energy=0.d0
        do iab=1,2

          if(iab.eq.1) then
            ish=0
            nel=nup
           else
            ish=nup
            nel=ndn
          endif

          if(io.ge.ivirt(iab)) then
            do i=1,nel
              dorb_psi=dorb_psi+zmat(io,i,iab)*orb(i+ish,jo,o)
              dorb_energy=dorb_energy+dzmat(io,i,iab)*orb(i+ish,jo,o)+zmat(io,i,iab)*b(jo,i+ish,x)
            enddo
          endif
          if(ideriv_ref(iterm,iab).gt.0) then
            irep=irepcol_ref(iterm,iab)

            dorb_psi_ref=dorb_psi_ref+aa(irep,jo,iab,o)
            dorb_energy_ref=dorb_energy_ref+tildem(irep,jo,iab,x)

            do i=1,nel
              dorb_psi=dorb_psi-aaz(irep,i,iab)*orb(i+ish,jo,o)
              dorb_energy=dorb_energy-emz(irep,i,iab)*orb(i+ish,jo,o)-aaz(irep,i,iab)*b(jo,i+ish,x)
            enddo
          endif

        enddo

        orbprim(iterm)=dorb_psi*detratio
        eorbprim(iterm)=dorb_energy*detratio+dorb_energy_ref-denergy*orbprim(iterm)

        orbprim(iterm)=orbprim(iterm)+dorb_psi_ref

      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine optorb_compute(psid,eloc,deloc)

      use contrl_file, only: ounit
      use csfs,    only: nstates
      use optorb_cblock, only: norbterm
      use optwf_control, only: ioptorb
      use orb_mat_001, only: orb_ho,orb_o,orb_oe
      use precision_kinds, only: dp
      use zcompact, only: aaz,dzmat,emz,zmat
      implicit none

      integer :: i, istate
      real(dp), dimension(*) :: psid
      real(dp), dimension(*) :: eloc
      real(dp), dimension(*) :: deloc

      if(ioptorb.eq.0) return

      do istate=1,nstates

        call optorb_deriv(psid(istate),deloc(istate)
     &   ,zmat(1,1,1,istate),dzmat(1,1,1,istate),emz(1,1,1,istate),aaz(1,1,1,istate)
     &   ,orb_o(1,istate),orb_ho(1,istate),istate)

        do i=1,norbterm
            orb_oe(i,istate)=orb_o(i,istate)*eloc(istate)
            orb_ho(i,istate)=orb_ho(i,istate)+eloc(istate)*orb_o(i,istate)
        enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine optorb_sum(wtg_new,wtg_old,enew,eold,iflag)

      use csfs,    only: nstates
      use optorb_cblock, only: isample_cmat,norbterm,nreduced
      use optwf_control, only: iapprox,ioptorb
      use orb_mat_001, only: orb_ho,orb_o,orb_oe
      use orb_mat_002, only: orb_ho_old,orb_o_old,orb_oe_old
      use orb_mat_003, only: orb_o_sum
      use orb_mat_004, only: orb_oe_sum
      use orb_mat_005, only: orb_ho_cum
      use orb_mat_006, only: orb_oo_cum
      use orb_mat_007, only: orb_oho_cum
      use orb_mat_030, only: orb_ecum,orb_wcum
      use precision_kinds, only: dp

      implicit none

      integer :: i, idiag_only, idx, ie, iflag
      integer :: istate, j, je
      real(dp) :: p, q
      real(dp), dimension(*) :: wtg_new
      real(dp), dimension(*) :: wtg_old
      real(dp), dimension(*) :: enew
      real(dp), dimension(*) :: eold

      if(ioptorb.eq.0) return

c     if(ns_current.ne.iorbsample) return
c ns_current reset
c     ns_current=0

      idiag_only=0
      if(iapprox.gt.0) idiag_only=1

      do istate=1,nstates

      p=wtg_new(istate)

      do i=1,norbterm
       orb_o_sum(i,istate)=orb_o_sum(i,istate)+p*orb_o(i,istate)
       orb_oe_sum(i,istate) =orb_oe_sum(i,istate)+p*orb_oe(i,istate)
       orb_ho_cum(i,istate) =orb_ho_cum(i,istate)+p*orb_ho(i,istate)
      enddo

      orb_wcum(istate)=orb_wcum(istate)+p
      orb_ecum(istate)=orb_ecum(istate)+p*enew(istate)

      if(isample_cmat.eq.0) go to 200

      if(idiag_only.eq.0) then
        idx=0
        do i=1,nreduced
         ie=i
         do j=1,i
          idx=idx+1
          je=j
          orb_oo_cum(idx,istate)=orb_oo_cum(idx,istate)+p*orb_o(ie,istate)*orb_o(je,istate)
         enddo
        enddo

        idx=0
        do i=1,nreduced
         ie=i
         do j=1,nreduced
          idx=idx+1
          je=j
          orb_oho_cum(idx,istate)=orb_oho_cum(idx,istate)+p*orb_o(je,istate)*orb_ho(ie,istate)
         enddo
        enddo
       else
        do i=1,nreduced
          ie=i
          orb_oo_cum(i,istate)=orb_oo_cum(i,istate)+p*orb_o(ie,istate)*orb_o(ie,istate)
          orb_oho_cum(i,istate)=orb_oho_cum(i,istate)+p*orb_o(ie,istate)*orb_ho(ie,istate)
        enddo
      endif

  200 continue
      enddo

      if(iflag.eq.0) return

      do istate=1,nstates

      q=wtg_old(istate)

      do i=1,norbterm
       orb_o_sum(i,istate)=orb_o_sum(i,istate)+q*orb_o_old(i,istate)
       orb_oe_sum(i,istate) =orb_oe_sum(i,istate)+q*orb_oe_old(i,istate)
      enddo

      orb_wcum(istate)=orb_wcum(istate)+q
      orb_ecum(istate)=orb_ecum(istate)+q*eold(istate)

      if(isample_cmat.eq.0) go to 300

      if(idiag_only.eq.0) then
        idx=0
        do i=1,nreduced
         ie=i
         do j=1,i
          idx=idx+1
          je=j
          orb_oo_cum(idx,istate)=orb_oo_cum(idx,istate)+q*orb_o_old(ie,istate)*orb_o_old(je,istate)
         enddo
        enddo

        idx=0
        do i=1,nreduced
         ie=i
         do j=1,nreduced
          idx=idx+1
          je=j
          orb_oho_cum(idx,istate)=orb_oho_cum(idx,istate)+q*orb_o_old(je,istate)*orb_ho_old(ie,istate)
         enddo
        enddo
       else
        do i=1,nreduced
          ie=i
          orb_oo_cum(i,istate)=orb_oo_cum(i,istate)+q*orb_o_old(ie,istate)*orb_o_old(ie,istate)
          orb_oho_cum(i,istate)=orb_oho_cum(i,istate)+q*orb_o_old(ie,istate)*orb_ho_old(ie,istate)
        enddo
      endif

  300 continue
      enddo
      end
c-----------------------------------------------------------------------
      subroutine optorb_cum(wsum,esum)

      use csfs,    only: nstates
      use optorb_cblock, only: idump_blockav,nb_current,nefp_blocks
      use optorb_cblock, only: norb_f_bcum,norbterm
      use optwf_control, only: ioptorb
      use orb_mat_003, only: orb_o_cum,orb_o_sum
      use orb_mat_004, only: orb_oe_cum,orb_oe_sum
      use orb_mat_024, only: orb_e_bsum,orb_f_bcm2,orb_f_bcum,orb_o_bsum
      use orb_mat_024, only: orb_oe_bsum,orb_w_bsum
      use precision_kinds, only: dp

      implicit none

      integer :: i, istate
      real(dp) :: eb, fnow
      real(dp), dimension(*) :: wsum
      real(dp), dimension(*) :: esum

      if(ioptorb.eq.0) return

      nb_current=nb_current+1

      do istate=1,nstates

      orb_e_bsum(istate)=orb_e_bsum(istate)+esum(istate)
      orb_w_bsum(istate)=orb_w_bsum(istate)+wsum(istate)
      do i=1,norbterm
       orb_o_bsum(i,istate)=orb_o_bsum(i,istate)+orb_o_sum(i,istate)
       orb_oe_bsum(i,istate)=orb_oe_bsum(i,istate)+orb_oe_sum(i,istate)
      enddo

      if(nb_current.eq.nefp_blocks)then
       eb=orb_e_bsum(istate)/orb_w_bsum(istate)

       do i=1,norbterm
         fnow=orb_oe_bsum(i,istate)/orb_w_bsum(istate)-orb_o_bsum(i,istate)/orb_w_bsum(istate)*eb
         orb_f_bcum(i,istate)=orb_f_bcum(i,istate)+fnow
         orb_f_bcm2(i,istate)=orb_f_bcm2(i,istate)+fnow**2
       enddo

       orb_e_bsum(istate)=0.d0
       orb_w_bsum(istate)=0.d0
       do i=1,norbterm
        orb_o_bsum(i,istate)=0.d0
        orb_oe_bsum(i,istate)=0.d0
       enddo
      endif

      do i=1,norbterm
       orb_o_cum(i,istate)=orb_o_cum(i,istate) + orb_o_sum(i,istate)
       orb_oe_cum(i,istate)=orb_oe_cum(i,istate) + orb_oe_sum(i,istate)
      enddo

  200 continue
      enddo

      if(nb_current.eq.nefp_blocks) then
        nb_current=0
        norb_f_bcum=norb_f_bcum+1
      endif

      if(idump_blockav.ne.0)then
       write(idump_blockav) esum(1)/wsum(1),(orb_o_sum(i,1)/wsum(1),orb_oe_sum(i,1)/wsum(1),i=1,norbterm)
      endif

      end

c-----------------------------------------------------------------------
      subroutine optorb_init(iflg)

      use csfs,    only: nstates
      use optorb_cblock, only: isample_cmat,nb_current,norb_f_bcum
      use optorb_cblock, only: norbterm,nreduced
      use optwf_control, only: iapprox,ioptorb
      use orb_mat_003, only: orb_o_cum,orb_o_sum
      use orb_mat_004, only: orb_oe_cum,orb_oe_sum
      use orb_mat_005, only: orb_ho_cum
      use orb_mat_006, only: orb_oo_cum
      use orb_mat_007, only: orb_oho_cum
      use orb_mat_024, only: orb_e_bsum,orb_f_bcm2,orb_f_bcum,orb_o_bsum
      use orb_mat_024, only: orb_oe_bsum,orb_w_bsum
      use orb_mat_030, only: orb_ecum,orb_wcum

      implicit none

      integer :: i, idiag_only, idx, iflg, istate
      integer :: j, ns_current

      if(ioptorb.eq.0) return

      idiag_only=0
      if(iapprox.gt.0) idiag_only=1

      do istate=1,nstates

      do i=1,norbterm
       orb_o_sum(i,istate)=0.d0
       orb_oe_sum(i,istate) =0.d0
       orb_o_bsum(i,istate)=0.d0
       orb_oe_bsum(i,istate)=0.d0
      enddo
      orb_e_bsum(istate)=0.d0
      orb_w_bsum(istate)=0.d0

      enddo
c     $ iflg = 0: init *cum, *cm2 as well
      if(iflg.gt.0) return

      ns_current=0
      nb_current=0
      norb_f_bcum=0

      do istate=1,nstates

      do i=1,norbterm
       orb_o_cum(i,istate)=0.d0
       orb_oe_cum(i,istate) =0.d0
       orb_ho_cum(i,istate) =0.d0
       orb_f_bcum(i,istate)=0.d0
       orb_f_bcm2(i,istate)=0.d0
      enddo
      orb_wcum(istate)=0.d0
      orb_ecum(istate)=0.d0

      if(isample_cmat.ne.0) then
       if(idiag_only.eq.0) then
         idx=0
         do i=1,nreduced
          do j=1,i
           idx=idx+1
           orb_oo_cum(idx,istate)=0.d0
          enddo
         enddo

         idx=0
         do i=1,nreduced
          do j=1,nreduced
           idx=idx+1
           orb_oho_cum(idx,istate)=0.d0
          enddo
         enddo
       else
         do i=1,nreduced
           orb_oo_cum(i,istate)=0.d0
           orb_oho_cum(i,istate)=0.d0
         enddo
       endif
      endif

  200 continue
      enddo

      end
c-----------------------------------------------------------------------
      subroutine optorb_save

      use csfs,    only: nstates
      use optorb_cblock, only: norbterm
      use optwf_control, only: ioptorb
      use orb_mat_001, only: orb_ho,orb_o,orb_oe
      use orb_mat_002, only: orb_ho_old,orb_o_old,orb_oe_old

      implicit none

      integer :: i, istate

      if(ioptorb.eq.0) return

      do istate=1,nstates

      do i=1,norbterm
       orb_o_old(i,istate)=orb_o(i,istate)
       orb_oe_old(i,istate)=orb_oe(i,istate)
       orb_ho_old(i,istate)=orb_ho(i,istate)
      enddo

  200 continue
      enddo

      end
c-----------------------------------------------------------------------
      subroutine optorb_restore

      use csfs,    only: nstates
      use optorb_cblock, only: norbterm
      use optwf_control, only: ioptorb
      use orb_mat_001, only: orb_ho,orb_o,orb_oe
      use orb_mat_002, only: orb_ho_old,orb_o_old,orb_oe_old

      implicit none

      integer :: i, istate

      if(ioptorb.eq.0) return

      do istate=1,nstates

      do i=1,norbterm
       orb_o(i,istate)=orb_o_old(i,istate)
       orb_oe(i,istate)=orb_oe_old(i,istate)
       orb_ho(i,istate)=orb_ho_old(i,istate)
      enddo

  200 continue
      enddo

      end
c-----------------------------------------------------------------------
      subroutine optorb_avrg(wcum,eave,oav,eoav,fo,foerr,istate)

      use contrl_file, only: ounit
      use optorb_cblock, only: norb_f_bcum,norbterm
      use optwf_control, only: ioptorb
      use orb_mat_003, only: orb_o_cum
      use orb_mat_004, only: orb_oe_cum
      use orb_mat_024, only: orb_f_bcm2,orb_f_bcum
      use precision_kinds, only: dp

      implicit none

      integer :: i, istate, n
      real(dp) :: dble, eave, errn, wcum
      real(dp) :: x, x2
      real(dp), dimension(*) :: oav
      real(dp), dimension(*) :: eoav
      real(dp), dimension(*) :: fo
      real(dp), dimension(*) :: foerr


      if(ioptorb.eq.0) return

      do i=1,norbterm
        oav(i)=orb_o_cum(i,istate)/wcum
        eoav(i)=orb_oe_cum(i,istate)/wcum
        fo(i)=eoav(i)-eave*oav(i)
        x = orb_f_bcum(i,istate)
        x2= orb_f_bcm2(i,istate)
        n = norb_f_bcum
        errn = dsqrt(dabs(x2/dble(n)-(x/dble(n))**2)/dble(n))
        foerr(i)=errn
      enddo

      write(ounit,'(''ORB-PT: forces collected'',i4)') norb_f_bcum
      end subroutine
c-----------------------------------------------------------------------
      subroutine optorb_dump(iu)

      use csfs,    only: nstates
      use optorb_cblock, only: nefp_blocks,norb_f_bcum,norbprim,norbterm
      use optorb_cblock, only: nreduced
      use optwf_control, only: iapprox,ioptorb
      use orb_mat_003, only: orb_o_cum
      use orb_mat_004, only: orb_oe_cum
      use orb_mat_005, only: orb_ho_cum
      use orb_mat_006, only: orb_oo_cum
      use orb_mat_007, only: orb_oho_cum
      use orb_mat_024, only: orb_f_bcm2,orb_f_bcum
      use orb_mat_030, only: orb_ecum,orb_wcum

      implicit none

      integer :: i, istate, iu, matdim


      if(ioptorb.eq.0) return

      matdim=nreduced*(nreduced+1)/2
      if(iapprox.gt.0) matdim=nreduced

      write(iu) norbprim,norbterm,nreduced
      write(iu) nefp_blocks,norb_f_bcum
      do istate=1,nstates
      write(iu) (orb_o_cum(i,istate),i=1,norbterm)
      write(iu) (orb_oe_cum(i,istate),i=1,norbterm)
      write(iu) (orb_ho_cum(i,istate),i=1,norbterm)
      write(iu) (orb_f_bcum(i,istate),orb_f_bcm2(i,istate),i=1,norbterm)
      write(iu) (orb_oo_cum(i,istate),i=1,matdim)
      write(iu) (orb_oho_cum(i,istate),i=1,nreduced*nreduced)
      write(iu) orb_wcum(istate),orb_ecum(istate)
  200 continue
      enddo

      end
c-----------------------------------------------------------------------
      subroutine optorb_rstrt(iu)

      use contrl_file, only: ounit
      use csfs,    only: nstates
      use optorb_cblock, only: nefp_blocks,norb_f_bcum,norbprim,norbterm
      use optorb_cblock, only: nreduced
      use optwf_control, only: iapprox,ioptorb
      use orb_mat_003, only: orb_o_cum
      use orb_mat_004, only: orb_oe_cum
      use orb_mat_005, only: orb_ho_cum
      use orb_mat_006, only: orb_oo_cum
      use orb_mat_007, only: orb_oho_cum
      use orb_mat_024, only: orb_f_bcm2,orb_f_bcum
      use orb_mat_030, only: orb_ecum,orb_wcum

      implicit none

      integer :: i, istate, iu, matdim, morbprim
      integer :: morbterm, mreduced


      if(ioptorb.eq.0) return
      read(iu) morbprim,morbterm,mreduced
      if(morbprim.ne.norbprim) then
       write(ounit,*) 'wrong number of primitive orb terms!'
       write(ounit,*) 'old ',morbprim,' new ',norbprim
       call fatal_error('OPTORB_RSTRT: Restart, inconsistent ORB information')
      endif
      if(morbterm.ne.norbterm) then
       write(ounit,*) 'wrong number of orb terms!'
       write(ounit,*) 'old ',morbterm,' new ',norbterm
       call fatal_error('OPTORB_RSTRT: Restart, inconsistent ORB information')
      endif

c nreduced has to be set since it will only be known for non-continuation runs
      nreduced=mreduced
      matdim=nreduced*(nreduced+1)/2
      if(iapprox.gt.0) matdim=nreduced

      read(iu) nefp_blocks,norb_f_bcum

      do istate=1,nstates
      read(iu) (orb_o_cum(i,istate),i=1,norbterm)
      read(iu) (orb_oe_cum(i,istate),i=1,norbterm)
      read(iu) (orb_ho_cum(i,istate),i=1,norbterm)
      read(iu) (orb_f_bcum(i,istate),orb_f_bcm2(i,istate),i=1,norbterm)
      read(iu) (orb_oo_cum(i,istate),i=1,matdim)
      read(iu) (orb_oho_cum(i,istate),i=1,nreduced*nreduced)
      read(iu) orb_wcum,orb_ecum
  200 continue
      enddo
      end
c-----------------------------------------------------------------------
      subroutine optorb_fin(wcum,ecum)

      use ci000,   only: nciterm
      use csfs,    only: nstates
      use gradhess_all, only: grad,h,s
      use optorb_cblock, only: idump_blockav,norbterm,nreduced
      use optwf_control, only: iapprox,ioptorb,iuse_orbeigv,method
      use optwf_parms, only: nparmd,nparmj
      use orb_mat_005, only: orb_ho_cum
      use orb_mat_006, only: orb_oo_cum
      use orb_mat_007, only: orb_oho_cum
      use precision_kinds, only: dp
      use sa_weights, only: weights

      implicit none

      integer :: i, i0, i1, idx, ish
      integer :: istate, j
      real(dp) :: eave, orb_oho, orb_oo, passes, passesi
      real(dp) :: wts
      real(dp), dimension(norbterm) :: oav
      real(dp), dimension(norbterm) :: eoav
      real(dp), dimension(norbterm) :: fo
      real(dp), dimension(norbterm) :: foerr
      real(dp), dimension(*) :: wcum
      real(dp), dimension(*) :: ecum

      if(ioptorb.eq.0.or.method.eq.'sr_n'.or.method.eq.'lin_d') return

      nparmd=max(nciterm-1,0)
      ish=nparmj+nparmd
      if(method.eq.'linear') ish=ish+1

      s(1,1)=0
      h(1,1)=0
      do j=1,nreduced
        grad(j+ish)=0
        s(j+ish,1)=0
        h(j+ish,1)=0
        s(1,j+ish)=0
        h(1,j+ish)=0
        do i=1,nreduced
          s(i+ish,j+ish)=0
          h(i+ish,j+ish)=0
        enddo
      enddo

      do istate=1,nstates

      wts=weights(istate)

      passes=wcum(istate)
      passesi=1/passes
      eave=ecum(istate)*passesi

c     if(iorbsample.ne.1) then
c       passes=orb_wcum(istate)
c       passesi=1/passes
c       eave=orb_ecum(istate)*passesi
c     endif

      call optorb_avrg(passes,eave,oav(1),eoav(1),fo(1),foerr(1),istate)

c Hessian method
      if(method.eq.'hessian') then

        if(iuse_orbeigv.eq.0) then
c Formulas for exact orbital hessian not implemented
          call fatal_error('OPTORB_FIN: formulas for exact hessian not implemented')
        endif

c Linear method
       elseif(method.eq.'linear') then

        s(1,1)=1
        h(1,1)=h(1,1)+wts*eave
c Exact Hamiltonian
        if(iuse_orbeigv.eq.0) then

c Hamiltonian on semi-orthogonal basis
        idx=0
        do i=1,nreduced
          s(i+ish,1)=0
          s(1,i+ish)=0
          h(i+ish,1)=h(i+ish,1)+wts*(eoav(i)-eave*oav(i))
          h(1,i+ish)=h(1,i+ish)+wts*(orb_ho_cum(i,istate)*passesi-eave*oav(i))
c         write(ounit,*) 'H',wts,eoav(i)-eave*oav(i),orb_ho_cum(i,istate)*passesi-eave*oav(i)
          i0=1
          if(iapprox.gt.0) i0=i
          do j=i0,i
            idx=idx+1
            orb_oo=orb_oo_cum(idx,istate)*passesi-oav(i)*oav(j)
            s(i+ish,j+ish)=s(i+ish,j+ish)+wts*orb_oo
            s(j+ish,i+ish)=s(i+ish,j+ish)
          enddo
        enddo

        i0=1
        i1=nreduced
        idx=0
        do i=1,nreduced
          if(iapprox.gt.0) then
            i0=i
            i1=i
          endif
          do j=i0,i1
            idx=idx+1
            orb_oho=(orb_oho_cum(idx,istate)-oav(j)*orb_ho_cum(i,istate))*passesi
     &             -oav(i)*eoav(j)+eave*oav(i)*oav(j)
            h(j+ish,i+ish)=h(j+ish,i+ish)+wts*orb_oho
          enddo
        enddo

       endif

c Perturbative method
       elseif(method.eq.'perturbative') then

        if(iuse_orbeigv.eq.0) then
c Formulas for exact orbital perturbative not implemented
          call fatal_error('OPTORB_FIN: formulas for exact perturbative not implemented')
         else
          do i=1,nreduced
            grad(i)=grad(i)+wts*fo(i)
          enddo
          idx=0
          do i=1,nreduced
            do j=1,i
              idx=idx+1
              s(i,j)=s(i,j)+wts*(orb_oo_cum(idx,istate)*passesi-oav(i)*oav(j))
              s(j,i)=s(i,j)
            enddo
          enddo
        endif
      endif

  200 continue
      enddo

c Approximations on matrix elements
      if(method.eq.'linear') then
        if(iapprox.gt.0) then
          do i=1,nreduced
            do j=1,i-1
              s(i+ish,j+ish)=0
              s(j+ish,i+ish)=0
              h(i+ish,j+ish)=0
              h(j+ish,i+ish)=0
            enddo
          enddo
          if(iapprox.eq.2) then
            do i=1,nreduced
              h(1,i+ish)=h(i+ish,1)
            enddo
          endif
         elseif(iapprox.lt.0) then
          if(iapprox.eq.-1) then
            do i=1,nreduced
              h(1,i+ish)=h(i+ish,1)
            enddo
           elseif(iapprox.eq.-2) then
            do i=1,nreduced
              h(1,i+ish)=h(i+ish,1)
              do j=1,i-1
                h(i+ish,j+ish)=0.5*(h(i+ish,j+ish)+h(j+ish,i+ish))
                h(j+ish,i+ish)=h(i+ish,j+ish)
              enddo
            enddo
           elseif(iapprox.eq.-3) then
            do i=1,nreduced
              h(1,i+ish)=0.5*(h(i+ish,1)+h(1,i+ish))
              h(i+ish,1)=h(1,i+ish)
              do j=1,i-1
                h(i+ish,j+ish)=0.5*(h(i+ish,j+ish)+h(j+ish,i+ish))
                h(j+ish,i+ish)=h(i+ish,j+ish)
              enddo
            enddo
          endif
        endif
       elseif(method.eq.'perturbative') then
c Approximation: diagonal perturbative approach
        if(iapprox.gt.0) then
          do i=1,nreduced
            do j=1,i-1
              s(j,i)=0
              s(i,j)=0
            enddo
          enddo
        endif
      endif

      if(idump_blockav.ne.0) close(idump_blockav)

      end
c-----------------------------------------------------------------------
      subroutine detratio_col(nel,orb,icol,sinvt,ratio,isltnew)

      use precision_kinds, only: dp

      implicit none

      integer :: icol, ie, isltnew, jcol, je
      integer :: nel
      real(dp) :: ratio, sum
      real(dp), dimension(nel) :: orb
      real(dp), dimension(nel, nel) :: sinvt

c values of new orbital
c inverse transposed slater matrix (first index electron, 2nd orbital)
c compute ratio of new and old determinant, if isltnew is
c not zero, update inverse slater matrix as well
c the new determinant differs from the old by replacing column icol
c with the orbital values in orb

      ratio=0.d0
      do ie=1,nel
       ratio=ratio+sinvt(icol,ie)*orb(ie)
      enddo
      if(isltnew.gt.0) then
c matrix except replaced column
       do jcol=1,nel
        if(jcol.ne.icol) then
         sum=0.d0
         do je=1,nel
          sum=sum+orb(je)*sinvt(jcol,je)
         enddo
         sum=sum/ratio
         do je=1,nel
          sinvt(jcol,je)=sinvt(jcol,je)-sum*sinvt(icol,je)
         enddo
        endif
       enddo
c replaced column
       do ie=1,nel
        sinvt(icol,ie)=sinvt(icol,ie)/ratio
       enddo
      endif

      end
c-----------------------------------------------------------------------
      subroutine optorb_define

      use coefs,   only: next_max
      use contrl_file, only: errunit,ounit
      use dorb_m,  only: iworbd
      use optorb,  only: irrep
      use optorb_cblock, only: norbterm,nreduced
      use optorb_mix, only: iwmix_virt,norbopt,norbvirt
      use optorb_mod, only: mxreduced
      use optwf_control, only: method,ncore,no_active
      use orb_mat_022, only: ideriv
      use orb_mat_033, only: ideriv_iab,ideriv_ref,irepcol_ref
      use orbval,  only: nadorb,ndetorb,orb
      use slater,  only: kref,ndet,norb
      use system,  only: ndn,nelec,nup
      use vmc_mod, only: norb_tot

      implicit none

      integer :: i, iab, icount_orbdef, ie, iesave
      integer :: io, iocc, iprt, iterm
      integer :: j, jo, k, n0
      integer :: n1, noporb
      integer, dimension(2, ndet) :: iodet
      integer, dimension(2, ndet) :: iopos
      integer, dimension(2, norb_tot) :: iflag
      integer, dimension(2) :: ne
      integer, dimension(2) :: m

      data icount_orbdef /1/

      save icount_orbdef

      iprt=3

      ndn=nelec-nup

      ne(1)=nup
      ne(2)=nelec
c orbital indices in determinants of trial wave function
      ndetorb=0

      do i=1,ndet
       do j=1,nelec
        if(iworbd(j,i).gt.norb) then
         write(ounit,1) i,j,iworbd(j,i),norb
         call fatal_error('VERIFY: orbital index out of range')
        endif
        if(iworbd(j,i).gt.ndetorb) then
         ndetorb=iworbd(j,i)
        endif
       enddo
      enddo
  1   format('Det ',i4,' column ',i4,' orb index ',i4,' norb ',i4)

c Number of external orbitals for orbital optimization
      next_max=norb-ndetorb
      if(nadorb.gt.next_max) nadorb=next_max
       write(ounit, '(a, t40, i0)' ) 'norb', norb
       write(ounit, '(a, t40, i0)') 'nadorb', nadorb
       write(ounit, '(a, t40, i0)') 'ndet_orb', ndetorb
       write(ounit, '(a, t40, i0)') 'next_max', next_max
      ! call fatal_error('optorb.f')

      if(iprt.gt.0) then
       write(ounit,'(''Determinantal orbitals in orbital optimization: '',i0)') ndetorb
       write(ounit,'(''External orbitals in orbital optimization:      '',i0)') nadorb
       write(ounit,'(''Total orbitals in orbital optimization:         '',i0)') nadorb+ndetorb-ncore
      endif
      norb=ndetorb


c Omit doubly occupied in all input determinants
      do i=1,ndetorb
        iflag(1,i)=0
        do k=1,ndet
          iocc=0
          do j=1,nelec
            if(iworbd(j,k).eq.i) iocc=iocc+1
          enddo
          if(iocc.ne.2) then
            iflag(1,i)=1
            goto 5
          endif
        enddo
   5  continue
      enddo

c Omit empty orbitals

      do i=1,ndetorb
       iflag(2,i)=0
       do k=1,ndet
        do j=1,nelec
          if(iworbd(j,k).eq.i) iflag(2,i)=1
        enddo
       enddo
      enddo
      do i=ndetorb+1,ndetorb+nadorb
       iflag(1,i)=1
       iflag(2,i)=0
      enddo

      if(norbopt.eq.0.or.norbvirt.eq.0) then
        do io=1,ndetorb
         do jo=ncore+1,ndetorb+nadorb
          iwmix_virt(io,jo)=jo
         enddo
        enddo
      elseif(norbopt.ne.ndetorb.or.norbvirt.lt.nadorb) then
       write(ounit,'(''OPTORB_DEFINE: norbopt,ndetorb'',2i6)') norbopt,ndetorb
       write(ounit,'(''OPTORB_DEFINE: noptvirt,nadorb'',2i6)') norbvirt,nadorb
       call fatal_error('OPTORB_DEFINE: Mixvirt block, inconsistent')
      endif



c Orbital variation io -> io+a*jo
c io: occupied orbitals in twf
c jo: all orbitals
c omitted if not same symmetry, or io empty, or both doubly occupied
      noporb=0
      iterm=0

      if(iprt.gt.2) then
       write(ounit,*) '(''=========== orbital pair list =========='')'
      endif

      do io=ncore+1,ndetorb
c Omit empty orbitals
       if(iflag(2,io).eq.0) goto 60
       do jo=ncore+1,ndetorb+nadorb
c Omit if io and jo are the same
        if(io.eq.jo) goto 50
c Omit if io and jo have different symmetry
        if(irrep(io).ne.irrep(jo)) goto 50
c Omit if io and jo are both doubly occupied in all determinants
        if((iflag(1,io).eq.0).and.(iflag(1,jo).eq.0)) goto 50
c Omit if io and jo are both active orbitals
        if(no_active.ne.0.and.iflag(1,io).ne.0.and.iflag(2,jo).ne.0) goto 50
c Omit if we only want to mix according to the table mixvirt
        if(iwmix_virt(io,jo).eq.0) goto 50
c Include: io is occupied in some determinant and jo not
        do iab=1,2
          n0=0
          n1=nup
          if(iab.eq.2) then
            n0=nup
            n1=ndn
          endif
          m(iab)=0
          do k=1,ndet
            do ie=1,n1
              if(iworbd(ie+n0,k).eq.io) then
                iesave=ie
                goto 20
              endif
            enddo
            goto 30
 20         continue
            do ie=1,n1
              if(iworbd(ie+n0,k).eq.jo) goto 30
            enddo
            m(iab)=m(iab)+1
            iodet(iab,m(iab))=k
            iopos(iab,m(iab))=iesave
 30       continue
          enddo
        enddo
        if(m(1)+m(2).eq.0) then
          if(iprt.gt.3) write(ounit,'(''no appropriate determinant for '',2i4)') io,jo
          goto 50
        endif

c Define new operator (new variation) and its terms
        noporb=noporb+1
        ! if(noporb.gt.norbterm) then
        !   write(ounit,'(''noporb,max_orb'',2i5)') noporb,norbterm
        !   call fatal_error('ORB_DEFINE: too many terms, increase norbterm')
        ! endif

        ideriv(1,noporb)=io
        ideriv(2,noporb)=jo
        ideriv_iab(noporb)=0
        if(m(1).gt.0) ideriv_iab(noporb)=1
        if(m(2).gt.0) ideriv_iab(noporb)=ideriv_iab(noporb)+2

        do iab=1,2
          n0=0
          n1=nup
          if(iab.eq.2) then
            n0=nup
            n1=ndn
          endif
          ideriv_ref(noporb,iab)=0
          do i=1,n1
            if(iworbd(i+n0,kref).eq.io) then
              ideriv_ref(noporb,iab)=1
              irepcol_ref(noporb,iab)=i
            endif
          enddo
        enddo
        if(iprt.gt.2) write(ounit,'(a16,i4,a8,i4,i5,a15,i4)') 'new variation: ',noporb,' pair ',io,jo,' spin ',ideriv_iab(noporb)

 50    continue
       enddo
 60   continue
      enddo

      norbterm=noporb
      write(ounit,'(''number of orbital variations: '',2i8)') norbterm

c if mix_n, optorb_define called mutiple times with method=sr_n or lin_d
      if(method.eq.'linear') then

        if(mxreduced.ne.norbterm) call fatal_error('READ_INPUT: mxreduced .ne. norbterm')
        nreduced=norbterm
       elseif(method.eq.'sr_n'.or.method.eq.'lin_d'.or.method.eq.'mix_n') then
        nreduced=1
      endif

      icount_orbdef=icount_orbdef+1

      write(ounit,'(''Done with  optorb_define'')')
      return
      end
c-----------------------------------------------------------------------
      subroutine check_orbitals
c Do not compute virtual orbitals during single-electron move
      use orbval,  only: nadorb
      use optorb_mod, only: nadorb_save
      implicit none

      nadorb_save=nadorb
      nadorb=0

      end subroutine

      subroutine check_orbitals_reset
      use orbval,  only: nadorb
      use optorb_mod, only: nadorb_save
      implicit none

      nadorb=nadorb_save

      end subroutine

      end module
