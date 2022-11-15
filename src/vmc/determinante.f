      module determinante_mod
      contains
      subroutine determinante(iel,x,rvec_en,r_en,iflag)

      use dorb_m,  only: iworbd
      use multislater, only: detiab
      use multislatern, only: detn,orbn
      use orbitals_mod, only: orbitalse
      use precision_kinds, only: dp
      use slater,  only: kref,slmi
      use slatn,   only: slmin
      use system,  only: ncent_tot,ndn,nelec,nup
      implicit none

      integer :: i, iab, iel, iflag, ik
      integer :: ikel, ish, j, nel
      real(dp) :: ratio_kref, sum
      real(dp), dimension(3, *) :: x
      real(dp), dimension(3, nelec, ncent_tot) :: rvec_en
      real(dp), dimension(nelec, ncent_tot) :: r_en
      real(dp), allocatable :: orbn_ordered(:)
      real(dp) :: ratio_kref_inv


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

      allocate(orbn_ordered(nel))
      do j=1,nel
        orbn_ordered(j) = orbn(iworbd(j+ish,kref))
      enddo

      ratio_kref=0.d0
      do j=1,nel
        ratio_kref=ratio_kref+slmi(j+ikel,iab)*orbn_ordered(j)
      enddo

      detn(kref)=detiab(kref,iab)*ratio_kref

      if(ratio_kref.eq.0.d0) return

      ratio_kref_inv = 1.d0/ratio_kref


      ik=nel*(iel-ish-1)
      sum=0.d0
      do j=1,nel
        sum=sum+slmi(j+ik,iab)*orbn_ordered(j)
      enddo
      sum=-sum*ratio_kref_inv
      do j=1,nel
       slmin(j+ik)=slmi(j+ik,iab)-slmi(j+ikel,iab)*sum
      enddo
      do i=1,nel
!        if(i+ish.ne.iel) then
          ik=nel*(i-1)
          sum=0.d0
          do j=1,nel
            sum=sum+slmi(j+ik,iab)*orbn_ordered(j)
          enddo
          sum=sum*ratio_kref_inv
          do j=1,nel
           slmin(j+ik)=slmi(j+ik,iab)-slmi(j+ikel,iab)*sum
          enddo
!       endif
      enddo
      do j=1,nel
        slmin(j+ikel)=slmi(j+ikel,iab)*ratio_kref_inv
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine compute_determinante_grad(iel,psig,psid,vd,iflag_move)

      use csfs,    only: nstates
      use mstates3, only: iweight_g,weights_g
      use mstates_ctrl, only: iguiding
      use multideterminant_mod, only: compute_ymat
      use multideterminante_mod, only: multideterminante_grad
      use multimat, only: aa,wfmat
      use multimatn, only: aan
      use multislater, only: detiab
      use multislatern, only: detn,dorbn
      use orbval,  only: dorb
      use precision_kinds, only: dp
      use slater,  only: kref,norb,slmi
      use slatn,   only: slmin
      use system,  only: nelec,nup
      use velocity_jastrow, only: vj,vjn
      use vmc_mod, only: norb_tot
      use ycompact, only: ymat
      use ycompactn, only: ymatn


      implicit none

      integer :: i, iab, iel, iflag_move, iorb
      integer :: istate, kk
      real(dp) :: detratio, psi2g, psi2gi, psig
      real(dp), dimension(*) :: psid
      real(dp), dimension(3) :: vd
      real(dp), dimension(3) :: vref
      real(dp), dimension(3) :: vd_s
      real(dp), dimension(norb,3) :: dorb_tmp

      ! NR : ymat_tmp was not saved ....
      ! it has the save keywoprd in the dev branch ...
      ! real(dp), dimension(norb_tot, nelec) :: ymat_tmp

      real(dp), allocatable, save :: ymat_tmp(:,:)
      if (.not. allocated(ymat_tmp)) then
        ! CF : ymat_tmp(norb_tot,nelec) max value of # orb
        allocate(ymat_tmp(norb_tot,nelec))
      endif

      ! save ymat_tmp

      if(iel.le.nup) then
        iab=1
       else
        iab=2
      endif

      psi2g=psig*psig
      psi2gi=1.d0/psi2g
c      print*,"norb",norb,"norb_tot",norb_tot
      
c All quantities saved (old) avaliable
      if(iflag_move.eq.1) then

        do kk=1,3
          do iorb=1,norb
            dorb_tmp(iorb,kk)=dorb(iorb,iel,kk)
          enddo
        enddo
        
        call determinante_ref_grad(iel,slmi(1,iab),dorb_tmp,norb,vref)

        if(iguiding.eq.0) then
          detratio=detiab(kref,1)*detiab(kref,2)/psid(1)

          call multideterminante_grad(iel,dorb_tmp,norb,detratio,slmi(1,iab),aa(1,1,iab),ymat(1,1,iab,1),vd)

          do kk=1,3
            vd(kk)=vd(kk)+vref(kk)
          enddo
         else
          do kk=1,3
            vd(kk)=0.d0
          enddo
          do i=1,nstates
            istate=iweight_g(i)

            detratio=detiab(kref,1)*detiab(kref,2)/psid(istate)

            call multideterminante_grad(iel,dorb_tmp,norb,detratio,slmi(1,iab),aa(1,1,iab),ymat(1,1,iab,istate),vd_s)

            do kk=1,3
              vd(kk)=vd(kk)+weights_g(i)*psid(istate)*psid(istate)*(vd_s(kk)+vref(kk))
            enddo
          enddo
          vd(1)=vd(1)*psi2gi
          vd(2)=vd(2)*psi2gi
          vd(3)=vd(3)*psi2gi
        endif

c       write(ounit,*) 'VJ',(vj(kk,iel),kk=1,3)
c       write(ounit,*) 'V0',(vref(kk),kk=1,3)
c       write(ounit,*) 'VD',(vd(kk),kk=1,3)

        vd(1)=vj(1,iel)+vd(1)
        vd(2)=vj(2,iel)+vd(2)
        vd(3)=vj(3,iel)+vd(3)

c Within single-electron move - quantities of electron iel not saved
       elseif(iflag_move.eq.0) then

        call determinante_ref_grad(iel,slmin,dorbn,norb_tot,vref)

        if(iguiding.eq.0) then

          if(iab.eq.1) then
            detratio=detn(kref)*detiab(kref,2)/psid(1)
           else
            detratio=detiab(kref,1)*detn(kref)/psid(1)
          endif

          call multideterminante_grad(iel,dorbn,norb_tot,detratio,slmin,aan,ymatn,vd)


          do kk=1,3
            vd(kk)=vd(kk)+vref(kk)
          enddo

         else

          do kk=1,3
            vd(kk)=0.d0
          enddo
          do i=1,nstates
            istate=iweight_g(i)

            if(iab.eq.1) then
              detratio=detn(kref)*detiab(kref,2)/psid(istate)
             else
              detratio=detiab(kref,1)*detn(kref)/psid(istate)
            endif

            call multideterminante_grad(iel,dorbn,norb_tot,detratio,slmin,aan,ymatn(1,1,istate),vd_s)

            do kk=1,3
              vd(kk)=vd(kk)+weights_g(i)*psid(istate)*psid(istate)*(vd_s(kk)+vref(kk))
            enddo
          enddo
          vd(1)=vd(1)*psi2gi
          vd(2)=vd(2)*psi2gi
          vd(3)=vd(3)*psi2gi
        endif

c       write(ounit,*) 'VJ',(vjn(kk,iel),kk=1,3)
c       write(ounit,*) 'V0',(vref(kk),kk=1,3)
c       write(ounit,*) 'VD',(vd(kk),kk=1,3)

        vd(1)=vjn(1,iel)+vd(1)
        vd(2)=vjn(2,iel)+vd(2)
        vd(3)=vjn(3,iel)+vd(3)

       else

c Within single-electron move - iel not equal to electron moved - quantities of electron iel not saved
        do kk=1,3
          do iorb=1,norb
            dorb_tmp(iorb,kk)=dorb(iorb,iel,kk)
          enddo
        enddo


c iel has same spin as electron moved
        if(iflag_move.eq.2) then

          if(iab.eq.1) then
            detratio=detn(kref)*detiab(kref,2)/psid(1)
           else
            detratio=detiab(kref,1)*detn(kref)/psid(1)
          endif

          call determinante_ref_grad(iel,slmin,dorb_tmp,norb,vref)

          call multideterminante_grad(iel,dorb_tmp,norb,detratio,slmin,aan,ymatn,vd)


c iel has different spin than the electron moved
         else
          if(iab.eq.1) then
            detratio=detiab(kref,1)*detn(kref)/psid(1)
           else
            detratio=detn(kref)*detiab(kref,2)/psid(1)
          endif

          call determinante_ref_grad(iel,slmi(1,iab),dorb_tmp,norb,vref)

          if(iel.eq.1) call compute_ymat(1,detiab(1,1),detn,wfmat(:,:,1),ymat_tmp,1)

          if(iel.eq.nup+1) call compute_ymat(2,detn,detiab(1,2),wfmat(:,:,2),ymat_tmp,1)

          call multideterminante_grad(iel,dorb_tmp,norb,detratio,slmi(1,iab),aa(1,1,iab),ymat_tmp(1,1),vd)

        endif

        vd(1)=vjn(1,iel)+vd(1)+vref(1)
        vd(2)=vjn(2,iel)+vd(2)+vref(2)
        vd(3)=vjn(3,iel)+vd(3)+vref(3)
      endif


      return
      end
c-----------------------------------------------------------------------
      subroutine determinante_ref_grad(iel,slmi,dorb,norbs,ddx_ref)

      use dorb_m,  only: iworbd
      use precision_kinds, only: dp
      use slater,  only: kref
      use system,  only: ndn,nup
      use vmc_mod, only: nmat_dim,norb_tot

      implicit none

      integer :: iel, ik, ish, j, jel, norbs
      integer :: nel

      real(dp), dimension(nmat_dim) :: slmi
      real(dp), dimension(norbs,3) :: dorb
      real(dp), dimension(3) :: ddx_ref


      
      ddx_ref(1)=0
      ddx_ref(2)=0
      ddx_ref(3)=0

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
        ddx_ref(1)=ddx_ref(1)+slmi(j+ik)*dorb(iworbd(j+ish,kref),1)
        ddx_ref(2)=ddx_ref(2)+slmi(j+ik)*dorb(iworbd(j+ish,kref),2)
        ddx_ref(3)=ddx_ref(3)+slmi(j+ik)*dorb(iworbd(j+ish,kref),3)
      enddo

      return
      end
c-----------------------------------------------------------------------
      end module
