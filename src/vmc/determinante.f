      module determinante_mod
      contains
      subroutine determinante(iel,x,rvec_en,r_en,iflag)

      use dorb_m, only: iworbd
      use multislater, only: detiab
      use multislatern, only: detn, orbn
      use orbitals_mod, only: orbitalse
      use precision_kinds, only: dp
      use slater, only: slmi, kref
      use slatn, only: slmin
      use system, only: ndn, nup, nelec, ncent_tot
      use contrl_file, only: ounit
      use vmc_mod, only: nwftypeorb
      implicit none

      integer :: i, iab, iel, iflag, ik
      integer :: ikel, ish, j, nel, k
      real(dp) :: ratio_kref, sum
      real(dp), dimension(3, *) :: x
      real(dp), dimension(3, nelec, ncent_tot) :: rvec_en
      real(dp), dimension(nelec, ncent_tot) :: r_en







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
      do k=1,nwftypeorb !STU check mapping, use nwftypemax? slmi/n, orbn, detiab
        ratio_kref=0
        do j=1,nel
          ratio_kref=ratio_kref+slmi(j+ikel,iab,k)*orbn(iworbd(j+ish,kref),k)
        enddo

        detn(kref,k)=detiab(kref,iab,k)*ratio_kref
        if(ratio_kref.eq.0.d0) return
        do i=1,nel
          if(i+ish.ne.iel) then
            ik=nel*(i-1)
            sum=0
            do j=1,nel
              sum=sum+slmi(j+ik,iab,k)*orbn(iworbd(j+ish,kref),k)
            enddo
            sum=sum/ratio_kref
            do j=1,nel
              slmin(j+ik,k)=slmi(j+ik,iab,k)-slmi(j+ikel,iab,k)*sum
            enddo
          endif
        enddo
        do j=1,nel
          slmin(j+ikel,k)=slmi(j+ikel,iab,k)/ratio_kref
        enddo

      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine compute_determinante_grad(iel,psig,psid,psij,vd,iflag_move)

      use precision_kinds, only: dp
      use vmc_mod, only: norb_tot, nwftypeorb
      use csfs, only: nstates, anormo
      use system, only: nup, nelec
      use slatn, only: slmin
      use ycompact, only: ymat
      use ycompactn, only: ymatn
      use slater, only: norb
      use multimat, only: aa, wfmat
      use multimatn, only: aan
      use velocity_jastrow, only: vj, vjn
      use mstates_ctrl, only: iguiding
      use mstates3, only: iweight_g, weights_g
      use multislatern, only: detn, dorbn

      use orbval, only: dorb
      use slater, only: slmi, kref
      use multislater, only: detiab
      use multideterminante_mod, only: multideterminante_grad
      use multideterminant_mod, only: compute_ymat

      implicit none

      integer :: i, iab, iel, iflag_move, iorb
      integer :: istate, kk, k
      real(dp) :: detratio, psi2g, psi2gi, psig
      real(dp), dimension(*) :: psid
      real(dp), dimension(*) :: psij
      real(dp), dimension(3) :: vd
      real(dp), dimension(3, nstates) :: vref
      real(dp), dimension(3, nstates) :: vd_s
      real(dp), dimension(norb, 3, nstates) :: dorb_tmp

      ! NR : ymat_tmp was not saved ....
      ! it has the save keywoprd in the dev branch ...
      ! real(dp), dimension(norb_tot, nelec) :: ymat_tmp

      real(dp), allocatable, save :: ymat_tmp(:,:,:)
      if (.not. allocated(ymat_tmp)) then
        ! CF : ymat_tmp(norb_tot,nelec) max value of # orb
        allocate(ymat_tmp(norb_tot,nelec,nwftypeorb))
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
        do istate=1,nstates !STU added, use mapping
          do kk=1,3
            do iorb=1,norb
              dorb_tmp(iorb,kk,istate)=dorb(iorb,iel,kk,istate)
            enddo
          enddo
        
          call determinante_ref_grad(iel,slmi(1,iab,istate),dorb_tmp(1,1,istate),norb,vref(1,istate))
        enddo

        if(iguiding.eq.0) then
          detratio=detiab(kref,1,1)*detiab(kref,2,1)/psid(1)

          call multideterminante_grad(iel,dorb_tmp(1,1,1),norb,detratio,slmi(1,iab,1),aa(1,1,iab,1),ymat(1,1,iab,1),vd)

          do kk=1,3
            vd(kk)=vd(kk)+vref(kk,1)+vj(kk,iel,1)
          enddo
        else
          do kk=1,3
            vd(kk)=0.d0
          enddo
          do i=1,nstates !STU was already a loop here, check mapping for detiab, slmi, aa, vref
            istate=iweight_g(i)

            detratio=detiab(kref,1,istate)*detiab(kref,2,istate)/psid(istate)

            call multideterminante_grad(iel,dorb_tmp,norb,detratio,slmi(1,iab,istate),aa(1,1,iab,istate),ymat(1,1,iab,istate),vd_s)

            do kk=1,3 !STU added vj to velocity term
              vd(kk)=vd(kk)+weights_g(i)*psid(istate)*psid(istate)*exp(2*psij(istate))
     &               *(vd_s(kk,istate)+vref(kk,istate)+vj(kk,iel,istate))/1.0d0!anormo(istate)
c missing a vj(kk,iel,istate) term, also no istate dim on vd_s and vref,
c also missing exp(2*psij(istate))
            enddo
          enddo
          vd(1)=vd(1)*psi2gi
          vd(2)=vd(2)*psi2gi
          vd(3)=vd(3)*psi2gi
        endif

c       write(ounit,*) 'VJ',(vj(kk,iel),kk=1,3)
c       write(ounit,*) 'V0',(vref(kk),kk=1,3)
c       write(ounit,*) 'VD',(vd(kk),kk=1,3)
c!STU these may need to be changed, im removing since theyre added earlier
c        vd(1)=vj(1,iel)+vd(1)
c        vd(2)=vj(2,iel)+vd(2)
c        vd(3)=vj(3,iel)+vd(3)

c Within single-electron move - quantities of electron iel not saved
      elseif(iflag_move.eq.0) then

        do istate=1,nstates !STU added istate loop check
          call determinante_ref_grad(iel,slmin(1,istate),dorbn(1,1,istate),norb_tot,vref(1,istate))
        enddo
        if(iguiding.eq.0) then

          if(iab.eq.1) then
            detratio=detn(kref,1)*detiab(kref,2,1)/psid(1)
           else
            detratio=detiab(kref,1,1)*detn(kref,1)/psid(1)
          endif

          call multideterminante_grad(iel,dorbn(1,1,1),norb_tot,detratio,slmin(1,1),aan(1,1,1),ymatn(1,1,1),vd)


          do kk=1,3
            vd(kk)=vd(kk)+vref(kk,1)+vjn(kk,iel,1)
          enddo

         else

          do kk=1,3
            vd(kk)=0.d0
          enddo
          do i=1,nstates !STU was already here, check detn, dorbn, slmin, aan ,vd_s, vjn
            istate=iweight_g(i)

            if(iab.eq.1) then
              detratio=detn(kref,istate)*detiab(kref,2,istate)/psid(istate)
             else
              detratio=detiab(kref,1,istate)*detn(kref,istate)/psid(istate)
            endif

            call multideterminante_grad(iel,dorbn(1,1,istate),norb_tot,
     &                                  detratio,slmin(1,istate),aan(1,1,istate),
     &                                  ymatn(1,1,istate),vd_s(1,istate))

            do kk=1,3
              vd(kk)=vd(kk)+weights_g(i)*psid(istate)*psid(istate)*exp(2*psij(istate))
     &               *(vd_s(kk,istate)+vref(kk,istate)+vjn(kk,iel,istate))/1.0d0!anormo(istate)
c missing a vj(kk,iel,istate) term, also no istate dim on vd_s and vref,
c also missing exp(2*psij(istate))
            enddo
          enddo
          vd(1)=vd(1)*psi2gi
          vd(2)=vd(2)*psi2gi
          vd(3)=vd(3)*psi2gi
        endif

c       write(ounit,*) 'VJ',(vjn(kk,iel),kk=1,3)
c       write(ounit,*) 'V0',(vref(kk),kk=1,3)
c       write(ounit,*) 'VD',(vd(kk),kk=1,3)

!STU this may need to be changed, already added these in. but before
!multiplies by psi2gi
c        vd(1)=vjn(1,iel)+vd(1)
c        vd(2)=vjn(2,iel)+vd(2)
c        vd(3)=vjn(3,iel)+vd(3)

      else

c Within single-electron move - iel not equal to electron moved - quantities of electron iel not saved
        do k=1,nwftypeorb !STU added so no mapping needed
          do kk=1,3
            do iorb=1,norb
              dorb_tmp(iorb,kk,k)=dorb(iorb,iel,kk,k)
            enddo
          enddo
        enddo


c iel has same spin as electron moved
        if(iflag_move.eq.2) then
         
          do istate=1,nstates !STU added, need mapping maybe
            if(iab.eq.1) then
              detratio=detn(kref,istate)*detiab(kref,2,istate)/psid(1)
            else
              detratio=detiab(kref,1,istate)*detn(kref,istate)/psid(1)
            endif

            call determinante_ref_grad(iel,slmin(1,istate),dorb_tmp(1,1,istate),norb,vref(1,istate))

            call multideterminante_grad(iel,dorb_tmp(1,1,istate),norb,detratio,slmin(1,istate),aan(1,1,istate),ymatn(1,1,istate),vd)
          enddo

c iel has different spin than the electron moved
         else
           do istate=1,nstates !STU added check, 
             if(iab.eq.1) then
               detratio=detiab(kref,1,istate)*detn(kref,istate)/psid(1)
             else
               detratio=detn(kref,istate)*detiab(kref,2,istate)/psid(1)
             endif

             call determinante_ref_grad(iel,slmi(1,iab,istate),dorb_tmp(1,1,istate),norb,vref(1,istate))

             if(iel.eq.1) call compute_ymat(1,detiab(1,1,istate),detn(1,istate),
     &                                      wfmat(1,1,1,istate),ymat_tmp(1,1,istate),1)

             if(iel.eq.nup+1) call compute_ymat(2,detn,detiab(1,2,istate),
     &                                          wfmat(1,1,2,istate),
     &                                          ymat_tmp(1,1,istate),1)

             call multideterminante_grad(iel,dorb_tmp(1,1,istate),norb,
     &                                   detratio,slmi(1,iab,istate),aa(1,1,iab,istate),
     &                                   ymat_tmp(1,1,istate),vd)
            enddo
          endif

        vd(1)=vjn(1,iel,1)+vd(1)+vref(1,1)
        vd(2)=vjn(2,iel,1)+vd(2)+vref(2,1)
        vd(3)=vjn(3,iel,1)+vd(3)+vref(3,1)
      endif


      return
      end
c-----------------------------------------------------------------------
      subroutine determinante_ref_grad(iel,slmi,dorb,norbs,ddx_ref)

      use precision_kinds, only: dp
      use vmc_mod, only: norb_tot
      use vmc_mod, only: nmat_dim
      use system, only: ndn, nup
      use slater, only: kref
      use dorb_m, only: iworbd

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
