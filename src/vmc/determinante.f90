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
      use vmc_mod, only: nwftypeorb, norb_tot
      use m_backflow, only: ibackflow, rvec_en_bf, r_en_bf, quasi_x
      use backflow_mod, only: backflow
      use m_backflow, only: slmin_bf, orbn_bf, detn_bf, dslm, dorbn_bf
      use orbitals_no_qmckl_mod, only: orbitalse_no_qmckl_bf
      use dorb_m,  only: iworbd
      use matinv_mod, only: matinv
      implicit none

      integer :: i, iab, iel, iflag, ik
      integer :: ikel, ish, j, nel, k, jk, jorb
      real(dp) :: ratio_kref, sum
      real(dp), dimension(3, *) :: x
      real(dp), dimension(3, nelec, ncent_tot) :: rvec_en
      real(dp), dimension(nelec, ncent_tot) :: r_en


      if (ibackflow.gt.0) then
        call backflow(x)
        call orbitalse_no_qmckl_bf(quasi_x, rvec_en_bf, r_en_bf, 0)

        do iab=1,2
          if(iab.eq.1) then
              ish=0
              nel=nup
          else
              ish=nup
              nel=ndn
          endif


         do k=1,nwftypeorb

            jk=-nel
            do j=1,nel
              jorb=iworbd(j+ish,kref)
              jk=jk+nel
              call dcopy(nel,orbn_bf(1+ish,jorb,k),1,slmin_bf(1+jk,iab,k),1)
              call dcopy(nel,dorbn_bf(jorb,1+ish,1,k),norb_tot,dslm(1,1+jk,iab,k),3)
              call dcopy(nel,dorbn_bf(jorb,1+ish,2,k),norb_tot,dslm(2,1+jk,iab,k),3)
              call dcopy(nel,dorbn_bf(jorb,1+ish,3,k),norb_tot,dslm(3,1+jk,iab,k),3)
            enddo

            if(nel.gt.0) call matinv(slmin_bf(1,iab,k),nel,detn_bf(iab,k))
          enddo
        enddo

      else
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

      do k=1,nwftypeorb

        ratio_kref=0.0d0
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

      endif

      return
      end
!-----------------------------------------------------------------------
      subroutine compute_determinante_grad(iel,psi2g,psid,psij,vd,iflag_move)

      use csfs,    only: nstates
      use mstates3, only: iweight_g,weights_g
      use mstates_ctrl, only: iguiding
      use multideterminant_mod, only: compute_ymat
      use multideterminante_mod, only: multideterminante_grad
      use multimat, only: aa,wfmat
      use multimatn, only: aan
      use multislater, only: detiab
      use multislatern, only: detn,dorbn
      use orbval,  only: dorb, orb
      use precision_kinds, only: dp
      use vmc_mod, only: norb_tot, nwftypeorb, stoo, stoj
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
      use contrl_file, only: ounit

      use orbval, only: dorb
      use slater, only: slmi, kref
      use multislater, only: detiab
      use multideterminante_mod, only: multideterminante_grad
      use multideterminant_mod, only: compute_ymat
      use system, only: ndn, nup
      use m_backflow, only: ibackflow, detn_bf, dslm, slmin_bf
      use backflow_mod, only: backflow
      use config,  only: xold

      implicit none

      integer :: i, iab, iel, iflag_move, iorb
      integer :: istate, kk, k, isorb, isjas, isjas1
      real(dp) :: detratio, psi2g, psi2gi, ratio
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

      psi2gi=1.d0/psi2g
      isjas1=stoj(iweight_g(1))

      ! print *, 'iel, iflag_move, ibackflow = ', iel, iflag_move, ibackflow
      ! stop



! All quantities saved (old) avaliable
      if(iflag_move.eq.1) then
        do istate=1,nstates
          isorb=stoo(istate)
          do kk=1,3
            do iorb=1,norb
              dorb_tmp(iorb,kk,isorb)=dorb(iorb,iel,kk,isorb)
            enddo
          enddo
          if (ibackflow.gt.0) then
            call backflow(xold)
            call determinante_ref_grad_bf(iel,slmi(1,1,isorb),dslm(1,1,1,isorb),vref(1,isorb))
          else
            call determinante_ref_grad(iel,slmi(1,iab,isorb),dorb_tmp(1,1,isorb),norb,vref(1,isorb))
          endif
        enddo

        if(iguiding.eq.0) then
          detratio=detiab(kref,1,1)*detiab(kref,2,1)/psid(1)

          call multideterminante_grad(iel,dorb_tmp(1,1,1),norb,detratio,slmi(1,iab,1),aa(1,1,iab,1),ymat(1,1,iab,1),vd)

          do kk=1,3
            vd(kk)=vd(kk)+vref(kk,1)!+vj(kk,iel,1)
          enddo
        else
          do kk=1,3
            vd(kk)=0.d0
          enddo
          do i=1,nstates 

            istate=iweight_g(i)
            isorb=stoo(istate)
            isjas=stoj(istate)

            detratio=detiab(kref,1,isorb)*detiab(kref,2,isorb)/psid(istate)

            call multideterminante_grad(iel,dorb_tmp(1,1,isorb),norb,detratio,slmi(1,iab,isorb), &
               aa(1,1,iab,isorb),ymat(1,1,iab,istate),vd_s(1,istate))

            do kk=1,3 
              ratio=psid(istate)/psid(1)*exp(psij(isjas)-psij(isjas1))
              vd(kk)=vd(kk)+weights_g(i)/anormo(istate)*ratio*ratio &
              *(vd_s(kk,istate)+vref(kk,isorb)+vj(kk,iel,isjas))

            enddo
          enddo
          vd(1)=vd(1)*psi2gi
          vd(2)=vd(2)*psi2gi
          vd(3)=vd(3)*psi2gi
        endif

! Within single-electron move - quantities of electron iel not saved
      elseif(iflag_move.eq.0) then

        do isorb=1,nwftypeorb
          if (ibackflow.gt.0) then
            call determinante_ref_grad_bf(iel,slmin_bf(1,1,isorb),dslm(1,1,1,isorb),vref(1,isorb))
          else
            call determinante_ref_grad(iel,slmin(1,isorb),dorbn(1,1,isorb),norb_tot,vref(1,isorb))
          endif
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

          if (ibackflow.gt.0) then
            detratio = detn_bf(1,1) * detn_bf(2,1) / psid(1)
          endif

         else

          do kk=1,3
            vd(kk)=0.d0
          enddo

          do i=1,nstates
            istate=iweight_g(i)
            isorb=stoo(istate)
            isjas=stoj(istate)

            if(iab.eq.1) then
              detratio=detn(kref,isorb)*detiab(kref,2,isorb)/psid(istate)
             else
              detratio=detiab(kref,1,isorb)*detn(kref,isorb)/psid(istate)
            endif

            call multideterminante_grad(iel,dorbn(1,1,isorb),norb_tot,detratio,slmin(1,isorb), &
                        aan(1,1,isorb),ymatn(1,1,istate),vd_s(1,istate))

            do kk=1,3
              ratio=psid(istate)/psid(1)*exp(psij(isjas)-psij(isjas1))
              vd(kk)=vd(kk)+weights_g(i)/anormo(istate)*ratio*ratio &
              *(vd_s(kk,istate)+vref(kk,isorb)+vjn(kk,iel,isjas))

            enddo
          enddo
          vd(1)=vd(1)*psi2gi
          vd(2)=vd(2)*psi2gi
          vd(3)=vd(3)*psi2gi
        endif

      else

! Within single-electron move - iel not equal to electron moved - quantities of electron iel not saved
        do kk=1,3
          do iorb=1,norb
            dorb_tmp(iorb,kk,1)=dorb(iorb,iel,kk,1)
          enddo
        enddo

! iel has same spin as electron moved
        if(iflag_move.eq.2) then
 
          if(iab.eq.1) then
            detratio=detn(kref,1)*detiab(kref,2,1)/psid(1)
           else
            detratio=detiab(kref,1,1)*detn(kref,1)/psid(1)
          endif


          if (ibackflow.gt.0) then
            call determinante_ref_grad_bf(iel,slmin_bf(1,1,isorb),dslm(1,1,1,isorb),vref(1,isorb))
          else
            call determinante_ref_grad(iel,slmin,dorb_tmp,norb,vref)
          endif

          call multideterminante_grad(iel,dorb_tmp,norb,detratio,slmin,aan,ymatn,vd)

! iel has different spin than the electron moved
         else
          if(iab.eq.1) then
            detratio=detiab(kref,1,1)*detn(kref,1)/psid(1)
           else
            detratio=detn(kref,1)*detiab(kref,2,1)/psid(1)
          endif

          if (ibackflow.gt.0) then
            call determinante_ref_grad_bf(iel,slmin_bf(1,1,isorb),dslm(1,1,1,isorb),vref(1,isorb))
          else
            call determinante_ref_grad(iel,slmi(1,iab,1),dorb_tmp,norb,vref)
          endif

          if(iel.eq.1) call compute_ymat(1,detiab(1,1,1),detn,wfmat(1,1,1,1),ymat_tmp,1)

          if(iel.eq.nup+1) call compute_ymat(2,detn,detiab(1,2,1),wfmat(1,1,2,1),ymat_tmp,1)

          call multideterminante_grad(iel,dorb_tmp,norb, &
                    detratio,slmi(1,iab,1),aa(1,1,iab,1),ymat_tmp,vd)
        endif

        vd(1)=vjn(1,iel,1)+vd(1)+vref(1,1)
        vd(2)=vjn(2,iel,1)+vd(2)+vref(2,1)
        vd(3)=vjn(3,iel,1)+vd(3)+vref(3,1)

        if (ibackflow.gt.0) then
          detratio = detn_bf(1,1) * detn_bf(2,1) / psid(1)
        endif
      endif

      return
      end
!-----------------------------------------------------------------------
      subroutine determinante_ref_grad(iel,slmi,dorb,norbs,ddx_ref)

      use dorb_m,  only: iworbd
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
        ddx_ref(1)=ddx_ref(1)+slmi(j+ik)*dorb(iworbd(j+ish,kref),1)
        ddx_ref(2)=ddx_ref(2)+slmi(j+ik)*dorb(iworbd(j+ish,kref),2)
        ddx_ref(3)=ddx_ref(3)+slmi(j+ik)*dorb(iworbd(j+ish,kref),3)
      enddo

      return
      end

      subroutine determinante_ref_grad_bf(iel,slmi,dslm,ddx_ref)

      use dorb_m,  only: iworbd
      use precision_kinds, only: dp
      use vmc_mod, only: norb_tot
      use vmc_mod, only: nmat_dim
      use system, only: ndn, nup
      use slater, only: kref
      use dorb_m, only: iworbd
      use system, only: nelec
      use m_backflow, only: dquasi_dx
      implicit none

      integer :: iel, ik, ish, j, jel, norbs
      integer :: nel, iab, i, l, kk, jj

      real(dp), dimension(nmat_dim,2) :: slmi
      real(dp), dimension(3, nmat_dim,2) :: dslm
      real(dp), dimension(3) :: ddx_ref


      ddx_ref=0.0d0
        
      do iab=1,2

         if(iab.eq.1) then
            ish=0
            nel=nup
         else
            ish=nup
            nel=ndn
         endif

          do j=1,nel
            do l=1,nel
              do kk=1,3
                do jj=1,3
                  ddx_ref(kk)=ddx_ref(kk)+slmi((j-1)*nel + l,iab)&
                  *dslm(jj,(l-1)*nel + j,iab) &
                  * dquasi_dx(kk,iel,jj,j+ish)

                enddo
              enddo
            enddo
          enddo
        enddo


      return
      end
!-----------------------------------------------------------------------
end module
