      subroutine determinante(iel,x,rvec_en,r_en,iflag)

      use const, only: pi, hb, etrial, delta, deltai, fbias, nelec, imetro, ipr
      use dets, only: cdet, ndet
      use elec, only: ndn, nup
      implicit real*8(a-h,o-z)



      include 'vmc.h'
      include 'force.h'
      include 'mstates.h'
      include '3dgrid_flags.h'

      parameter(one=1.d0)

      common /contrl_per/ iperiodic,ibasis
      common /phifun/ phin(MBASIS,MELEC),dphin(3,MBASIS,MELEC)
     &,d2phin(MBASIS,MELEC),d2phin_all(3,3,MBASIS,MELEC),d3phin(3,MBASIS,MELEC)
     &,n0_nbasis(MELEC),n0_ibasis(MBASIS,MELEC),n0_ic(MBASIS,MELEC)
      common /coefs/ coef(MBASIS,MORB,MWF),nbasis,norb
c     common /kinet/ dtdx2o(MELEC),dtdx2n(MELEC)
      common /dorb/ iworbd(MELEC,MDET)

      common /multidet/ kref,numrep_det(MDET,2),irepcol_det(MELEC,MDET,2),ireporb_det(MELEC,MDET,2)
     & ,iwundet(MDET,2),iactv(2),ivirt(2)

      common /slater/ slmui(MMAT_DIM),slmdi(MMAT_DIM)
     &,fpu(3,MMAT_DIM),fpd(3,MMAT_DIM)
     &,fppu(MMAT_DIM),fppd(MMAT_DIM)
     &,ddx(3,MELEC),d2dx2(MELEC)

      common /multislater/ detu(MDET),detd(MDET)

      common /slatn/ slmin(MMAT_DIM)
      common /multislatern/ detn(MDET)
     &,orb(MORB),dorb(3,MORB),ddorb(MORB)

      common /wfsec/ iwftype(MFORCE),iwf,nwftype

      dimension x(3,*),rvec_en(3,MELEC,MCENT),r_en(MELEC,MCENT)

      call orbitalse(iel,x,rvec_en,r_en,iflag)

      if(iel.le.nup) then

      ikel=nup*(iel-1)

      ratio_kref=0
      do 30 j=1,nup
   30   ratio_kref=ratio_kref+slmui(j+ikel)*orb(iworbd(j,kref))

      detn(kref)=detu(kref)*ratio_kref

      if(ratio_kref.eq.0.d0) return

      do 45 i=1,nup
        if(i.ne.iel) then
          ik=nup*(i-1)
          sum=0
          do 35 j=1,nup
   35       sum=sum+slmui(j+ik)*orb(iworbd(j,kref))
          sum=sum/ratio_kref
          do 40 j=1,nup
   40       slmin(j+ik)=slmui(j+ik)-slmui(j+ikel)*sum
          endif
   45   continue
        do 50 j=1,nup
   50     slmin(j+ikel)=slmui(j+ikel)/ratio_kref

      else

      ikel=ndn*(iel-nup-1)

      ratio_kref=0
      do 55 j=1,ndn
   55   ratio_kref=ratio_kref+slmdi(j+ikel)*orb(iworbd(j+nup,kref))

      detn(kref)=detd(kref)*ratio_kref
      do 70 i=1,ndn
        if(i+nup.ne.iel) then
          ik=ndn*(i-1)
          sum=0
          do 60 j=1,ndn
   60       sum=sum+slmdi(j+ik)*orb(iworbd(j+nup,kref))
          sum=sum/ratio_kref
          do 65 j=1,ndn
   65      slmin(j+ik)=slmdi(j+ik)-slmdi(j+ikel)*sum
        endif
   70 continue
      do 75 j=1,ndn
   75   slmin(j+ikel)=slmdi(j+ikel)/ratio_kref

      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine compute_determinante_grad(iel,psig,psid,vd,iflag_move)

      use const, only: pi, hb, etrial, delta, deltai, fbias, nelec, imetro, ipr
      use csfs, only: ccsf, cxdet, iadet, ibdet, icxdet, ncsf, nstates

      use dets, only: cdet, ndet
      use elec, only: ndn, nup
      implicit real*8(a-h,o-z)




      include 'vmc.h'
      include 'force.h'
      include 'mstates.h'

      parameter (MEXCIT=10)
      parameter (one=1.d0,half=0.5d0)

      common /dorb/ iworbd(MELEC,MDET)

      common /coefs/ coef(MBASIS,MORB,MWF),nbasis,norb


      common /slater/ slmi(MMAT_DIM,2)
     &,fpu(3,MMAT_DIM),fpd(3,MMAT_DIM)
     &,fppu(MMAT_DIM),fppd(MMAT_DIM)
      common /slatn/ slmin(MMAT_DIM)

      common /multislater/ detu(MDET),detd(MDET)
      common /multislatern/ detn(MDET)
     &,orbn(MORB),dorbn(3,MORB),ddorbn(MORB)

      common /multimat/ aa(MELEC,MORB,2),wfmat(MEXCIT**2,MDET,2)
      common /multimatn/ aan(MELEC,MORB),wfmatn(MEXCIT**2,MDET)

      common /ycompact/ ymat(MORB,MELEC,2,MSTATES),dymat(MORB,MELEC,2,MSTATES)
      common /ycompactn/ ymatn(MORB,MELEC,MSTATES)

      common /orbval/ orb(MELEC,MORB),dorb(3,MELEC,MORB),ddorb(MELEC,MORB),ndetorb,nadorb

      common /multidet/ kref,numrep_det(MDET,2),irepcol_det(MELEC,MDET,2),ireporb_det(MELEC,MDET,2)
     & ,iwundet(MDET,2),iactv(2),ivirt(2)

      common /velocity_jastrow/vj(3,MELEC),vjn(3,MELEC)

      dimension psid(*),vd(3),vref(3),vd_s(3),dorb_tmp(3,MORB)
      dimension ymat_tmp(MORB,MELEC)

      save ymat_tmp

      if(iel.le.nup) then
        iab=1
       else
        iab=2
      endif

      psi2g=psig*psig
      psi2gi=1.d0/psi2g

c All quantities saved (old) avaliable
      if(iflag_move.eq.1) then

        do kk=1,3
          do iorb=1,norb
            dorb_tmp(kk,iorb)=dorb(kk,iel,iorb)
          enddo
        enddo

        call determinante_ref_grad(iel,slmi(1,iab),dorb_tmp,vref)

        if(iguiding.eq.0) then
          detratio=detu(kref)*detd(kref)/psid(1)
          call multideterminante_grad(iel,dorb_tmp,detratio,slmi(1,iab),aa(1,1,iab),wfmat(1,1,iab),ymat(1,1,iab,1),vd)

          do kk=1,3
            vd(kk)=vd(kk)+vref(kk)
          enddo
         else
          do kk=1,3
            vd(kk)=0.d0
          enddo
          do i=1,nstates
            istate=iweight_g(i)

            detratio=detu(kref)*detd(kref)/psid(istate)
            call multideterminante_grad(iel,dorb_tmp,detratio,slmi(1,iab),aa(1,1,iab),wfmat(1,1,iab),ymat(1,1,iab,istate),vd_s)

            do kk=1,3
              vd(kk)=vd(kk)+weights_g(i)*psid(istate)*psid(istate)*(vd_s(kk)+vref(kk))
            enddo
          enddo
          vd(1)=vd(1)*psi2gi
          vd(2)=vd(2)*psi2gi
          vd(3)=vd(3)*psi2gi
        endif

c       write(6,*) 'VJ',(vj(kk,iel),kk=1,3)
c       write(6,*) 'V0',(vref(kk),kk=1,3)
c       write(6,*) 'VD',(vd(kk),kk=1,3)

        vd(1)=vj(1,iel)+vd(1)
        vd(2)=vj(2,iel)+vd(2)
        vd(3)=vj(3,iel)+vd(3)

c Within single-electron move - quantities of electron iel not saved 
       elseif(iflag_move.eq.0) then
       
        call determinante_ref_grad(iel,slmin,dorbn,vref)

        if(iguiding.eq.0) then

          if(iab.eq.1) then
            detratio=detn(kref)*detd(kref)/psid(1)
           else
            detratio=detu(kref)*detn(kref)/psid(1)
          endif
          call multideterminante_grad(iel,dorbn,detratio,slmin,aan,wfmatn,ymatn,vd)

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
              detratio=detn(kref)*detd(kref)/psid(istate)
             else
              detratio=detu(kref)*detn(kref)/psid(istate)
            endif
            call multideterminante_grad(iel,dorbn,detratio,slmin,aan,wfmatn,ymatn(1,1,istate),vd_s)

            do kk=1,3
              vd(kk)=vd(kk)+weights_g(i)*psid(istate)*psid(istate)*(vd_s(kk)+vref(kk))
            enddo
          enddo
          vd(1)=vd(1)*psi2gi
          vd(2)=vd(2)*psi2gi
          vd(3)=vd(3)*psi2gi
        endif

c       write(6,*) 'VJ',(vjn(kk,iel),kk=1,3)
c       write(6,*) 'V0',(vref(kk),kk=1,3)
c       write(6,*) 'VD',(vd(kk),kk=1,3)

        vd(1)=vjn(1,iel)+vd(1)
        vd(2)=vjn(2,iel)+vd(2)
        vd(3)=vjn(3,iel)+vd(3)

       else

c Within single-electron move - iel not equal to electron moved - quantities of electron iel not saved 
        do kk=1,3
          do iorb=1,norb
            dorb_tmp(kk,iorb)=dorb(kk,iel,iorb)
          enddo
        enddo


c iel has same spin as electron moved
        if(iflag_move.eq.2) then

          if(iab.eq.1) then
            detratio=detn(kref)*detd(kref)/psid(1)
           else
            detratio=detu(kref)*detn(kref)/psid(1)
          endif

          call determinante_ref_grad(iel,slmin,dorb_tmp,vref)

          call multideterminante_grad(iel,dorb_tmp,detratio,slmin,aan,wfmatn,ymatn,vd)

c iel has different spin than the electron moved
         else
          if(iab.eq.1) then
            detratio=detu(kref)*detn(kref)/psid(1)
           else
            detratio=detn(kref)*detd(kref)/psid(1)
          endif

          call determinante_ref_grad(iel,slmi(1,iab),dorb_tmp,vref)

          if(iel.eq.1) call compute_ymat(1,detu,detn,wfmat(1,1,1),ymat_tmp,1)

          if(iel.eq.nup+1) call compute_ymat(2,detn,detd,wfmat(1,1,2),ymat_tmp,1)

          call multideterminante_grad(iel,dorb_tmp,detratio,slmi(1,iab),aa(1,1,iab),wfmat(1,1,iab),ymat_tmp(1,1),vd)
        endif

        vd(1)=vjn(1,iel)+vd(1)+vref(1)
        vd(2)=vjn(2,iel)+vd(2)+vref(2)
        vd(3)=vjn(3,iel)+vd(3)+vref(3)
      endif


      return 
      end
c-----------------------------------------------------------------------
      subroutine determinante_ref_grad(iel,slmi,dorb,ddx_ref)

      use const, only: pi, hb, etrial, delta, deltai, fbias, nelec, imetro, ipr
      use dets, only: cdet, ndet
      use elec, only: ndn, nup
      implicit real*8(a-h,o-z)



      include 'vmc.h'
      include 'force.h'
      include 'mstates.h'

      parameter (MEXCIT=10)
      parameter (one=1.d0,half=0.5d0)

      common /dorb/ iworbd(MELEC,MDET)

      common /multidet/ kref,numrep_det(MDET,2),irepcol_det(MELEC,MDET,2),ireporb_det(MELEC,MDET,2)
     & ,iwundet(MDET,2),iactv(2),ivirt(2)

      dimension slmi(MMAT_DIM),dorb(3,MORB)
      dimension ddx_ref(3)

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
      do 84 j=1,nel
        ddx_ref(1)=ddx_ref(1)+slmi(j+ik)*dorb(1,iworbd(j+ish,kref))
        ddx_ref(2)=ddx_ref(2)+slmi(j+ik)*dorb(2,iworbd(j+ish,kref))
   84   ddx_ref(3)=ddx_ref(3)+slmi(j+ik)*dorb(3,iworbd(j+ish,kref))

      return
      end
c-----------------------------------------------------------------------
