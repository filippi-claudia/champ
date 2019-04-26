      subroutine multideterminante(iel)

      implicit real*8(a-h,o-z)
      include 'vmc.h'
      include 'force.h'
      include 'mstates.h'

      parameter (one=1.d0,half=0.5d0)
      parameter (MEXCIT=10)

      common /const/ pi,hb,etrial,delta,deltai,fbias,nelec,imetro,ipr
      common /elec/ nup,ndn
      common /dets/ cdet(MDET,MSTATES,MWF),ndet
      common /csfs/ ccsf(MDET,MSTATES,MWF),cxdet(MDET*MDETCSFX)
     &,icxdet(MDET*MDETCSFX),iadet(MDET),ibdet(MDET),ncsf,nstates

      common /dorb/ iworbd(MELEC,MDET)

      common /multislater/ detu(MDET),detd(MDET)

      common /slatn/ slmin(MMAT_DIM)
      common /multislatern/ detn(MDET)
     &,orbn(MORB),dorbn(3,MORB),ddorbn(MORB)

      common /orbval/ orb(MELEC,MORB),dorb(3,MELEC,MORB),ddorb(MELEC,MORB),ndetorb,nadorb
      common /coefs/ coef(MBASIS,MORB,MWF),nbasis,norb

      common /multidet/ kref,numrep_det(MDET,2),irepcol_det(MELEC,MDET,2),ireporb_det(MELEC,MDET,2)
     & ,iwundet(MDET,2),iactv(2),ivirt(2)

      common /multimat/ aa(MELEC,MORB,2),wfmat(MEXCIT**2,MDET,2)

      common /multimatn/ aan(MELEC,MORB),wfmatn(MEXCIT**2,MDET)

      common /ycompactn/ ymatn(MORB,MELEC,MSTATES)

      dimension gmat(MELEC,MORB,3),gmatn(MEXCIT**2,3)
      dimension b(MORB,3),ddx_mdet(3)

      dimension orb_sav(MORB)

      if(ndet.eq.1) return

      iab=1
      nel=nup
      ish=0
      if(iel.gt.nup) then
        iab=2
        nel=ndn
        ish=nup
      endif

c temporarely copy orbn to orb 
      do iorb=1,norb
        orb_sav(iorb)=orb(iel,iorb)
        orb(iel,iorb)=orbn(iorb)
      enddo

      do jrep=ivirt(iab),norb
        do irep=1,nel

          dum1=0.d0
          do i=1,nel
           dum1=dum1+slmin(irep+(i-1)*nup)*orb(i+ish,jrep)
          enddo
          aan(irep,jrep)=dum1

        enddo
      enddo

c compute wave function 
      do 200 k=1,ndet

        if(k.ne.kref) then

        if(iwundet(k,iab).eq.k) then

          ndim=numrep_det(k,iab)

          jj=0
          do jrep=1,ndim
            jorb=ireporb_det(jrep,k,iab)
            do irep=1,ndim
              iorb=irepcol_det(irep,k,iab)
              jj=jj+1

              wfmatn(jj,k)=aan(iorb,jorb)
            enddo
          enddo

          call matinv(wfmatn(1,k),ndim,det)

          detn(k)=det

         else
          index_det=iwundet(k,iab)
          detn(k)=detn(index_det)

        endif

        endif

 200  continue

      do 400 k=1,ndet
        if(k.ne.kref.and.iwundet(k,iab).ne.kref) then
          detn(k)=detn(k)*detn(kref)
        endif
 400  continue

      do 800 istate=1,nstates
        if(iab.eq.1) call compute_ymat(iab,detn,detd,wfmatn,ymatn(1,1,istate),istate)
        if(iab.eq.2) call compute_ymat(iab,detu,detn,wfmatn,ymatn(1,1,istate),istate)
 800  continue

      do iorb=1,norb
        orb(iel,iorb)=orb_sav(iorb)
      enddo

      return
      end

c-----------------------------------------------------------------------
      subroutine multideterminante_grad(iel,dorb,detratio,slmi,aa,wfmat,ymat,velocity)

      implicit real*8(a-h,o-z)
      include 'vmc.h'
      include 'force.h'
      include 'mstates.h'

      parameter (one=1.d0,half=0.5d0)
      parameter (MEXCIT=10)

      common /const/ pi,hb,etrial,delta,deltai,fbias,nelec,imetro,ipr
      common /elec/ nup,ndn
      common /dets/ cdet(MDET,MSTATES,MWF),ndet
      common /dorb/ iworbd(MELEC,MDET)


      common /coefs/ coef(MBASIS,MORB,MWF),nbasis,norb

      common /multidet/ kref,numrep_det(MDET,2),irepcol_det(MELEC,MDET,2),ireporb_det(MELEC,MDET,2)
     & ,iwundet(MDET,2),iactv(2),ivirt(2)


      common /wfsec/ iwftype(MFORCE),iwf,nwftype

      dimension slmi(MMAT_DIM)
      dimension aa(MELEC,MORB),wfmat(MEXCIT**2,MDET),ymat(MORB,MELEC)

      dimension b(MORB,3),dorb(3,MORB)
      dimension gmat(MELEC,MORB,3)
      dimension velocity(3)

      do k=1,3
        velocity(k)=0.d0
      enddo
      if(ndet.eq.1) return

      if(iel.le.nup) then
        iab=1
        nel=nup
        ish=0
       else
        iab=2
        nel=ndn
        ish=nup
      endif

      jel=iel-ish

c TMP to fix
      do kk=1,3
        do iorb=1,norb
          b(iorb,kk)=dorb(kk,iorb)
        enddo
      enddo

      do 50 kk=1,3

        do jrep=ivirt(iab),norb
          dum=0
          do j=1,nel
            dum=dum+b(iworbd(j+ish,kref),kk)*aa(j,jrep)
          enddo
          dum=b(jrep,kk)-dum

          do irep=iactv(iab),nel
            gmat(irep,jrep,kk)=dum*slmi(irep+(jel-1)*nel)
          enddo
        enddo

 50   continue

c     if(iab.eq.2) write(6,*) 'gmat ',(((gmat(irep,jrep,kk),irep=iactv(iab),nel),jrep=ivirt(iab),norb),kk=1,3)

      do kk=1,3
        dum=0
        do jrep=ivirt(iab),norb
          do irep=iactv(iab),nel
            dum=dum+ymat(jrep,irep)*gmat(irep,jrep,kk)
          enddo
        enddo
        velocity(kk)=dum*detratio
      enddo

c     if(iab.eq.2) write(6,*) 'ymat ',((ymat(jrep,irep),irep=iactv(iab),nel),jrep=ivirt(iab),norb)

      return
      end
c-----------------------------------------------------------------------
