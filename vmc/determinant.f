      subroutine determinant(ipass,x,rvec_en,r_en)
c Written by Cyrus Umrigar starting from Kevin Schmidt's routine
c Modified by A. Scemama

      implicit real*8(a-h,o-z)
      include 'vmc.h'
      include 'force.h'
      include 'optjas.h'
      include 'mstates.h'
      include '3dgrid_flags.h'

      parameter (one=1.d0,half=0.5d0)

      common /const/ pi,hb,etrial,delta,deltai,fbias,nelec,imetro,ipr
      common /elec/ nup,ndn
      common /dorb/ iworbd(MELEC,MDET)

      common /slater/ slmi(MMAT_DIM,2)
     &,fp(3,MMAT_DIM,2)
     &,fpp(MMAT_DIM,2)
     &,ddx(3,MELEC),d2dx2(MELEC)

      common /multidet/ kref,numrep_det(MDET,2),irepcol_det(MELEC,MDET,2),ireporb_det(MELEC,MDET,2)
     & ,iwundet(MDET,2),iactv(2),ivirt(2)

      common /multislater/detiab(MDET,2)

      common /orbval/ orb(MELEC,MORB),dorb(3,MELEC,MORB),ddorb(MELEC,MORB),ndetorb,nadorb

      dimension x(3,*),rvec_en(3,MELEC,MCENT),r_en(MELEC,MCENT)

c compute orbitals
      call orbitals(x,rvec_en,r_en)

      icheck=0
  10  continue

      do 400 iab=1,2

      if(iab.eq.1) then
        ish=0
        nel=nup
       else
        ish=nup
        nel=ndn
      endif

      detiab(kref,iab)=one

      jk=-nel
      do j=1,nel
        jorb=iworbd(j+ish,kref)
        jk=jk+nel
        call dcopy(nel,orb   (1+ish,jorb),1,slmi(1+jk,iab),1)
        call dcopy(nel,dorb(1,1+ish,jorb),3,fp(1,j,iab),nel*3)
        call dcopy(nel,dorb(2,1+ish,jorb),3,fp(2,j,iab),nel*3)
        call dcopy(nel,dorb(3,1+ish,jorb),3,fp(3,j,iab),nel*3)
        call dcopy(nel,ddorb (1+ish,jorb),1,fpp (j,iab),nel)
      enddo

c calculate the inverse transpose matrix and itsdeterminant
      if(nel.gt.0) call matinv(slmi(1,iab),nel,detiab(kref,iab))

c loop through up spin electrons
c take inner product of transpose inverse with derivative
c vectors to get (1/detup)*d(detup)/dx and (1/detup)*d2(detup)/dx**2
      ik=-nel
      do i=1,nel
        ik=ik+nel
        ddx(1,i+ish)=ddot(nel,slmi(1+ik,iab),1,fp(1,1+ik,iab),3)
        ddx(2,i+ish)=ddot(nel,slmi(1+ik,iab),1,fp(2,1+ik,iab),3)
        ddx(3,i+ish)=ddot(nel,slmi(1+ik,iab),1,fp(3,1+ik,iab),3)
        d2dx2(i+ish)=ddot(nel,slmi(1+ik,iab),1,fpp( 1+ik,iab),1)
      enddo

       if(ipr.ge.4) then
          ik=-nel
          do i=1,nel
            ik=ik+nel
            write(6,*) 'slmi',iab,'M',(slmi(ii+ik,iab),ii=1,nel)
          enddo
        endif
 400  continue

      if(ipr.ge.4) write(6,'(''detu,detd'',9d12.5)') detiab(kref,1),detiab(kref,2)

      icheck=icheck+1
      if(icheck.le.10) then
        call check_detref(ipass,icheck,newref)
        if(newref.gt.0) goto 10
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine check_detref(ipass,icheck,iflag)

      implicit real*8(a-h,o-z)
      include 'vmc.h'
      include 'force.h'
      include 'optjas.h'
      include 'mstates.h'

      common /const/ pi,hb,etrial,delta,deltai,fbias,nelec,imetro,ipr

      common /optwf_contrl/ ioptjas,ioptorb,ioptci,nparm

      common /coefs/ coef(MBASIS,MORB,MWF),nbasis,norb

      common /orbval/ orb(MELEC,MORB),dorb(3,MELEC,MORB),ddorb(MELEC,MORB),ndetorb,nadorb

      common /multidet/ kref,numrep_det(MDET,2),irepcol_det(MELEC,MDET,2),ireporb_det(MELEC,MDET,2)
     & ,iwundet(MDET,2),iactv(2),ivirt(2)

      common /multislater/detiab(MDET,2)

      common /estpsi/ detref(2),apsi(MSTATES),aref

      parameter (one=1.d0,half=0.5d0)

      iflag=0
      if(ipass.le.2) return

      do iab=1,2
        dlogdet=dlog10(dabs(detiab(kref,iab)))
c       dcheck=dabs(dlogdet-detref(iab)/ipass)
c       if(iab.eq.1.and.dcheck.gt.6) iflag=1
c       if(iab.eq.2.and.dcheck.gt.6) iflag=2
        dcheck=detref(iab)/ipass-dlogdet
        if(iab.eq.1.and.dcheck.gt.6) iflag=1
        if(iab.eq.2.and.dcheck.gt.6) iflag=2
        if(ipr.ge.2) write(6,*) 'check',dlogdet,detref(iab)/ipass
      enddo

      if(ipr.ge.2) write(6,*) 'check detref',iflag
      if(iflag.gt.0) then
        call multideterminants_define(iflag,icheck)
        if (ioptorb.ne.0) then
          norb=norb+nadorb
          call optorb_define
        endif
      endif
    
      return
      end
c-----------------------------------------------------------------------
      subroutine compute_bmatrices_kin

      implicit real*8(a-h,o-z)
      include 'vmc.h'
      include 'force.h'
      include 'optjas.h'
      include 'mstates.h'

      parameter (one=1.d0,half=0.5d0)

      common /const/ pi,hb,etrial,delta,deltai,fbias,nelec,imetro,ipr

      common /atom/ znuc(MCTYPE),cent(3,MCENT),pecent
     &,iwctype(MCENT),nctype,ncent

      common /coefs/ coef(MBASIS,MORB,MWF),nbasis,norb

      common /optwf_contrl/ ioptjas,ioptorb,ioptci,nparm
      common /optwf_parms/ nparml,nparme,nparmd,nparms,nparmg,nparmj

      common /orbval/ orb(MELEC,MORB),dorb(3,MELEC,MORB),ddorb(MELEC,MORB),ndetorb,nadorb

      common /velocity_jastrow/vj(3,MELEC),vjn(3,MELEC)
      common /derivjas/ gvalue(MPARMJ),g(3,MELEC,MPARMJ)
     &,d2g(MPARMJ),go(MELEC,MELEC,MPARMJ)

      common /force_analy/ iforce_analy

      common /da_orbval/ da_orb(3,MELEC,MORB,MCENT)
     &,da_d2orb(3,MELEC,MORB,MCENT),da_dorb(3,3,MELEC,MORB,MCENT)
      common /da_jastrow4val/ da_j(3,MELEC,MCENT),da_d2j(3,MELEC,MCENT),da_vj(3,3,MELEC,MCENT)

      common /Bloc/ b(MORB,MELEC),xmatu(MELEC**2),xmatd(MELEC**2)
     & ,tildem(MELEC,MORB,2)
      common /Bloc_dj/ b_dj(MORB,MELEC,MPARMJ)
      common /Bloc_da/ b_da(3,MELEC,MORB,MCENT)

c compute kinetic contribution of B+Btilde to compute Eloc
      do i=1,nelec
        do iorb=1,norb+nadorb
          b(iorb,i)=-hb*(ddorb(i,iorb)+2*(vj(1,i)*dorb(1,i,iorb)+vj(2,i)*dorb(2,i,iorb)+vj(3,i)*dorb(3,i,iorb)))
        enddo
      enddo

c compute derivative of kinetic contribution of B+Btilde wrt jastrow parameters
      if(ioptjas.gt.0) then
        do iparm=1,nparmj
          do i=1,nelec
            do iorb=1,norb
              b_dj(iorb,i,iparm)=-2*hb*(g(1,i,iparm)*dorb(1,i,iorb)+g(2,i,iparm)*dorb(2,i,iorb)+g(3,i,iparm)*dorb(3,i,iorb))
            enddo
          enddo
        enddo
      endif

c compute derivative of kinetic contribution of B+Btilde wrt nuclear coordinates
      if(iforce_analy.eq.1) then
        do ic=1,ncent
          do iorb=1,norb
            call dcopy(3*nelec,da_d2orb(1,1,iorb,ic),1,b_da(1,1,iorb,ic),1)
          enddo
        enddo
        do ic=1,ncent
          do i=1,nelec
            do l=1,3
              call daxpy(norb,2*vj(1,i),da_dorb(l,1,i,1,ic),9*melec,b_da(l,i,1,ic),3*melec)
              call daxpy(norb,2*vj(2,i),da_dorb(l,2,i,1,ic),9*melec,b_da(l,i,1,ic),3*melec)
              call daxpy(norb,2*vj(3,i),da_dorb(l,3,i,1,ic),9*melec,b_da(l,i,1,ic),3*melec)
              call daxpy(norb,2*da_vj(l,1,i,ic),dorb(1,i,1),3*melec,b_da(l,i,1,ic),3*melec)
              call daxpy(norb,2*da_vj(l,2,i,ic),dorb(2,i,1),3*melec,b_da(l,i,1,ic),3*melec)
              call daxpy(norb,2*da_vj(l,3,i,ic),dorb(3,i,1),3*melec,b_da(l,i,1,ic),3*melec)
              do iorb=1,norb
                b_da(l,i,iorb,ic)=-hb*b_da(l,i,iorb,ic)
              enddo
            enddo
          enddo
        enddo
      endif
c     do 10 ic=1,ncent
c       do 10 iorb=1,norb
c         do 10 i=1,nelec
c           do 10 l=1,3
c 10          db(l,i,iorb,ic)=da_d2orb(l,i,iorb,ic)+two*(
c    &           vj(1,i)*da_dorb(l,1,i,iorb,ic)
c    &          +vj(2,i)*da_dorb(l,2,i,iorb,ic)
c    &          +vj(3,i)*da_dorb(l,3,i,iorb,ic)
c    &          +da_vj(l,1,i,ic)*dorb(1,i,iorb)
c    &          +da_vj(l,2,i,ic)*dorb(2,i,iorb)
c    &          +da_vj(l,3,i,ic)*dorb(3,i,iorb))

      return
      end
