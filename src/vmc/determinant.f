      subroutine determinant(ipass,x,rvec_en,r_en)
c     Written by Cyrus Umrigar starting from Kevin Schmidt's routine
c     Modified by A. Scemama

      use vmc_mod, only: MELEC, MORB, MDET, MCENT
      use vmc_mod, only: MMAT_DIM
      use const, only: ipr
      use dets, only: ndet
      use elec, only: ndn, nup
      use multidet, only: kref
      use dorb_m, only: iworbd
      use contr3, only: mode
      use orbval, only: ddorb, dorb, nadorb, ndetorb, orb
      use slater, only: d2dx2, ddx, fp, fpp, slmi
      use multislater, only: detiab
      use csfs, only: nstates

      implicit real*8(a-h,o-z)

      parameter (one=1.d0,half=0.5d0)
      dimension x(3,*),rvec_en(3,MELEC,MCENT),r_en(MELEC,MCENT)

c     compute orbitals
      call orbitals(x,rvec_en,r_en)

      icheck=0
 10   continue
      do iab=1,2
         if(iab.eq.1) then
            ish=0
            nel=nup
         else
            ish=nup
            nel=ndn
         endif

         do istate=1,nstates
            detiab(kref,istate,iab)=1.0d0
            jk=-nel
            do j=1,nel
               jorb=iworbd(j+ish,kref)
               jk=jk+nel
               call dcopy(nel,orb(1+ish,jorb,istate),1,slmi(1+jk,istate,iab),1)
               call dcopy(nel,dorb(1,1+ish,jorb,istate),3,fp(1,j,istate,iab),nel*3)
               call dcopy(nel,dorb(2,1+ish,jorb,istate),3,fp(2,j,istate,iab),nel*3)
               call dcopy(nel,dorb(3,1+ish,jorb,istate),3,fp(3,j,istate,iab),nel*3)
               call dcopy(nel,ddorb(1+ish,jorb,istate),1,fpp(j,istate,iab),nel)
            enddo

c     calculate the inverse transpose matrix and itsdeterminant
            if(nel.gt.0) call matinv(slmi(1,istate,iab),nel,detiab(kref,istate,iab))

c     loop through up spin electrons
c     take inner product of transpose inverse with derivative
c     vectors to get (1/detup)*d(detup)/dx and (1/detup)*d2(detup)/dx**2
            ik=-nel
            do i=1,nel
               ik=ik+nel
               ddx(1,i+ish,istate)=ddot(nel,slmi(1+ik,istate,iab),1,fp(1,1+ik,istate,iab),3)
               ddx(2,i+ish,istate)=ddot(nel,slmi(1+ik,istate,iab),1,fp(2,1+ik,istate,iab),3)
               ddx(3,i+ish,istate)=ddot(nel,slmi(1+ik,istate,iab),1,fp(3,1+ik,istate,iab),3)
               d2dx2(i+ish,istate)=ddot(nel,slmi(1+ik,istate,iab),1,fpp( 1+ik,istate,iab),1)
            enddo

            if(ipr.ge.4) then
               write(6,*) "State ", istate
               ik=-nel
               do i=1,nel
                  ik=ik+nel
                  write(6,*) 'slmi',iab,'M',(slmi(ii+ik,istate,iab),ii=1,nel)
               enddo
            endif
         enddo

      enddo

      if(ipr.ge.4) then
         do istate=1,nstates
            write(6,*) "STATE", istate
            write(6,'(''detu,detd'',9d12.5)') detiab(kref,istate,1),detiab(kref,istate,2)
         enddo
      endif

c     for dmc must be implemented: for each iw, must save not only kref,kref_old but also cdet etc.
c     if(index(mode,'dmc').eq.0) then
c        icheck=icheck+1
c        if(ndet.gt.1.and.kref.lt.ndet.and.icheck.le.10) then
c           call check_detref(ipass,icheck,newref)
c           if(newref.gt.0) goto 10
c        endif
c     endif

      end subroutine

c-----------------------------------------------------------------------

c     subroutine check_detref(ipass,icheck,iflag)

c     use vmc_mod, only: MELEC, MORB, MDET
c     use const, only: ipr
c     use estpsi, only: detref
c     use multidet, only: kref
c     use optwf_contrl, only: ioptorb
c     use coefs, only: norb
c     use orbval, only: ddorb, dorb, nadorb, ndetorb, orb
c     use multislater, only: detiab

c     implicit real*8(a-h,o-z)

c     iflag=0
c     if(ipass.le.2) return

c     do iab=1,2
c        dlogdet=dlog10(dabs(detiab(kref,iab)))
c        dcheck=detref(iab)/ipass-dlogdet
c        if(iab.eq.1.and.dcheck.gt.6) iflag=1
c        if(iab.eq.2.and.dcheck.gt.6) iflag=2
c        if(ipr.ge.2) write(6,*) 'check',dlogdet,detref(iab)/ipass
c     enddo

c     if(ipr.ge.2) write(6,*) 'check detref',iflag
c     if(iflag.gt.0) then
c        call multideterminants_define(iflag,icheck)
c        if (ioptorb.ne.0) then
c           norb=norb+nadorb
c           call optorb_define
c        endif
c     endif
c     
c     end subroutine

c-----------------------------------------------------------------------

      subroutine compute_bmatrices_kin

      use vmc_mod, only: MELEC, MORB
      use atom, only: ncent
      use const, only: hb, nelec
      use da_jastrow4val, only: da_vj
      use da_orbval, only: da_d2orb, da_dorb
      use derivjas, only: g
      use optwf_contrl, only: ioptjas
      use optwf_parms, only: nparmj
      use Bloc, only: b, b_dj, b_da
      use coefs, only: norb
      use force_analy, only: iforce_analy
      use velocity_jastrow, only: vj
      use orbval, only: ddorb, dorb, nadorb, ndetorb, orb
      use csfs, only: nstates

      implicit real*8(a-h,o-z)

      parameter(one=1.d0,half=0.5d0)

      do istate=1,nstates

c     compute kinetic contribution of B+Btilde to compute Eloc
         do i=1,nelec
            do iorb=1,norb+nadorb
               b(iorb,i,istate)=-hb*(ddorb(i,iorb,istate)+2*(vj(1,i)*dorb(1,i,iorb,istate)
     &              +vj(2,i)*dorb(2,i,iorb,istate)+vj(3,i)*dorb(3,i,iorb,istate)))
            enddo
         enddo

c     compute derivative of kinetic contribution of B+Btilde wrt jastrow parameters

         if(ioptjas.gt.0) then
            do iparm=1,nparmj
               do i=1,nelec
                  do iorb=1,norb
                     b_dj(iorb,i,iparm,istate)=-2*hb*(g(1,i,iparm)*dorb(1,i,iorb,istate)
     &                    +g(2,i,iparm)*dorb(2,i,iorb,istate)+g(3,i,iparm)*dorb(3,i,iorb,istate))
                  enddo
               enddo
            enddo
         endif

c     compute derivative of kinetic contribution of B+Btilde wrt nuclear coordinates

         if(iforce_analy.eq.1) then
            do ic=1,ncent
               do iorb=1,norb
                  call dcopy(3*nelec,da_d2orb(1,1,iorb,ic,istate),1,b_da(1,1,iorb,ic,istate),1)
               enddo
            enddo
            do ic=1,ncent
               do i=1,nelec
                  do l=1,3
                     call daxpy(norb,2*vj(1,i),da_dorb(l,1,i,1,ic,istate),
     &                    9*MELEC,b_da(l,i,1,ic,istate),3*MELEC)
                     call daxpy(norb,2*vj(2,i),da_dorb(l,2,i,1,ic,istate),
     &                    9*MELEC,b_da(l,i,1,ic,istate),3*MELEC)
                     call daxpy(norb,2*vj(3,i),da_dorb(l,3,i,1,ic,istate),
     &                    9*MELEC,b_da(l,i,1,ic,istate),3*MELEC)
                     call daxpy(norb,2*da_vj(l,1,i,ic),dorb(1,i,1,istate),
     &                    3*MELEC,b_da(l,i,1,ic,istate),3*MELEC)
                     call daxpy(norb,2*da_vj(l,2,i,ic),dorb(2,i,1,istate),
     &                    3*MELEC,b_da(l,i,1,ic,istate),3*MELEC)
                     call daxpy(norb,2*da_vj(l,3,i,ic),dorb(3,i,1,istate),
     &                    3*MELEC,b_da(l,i,1,ic,istate),3*MELEC)
                     do iorb=1,norb
                        b_da(l,i,iorb,ic,istate)=-hb*b_da(l,i,iorb,ic,istate)
                     enddo
                  enddo
               enddo
            enddo
         endif
      enddo

      end subroutine
