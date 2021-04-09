      subroutine orbitals(x,rvec_en,r_en)
c     Written by Cyrus Umrigar starting from Kevin Schmidt's routine
c     Modified by A. Scemama

      use vmc_mod, only: MELEC, MORB, MBASIS, MCENT
      use const, only: nelec, ipr
      use phifun, only: d2phin, dphin, n0_ibasis, n0_nbasis
      use phifun, only: phin
      use wfsec, only: iwf
      use coefs, only: coef, nbasis, norb
      use contrl_per, only: iperiodic
      use force_analy, only: iforce_analy
      use grid3dflag, only: i3dlagorb, i3dsplorb
      use csfs, only: nstates
      use orbval, only: ddorb, dorb, nadorb, ndetorb, orb

      implicit real*8(a-h,o-z)

      dimension x(3,*),rvec_en(3,MELEC,MCENT),r_en(MELEC,MCENT)
      dimension bhin(MELEC,MBASIS),dbhin(3*MELEC,MBASIS),d2bhin(MELEC,MBASIS)

      ier=1
      if(iperiodic.eq.0) then
c     Spline interpolation
         if(i3dsplorb.eq.2) then
            do istate=1,nstates
               do i=1,nelec
                  ier=0
                  do iorb=1,norb+nadorb
                     ddorb(i,iorb,istate)=1.0d0    ! compute the laplacian
                     dorb(1:3,i,iorb,istate)=1.0d0 ! compute the gradients
                     call spline_mo(x(1,i),iorb,orb(i,iorb,istate), 
     &                             dorb(1,i,iorb,istate),ddorb(i,iorb,istate),ier)
                  enddo
                  if(ier.eq.1) then
                     call basis_fnse_vgl(i,rvec_en,r_en)
                     do iorb=1,norb+nadorb
                        orb(i,iorb,istate)=0.0d0
                        dorb(1:3,i,iorb,istate)=0.0d0
                        ddorb(i,iorb,istate)=0.d0
                        do m0=1,n0_nbasis(i)
                           m=n0_ibasis(m0,i)
                           orb(i,iorb,istate)=orb(i,iorb,istate)+coef(m,iorb,istate,iwf)*phin(m,i)
                           dorb(1,i,iorb,istate)=dorb(1,i,iorb,istate)+coef(m,iorb,istate,iwf)*dphin(1,m,i)
                           dorb(2,i,iorb,istate)=dorb(2,i,iorb,istate)+coef(m,iorb,istate,iwf)*dphin(2,m,i)
                           dorb(3,i,iorb,istate)=dorb(3,i,iorb,istate)+coef(m,iorb,istate,iwf)*dphin(3,m,i)
                           ddorb(i,iorb,istate)=ddorb(i,iorb,istate)+coef(m,iorb,istate,iwf)*d2phin(m,i)
                        enddo
                     enddo
                  endif
               enddo
            enddo
c     Lagrange interpolation
         elseif(i3dlagorb.eq.2) then
            do istate=1,nstates
               do i=1,nelec
                  ier=0
                  call lagrange_mos(1,x(1,i),orb(1,1,istate),i,ier)
                  call lagrange_mos_grad(2,x(1,i),dorb(1,1,1,istate),i,ier)
                  call lagrange_mos_grad(3,x(1,i),dorb(1,1,1,istate),i,ier)
                  call lagrange_mos_grad(4,x(1,i),dorb(1,1,1,istate),i,ier)
                  call lagrange_mos(5,x(1,i),ddorb(1,1,istate),i,ier)
                  if(ier.eq.1) then
                     call basis_fnse_vgl(i,rvec_en,r_en)
                     do iorb=1,norb+nadorb
                        orb(i,iorb,istate)=0.0d0
                        dorb(1:3,i,iorb,istate)=0.0d0
                        ddorb(i,iorb,istate)=0.0d0
                        do m0=1,n0_nbasis(i)
                           m=n0_ibasis(m0,i)
                           orb(i,iorb,istate)=orb(i,iorb,istate)+coef(m,iorb,istate,iwf)*phin(m,i)
                           dorb(1,i,iorb,istate)=dorb(1,i,iorb,istate)+coef(m,iorb,istate,iwf)*dphin(1,m,i)
                           dorb(2,i,iorb,istate)=dorb(2,i,iorb,istate)+coef(m,iorb,istate,iwf)*dphin(2,m,i)
                           dorb(3,i,iorb,istate)=dorb(3,i,iorb,istate)+coef(m,iorb,istate,iwf)*dphin(3,m,i)
                           ddorb(i,iorb,istate)=ddorb(i,iorb,istate)+coef(m,iorb,istate,iwf)*d2phin(m,i)
                        enddo
                     enddo    
                  endif
               enddo
            enddo
         else 
c     no 3d interpolation
c     get basis functions for all electrons
            call basis_fns_vgl(x,rvec_en,r_en)
            do istate=1,nstates
               do i=1,nelec
                  do iorb=1,norb+nadorb
                     orb(i,iorb,istate)=0.0d0
                     dorb(1:3,i,iorb,istate)=0.0d0
                     ddorb(i,iorb,istate)=0.0d0
                     do m0=1,n0_nbasis(i)
                        m=n0_ibasis(m0,i)
                        orb(i,iorb,istate)=orb(i,iorb,istate)+coef(m,iorb,istate,iwf)*phin(m,i)
                        dorb(1,i,iorb,istate)=dorb(1,i,iorb,istate)+coef(m,iorb,istate,iwf)*dphin(1,m,i)
                        dorb(2,i,iorb,istate)=dorb(2,i,iorb,istate)+coef(m,iorb,istate,iwf)*dphin(2,m,i)
                        dorb(3,i,iorb,istate)=dorb(3,i,iorb,istate)+coef(m,iorb,istate,iwf)*dphin(3,m,i)
                        ddorb(i,iorb,istate)=ddorb(i,iorb,istate)+coef(m,iorb,istate,iwf)*d2phin(m,i)
                     enddo
                  enddo     
               enddo
            enddo
         endif

         if(iforce_analy.eq.1) call da_orbitals

      else
         call orbitals_pw(x,orb,dorb,ddorb)
      endif

      if(ipr.ge.4) then
         do istate=1,nstates
            print *, "ORBITALS STATE", istate
            do iorb=1,norb+nadorb
               write(6,'(''iorb,orb='',i4,1000f15.11)')iorb,(orb(i,iorb,istate),i=1,nelec)
            enddo
            do iorb=1,norb+nadorb
               write(6,'(''iorb,d2orb='',i4,1000f15.11)')iorb,(ddorb(i,iorb,istate),i=1,nelec)
            enddo
            do k=1,3
               do iorb=1,norb+nadorb
                  write(6,'(''iorb,dorb='',2i4,1000f12.8)')k,iorb,(dorb(k,i,iorb,istate),i=1,nelec)
               enddo
            enddo
         enddo
      endif

      end subroutine

c------------------------------------------------------------------------------------

      subroutine virtual_orbitals

      use vmc_mod, only: MELEC, MORB, MBASIS
      use const, only: nelec
      use optwf_contrl, only: ioptci, ioptorb
      use phifun, only: d2phin, dphin
      use phifun, only: phin
      use coefs, only: coef, nbasis, norb
      use orbval, only: ddorb, dorb, nadorb, ndetorb, orb
      use csfs, only: nstates

      implicit real*8(a-h,o-z)

c     compute values of extra ('virtual') orbitals needed for optorb operators
c     assuming that basis function values in phin are up to date

      dimension bhin(MELEC,MBASIS),dbhin(3,MELEC,MBASIS),d2bhin(MELEC,MBASIS)

      if (nadorb.eq.0.or.(ioptorb.eq.0.and.ioptci.eq.0)) return

c     primary geometry only:

      iwf=1
      if(norb+nadorb.gt.MORB) then
         write(6,'(''VIRTUAL_ORB: Too many orbitals, norb + nadorb='',
     &i3,'' > MORB='',i4)') norb+nadorb,MORB
         call fatal_error('Aborted')
      endif

      do i=1,nelec
         call dcopy(nbasis,phin(1,i),1,bhin(i,1),MELEC)
         call dcopy(nbasis,dphin(1,1,i),3,dbhin(1,i,1),3*MELEC)
         call dcopy(nbasis,dphin(2,1,i),3,dbhin(2,i,1),3*MELEC)
         call dcopy(nbasis,dphin(3,1,i),3,dbhin(3,i,1),3*MELEC)
         call dcopy(nbasis,d2phin(1,i),1,d2bhin(i,1),MELEC)
      enddo

      do istate=1,nstates
         call dgemm('n','n',nelec,nadorb,nbasis,1.d0,bhin,
     &        MELEC,coef(1,norb+1,istate,iwf),MBASIS,0.d0,orb(1,norb+1,istate),MELEC)
         call dgemm('n','n',3*nelec,nadorb,nbasis,1.d0,dbhin,
     &        3*MELEC,coef(1,norb+1,istate,iwf),MBASIS,0.d0,dorb(1,1,norb+1,istate),3*MELEC)
         call dgemm('n','n',nelec,nadorb,nbasis,1.d0,d2bhin,
     &        MELEC,coef(1,norb+1,istate,iwf),MBASIS,0.d0,ddorb(1,norb+1,istate),MELEC)
      enddo

      end subroutine

c------------------------------------------------------------------------------------

      subroutine da_orbitals

      use vmc_mod, only: MELEC, MORB, MBASIS
      use atom, only: ncent
      use const, only: nelec
      use da_orbval, only: da_d2orb, da_dorb, da_orb
      use numbas2, only: ibas0, ibas1
      use phifun, only: d2phin_all, d3phin, dphin
      use wfsec, only: iwf
      use coefs, only: coef, nbasis, norb
      use contrl_per, only: ibasis
      use orbval, only: ddorb, dorb, nadorb, ndetorb, orb
      use csfs, only: nstates

      implicit real*8(a-h,o-z)

      dimension tphin(3*MELEC,MBASIS),t2phin_all(3*3*MELEC,MBASIS),t3phin(3*MELEC,MBASIS)

      do ibasis=1,nbasis
         i=0
         j=0
         do ielec=1,nelec
            do l=1,3
               i=i+1
               tphin(i,ibasis)=dphin(l,ibasis,ielec)
               t3phin(i,ibasis)=d3phin(l,ibasis,ielec)
               do k=1,3
                  j=j+1
                  t2phin_all(j,ibasis)=d2phin_all(k,l,ibasis,ielec)
               enddo
            enddo
         enddo
      enddo

      n=3*nelec
      m=3*MELEC

      do istate=1,nstates
         do ic=1,ncent
            k=ibas1(ic)-ibas0(ic)+1
            j=ibas0(ic)
            call dgemm('n','n',n,norb,k,-1.d0,tphin(1,j),
     &           m,coef(j,1,istate,iwf),MBASIS,0.d0,da_orb(1,1,1,ic,istate),m)
            call dgemm('n','n',n,norb,k,-1.d0,t3phin(1,j),
     &           m,coef(j,1,istate,iwf),MBASIS,0.d0,da_d2orb(1,1,1,ic,istate),m)
            call dgemm('n','n',3*n,norb,k,-1.d0,t2phin_all(1,j),
     &           3*m,coef(j,1,istate,iwf),MBASIS,0.d0,da_dorb(1,1,1,1,ic,istate),3*m)
         enddo
      enddo

      end subroutine

c------------------------------------------------------------------------------------

      subroutine orbitalse(iel,x,rvec_en,r_en,iflag)

      use vmc_mod, only: MELEC, MORB, MDET, MCENT
      use phifun, only: d2phin, dphin, n0_ibasis, n0_nbasis
      use phifun, only: phin
      use wfsec, only: iwf
      use coefs, only: coef, nbasis, norb
      use contrl_per, only: iperiodic
      use grid3dflag, only: i3dlagorb, i3dsplorb
      use multislatern, only: ddorbn, detn, dorbn, orbn
      use csfs, only: nstates

      implicit real*8(a-h,o-z)

      dimension x(3,*),rvec_en(3,MELEC,MCENT),r_en(MELEC,MCENT)

      if(iperiodic.eq.0) then
c     get the value and gradients from the 3d-interpolated orbitals
         ier=0
c     spline interplolation
         if(i3dsplorb.ge.1) then
            do istate=1,nstates
               do iorb=1,norb
                  ddorbn(iorb,istate)=0.0d0     ! Don't compute the laplacian
                  dorbn(1:3,iorb,istate)=1.0d0  ! compute the gradients
                  call spline_mo(x(1,iel),iorb,orbn(iorb,istate),
     &                          dorbn(1,iorb,istate),ddorbn(iorb,istate),ier)
               enddo
            enddo
c     Lagrange interpolation
         elseif(i3dlagorb.ge.1) then
            do istate=1,nstates
               call lagrange_mose(1,x(1,iel),orbn(1,istate),ier)
               call lagrange_mos_grade(2,x(1,iel),dorbn(1,1,istate),ier)
               call lagrange_mos_grade(3,x(1,iel),dorbn(1,1,istate),ier)
               call lagrange_mos_grade(4,x(1,iel),dorbn(1,1,istate),ier)
            enddo
         else
            ier=1
         endif 

         if(ier.eq.1) then
c     get basis functions for electron iel
            if(iflag.eq.0) then
               call basis_fnse_vg(iel,rvec_en,r_en)
            else
               call basis_fnse_vgl(iel,rvec_en,r_en)
            endif
            do istate=1,nstates
               do iorb=1,norb
                  orbn(iorb,istate)=0.0d0
                  dorbn(1:3,iorb,istate)=0.0d0
                  ddorbn(iorb,istate)=0.0d0
                  do m0=1,n0_nbasis(iel)
                     m=n0_ibasis(m0,iel)
                     orbn(iorb,istate)=orbn(iorb,istate)+coef(m,iorb,istate,iwf)*phin(m,iel)
                     dorbn(1,iorb,istate)=dorbn(1,iorb,istate)+coef(m,iorb,istate,iwf)*dphin(1,m,iel)
                     dorbn(2,iorb,istate)=dorbn(2,iorb,istate)+coef(m,iorb,istate,iwf)*dphin(2,m,iel)
                     dorbn(3,iorb,istate)=dorbn(3,iorb,istate)+coef(m,iorb,istate,iwf)*dphin(3,m,iel)
                     if(iflag.gt.0) then
                        ddorbn(iorb,istate)=ddorbn(iorb,istate)+coef(m,iorb,istate,iwf)*d2phin(m,iel)
	             endif
                  enddo
               enddo     
            enddo     
         endif
      else
         do istate=1,nstates
            call orbitals_pw_grade(iel,x(1,iel),orbn(1,istate),dorbn(1,1,istate),ddorbn(1,istate))
	 enddo
      endif

      end subroutine
