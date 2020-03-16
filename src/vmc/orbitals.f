      subroutine orbitals(x,rvec_en,r_en)
c Written by Cyrus Umrigar starting from Kevin Schmidt's routine
c Modified by A. Scemama

      use const, only: pi, hb, etrial, delta, deltai, fbias, nelec, imetro, ipr
      use da_orbval, only: da_d2orb, da_dorb, da_orb

      use elec, only: ndn, nup
      implicit real*8(a-h,o-z)



      include 'vmc.h'
      include 'force.h'
      include 'mstates.h'
      include '3dgrid_flags.h'

      parameter (one=1.d0,half=0.5d0)

c note that the dimension of the slater matrices is assumed
c to be given by MMAT_DIM = (MELEC/2)**2, that is there are
c as many ups as downs. If this is not true then be careful if
c nelec is close to MELEC. The Slater matrices must be
c dimensioned at least max(nup**2,ndn**2)

      common /contrl_per/ iperiodic,ibasis
      common /wfsec/ iwftype(MFORCE),iwf,nwftype
      common /force_analy/ iforce_analy

      common /phifun/ phin(MBASIS,MELEC),dphin(3,MBASIS,MELEC)
     &,d2phin(MBASIS,MELEC),d2phin_all(3,3,MBASIS,MELEC),d3phin(3,MBASIS,MELEC)
     &,n0_nbasis(MELEC),n0_ibasis(MBASIS,MELEC),n0_ic(MBASIS,MELEC)
      common /coefs/ coef(MBASIS,MORB,MWF),nbasis,norb
      common /orbval/ orb(MELEC,MORB),dorb(3,MELEC,MORB),ddorb(MELEC,MORB),ndetorb,nadorb

      dimension x(3,*),rvec_en(3,MELEC,MCENT),r_en(MELEC,MCENT)
      dimension bhin(melec,mbasis),dbhin(3*melec,mbasis),d2bhin(melec,mbasis)

      ier=1
      if(iperiodic.eq.0) then

c spline interpolation
        if(i3dsplorb.eq.2) then
          do 23 i=1,nelec
            ier = 0.d0
            do 21 iorb=1,norb+nadorb
              ddorb(i,iorb)=1.d0    ! compute the laplacian
              dorb(1,i,iorb)=1.d0   ! compute the gradients
              dorb(2,i,iorb)=1.d0   ! compute the gradients
              dorb(3,i,iorb)=1.d0   ! compute the gradients
   21         call spline_mo (x(1,i),iorb,orb(i,iorb),dorb(1,i,iorb),ddorb(i,iorb),ier)

            if(ier.eq.1) then
              call basis_fnse_vgl(i,rvec_en,r_en)
              do 22 iorb=1,norb+nadorb
                orb(i,iorb)=0.d0
                dorb(1,i,iorb)=0.d0
                dorb(2,i,iorb)=0.d0
                dorb(3,i,iorb)=0.d0
                ddorb(i,iorb)=0.d0
c               do 22 m=1,nbasis
                do 22 m0=1,n0_nbasis(i)
                  m=n0_ibasis(m0,i)
                  orb(i,iorb)=orb(i,iorb)+coef(m,iorb,iwf)*phin(m,i)
                  dorb(1,i,iorb)=dorb(1,i,iorb)+coef(m,iorb,iwf)*dphin(1,m,i)
                  dorb(2,i,iorb)=dorb(2,i,iorb)+coef(m,iorb,iwf)*dphin(2,m,i)
                  dorb(3,i,iorb)=dorb(3,i,iorb)+coef(m,iorb,iwf)*dphin(3,m,i)
                  ddorb(i,iorb)=ddorb(i,iorb)+coef(m,iorb,iwf)*d2phin(m,i)
   22         continue
            endif
   23     continue

c Lagrange interpolation
        elseif(i3dlagorb.eq.2) then
         do 25 i=1,nelec
           ier=0
           call lagrange_mos(1,x(1,i),orb,i,ier)
           call lagrange_mos_grad(2,x(1,i),dorb,i,ier)
           call lagrange_mos_grad(3,x(1,i),dorb,i,ier)
           call lagrange_mos_grad(4,x(1,i),dorb,i,ier)
           call lagrange_mos(5,x(1,i),ddorb,i,ier)

           if(ier.eq.1) then
             call basis_fnse_vgl(i,rvec_en,r_en)
             do 24 iorb=1,norb+nadorb
               orb(i,iorb)=0.d0
               dorb(1,i,iorb)=0.d0
               dorb(2,i,iorb)=0.d0
               dorb(3,i,iorb)=0.d0
               ddorb(i,iorb)=0.d0
c              do 24 m=1,nbasis
               do 24 m0=1,n0_nbasis(i)
                 m=n0_ibasis(m0,i)
                 orb(i,iorb)=orb(i,iorb)+coef(m,iorb,iwf)*phin(m,i)
                 dorb(1,i,iorb)=dorb(1,i,iorb)+coef(m,iorb,iwf)*dphin(1,m,i)
                 dorb(2,i,iorb)=dorb(2,i,iorb)+coef(m,iorb,iwf)*dphin(2,m,i)
                 dorb(3,i,iorb)=dorb(3,i,iorb)+coef(m,iorb,iwf)*dphin(3,m,i)
   24            ddorb(i,iorb)=ddorb(i,iorb)+coef(m,iorb,iwf)*d2phin(m,i)
           endif
   25    continue

c no 3d interpolation
        else 

c get basis functions for all electrons
         call basis_fns_vgl(x,rvec_en,r_en)

c in alternativa al loop 26
c        do jbasis=1,nbasis
c         i=0
c         do ielec=1,nelec
c          bhin(ielec,jbasis)=phin(jbasis,ielec)
c          do l=1,3
c           i=i+1
c           dbhin(i,jbasis)=dphin(l,jbasis,ielec)
c          enddo
c          d2bhin(ielec,jbasis)=d2phin(jbasis,ielec)
c         enddo
c        enddo
c        call dgemm('n','n',  nelec,norb,nbasis,1.d0,bhin,   melec,  coef(1,1,iwf),mbasis,0.d0,orb,   melec)
c        call dgemm('n','n',3*nelec,norb,nbasis,1.d0,dbhin,3*melec,  coef(1,1,iwf),mbasis,0.d0,dorb,3*melec)
c        call dgemm('n','n',  nelec,norb,nbasis,1.d0,d2bhin, melec,  coef(1,1,iwf),mbasis,0.d0,ddorb, melec)
         do 26 iorb=1,norb+nadorb
           do 26 i=1,nelec
            orb(i,iorb)=0
            dorb(1,i,iorb)=0
            dorb(2,i,iorb)=0
            dorb(3,i,iorb)=0
            ddorb(i,iorb)=0
c           do 26 m=1,nbasis
            do 26 m0=1,n0_nbasis(i)
             m=n0_ibasis(m0,i)
             orb  (  i,iorb)=orb  (  i,iorb)+coef(m,iorb,iwf)*phin  ( m,i)
             dorb (1,i,iorb)=dorb (1,i,iorb)+coef(m,iorb,iwf)*dphin (1,m,i)
             dorb (2,i,iorb)=dorb (2,i,iorb)+coef(m,iorb,iwf)*dphin (2,m,i)
             dorb (3,i,iorb)=dorb (3,i,iorb)+coef(m,iorb,iwf)*dphin (3,m,i)
   26        ddorb(  i,iorb)=ddorb(  i,iorb)+coef(m,iorb,iwf)*d2phin( m,i)
       endif

       if(iforce_analy.eq.1) call da_orbitals

       else
        call orbitals_pw(x,orb,dorb,ddorb)
      endif

      if(ipr.ge.4) then
        do 260 iorb=1,norb+nadorb
  260     write(6,'(''iorb,orb='',i4,1000f15.11)') iorb,(orb(i,iorb),i=1,nelec)
        do 270 iorb=1,norb+nadorb
  270     write(6,'(''iorb,d2orb='',i4,1000f15.11)') iorb,(ddorb(i,iorb),i=1,nelec)
        do 280 k=1,3
          do 280 iorb=1,norb+nadorb
  280       write(6,'(''iorb,dorb='',2i4,1000f12.8)') k,iorb,(dorb(k,i,iorb),i=1,nelec)
      endif

      return
      end
c------------------------------------------------------------------------------------
      subroutine virtual_orbitals
      use const, only: pi, hb, etrial, delta, deltai, fbias, nelec, imetro, ipr
      implicit real*8(a-h,o-z)

c compute values of extra ('virtual') orbitals needed for optorb operators
c assuming that basis function values in phin are up to date

      include 'vmc.h'
      include 'force.h'
      include 'optorb.h'

      common /phifun/ phin(MBASIS,MELEC),dphin(3,MBASIS,MELEC)
     &,d2phin(MBASIS,MELEC),d2phin_all(3,3,MBASIS,MELEC),d3phin(3,MBASIS,MELEC)
     &,n0_nbasis(MELEC),n0_ibasis(MBASIS,MELEC),n0_ic(MBASIS,MELEC)
      common /coefs/ coef(MBASIS,MORB,MWF),nbasis,norb
      common /orbval/ orb(MELEC,MORB),dorb(3,MELEC,MORB),ddorb(MELEC,MORB),ndetorb,nadorb

      common /optwf_contrl/ ioptjas,ioptorb,ioptci,nparm

      dimension bhin(melec,mbasis),dbhin(3,melec,mbasis),d2bhin(melec,mbasis)

      if (nadorb.eq.0.or.(ioptorb.eq.0.and.ioptci.eq.0)) return

c primary geometry only:
      iwf=1
      if(norb+nadorb.gt.MORB) then
        write(6,'(''VIRTUAL_ORB: Too many orbitals, norb + nadorb='',
     &  i4,'' > MORB='',i4)') norb+nadorb,MORB
        call fatal_error('Aborted')
      endif

      do i=1,nelec
        call dcopy(nbasis,phin(1,i),1,bhin(i,1),melec)
        call dcopy(nbasis,dphin(1,1,i),3,dbhin(1,i,1),3*melec)
        call dcopy(nbasis,dphin(2,1,i),3,dbhin(2,i,1),3*melec)
        call dcopy(nbasis,dphin(3,1,i),3,dbhin(3,i,1),3*melec)
        call dcopy(nbasis,d2phin(1,i),1,d2bhin(i,1),melec)
      enddo
      call dgemm('n','n',  nelec,nadorb,nbasis,1.d0,  bhin,  melec,coef(1,norb+1,iwf),mbasis,0.d0,  orb(  1,norb+1),  melec)
      call dgemm('n','n',3*nelec,nadorb,nbasis,1.d0, dbhin,3*melec,coef(1,norb+1,iwf),mbasis,0.d0, dorb(1,1,norb+1),3*melec)
      call dgemm('n','n',  nelec,nadorb,nbasis,1.d0,d2bhin,  melec,coef(1,norb+1,iwf),mbasis,0.d0,ddorb(  1,norb+1),  melec)

c     do 25 iorb=norb+1,norb+nadorb
c       do 25 i=1,nelec
c         orb(i,iorb)=0.d0
c         dorb(1,i,iorb)=0.d0
c         dorb(2,i,iorb)=0.d0
c         dorb(3,i,iorb)=0.d0
c         ddorb(i,iorb)=0.d0
c         do 25 m=1,nbasis
c           orb(i,iorb)=orb(i,iorb)+coef(m,iorb,iwf)*phin(m,i)
c           dorb(1,i,iorb)=dorb(1,i,iorb)+coef(m,iorb,iwf)*dphin(1,m,i)
c           dorb(2,i,iorb)=dorb(2,i,iorb)+coef(m,iorb,iwf)*dphin(2,m,i)
c           dorb(3,i,iorb)=dorb(3,i,iorb)+coef(m,iorb,iwf)*dphin(3,m,i)
c           ddorb(i,iorb)=ddorb(i,iorb)+coef(m,iorb,iwf)*d2phin(m,i)
c25   continue

      return
      end
c------------------------------------------------------------------------------------
      subroutine da_orbitals
      use atom, only: znuc, cent, pecent, iwctype, nctype, ncent

      use const, only: pi, hb, etrial, delta, deltai, fbias, nelec, imetro, ipr
      use da_orbval, only: da_d2orb, da_dorb, da_orb

      use numbas1, only: iwlbas, nbastyp
      implicit real*8(a-h,o-z)




      include 'vmc.h'
      include 'numbas.h'
      include 'force.h'

      common /phifun/ phin(MBASIS,MELEC),dphin(3,MBASIS,MELEC)
     &,d2phin(MBASIS,MELEC),d2phin_all(3,3,MBASIS,MELEC),d3phin(3,MBASIS,MELEC)
     &,n0_nbasis(MELEC),n0_ibasis(MBASIS,MELEC),n0_ic(MBASIS,MELEC)
      common /coefs/ coef(MBASIS,MORB,MWF),nbasis,norb


      common /numbas2/ ibas0(MCENT),ibas1(MCENT)


      common /orbval/ orb(MELEC,MORB),dorb(3,MELEC,MORB),ddorb(MELEC,MORB),ndetorb,nadorb

      common /wfsec/ iwftype(MFORCE),iwf,nwftype

      dimension dum(MELEC,MORB,MCENT),aux(MELEC)
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
      m=3*melec
      do 50 ic=1,ncent
        k=ibas1(ic)-ibas0(ic)+1
        j=ibas0(ic)
      call dgemm('n','n',  n,norb,k,-1.d0,tphin(1,j)     ,  m,coef(j,1,iwf),mbasis,0.d0,da_orb(1,1,1,ic)   ,  m)
      call dgemm('n','n',  n,norb,k,-1.d0,t3phin(1,j)    ,  m,coef(j,1,iwf),mbasis,0.d0,da_d2orb(1,1,1,ic) ,  m)
  50  call dgemm('n','n',3*n,norb,k,-1.d0,t2phin_all(1,j),3*m,coef(j,1,iwf),mbasis,0.d0,da_dorb(1,1,1,1,ic),3*m)

      return
      end
c------------------------------------------------------------------------------------
      subroutine orbitalse(iel,x,rvec_en,r_en,iflag)

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


      if(iperiodic.eq.0) then

c get the value and gradients from the 3d-interpolated orbitals
        ier=0
c spline interplolation
        if(i3dsplorb.ge.1) then
          do 15 iorb=1,norb
            ddorb(iorb)=0    ! Don't compute the laplacian
            dorb(1,iorb)=1   ! compute the gradients
            dorb(2,iorb)=1   ! compute the gradients
            dorb(3,iorb)=1   ! compute the gradients
   15       call spline_mo (x(1,iel),iorb,orb(iorb),dorb(1,iorb),ddorb(iorb),ier)

c Lagrange interpolation
         elseif(i3dlagorb.ge.1) then
          call lagrange_mose(1,x(1,iel),orb,ier)
          call lagrange_mos_grade(2,x(1,iel),dorb,ier)
          call lagrange_mos_grade(3,x(1,iel),dorb,ier)
          call lagrange_mos_grade(4,x(1,iel),dorb,ier)
         else
          ier=1
        endif 

        if(ier.eq.1) then
c get basis functions for electron iel
          if(iflag.eq.0) then
            call basis_fnse_vg(iel,rvec_en,r_en)
           else
            call basis_fnse_vgl(iel,rvec_en,r_en)
          endif

          do 25 iorb=1,norb
            orb(iorb)=0
            dorb(1,iorb)=0
            dorb(2,iorb)=0
            dorb(3,iorb)=0
            ddorb(iorb)=0
c           do 25 m=1,nbasis
            do 25 m0=1,n0_nbasis(iel)
             m=n0_ibasis(m0,iel)
             orb(iorb)=orb(iorb)+coef(m,iorb,iwf)*phin(m,iel)
             dorb(1,iorb)=dorb(1,iorb)+coef(m,iorb,iwf)*dphin(1,m,iel)
             dorb(2,iorb)=dorb(2,iorb)+coef(m,iorb,iwf)*dphin(2,m,iel)
             dorb(3,iorb)=dorb(3,iorb)+coef(m,iorb,iwf)*dphin(3,m,iel)
             if(iflag.gt.0) ddorb(iorb)=ddorb(iorb)+coef(m,iorb,iwf)*d2phin(m,iel)
   25     continue
        endif
       else
        call orbitals_pw_grade(iel,x(1,iel),orb,dorb,ddorb)
      endif

      return
      end
c------------------------------------------------------------------------------------
