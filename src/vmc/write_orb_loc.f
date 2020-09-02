      subroutine write_orb_loc
c Written by Cyrus Umrigar and Claudia Filippi, starting from Kevin Schmidt's routine
c Reads in localized orbitals, in either
c Modified by A. Scemama (printing in a GAMESS-like format)
c 1) a slater basis
c 2) a gaussian basis
      use atom, only: znuc, iwctype, nctype, ncent
      use ghostatom, only: newghostype
      use const, only: nelec
      use numbas, only: numr
      use coefs, only: coef, nbasis, norb
      use basis, only: zex, betaq, n1s, n2s, n2p, n3s, n3p, n3dzr, n3dx2, n3dxy, n3dxz, n3dyz
      use basis, only: n4s, n4p, n4fxxx, n4fyyy, n4fzzz, n4fxxy, n4fxxz, n4fyyx, n4fyyz
      use basis, only: n4fzzx, n4fzzy, n4fxyz, nsa, npa, ndzra, ndxya, ndxza, ndyza, ndx2a

      implicit real*8(a-h,o-z)

      parameter (zero=0.d0,one=1.d0)
      parameter(nprime=10)

      include 'vmc.h'
      include 'force.h'
      include 'ewald.h'
      include 'pseudo.h'
      include 'numbas.h'


c Check that nbasis in lcao matches specified basis on all centers
      ncheck=0
      do 20 ic=1,ncent
        i=iwctype(ic)
   20   ncheck=ncheck+ iabs(n1s(i))+iabs(n2s(i))
     &  +iabs(n2p(1,i))+iabs(n2p(2,i))+iabs(n2p(3,i))
     &  +iabs(n3s(i))+iabs(n3p(1,i))+iabs(n3p(2,i))+iabs(n3p(3,i))
     &  +iabs(n3dzr(i))+iabs(n3dx2(i))
     &  +iabs(n3dxy(i))+iabs(n3dxz(i))+iabs(n3dyz(i))
     &  +iabs(n4s(i))+iabs(n4p(1,i))+iabs(n4p(2,i))+iabs(n4p(3,i))
     &  +iabs(n4fxxx(i))+iabs(n4fyyy(i))+iabs(n4fzzz(i))
     &  +iabs(n4fxxy(i))+iabs(n4fxxz(i))+iabs(n4fyyx(i))
     &  +iabs(n4fyyz(i))+iabs(n4fzzx(i))+iabs(n4fzzy(i))
     &  +iabs(n4fxyz(i))
     &  +iabs(nsa(i))+iabs(npa(1,i))+iabs(npa(2,i))+iabs(npa(3,i))
     &  +iabs(ndzra(i))+iabs(ndx2a(i))
     &  +iabs(ndxya(i))+iabs(ndxza(i))+iabs(ndyza(i))


      if(ncheck.ne.nbasis) then
        write(6,'(''ncheck = '',i6)') ncheck
        write(6,'(''nbasis = '',i6)') nbasis
        call fatal_error('WRITE_BAS: ncheck not equal to nbasis')
      endif

c Exponent for asymptotic basis
      betaq=zero
      do 30 ic=1,ncent
   30   betaq=betaq+znuc(iwctype(ic))
      betaq=betaq-nelec+one

      write(6,'(/,''center type'',(12i4))') (i,i=1,nctype)
      write(6,'(/,''1s'',t11,(12i5))') (n1s(i),i=1,nctype)
      write(6,'(''2s'',t11,(12i5))') (n2s(i),i=1,nctype)
      write(6,'(''2px'',t11,(12i5))') (n2p(1,i),i=1,nctype)
      write(6,'(''2py'',t11,(12i5))') (n2p(2,i),i=1,nctype)
      write(6,'(''2pz'',t11,(12i5))') (n2p(3,i),i=1,nctype)
      write(6,'(''3s'',t11,(12i5))') (n3s(i),i=1,nctype)
      write(6,'(''3px'',t11,(12i5))') (n3p(1,i),i=1,nctype)
      write(6,'(''3py'',t11,(12i5))') (n3p(2,i),i=1,nctype)
      write(6,'(''3pz'',t11,(12i5))') (n3p(3,i),i=1,nctype)
      write(6,'(''3dzr'',t11,(12i5))') (n3dzr(i),i=1,nctype)
      write(6,'(''3dx2'',t11,(12i5))') (n3dx2(i),i=1,nctype)
      write(6,'(''3dxy'',t11,(12i5))') (n3dxy(i),i=1,nctype)
      write(6,'(''3dxz'',t11,(12i5))') (n3dxz(i),i=1,nctype)
      write(6,'(''3dyz'',t11,(12i5))') (n3dyz(i),i=1,nctype)
      write(6,'(''4s'',t11,(12i5))') (n4s(i),i=1,nctype)
      write(6,'(''4px'',t11,(12i5))') (n4p(1,i),i=1,nctype)
      write(6,'(''4py'',t11,(12i5))') (n4p(2,i),i=1,nctype)
      write(6,'(''4pz'',t11,(12i5))') (n4p(3,i),i=1,nctype)
      write(6,'(''4fxxx'',t11,(12i5))') (n4fxxx(i),i=1,nctype)
      write(6,'(''4fyyy'',t11,(12i5))') (n4fyyy(i),i=1,nctype)
      write(6,'(''4fzzz'',t11,(12i5))') (n4fzzz(i),i=1,nctype)
      write(6,'(''4fxxy'',t11,(12i5))') (n4fxxy(i),i=1,nctype)
      write(6,'(''4fxxz'',t11,(12i5))') (n4fxxz(i),i=1,nctype)
      write(6,'(''4fyyx'',t11,(12i5))') (n4fyyx(i),i=1,nctype)
      write(6,'(''4fyyz'',t11,(12i5))') (n4fyyz(i),i=1,nctype)
      write(6,'(''4fzzx'',t11,(12i5))') (n4fzzx(i),i=1,nctype)
      write(6,'(''4fzzy'',t11,(12i5))') (n4fzzy(i),i=1,nctype)
      write(6,'(''4fxyz'',t11,(12i5))') (n4fxyz(i),i=1,nctype)
      write(6,'(''sa'',t11,(12i5))') (nsa(i),i=1,nctype)
      write(6,'(''pxa'',t11,(12i5))') (npa(1,i),i=1,nctype)
      write(6,'(''pya'',t11,(12i5))') (npa(2,i),i=1,nctype)
      write(6,'(''pza'',t11,(12i5))') (npa(3,i),i=1,nctype)
      write(6,'(''dzra'',t11,(12i5))') (ndzra(i),i=1,nctype)
      write(6,'(''dx2a'',t11,(12i5))') (ndx2a(i),i=1,nctype)
      write(6,'(''dxya'',t11,(12i5))') (ndxya(i),i=1,nctype)
      write(6,'(''dxza'',t11,(12i5))') (ndxza(i),i=1,nctype)
      write(6,'(''dyza'',t11,(12i5))') (ndyza(i),i=1,nctype)
      write(6,'(/,''charge'',t12,(12f5.0))')(znuc(i),i=1,nctype)
      write(6,*)

      if(newghostype.gt.0) then
        write(6,'(/,''ghost type'',(12i4))') (i,i=nctype+1,nctype+newghostype)
        write(6,'(/,''1s'',t11,(12i5))') (n1s(i),i=nctype+1,nctype+newghostype)
        write(6,'(''2s'',t11,(12i5))')   (n2s(i),i=nctype+1,nctype+newghostype)
        write(6,'(''2px'',t11,(12i5))')  (n2p(1,i),i=nctype+1,nctype+newghostype)
        write(6,'(''2py'',t11,(12i5))')  (n2p(2,i),i=nctype+1,nctype+newghostype)
        write(6,'(''2pz'',t11,(12i5))')  (n2p(3,i),i=nctype+1,nctype+newghostype)
        write(6,'(''3s'',t11,(12i5))')   (n3s(i),i=nctype+1,nctype+newghostype)
        write(6,'(''3px'',t11,(12i5))')  (n3p(1,i),i=nctype+1,nctype+newghostype)
        write(6,'(''3py'',t11,(12i5))')  (n3p(2,i),i=nctype+1,nctype+newghostype)
        write(6,'(''3pz'',t11,(12i5))')  (n3p(3,i),i=nctype+1,nctype+newghostype)
        write(6,'(''3dzr'',t11,(12i5))') (n3dzr(i),i=nctype+1,nctype+newghostype)
        write(6,'(''3dx2'',t11,(12i5))') (n3dx2(i),i=nctype+1,nctype+newghostype)
        write(6,'(''3dxy'',t11,(12i5))') (n3dxy(i),i=nctype+1,nctype+newghostype)
        write(6,'(''3dxz'',t11,(12i5))') (n3dxz(i),i=nctype+1,nctype+newghostype)
        write(6,'(''3dyz'',t11,(12i5))') (n3dyz(i),i=nctype+1,nctype+newghostype)
        write(6,'(''4s'',t11,(12i5))')   (n4s(i),i=nctype+1,nctype+newghostype)
        write(6,'(''4px'',t11,(12i5))')  (n4p(1,i),i=nctype+1,nctype+newghostype)
        write(6,'(''4py'',t11,(12i5))')  (n4p(2,i),i=nctype+1,nctype+newghostype)
        write(6,'(''4pz'',t11,(12i5))')  (n4p(3,i),i=nctype+1,nctype+newghostype)
        write(6,'(''4fxxx'',t11,(12i5))') (n4fxxx(i),i=nctype+1,nctype+newghostype)
        write(6,'(''4fyyy'',t11,(12i5))') (n4fyyy(i),i=nctype+1,nctype+newghostype)
        write(6,'(''4fzzz'',t11,(12i5))') (n4fzzz(i),i=nctype+1,nctype+newghostype)
        write(6,'(''4fxxy'',t11,(12i5))') (n4fxxy(i),i=nctype+1,nctype+newghostype)
        write(6,'(''4fxxz'',t11,(12i5))') (n4fxxz(i),i=nctype+1,nctype+newghostype)
        write(6,'(''4fyyx'',t11,(12i5))') (n4fyyx(i),i=nctype+1,nctype+newghostype)
        write(6,'(''4fyyz'',t11,(12i5))') (n4fyyz(i),i=nctype+1,nctype+newghostype)
        write(6,'(''4fzzx'',t11,(12i5))') (n4fzzx(i),i=nctype+1,nctype+newghostype)
        write(6,'(''4fzzy'',t11,(12i5))') (n4fzzy(i),i=nctype+1,nctype+newghostype)
        write(6,'(''4fxyz'',t11,(12i5))') (n4fxyz(i),i=nctype+1,nctype+newghostype)
        write(6,'(''sa'',t11,(12i5))')   (nsa(i),i=nctype+1,nctype+newghostype)
        write(6,'(''pxa'',t11,(12i5))')  (npa(1,i),i=nctype+1,nctype+newghostype)
        write(6,'(''pya'',t11,(12i5))')  (npa(2,i),i=nctype+1,nctype+newghostype)
        write(6,'(''pza'',t11,(12i5))')  (npa(3,i),i=nctype+1,nctype+newghostype)
        write(6,'(''dzra'',t11,(12i5))') (ndzra(i),i=nctype+1,nctype+newghostype)
        write(6,'(''dx2a'',t11,(12i5))') (ndx2a(i),i=nctype+1,nctype+newghostype)
        write(6,'(''dxya'',t11,(12i5))') (ndxya(i),i=nctype+1,nctype+newghostype)
        write(6,'(''dxza'',t11,(12i5))') (ndxza(i),i=nctype+1,nctype+newghostype)
        write(6,'(''dyza'',t11,(12i5))') (ndyza(i),i=nctype+1,nctype+newghostype)
        write(6,'(/,''charge'',t12,(12f5.0))')(znuc(i),i=nctype+1,nctype+newghostype)
        write(6,*)
      endif

      write(6,'(/,''no. of orbitals ='',t30,i10)') norb

c the following sets up strings to identify basis functions on centers for
c the printout of coefficients and screening constants.


  210 continue

      iu=6
      if(nbasis.gt.15) then
        iu=45
        write(iu,'(/,''Orbitals'')')
      endif

      i=1
      do while (i.le.norb)
       imax=i+4
       if (i+4.gt.norb) then
        imax=norb
       endif
       write (iu, 211) (l, l=i,imax)
C      write (iu, 212) TODO : print out the symmetry

       j=1
       do ic=1,ncent 

        do k=1,abs(n1s(iwctype(ic)))
         write (iu, 213) j, iwctype(ic), ic, '1s', (coef(j,l,1), l=i, imax)
         j=j+1
        enddo


        do k=1,abs(n2s(iwctype(ic)))
         write (iu, 213) j, iwctype(ic), ic, '2s', (coef(j,l,1), l=i, imax)
         j=j+1
        enddo

        do k=1,abs(n2p(1,iwctype(ic)))
         write (iu, 213) j, iwctype(ic), ic, '2px', (coef(j,l,1), l=i, imax)
         j=j+1
        enddo

        do k=1,abs(n2p(2,iwctype(ic)))
         write (iu, 213) j, iwctype(ic), ic, '2py', (coef(j,l,1), l=i, imax)
         j=j+1
        enddo

        do k=1,abs(n2p(3,iwctype(ic)))
         write (iu, 213) j, iwctype(ic), ic, '2pz', (coef(j,l,1), l=i, imax)
         j=j+1
        enddo


        do k=1,abs(n3s(iwctype(ic)))
         write (iu, 213) j, iwctype(ic), ic, '3s', (coef(j,l,1), l=i, imax)
         j=j+1
        enddo

        do k=1,abs(n3p(1,iwctype(ic)))
         write (iu, 213) j, iwctype(ic), ic, '3px', (coef(j,l,1), l=i, imax)
         j=j+1
        enddo

        do k=1,abs(n3p(2,iwctype(ic)))
         write (iu, 213) j, iwctype(ic), ic, '3py', (coef(j,l,1), l=i, imax)
         j=j+1
        enddo

        do k=1,abs(n3p(3,iwctype(ic)))
         write (iu, 213) j, iwctype(ic), ic, '3pz', (coef(j,l,1), l=i, imax)
         j=j+1
        enddo

        do k=1,abs(n3dzr(iwctype(ic)))
         write (iu, 213) j, iwctype(ic), ic, '3dzr', (coef(j,l,1), l=i, imax)
         j=j+1
        enddo

        do k=1,abs(n3dx2(iwctype(ic)))
         write (iu, 213) j, iwctype(ic), ic, '3dx2', (coef(j,l,1), l=i, imax)
         j=j+1
        enddo

        do k=1,abs(n3dxy(iwctype(ic)))
         write (iu, 213) j, iwctype(ic), ic, '3dxy', (coef(j,l,1), l=i, imax)
         j=j+1
        enddo

        do k=1,abs(n3dxz(iwctype(ic)))
         write (iu, 213) j, iwctype(ic), ic, '3dxz', (coef(j,l,1), l=i, imax)
         j=j+1
        enddo

        do k=1,abs(n3dyz(iwctype(ic)))
         write (iu, 213) j, iwctype(ic), ic, '3dyz', (coef(j,l,1), l=i, imax)
         j=j+1
        enddo


        do k=1,abs(n4s(iwctype(ic)))
         write (iu, 213) j, iwctype(ic), ic, '4s', (coef(j,l,1), l=i, imax)
         j=j+1
        enddo

        do k=1,abs(n4p(1,iwctype(ic)))
         write (iu, 213) j, iwctype(ic), ic, '4px', (coef(j,l,1), l=i, imax)
         j=j+1
        enddo

        do k=1,abs(n4p(2,iwctype(ic)))
         write (iu, 213) j, iwctype(ic), ic, '4py', (coef(j,l,1), l=i, imax)
         j=j+1
        enddo

        do k=1,abs(n4p(3,iwctype(ic)))
         write (iu, 213) j, iwctype(ic), ic, '4pz', (coef(j,l,1), l=i, imax)
         j=j+1
        enddo

        do k=1,abs(n4fxxx(iwctype(ic)))
         write (iu, 213) j, iwctype(ic), ic, '4fxxx', (coef(j,l,1), l=i,imax)
         j=j+1
        enddo

        do k=1,abs(n4fyyy(iwctype(ic)))
         write (iu, 213) j, iwctype(ic), ic, '4fyyy', (coef(j,l,1), l=i,imax)
         j=j+1
        enddo

        do k=1,abs(n4fzzz(iwctype(ic)))
         write (iu, 213) j, iwctype(ic), ic, '4fzzz', (coef(j,l,1), l=i,imax)
         j=j+1
        enddo

        do k=1,abs(n4fxxy(iwctype(ic)))
         write (iu, 213) j, iwctype(ic), ic, '4fxxy', (coef(j,l,1), l=i,imax)
         j=j+1
        enddo

        do k=1,abs(n4fxxz(iwctype(ic)))
         write (iu, 213) j, iwctype(ic), ic, '4fxxz', (coef(j,l,1), l=i,imax)
         j=j+1
        enddo

        do k=1,abs(n4fyyx(iwctype(ic)))
         write (iu, 213) j, iwctype(ic), ic, '4fyyx', (coef(j,l,1), l=i,imax)
         j=j+1
        enddo

        do k=1,abs(n4fyyz(iwctype(ic)))
         write (iu, 213) j, iwctype(ic), ic, '4fyyz', (coef(j,l,1), l=i,imax)
         j=j+1
        enddo

        do k=1,abs(n4fzzx(iwctype(ic)))
         write (iu, 213) j, iwctype(ic), ic, '4fzzx', (coef(j,l,1), l=i,imax)
         j=j+1
        enddo

        do k=1,abs(n4fzzy(iwctype(ic)))
         write (iu, 213) j, iwctype(ic), ic, '4fzzy', (coef(j,l,1), l=i,imax)
         j=j+1
        enddo

        do k=1,abs(n4fxyz(iwctype(ic)))
         write (iu, 213) j, iwctype(ic), ic, '4fxyz', (coef(j,l,1), l=i,imax)
         j=j+1
        enddo


        do k=1,abs(nsa(iwctype(ic)))
         write (iu, 213) j, iwctype(ic), ic, 's', (coef(j,l,1), l=i, imax)
         j=j+1
        enddo

        do k=1,abs(npa(1,iwctype(ic)))
         write (iu, 213) j, iwctype(ic), ic, 'px', (coef(j,l,1), l=i, imax)
         j=j+1
        enddo

        do k=1,abs(npa(2,iwctype(ic)))
         write (iu, 213) j, iwctype(ic), ic, 'py', (coef(j,l,1), l=i, imax)
         j=j+1
        enddo

        do k=1,abs(npa(3,iwctype(ic)))
         write (iu, 213) j, iwctype(ic), ic, 'pz', (coef(j,l,1), l=i, imax)
         j=j+1
        enddo

        do k=1,abs(ndzra(iwctype(ic)))
         write (iu, 213) j, iwctype(ic), ic, 'dzr', (coef(j,l,1), l=i, imax)
         j=j+1
        enddo

        do k=1,abs(ndx2a(iwctype(ic)))
         write (iu, 213) j, iwctype(ic), ic, 'dx2', (coef(j,l,1), l=i, imax)
         j=j+1
        enddo

        do k=1,abs(ndxya(iwctype(ic)))
         write (iu, 213) j, iwctype(ic), ic, 'dxy', (coef(j,l,1), l=i, imax)
         j=j+1
        enddo

        do k=1,abs(ndxza(iwctype(ic)))
         write (iu, 213) j, iwctype(ic), ic, 'dxz', (coef(j,l,1), l=i, imax)
         j=j+1
        enddo

        do k=1,abs(ndyza(iwctype(ic)))
         write (iu, 213) j, iwctype(ic), ic, 'dyz', (coef(j,l,1), l=i, imax)
         j=j+1
        enddo

       enddo

       i=i+5
       write (iu, *) 
      enddo ! do while (i.le.norb)

  211 format (5X,8X,5(1X,I10))
  212 format (5X,8X,5(1X,A10))
  213 format (I5,I3,I3,A5,5(1X,F10.6))

      if(numr.eq.0) then
        write(iu,'(''screening constants'')')
        write(iu,'(12f10.6)') (zex(i,1),i=1,nbasis)
      endif

      return
      end
