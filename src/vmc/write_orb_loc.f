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
      use basis, only: zex, betaq
      use basis, only: ns, npx, npy, npz, ndxx, ndxy, ndxz, ndyy, ndyz, ndzz
      use basis, only: nfxxx, nfxxy, nfxxz, nfxyy, nfxyz, nfxzz, nfyyy, nfyyz, nfyzz, nfzzz
      use basis, only: ngxxxx, ngxxxy, ngxxxz, ngxxyy, ngxxyz, ngxxzz, ngxyyy, ngxyyz
      use basis, only: ngxyzz, ngxzzz, ngyyyy, ngyyyz, ngyyzz, ngyzzz, ngzzzz
      use contrl_file,    only: ounit, errunit

      use precision_kinds, only: dp
      implicit none

      integer :: i, iabs, ic, imax, iu
      integer :: j, k, l, ncheck
      integer, parameter :: nprime = 10
      real(dp), parameter :: zero = 0.d0
      real(dp), parameter :: one = 1.d0




c Check that nbasis in lcao matches specified basis on all centers
      ncheck=0
      do ic=1,ncent
        i=iwctype(ic)
        ncheck=ncheck+ ns(i)
     &  + npx(i) + npy(i) + npz(i)
     &  + ndxx(i) + ndxy(i) + ndxz(i) + ndyy(i) + ndyz(i) + ndzz(i)
     &  + nfxxx(i) + nfxxy(i) + nfxxz(i) + nfxyy(i) + nfxyz(i)
     &  + nfxzz(i) + nfyyy(i) + nfyyz(i) + nfyzz(i) + nfzzz(i)
     &  + ngxxxx(i) + ngxxxy(i) + ngxxxz(i) + ngxxyy(i) + ngxxyz(i)
     &  + ngxxzz(i) + ngxyyy(i) + ngxyyz(i) + ngxyzz(i) + ngxzzz(i)
     &  + ngyyyy(i) + ngyyyz(i) + ngyyzz(i) + ngyzzz(i) + ngzzzz(i)
      enddo


      if(ncheck.ne.nbasis) then
        write(ounit,'(''ncheck = '',i6)') ncheck
        write(ounit,'(''nbasis = '',i6)') nbasis
        call fatal_error('WRITE_BAS: ncheck not equal to nbasis')
      endif

c Exponent for asymptotic basis
      betaq=zero
      do ic=1,ncent
        betaq=betaq+znuc(iwctype(ic))
      enddo
      betaq=betaq-nelec+one

      write(ounit,'(/,''center type'',(12i4))') (i,i=1,nctype)
      write(ounit,'(/,''s'',t11,(12i5))') (ns(i),i=1,nctype)
      write(ounit,'(''px'',t11,(12i5))') (npx(i),i=1,nctype)
      write(ounit,'(''py'',t11,(12i5))') (npy(i),i=1,nctype)
      write(ounit,'(''pz'',t11,(12i5))') (npz(i),i=1,nctype)
      write(ounit,'(''dxx'',t11,(12i5))') (ndxx(i),i=1,nctype)
      write(ounit,'(''dxy'',t11,(12i5))') (ndxy(i),i=1,nctype)
      write(ounit,'(''dxz'',t11,(12i5))') (ndxz(i),i=1,nctype)
      write(ounit,'(''dyy'',t11,(12i5))') (ndyy(i),i=1,nctype)
      write(ounit,'(''dyz'',t11,(12i5))') (ndyz(i),i=1,nctype)
      write(ounit,'(''dzz'',t11,(12i5))') (ndzz(i),i=1,nctype)
      write(ounit,'(''fxxx'',t11,(12i5))') (nfxxx(i),i=1,nctype)
      write(ounit,'(''fxxy'',t11,(12i5))') (nfxxy(i),i=1,nctype)
      write(ounit,'(''fxxz'',t11,(12i5))') (nfxxz(i),i=1,nctype)
      write(ounit,'(''fxyy'',t11,(12i5))') (nfxyy(i),i=1,nctype)
      write(ounit,'(''fxyz'',t11,(12i5))') (nfxyz(i),i=1,nctype)
      write(ounit,'(''fxzz'',t11,(12i5))') (nfxzz(i),i=1,nctype)
      write(ounit,'(''fyyy'',t11,(12i5))') (nfyyy(i),i=1,nctype)
      write(ounit,'(''fyyz'',t11,(12i5))') (nfyyz(i),i=1,nctype)
      write(ounit,'(''fyzz'',t11,(12i5))') (nfyzz(i),i=1,nctype)
      write(ounit,'(''fzzz'',t11,(12i5))') (nfzzz(i),i=1,nctype)
      write(ounit,'(''gxxxx'',t11,(12i5))') (ngxxxx(i),i=1,nctype)
      write(ounit,'(''gxxxy'',t11,(12i5))') (ngxxxy(i),i=1,nctype)
      write(ounit,'(''gxxxz'',t11,(12i5))') (ngxxxz(i),i=1,nctype)
      write(ounit,'(''gxxyy'',t11,(12i5))') (ngxxyy(i),i=1,nctype)
      write(ounit,'(''gxxyz'',t11,(12i5))') (ngxxyz(i),i=1,nctype)
      write(ounit,'(''gxxzz'',t11,(12i5))') (ngxxzz(i),i=1,nctype)
      write(ounit,'(''gxyyy'',t11,(12i5))') (ngxyyy(i),i=1,nctype)
      write(ounit,'(''gxyyz'',t11,(12i5))') (ngxyyz(i),i=1,nctype)
      write(ounit,'(''gxyzz'',t11,(12i5))') (ngxyzz(i),i=1,nctype)
      write(ounit,'(''gxzzz'',t11,(12i5))') (ngxzzz(i),i=1,nctype)
      write(ounit,'(''gyyyy'',t11,(12i5))') (ngyyyy(i),i=1,nctype)
      write(ounit,'(''gyyyz'',t11,(12i5))') (ngyyyz(i),i=1,nctype)
      write(ounit,'(''gyyzz'',t11,(12i5))') (ngyyzz(i),i=1,nctype)
      write(ounit,'(''gyzzz'',t11,(12i5))') (ngyzzz(i),i=1,nctype)
      write(ounit,'(''gzzzz'',t11,(12i5))') (ngzzzz(i),i=1,nctype)
      write(ounit,'(/,''charge'',t12,(12f5.0))')(znuc(i),i=1,nctype)
      write(ounit,*)

      if(newghostype.gt.0) then
        write(ounit,'(/,''ghost type'',(12i4))') (i,i=nctype+1,nctype+newghostype)
        write(ounit,'(/,''s'',t11,(12i5))') (ns(i),i=nctype+1,nctype+newghostype)
        write(ounit,'(''px'',t11,(12i5))')  (npx(i),i=nctype+1,nctype+newghostype)
        write(ounit,'(''py'',t11,(12i5))')  (npx(i),i=nctype+1,nctype+newghostype)
        write(ounit,'(''pz'',t11,(12i5))')  (npy(i),i=nctype+1,nctype+newghostype)
        write(ounit,'(''dxx'',t11,(12i5))') (ndxx(i),i=nctype+1,nctype+newghostype)
        write(ounit,'(''dxy'',t11,(12i5))') (ndxy(i),i=nctype+1,nctype+newghostype)
        write(ounit,'(''dxz'',t11,(12i5))') (ndxz(i),i=nctype+1,nctype+newghostype)
        write(ounit,'(''dyy'',t11,(12i5))') (ndyy(i),i=nctype+1,nctype+newghostype)
        write(ounit,'(''dyz'',t11,(12i5))') (ndyz(i),i=nctype+1,nctype+newghostype)
        write(ounit,'(''dzz'',t11,(12i5))') (ndzz(i),i=nctype+1,nctype+newghostype)
        write(ounit,'(''fxxx'',t11,(12i5))') (nfxxx(i),i=nctype+1,nctype+newghostype)
        write(ounit,'(''fxxy'',t11,(12i5))') (nfxxy(i),i=nctype+1,nctype+newghostype)
        write(ounit,'(''fxxz'',t11,(12i5))') (nfxxz(i),i=nctype+1,nctype+newghostype)
        write(ounit,'(''fxyy'',t11,(12i5))') (nfxyy(i),i=nctype+1,nctype+newghostype)
        write(ounit,'(''fxyz'',t11,(12i5))') (nfxyz(i),i=nctype+1,nctype+newghostype)
        write(ounit,'(''fxzz'',t11,(12i5))') (nfxzz(i),i=nctype+1,nctype+newghostype)
        write(ounit,'(''fyyy'',t11,(12i5))') (nfyyy(i),i=nctype+1,nctype+newghostype)
        write(ounit,'(''fyyz'',t11,(12i5))') (nfyyz(i),i=nctype+1,nctype+newghostype)
        write(ounit,'(''fyzz'',t11,(12i5))') (nfyzz(i),i=nctype+1,nctype+newghostype)
        write(ounit,'(''fzzz'',t11,(12i5))') (nfzzz(i),i=nctype+1,nctype+newghostype)
        write(ounit,'(''gxxxx'',t11,(12i5))') (ngxxxx(i),i=nctype+1,nctype+newghostype)
        write(ounit,'(''gxxxy'',t11,(12i5))') (ngxxxy(i),i=nctype+1,nctype+newghostype)
        write(ounit,'(''gxxxz'',t11,(12i5))') (ngxxxz(i),i=nctype+1,nctype+newghostype)
        write(ounit,'(''gxxyy'',t11,(12i5))') (ngxxyy(i),i=nctype+1,nctype+newghostype)
        write(ounit,'(''gxxyz'',t11,(12i5))') (ngxxyz(i),i=nctype+1,nctype+newghostype)
        write(ounit,'(''gxxzz'',t11,(12i5))') (ngxxzz(i),i=nctype+1,nctype+newghostype)
        write(ounit,'(''gxyyy'',t11,(12i5))') (ngxyyy(i),i=nctype+1,nctype+newghostype)
        write(ounit,'(''gxyyz'',t11,(12i5))') (ngxyyz(i),i=nctype+1,nctype+newghostype)
        write(ounit,'(''gxyzz'',t11,(12i5))') (ngxyzz(i),i=nctype+1,nctype+newghostype)
        write(ounit,'(''gxzzz'',t11,(12i5))') (ngxzzz(i),i=nctype+1,nctype+newghostype)
        write(ounit,'(''gyyyy'',t11,(12i5))') (ngyyyy(i),i=nctype+1,nctype+newghostype)
        write(ounit,'(''gyyyz'',t11,(12i5))') (ngyyyz(i),i=nctype+1,nctype+newghostype)
        write(ounit,'(''gyyzz'',t11,(12i5))') (ngyyzz(i),i=nctype+1,nctype+newghostype)
        write(ounit,'(''gyzzz'',t11,(12i5))') (ngyzzz(i),i=nctype+1,nctype+newghostype)
        write(ounit,'(''gzzzz'',t11,(12i5))') (ngzzzz(i),i=nctype+1,nctype+newghostype)
        write(ounit,'(/,''charge'',t12,(12f5.0))')(znuc(i),i=nctype+1,nctype+newghostype)
        write(ounit,*)
      endif

      write(ounit,'(/,''no. of orbitals ='',t30,i10)') norb

c the following sets up strings to identify basis functions on centers for
c the printout of coefficients and screening constants.


  210 continue

      iu=ounit
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

        do k=1,ns(iwctype(ic))
         write (iu, 213) j, iwctype(ic), ic, 's', (coef(j,l,1), l=i, imax)
         j=j+1
        enddo


        do k=1,npx(iwctype(ic))
         write (iu, 213) j, iwctype(ic), ic, 'px', (coef(j,l,1), l=i, imax)
         j=j+1
        enddo

        do k=1,npy(iwctype(ic))
         write (iu, 213) j, iwctype(ic), ic, 'py', (coef(j,l,1), l=i, imax)
         j=j+1
        enddo

        do k=1,npz(iwctype(ic))
         write (iu, 213) j, iwctype(ic), ic, 'pz', (coef(j,l,1), l=i, imax)
         j=j+1
        enddo

        do k=1,ndxx(iwctype(ic))
         write (iu, 213) j, iwctype(ic), ic, 'dxx', (coef(j,l,1), l=i, imax)
         j=j+1
        enddo

        do k=1,ndxy(iwctype(ic))
         write (iu, 213) j, iwctype(ic), ic, 'dxy', (coef(j,l,1), l=i, imax)
         j=j+1
        enddo

        do k=1,ndxz(iwctype(ic))
         write (iu, 213) j, iwctype(ic), ic, 'dxz', (coef(j,l,1), l=i, imax)
         j=j+1
        enddo

        do k=1,ndyy(iwctype(ic))
         write (iu, 213) j, iwctype(ic), ic, 'dyy', (coef(j,l,1), l=i, imax)
         j=j+1
        enddo

        do k=1,ndyz(iwctype(ic))
         write (iu, 213) j, iwctype(ic), ic, 'dyz', (coef(j,l,1), l=i, imax)
         j=j+1
        enddo

        do k=1,ndzz(iwctype(ic))
         write (iu, 213) j, iwctype(ic), ic, 'dzz', (coef(j,l,1), l=i, imax)
         j=j+1
        enddo

        do k=1,nfxxx(iwctype(ic))
         write (iu, 213) j, iwctype(ic), ic, 'fxxx', (coef(j,l,1), l=i, imax)
         j=j+1
        enddo

        do k=1,nfxxy(iwctype(ic))
         write (iu, 213) j, iwctype(ic), ic, 'fxxy', (coef(j,l,1), l=i, imax)
         j=j+1
        enddo

        do k=1,nfxxz(iwctype(ic))
         write (iu, 213) j, iwctype(ic), ic, 'fxxz', (coef(j,l,1), l=i, imax)
         j=j+1
        enddo

        do k=1,nfxyy(iwctype(ic))
         write (iu, 213) j, iwctype(ic), ic, 'fxyy', (coef(j,l,1), l=i, imax)
         j=j+1
        enddo

        do k=1,nfxyz(iwctype(ic))
         write (iu, 213) j, iwctype(ic), ic, 'fxyz', (coef(j,l,1), l=i, imax)
         j=j+1
        enddo

        do k=1,nfxzz(iwctype(ic))
         write (iu, 213) j, iwctype(ic), ic, 'fxzz', (coef(j,l,1), l=i, imax)
         j=j+1
        enddo

        do k=1,nfyyy(iwctype(ic))
         write (iu, 213) j, iwctype(ic), ic, 'fyyy', (coef(j,l,1), l=i, imax)
         j=j+1
        enddo

        do k=1,nfyyz(iwctype(ic))
         write (iu, 213) j, iwctype(ic), ic, 'fyyz', (coef(j,l,1), l=i, imax)
         j=j+1
        enddo

        do k=1,nfyzz(iwctype(ic))
         write (iu, 213) j, iwctype(ic), ic, 'fyzz', (coef(j,l,1), l=i, imax)
         j=j+1
        enddo

        do k=1,nfzzz(iwctype(ic))
         write (iu, 213) j, iwctype(ic), ic, 'fzzz', (coef(j,l,1), l=i, imax)
         j=j+1
        enddo

        do k=1,ngxxxx(iwctype(ic))
         write (iu, 213) j, iwctype(ic), ic, 'gxxxx', (coef(j,l,1), l=i, imax)
         j=j+1
        enddo

        do k=1,ngxxxy(iwctype(ic))
         write (iu, 213) j, iwctype(ic), ic, 'gxxxy', (coef(j,l,1), l=i, imax)
         j=j+1
        enddo

        do k=1,ngxxxz(iwctype(ic))
         write (iu, 213) j, iwctype(ic), ic, 'gxxxz', (coef(j,l,1), l=i, imax)
         j=j+1
        enddo

        do k=1,ngxxyy(iwctype(ic))
         write (iu, 213) j, iwctype(ic), ic, 'gxxyy', (coef(j,l,1), l=i, imax)
         j=j+1
        enddo

        do k=1,ngxxyz(iwctype(ic))
         write (iu, 213) j, iwctype(ic), ic, 'gxxyz', (coef(j,l,1), l=i, imax)
         j=j+1
        enddo

        do k=1,ngxxzz(iwctype(ic))
         write (iu, 213) j, iwctype(ic), ic, 'gxxzz', (coef(j,l,1), l=i, imax)
         j=j+1
        enddo

        do k=1,ngxyyy(iwctype(ic))
         write (iu, 213) j, iwctype(ic), ic, 'gxyyy', (coef(j,l,1), l=i, imax)
         j=j+1
        enddo

        do k=1,ngxyyz(iwctype(ic))
         write (iu, 213) j, iwctype(ic), ic, 'gxyyz', (coef(j,l,1), l=i, imax)
         j=j+1
        enddo

        do k=1,ngxyzz(iwctype(ic))
         write (iu, 213) j, iwctype(ic), ic, 'gxyzz', (coef(j,l,1), l=i, imax)
         j=j+1
        enddo

        do k=1,ngxzzz(iwctype(ic))
         write (iu, 213) j, iwctype(ic), ic, 'gxzzz', (coef(j,l,1), l=i, imax)
         j=j+1
        enddo

        do k=1,ngyyyy(iwctype(ic))
         write (iu, 213) j, iwctype(ic), ic, 'gyyyy', (coef(j,l,1), l=i, imax)
         j=j+1
        enddo

        do k=1,ngyyyz(iwctype(ic))
         write (iu, 213) j, iwctype(ic), ic, 'gyyyz', (coef(j,l,1), l=i, imax)
         j=j+1
        enddo

        do k=1,ngyyzz(iwctype(ic))
         write (iu, 213) j, iwctype(ic), ic, 'gyyzz', (coef(j,l,1), l=i, imax)
         j=j+1
        enddo

        do k=1,ngyzzz(iwctype(ic))
         write (iu, 213) j, iwctype(ic), ic, 'gyzzz', (coef(j,l,1), l=i, imax)
         j=j+1
        enddo

        do k=1,ngzzzz(iwctype(ic))
         write (iu, 213) j, iwctype(ic), ic, 'gzzzz', (coef(j,l,1), l=i, imax)
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
