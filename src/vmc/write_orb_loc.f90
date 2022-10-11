      module write_orb_loc_mod
      contains
      subroutine write_orb_loc
! Written by Cyrus Umrigar and Claudia Filippi, starting from Kevin Schmidt's routine
! Reads in localized orbitals, in either
! Modified by A. Scemama (printing in a GAMESS-like format)
! 1) a slater basis
! 2) a gaussian basis
      use atom, only: znuc, iwctype, nctype, ncent
      use ghostatom, only: newghostype
      use const, only: nelec
      use numbas, only: numr
      use coefs, only: coef, nbasis, norb
      use basis, only: zex, betaq
      use basis, only: ns, np, nd, nf, ng
      use contrl_file,    only: ounit, errunit
      use error, only: fatal_error

      use precision_kinds, only: dp
      implicit none
      integer :: iu
      integer :: i, iabs, ic, imax
      integer :: j, k, l, m, ncheck
      character(len=2), dimension(3)  :: ao_label_p
      character(len=3), dimension(6)  :: ao_label_d
      character(len=4), dimension(10) :: ao_label_f
      character(len=5), dimension(15) :: ao_label_g

      real(dp), parameter :: zero = 0.d0
      real(dp), parameter :: one = 1.d0
      ao_label_p = (/ 'px', 'py', 'pz' /)
      ao_label_d = (/ 'dxx', 'dxy', 'dxz', 'dyy', 'dyz', 'dzz' /)
      ao_label_f = (/ 'fxxx', 'fxxy', 'fxxz', 'fxyy', 'fxyz', &
                      'fxzz', 'fyyy', 'fyyz', 'fyzz', 'fzzz' /)
      ao_label_g = (/ 'gxxxx', 'gxxxy', 'gxxxz', 'gxxyy', 'gxxyz', &
                      'gxxzz', 'gxyyy', 'gxyyz', 'gxyzz', 'gxzzz', &
                      'gyyyy', 'gyyyz', 'gyyzz', 'gyzzz', 'gzzzz' /)



      ! Check that nbasis in lcao matches specified basis on all centers
      ncheck=0
      do ic=1,ncent
        i = iwctype(ic)
        ncheck = ncheck + ns(i) + 3*np(i) + 6*nd(i) + 10*nf(i) + 15*ng(i)
      enddo


      if(ncheck.ne.nbasis) then
        write(ounit,'(''ncheck = '',i6)') ncheck
        write(ounit,'(''nbasis = '',i6)') nbasis
        call fatal_error('WRITE_BAS: ncheck not equal to nbasis')
      endif

! Exponent for asymptotic basis
      betaq=zero
      do ic=1,ncent
        betaq=betaq+znuc(iwctype(ic))
      enddo
      betaq=betaq-nelec+one

      write(ounit,'(/,''center type'',(12i4))') (i,i=1,nctype)
      write(ounit,'(/,''s'',t11,(12i5))') (ns(i),i=1,nctype)
      write(ounit,'(''p'',t11,(12i5))') (np(i),i=1,nctype)
      write(ounit,'(''d'',t11,(12i5))') (nd(i),i=1,nctype)
      write(ounit,'(''f'',t11,(12i5))') (nf(i),i=1,nctype)
      write(ounit,'(''g'',t11,(12i5))') (ng(i),i=1,nctype)
      write(ounit,'(/,''charge'',t12,(12f5.0))')(znuc(i),i=1,nctype)
      write(ounit,*)

      if(newghostype.gt.0) then
        write(ounit,'(/,''ghost type'',(12i4))') (i,i=nctype+1,nctype+newghostype)
        write(ounit,'(/,''s'',t11,(12i5))') (ns(i),i=nctype+1,nctype+newghostype)
        write(ounit,'(''p'',t11,(12i5))')  (np(i),i=nctype+1,nctype+newghostype)
        write(ounit,'(''d'',t11,(12i5))')  (nd(i),i=nctype+1,nctype+newghostype)
        write(ounit,'(''f'',t11,(12i5))')  (nf(i),i=nctype+1,nctype+newghostype)
        write(ounit,'(''g'',t11,(12i5))')  (ng(i),i=nctype+1,nctype+newghostype)
        write(ounit,'(/,''charge'',t12,(12f5.0))')(znuc(i),i=nctype+1,nctype+newghostype)
        write(ounit,*)
      endif

      write(ounit,'(/,''no. of orbitals ='',t30,i10)') norb

! the following sets up strings to identify basis functions on centers for
! the printout of coefficients and screening constants.


  210 continue

      iu=ounit
      if(nbasis.gt.15) then
        iu=45
        write(iu,'(/,''Orbitals'')')
      endif

      i=1
      do while (i.le.norb)
        imax=i+9
        if (i+9.gt.norb) then
          imax=norb
        endif
        write (iu, 211) (l, l=i,imax)
        ! TODO print symmetry

        j=1
        do ic=1,ncent

          do k=1,ns(iwctype(ic))
          write (iu, 213) j, iwctype(ic), ic, 's', (coef(j,l,1), l=i, imax)
          j=j+1
          enddo


          do k=1,np(iwctype(ic))
            do m = 1, 3           ! loop over px, py, pz
              write (iu, 213) j, iwctype(ic), ic, ao_label_p(m), (coef(j,l,1), l=i, imax)
              j=j+1
            enddo
          enddo

          do k=1,nd(iwctype(ic))
            do m = 1, 6           ! loop over dxx, dxy ...
              write (iu, 213) j, iwctype(ic), ic, ao_label_d(m), (coef(j,l,1), l=i, imax)
              j=j+1
            enddo
          enddo


          do k=1,nf(iwctype(ic))
            do m = 1, 10          ! loop over fxxx, fxxy ...
              write (iu, 213) j, iwctype(ic), ic, ao_label_f(m), (coef(j,l,1), l=i, imax)
              j=j+1
            enddo
          enddo

          do k=1,ng(iwctype(ic))
            do m = 1, 10          ! loop over gxxxx, gxxxy ...
              write (iu, 213) j, iwctype(ic), ic, ao_label_g(m), (coef(j,l,1), l=i, imax)
              j=j+1
            enddo
          enddo

        enddo ! loop on cent

        i=i+5
        write (iu, *)
      enddo ! do while (i.le.norb)

  211 format (5X,8X,10(1X,I10))
  212 format (5X,8X,10(1X,A10))
  213 format (I5,I3,I3,A5,10(1X,F10.6))

      if(numr.eq.0) then
        write(iu,'(''screening constants'')')
        write(iu,'(12f10.6)') (zex(i,1),i=1,nbasis)
      endif

      return
      end
      end module
