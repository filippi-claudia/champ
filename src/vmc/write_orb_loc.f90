      module write_orb_loc_mod
      contains
      subroutine write_orb_loc
! Written by Cyrus Umrigar and Claudia Filippi, starting from Kevin Schmidt's routine
! Reads in localized orbitals, in either
! Modified by A. Scemama (printing in a GAMESS-like format)
! 1) a slater basis
! 2) a gaussian basis
      use system, only: znuc, iwctype, nctype, ncent, newghostype, nelec
      use numbas, only: numr
      use coefs, only: nbasis
      use slater, only: coef
      use basis, only: zex, betaq
      use contrl_file,    only: ounit, errunit
      use error, only: fatal_error

      use precision_kinds, only: dp
      use slater,  only: coef
      use system,  only: iwctype,ncent,nctype,nelec,newghostype,znuc
      use slater, only: norb

      implicit none
      integer :: iu
      integer :: i, iabs, ic, imax
      integer :: j, k, l, ncheck
      integer, parameter :: nprime = 10
      real(dp), parameter :: zero = 0.d0
      real(dp), parameter :: one = 1.d0




! Exponent for asymptotic basis
      betaq=zero
      do ic=1,ncent
        betaq=betaq+znuc(iwctype(ic))
      enddo
      betaq=betaq-nelec+one

      write(ounit,'(/,''center type'',(12i4))') (i,i=1,nctype)

      write(ounit,'(/,''charge'',t12,(12f5.0))')(znuc(i),i=1,nctype)
      write(ounit,*)

      if(newghostype.gt.0) then
        write(ounit,'(/,''ghost type'',(12i4))') (i,i=nctype+1,nctype+newghostype)
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

      if(numr.eq.0) then
        write(iu,'(''screening constants'')')
        write(iu,'(12f10.6)') (zex(i,1),i=1,nbasis)
      endif

      return
      end
      end module
