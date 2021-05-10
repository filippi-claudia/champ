      subroutine read_bas_num(iwf)
      use numbas_mod, only: MRWF, MRWF_PTS
      use vmc_mod, only: MBASIS, MCTYPE
      use vmc_mod, only: NCOEF
      use atom, only: znuc, nctype, nctype_tot
! c Written by Claudia Filippi
! c Modified by F. Schautz to use fancy file names
! c Reads in localized orbitals on a radial grid

      use numbas_mod, only: MRWF, MRWF_PTS
      use vmc_mod, only: MBASIS, MCTYPE
      use vmc_mod, only: NCOEF
      use atom, only: znuc, nctype
      use ghostatom, only: newghostype
      use const, only: ipr
      use numbas, only: arg, d2rwf, igrid, nr, nrbas, r0, rwf
      use numbas, only: allocate_numbas
      use coefs, only: nbasis
      use numexp, only: ae, ce, allocate_numexp
      use pseudo, only: nloc
      use general, only: filename, filenames_bas_num, wforce

      use atom, 			        only: atomtyp   
      use general, 			      only: pooldir, bas_id
      use contrl_file,        only: ounit, errunit
      use precision_kinds,    only: dp

      implicit none

      integer         :: ic, ir, irb, ii, jj, ll, icoef, iff
      integer         ::  iwf, info
      real(dp)        :: val, dwf1, wfm, dwfn, dwfm
      integer         :: iunit, iostat, counter = 0
      logical         :: exist, skip = .true.

      real(dp), dimension(MRWF_PTS)       ::  x, work
      real(dp), dimension(ncoef)          ::  y
      real(dp), dimension(ncoef*ncoef)    ::  dmatr
      real(dp), dimension(nbasis)         ::  l
      integer, dimension(ncoef)           :: ipiv 
      integer, dimension(nctype)          :: icusp

! c nrbas = number of numerical orbitals for each center
! c igrid = 1 linear r(i+1)=arg+r(i), r(1)=r0
! c         2 exponential r(i+1)=arg*r(i), r(1)=r0
! c         3 shifted exponential r(i+1)=r0*(arg**(i-1)-1)

      if (.not. allocated(filenames_bas_num)) allocate(filenames_bas_num(nctype))

      call allocate_numbas()
      call allocate_numexp()


      do 100 ic=1,nctype+newghostype
        if (ic .gt. 999) call fatal_error('READ_BAS_NUM: atomtyp > 999')
        filename =  trim(pooldir) //  trim(bas_id) // ".basis." // atomtyp(ic)


        inquire(file=filename, exist=exist)
        if (exist) then
          open (newunit=iunit,file=filename, iostat=iostat, action='read', status='old')
          if (iostat .ne. 0) error stop "Problem in opening the numerical basis file"
        else
          error stop " Numerical Basis file "// filename // " does not exist."
        endif

        read(iunit,*, iostat=iostat) nrbas(ic),igrid(ic),nr(ic),arg(ic),r0(ic),icusp(ic)
        write(ounit,*) "reading the content ", nrbas(ic),igrid(ic),nr(ic),arg(ic),r0(ic),icusp(ic)        
        write(ounit,*) "(Reading basis grid file = [ ",  trim(filename), " ] )"
        write(ounit,'(''center type '',i4,'' nrbas,igrid,nr,arg,r0 ='',2i4,i5,2f10.5)') &
        ic,nrbas(ic),igrid(ic),nr(ic),arg(ic),r0(ic)

        if(nrbas(ic).gt.MRWF) call fatal_error('READ_BAS_NUM: nrbas gt MRWF')
        if(nr(ic).gt.MRWF_PTS) call fatal_error('READ_BAS_NUM: nr gt MRWF_PTS')
        if(igrid(ic).ne.1.and.igrid(ic).ne.2.and.igrid(ic).ne.3) &
        call fatal_error('READ_BAS_NUM: grid not implemented')

        if(nloc.eq.0) read(iunit,*) (l(irb),irb=1,nrbas(ic))

        do ir=1,nr(ic)
          read(iunit,*,iostat=iostat) x(ir),(rwf(ir,irb,ic,iwf),irb=1,nrbas(ic))
        enddo


        if(igrid(ic).eq.2.and.arg(ic).le.1.d0) arg(ic)=x(2)/x(1)
        if(igrid(ic).eq.3) r0(ic)=r0(ic)/(arg(ic)**(nr(ic)-1)-1.d0)



        do 100 irb=1,nrbas(ic)

        if(nloc.eq.0.and.l(irb).eq.0.and.icusp(ic).eq.1) then

! c small radii wf(r)=ce1-znuc*ce1*r+ce3*r**2+ce4*r**3+ce5*r**4
          do 15 ii=1,NCOEF-1
  15        dmatr(ii)=1.d0-znuc(ic)*x(ii)
          y(1)=rwf(1,irb,ic,iwf)
          ll=NCOEF-1
          do 16 jj=2,NCOEF-1
            y(jj)=rwf(jj,irb,ic,iwf)
            do 16 ii=2,NCOEF-1
              ll=ll+1
  16          dmatr(ll)=x(ii)**jj

          call dgesv(NCOEF-1,1,dmatr,NCOEF-1,ipiv,y,NCOEF,info)
          ce(1,irb,ic,iwf)=y(1)
          ce(2,irb,ic,iwf)=-znuc(ic)*ce(1,irb,ic,iwf)
          ce(3,irb,ic,iwf)=y(2)
          ce(4,irb,ic,iwf)=y(3)
          ce(5,irb,ic,iwf)=y(4)

         else

! c small radii wf(r)=ce1+ce2*r+ce3*r**2+ce4*r**3+ce5*r**4
          ll=0
          do 25 jj=1,NCOEF
            y(jj)=rwf(jj,irb,ic,iwf)
            do 25 ii=1,NCOEF
              ll=ll+1
  25          dmatr(ll)=x(ii)**(jj-1)
          call dgesv(NCOEF,1,dmatr,NCOEF,ipiv,y,NCOEF,info)

          do 26 icoef=1,NCOEF
  26        ce(icoef,irb,ic,iwf)=y(icoef)

        endif

! c       if(ipr.gt.1) then
          write(ounit,'(''basis = '',i4)') irb
          write(ounit,'(''check the small radius expansion'')')
          write(ounit,'(''coefficients'',1p10e22.10)') &
                    (ce(iff,irb,ic,iwf),iff=1,NCOEF)
          write(ounit,'(''check the small radius expansion'')')
          write(ounit,'(''irad, rad, extrapolated value, correct value'')')
          do 30 ir=1,10
            val=ce(1,irb,ic,iwf)
            do 28 icoef=2,NCOEF
  28        val=val+ce(icoef,irb,ic,iwf)*x(ir)**(icoef-1)
  30        write(ounit,'(i2,1p3e22.14)')ir,x(ir),val,rwf(ir,irb,ic,iwf)
! c       endif

        dwf1=0.d0
        do 32 icoef=2,NCOEF
  32      dwf1=dwf1+(icoef-1)*ce(icoef,irb,ic,iwf)*x(1)**(icoef-2)

! c large radii wf(r)=a0*exp(-ak*r)
! c       xm=0.5d0*(x(nr(ic))+x(nr(ic)-1))
        wfm=0.5d0*(rwf(nr(ic),irb,ic,iwf)+rwf(nr(ic)-1,irb,ic,iwf))
        dwfm=(rwf(nr(ic),irb,ic,iwf)-rwf(nr(ic)-1,irb,ic,iwf))/  &
        (x(nr(ic))-x(nr(ic)-1))
        if(dabs(wfm).gt.1.d-99) then
          ae(2,irb,ic,iwf)=-dwfm/wfm
          ae(1,irb,ic,iwf)=rwf(nr(ic),irb,ic,iwf)*    &
                           dexp(ae(2,irb,ic,iwf)*x(nr(ic)))
          dwfn=-ae(2,irb,ic,iwf)*rwf(nr(ic),irb,ic,iwf)
         else
          ae(1,irb,ic,iwf)=0.d0
          ae(2,irb,ic,iwf)=0.d0
          dwfn=0.d0
        endif
! c       if(ipr.gt.1) then
          write(ounit,'(''check the large radius expansion'')')
          write(ounit,'(''a0,ak'',1p2e22.10)')     &
                            ae(1,irb,ic,iwf),ae(2,irb,ic,iwf)
          write(ounit,'(''irad, rad, extrapolated value, correct value'')')
          do 40 ir=1,10
            val=ae(1,irb,ic,iwf)*dexp(-ae(2,irb,ic,iwf)*x(nr(ic)-ir))
  40        write(ounit,'(i2,1p3e22.14)')      &
            ir,x(nr(ic)-ir),val,rwf(nr(ic)-ir,irb,ic,iwf)
          write(ounit,*) 'dwf1,dwfn',dwf1,dwfn
! c       endif
        print*, "ae(2,irb,ic,iwf) ", irb, ic, iwf, ae(2,irb,ic,iwf)
        if(ae(2,irb,ic,iwf).lt.0) call fatal_error ('BASIS_READ_NUM: ak<0')

        call spline2(x,rwf(1,irb,ic,iwf),nr(ic),dwf1,dwfn,  &
        d2rwf(1,irb,ic,iwf),work)

      close(iunit)
 100  continue


      return
      end subroutine read_bas_num
