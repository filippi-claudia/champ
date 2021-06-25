      subroutine read_bas_num(iwf)
      use numbas_mod, only: MRWF, MRWF_PTS
      use vmc_mod, only: MBASIS, MCTYPE
      use vmc_mod, only: NCOEF
      use atom, only: znuc, nctype, nctype_tot
c Written by Claudia Filippi
c Modified by F. Schautz to use fancy file names
c Reads in localized orbitals on a radial grid

      use numbas_mod, only: MRWF, MRWF_PTS
      use vmc_mod, only: MBASIS, MCTYPE
      use vmc_mod, only: NCOEF
      use atom, only: znuc, nctype
      use ghostatom, only: newghostype
      use const, only: ipr
      use numbas, only: arg, d2rwf, igrid, nr, nrbas, r0, rwf
      use coefs, only: nbasis
      use numexp, only: ae, ce
      use pseudo, only: nloc
      use general, only: filename, filenames_bas_num, wforce

      use precision_kinds, only: dp
      implicit none

      integer :: ic, icoef, iff, ii, index
      integer :: info, ir, irb, iwf
      integer :: j, jj, ll
      integer, dimension(NCOEF) :: ipiv
      integer, dimension(nbasis) :: l
      integer, dimension(nctype_tot) :: icusp
      real(dp) :: c110, c120, c130, dabs, dwf1
      real(dp) :: dwfm, dwfn, val, wfm
      real(dp), dimension(MRWF_PTS) :: x
      real(dp), dimension(MRWF_PTS) :: work
      real(dp), dimension(NCOEF) :: y
      real(dp), dimension(NCOEF*NCOEF) :: dmatr


c nrbas = number of numerical orbitals for each center
c igrid = 1 linear r(i+1)=arg+r(i), r(1)=r0
c         2 exponential r(i+1)=arg*r(i), r(1)=r0
c         3 shifted exponential r(i+1)=r0*(arg**(i-1)-1)

      do 100 ic=1,nctype+newghostype
        if (ic .gt. 999) call fatal_error('READ_BAS_NUM: atomtyp > 999')
        filename=filenames_bas_num(ic)
        open(21,file=filename(1:index(filename,' ')-1),status='old')

        read(21,*) nrbas(ic),igrid(ic),nr(ic),arg(ic),r0(ic),icusp(ic)
        write(45,'(''Reading basis grid file = ['',a,'']'')') filename(1:index(filename,' ')-1)
        write(45,'(''center type '',i4,'' nrbas,igrid,nr,arg,r0 ='',2i4,i5,2f10.5)') 
     &  ic,nrbas(ic),igrid(ic),nr(ic),arg(ic),r0(ic)

        if(nrbas(ic).gt.MRWF) call fatal_error('READ_BAS_NUM: nrbas gt MRWF')
        if(nr(ic).gt.MRWF_PTS) call fatal_error('READ_BAS_NUM: nr gt MRWF_PTS')
        if(igrid(ic).ne.1.and.igrid(ic).ne.2.and.igrid(ic).ne.3)
     &  call fatal_error('READ_BAS_NUM: grid not implemented')

        if(nloc.eq.0) read(21,*) (l(irb),irb=1,nrbas(ic))

        do 10 ir=1,nr(ic)
          read(21,*) x(ir),(rwf(ir,irb,ic,iwf),irb=1,nrbas(ic))
   10     continue

        if(igrid(ic).eq.2.and.arg(ic).le.1.d0) arg(ic)=x(2)/x(1)
        if(igrid(ic).eq.3) r0(ic)=r0(ic)/(arg(ic)**(nr(ic)-1)-1.d0)

        do 100 irb=1,nrbas(ic)

        if(nloc.eq.0.and.l(irb).eq.0.and.icusp(ic).eq.1) then

c small radii wf(r)=ce1-znuc*ce1*r+ce3*r**2+ce4*r**3+ce5*r**4
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

c small radii wf(r)=ce1+ce2*r+ce3*r**2+ce4*r**3+ce5*r**4
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

c       if(ipr.gt.1) then
          write(45,'(''basis = '',i4)') irb
          write(45,'(''check the small radius expansion'')')
          write(45,'(''coefficients'',1p10e22.10)')
     &               (ce(iff,irb,ic,iwf),iff=1,NCOEF)
          write(45,'(''check the small radius expansion'')')
          write(45,'(''irad, rad, extrapolated value, correct value'')')
          do 30 ir=1,10
            val=ce(1,irb,ic,iwf)
            do 28 icoef=2,NCOEF
  28        val=val+ce(icoef,irb,ic,iwf)*x(ir)**(icoef-1)
  30        write(45,'(i2,1p3e22.14)')ir,x(ir),val,rwf(ir,irb,ic,iwf)
c       endif

        dwf1=0.d0
        do 32 icoef=2,NCOEF
  32      dwf1=dwf1+(icoef-1)*ce(icoef,irb,ic,iwf)*x(1)**(icoef-2)

c large radii wf(r)=a0*exp(-ak*r)
c       xm=0.5d0*(x(nr(ic))+x(nr(ic)-1))
        wfm=0.5d0*(rwf(nr(ic),irb,ic,iwf)+rwf(nr(ic)-1,irb,ic,iwf))
        dwfm=(rwf(nr(ic),irb,ic,iwf)-rwf(nr(ic)-1,irb,ic,iwf))/
     &  (x(nr(ic))-x(nr(ic)-1))
        if(dabs(wfm).gt.1.d-99) then
          ae(2,irb,ic,iwf)=-dwfm/wfm
          ae(1,irb,ic,iwf)=rwf(nr(ic),irb,ic,iwf)*
     &                     dexp(ae(2,irb,ic,iwf)*x(nr(ic)))
          dwfn=-ae(2,irb,ic,iwf)*rwf(nr(ic),irb,ic,iwf)
         else
          ae(1,irb,ic,iwf)=0.d0
          ae(2,irb,ic,iwf)=0.d0
          dwfn=0.d0
        endif
c       if(ipr.gt.1) then
          write(45,'(''check the large radius expansion'')')
          write(45,'(''a0,ak'',1p2e22.10)')
     &                       ae(1,irb,ic,iwf),ae(2,irb,ic,iwf)
          write(45,'(''irad, rad, extrapolated value, correct value'')')
          do 40 ir=1,10
            val=ae(1,irb,ic,iwf)*dexp(-ae(2,irb,ic,iwf)*x(nr(ic)-ir))
  40        write(45,'(i2,1p3e22.14)')
     &      ir,x(nr(ic)-ir),val,rwf(nr(ic)-ir,irb,ic,iwf)
          write(45,*) 'dwf1,dwfn',dwf1,dwfn
c       endif
        if(ae(2,irb,ic,iwf).lt.0) call fatal_error ('BASIS_READ_NUM: ak<0')

        call spline2(x,rwf(1,irb,ic,iwf),nr(ic),dwf1,dwfn,
     &  d2rwf(1,irb,ic,iwf),work)

      close(21)
 100  continue

c TEMPORARY debug
c     do 130 jwf=1,nforce
c       do 130 ic=1,nctype
c         if(ic.lt.10) then
c           write(lcent,'(i1)') ic
c          elseif(ic.lt.100) then
c           write(lcent,'(i2)') ic
c         else
c           call fatal_error ('problem with spline.test')
c         endif
c         if(jwf.lt.10) then
c           write(wforce,'(i1)') jwf
c          elseif(jwf.lt.100) then
c           write(wforce,'(i2)') jwf
c         endif
c         filename='spline.chk.'//lcent(1:index(lcent,' ')-1)//wforce
c         open(22,file=filename,status='unknown')
c         do 120 ir=1,nr(ic)-1
c           ii=1
c           do 110 irb=1,nrbas(ic)
c             call splfit(x(ir),irb,ic,jwf,work(ii),1)
c110          ii=ii+3
c120        write(22,'(1p40e20.12)') x(ir),(work(j),j=1,3*nrbas(ic))
c130      close(22)

      return
      end
