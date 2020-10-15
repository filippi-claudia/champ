      subroutine readps_tm
c Written by Claudia Filippi
c read Troullier-Martins pseudopotentials
c reads in r*v in ryd.
c does 3 conversions: a) ryd -> Har, b) r*v -> v and
c c) subtracts out local part from all except highest l component.
c Also eval pot. at 0 and initializes quadrature pts.
c 
c Modified by F. Schautz to use fancy file names
      use pseudo_mod, only: MPS_GRID
      use vmc_mod, only: NCOEF
      use atom, only: znuc, nctype
      use const, only: ipr
      use pseudo_tm, only: arg, d2pot, nr_ps, r0, rmax, vpseudo

      use pseudo, only: lpot, nloc, vps

      use qua, only: nquad, wq, xq, xq0, yq, yq0, zq, zq0

      implicit real*8(a-h,o-z)




      character*2 icorr,nameat
      character*3 irel
      character*4 nicore
      character*10 ititle(7),iray(6)
      character*20 atomtyp,atomsymbol
      character*256 filename,pooldir,pp_id



      dimension r(MPS_GRID),y(NCOEF),ce(NCOEF),dmatr(NCOEF*NCOEF),ipiv(NCOEF)
      dimension work(MPS_GRID)

c pool directory for pseudopotentials
      call p2gtad('general:pool',pooldir,'.',1)
      call stripquotes(pooldir)
      call p2gtad('general:pseudopot',pp_id,'none',1)
      call stripquotes(pp_id)
CVARDOC String to identify pseudopotential. If set, fancy names for 
CVARDOC the pseudopotential files will be used.  

      do 200 ic=1,nctype

      if(ic.lt.10) then
        write(atomtyp,'(i1)') ic
       elseif(ic.lt.100) then
        write(atomtyp,'(i2)') ic
       else
        call fatal_error('READPS_TM: nctype>100')
      endif


      if(pp_id.eq.'none') then
c old naming convention
        if(nloc.eq.2) then
          filename=pooldir(1:index(pooldir,' ')-1)//'/'//
     &             'pseudo.dat.'//atomtyp(1:index(atomtyp,' ')-1)
         elseif(nloc.eq.3) then
          filename=pooldir(1:index(pooldir,' ')-1)//'/'//
     &             'pseudopot'//atomtyp(1:index(atomtyp,' ')-1)
        endif
       else
c new naming convention
        call p2gtad('atom_types:'//atomtyp(1:index(atomtyp,' ')-1)
     &       ,atomsymbol,'X',1)
        if(nloc.eq.2) then
          filename=pooldir(1:index(pooldir,' ')-1)//
     &             '/'//
     &             pp_id(1:index(pp_id,' ')-1)//
     &             '.pseudo.dat.'//
     &             atomsymbol(1:index(atomsymbol,' ')-1)
         elseif(nloc.eq.3) then
          filename=pooldir(1:index(pooldir,' ')-1)//
     &             '/'//
     &             pp_id(1:index(pp_id,' ')-1)//
     &             '.pseudopot.'//
     &             atomsymbol(1:index(atomsymbol,' ')-1)
        endif
      endif

      if(nloc.eq.2) then
        open(1,file=filename(1:index(filename,' ')),status='old',form='unformatted')
       elseif(nloc.eq.3) then
        open(1,file=filename(1:index(filename,' ')),status='old')
      endif
c
      write(45,'(''Reading pseudopotential file '',a)')
     &     filename(1:index(filename,' '))

      if(nloc.eq.2) then
        read(1) nameat,icorr,irel,nicore,(iray(i),i=1,6),
     &  (ititle(i),i=1,7),npotd,npotu,nrm,r0(ic),arg(ic),zion
       elseif(nloc.eq.3) then
        read(1,'(a2,x,a2,x,a3,x,a4,x,13(a10,x))') nameat,icorr,irel,nicore
     &  ,(iray(i),i=1,6),(ititle(i),i=1,7)
        read(1,*) npotd,npotu,nrm,r0(ic),arg(ic),zion
      endif

      nr=nrm+1
      nr_ps(ic)=nr

      if(znuc(ic).ne.zion) then
        write(6,'(''znuc(ic) != zion in readps_tm'',2f6.1)') znuc(ic),zion
        call fatal_error('READPS_TM: znuc(ic) != zion')
      endif
      if(nr.gt.MPS_GRID) then
        write(6,'(''nr > MPS_GRID'',2i6)') nr,MPS_GRID
        call fatal_error('READPS_TM: nr > MPS_GRID')
      endif

      r(1)=0.d0
      if(nloc.eq.2) then
        read(1) (r(i),i=2,nr)
       elseif(nloc.eq.3) then
        read(1,*) (r(i),i=2,nr)
      endif

      rmax(ic)=0.d0
      jmax=0
      do 100 i=1,npotd
        if(nloc.eq.2) then
          read(1) ii,(vpseudo(j,ic,i),j=2,nr)
         elseif(nloc.eq.3) then
          read(1,*) ii,(vpseudo(j,ic,i),j=2,nr)
        endif
        vpseudo(1,ic,i)=0.d0
        do 50 j=nr,1,-1
c         write(6,'(''vps'',f8.4,f12.6,d12.4)') r(j),vpseudo(j,ic,i)/2,vpseudo(j,ic,i)/2+zion
          if(dabs(vpseudo(j,ic,i)+2.d0*zion).gt.1.d-6) then
            rmax(ic)=max(rmax(ic),r(j))
            jmax=max(jmax,j)
            goto 60
          endif
   50   continue
   60   do 100 j=2,nr
          vpseudo(j,ic,i)=0.5d0*vpseudo(j,ic,i)/r(j)
  100 continue
      jmax=jmax+5
      write(45,'(''center '',i3,'' pseudopot rmax,jmax= '',f6.2,i4)')
     & ic,rmax(ic),jmax

      close(1)

      lpot(ic)=npotd
      arg(ic)=dexp(arg(ic))
      write(45,'(''potential grid param: r0, alpha, nr = '',2f10.4,i5)') r0(ic),arg(ic),nr_ps(ic)

      do 110 i=1,npotd-1
        do 110 j=2,nr
          if(ipr.ge.1) then
            if(i.eq.1) write(38,*) r(j),vpseudo(j,ic,i),vpseudo(j,ic,npotd),-znuc(ic)/r(j)
            if(i.eq.2) write(39,*) r(j),vpseudo(j,ic,i),vpseudo(j,ic,npotd),-znuc(ic)/r(j)
          endif
          vpseudo(j,ic,i)=vpseudo(j,ic,i)-vpseudo(j,ic,npotd)
  110 continue

      if(rmax(ic).eq.0.d0) goto 200
      do 190 i=1,npotd
c small radii pot(r)=ce1+ce2*r+ce3*r**2+ce4*r**3+ce5*r**4
        ll=0
        do 120 jj=1,NCOEF
          y(jj)=vpseudo(jj+1,ic,i)
          do 120 ii=1,NCOEF
            ll=ll+1
  120       dmatr(ll)=r(ii+1)**(jj-1)
        call dgesv(NCOEF,1,dmatr,NCOEF,ipiv,y,NCOEF,info)

        do 125 icoef=1,NCOEF
  125     ce(icoef)=y(icoef)

        if(ipr.gt.1) then
c         write(45,'(''coefficients'',1p10e22.10)') (ce(iff),iff=1,NCOEF)
          write(45,'(''check the small radius expansion for l= '',i3)')i-1
          write(45,'(''irad, rad, extrapolated, correct value, diff'')')
          do 130 ir=1,10
            val=ce(1)
            do 128 icoef=2,NCOEF
  128         val=val+ce(icoef)*r(ir)**(icoef-1)
c             if(ir.eq.1) vpseudo(ir,ic,i)=val
  130       write(45,'(i2,1p3e16.8,1p1e11.2)')
     &      ir,r(ir),val,vpseudo(ir,ic,i),val-vpseudo(ir,ic,i)
        endif

        vpseudo(1,ic,i)=ce(1)

        dpot1=ce(2)
        do 135 icoef=3,NCOEF
  135     dpot1=dpot1+(icoef-1)*ce(icoef)*r(1)**(icoef-2)

        dpotn=0.d0
        if(i.eq.npotd) dpotn=zion/r(jmax)**2

        if(ipr.gt.1) write(45,'(''dpot1,dpotn'',1p2e15.5)') dpot1,dpotn

c get second derivative for spline fit
        call spline2(r,vpseudo(1,ic,i),jmax,dpot1,dpotn,
     &  d2pot(1,ic,i),work)

        if(ipr.ge.1) then
          if(i.eq.1) ipot_io=35
          if(i.eq.2) ipot_io=36
          if(i.eq.3) ipot_io=37
          do 140 j=1,jmax
  140       write(ipot_io,*)
c    &      r(j),vpseudo(j,ic,i),d2pot(j,ic,i),-znuc(ic)/r(j),-2*znuc(ic)/r(j)**3
     &      r(j),vpseudo(j,ic,i),d2pot(j,ic,i)
        endif
  190 continue

c Warning: Temporary print of psp
c     write(20,'(2f9.5,'' znuc,rpotc'')') znuc(1)
c     write(20,'(i6,f18.14,'' nr_ps(ict),arg(ict)'')') nr,arg(ic)
c     do 195 ir=1,nr
c 195   write(20,'(2d20.12)') r(ir),vpseudo(ir,ic,lpot(ic))

  200 continue

      call gesqua(nquad,xq0,yq0,zq0,wq)
c     call gesqua(nquad,xq,yq,zq,wq)

      write(45,'(''Quadrature points'',i3)') nquad
      do 210 i=1,nquad
  210   write(45,'(''xyz,w'',4f10.5)') xq0(i),yq0(i),zq0(i),wq(i)

      return
      end
c-----------------------------------------------------------------------

c compute tm-pseudopotential for electron iel
      subroutine getvps_tm(r_en,iel)

      use vmc_mod, only: MELEC, MCENT
      use atom, only: znuc, iwctype, ncent, ncent_tot
      use pseudo_tm, only: rmax
      use const, only: nelec
      use pseudo, only: lpot, vps

      implicit real*8(a-h,o-z)





      dimension r_en(nelec,ncent_tot)

      do 10 ic=1,ncent
        ict=iwctype(ic)

        r=r_en(iel,ic)
c local potential
        if(r.lt.rmax(ict)) then
          call splfit_tm(r,lpot(ict),ict,vpot)
          vps(iel,ic,lpot(ict))=vpot
         else
          vps(iel,ic,lpot(ict))=-znuc(ict)/r
        endif
c non-local pseudopotential
        do 10 l=1,lpot(ict)-1
          if(r.lt.rmax(ict)) then
            call splfit_tm(r,l,ict,vpot)
            vps(iel,ic,l)=vpot
           else
            vps(iel,ic,l)=0.d0
          endif
   10 continue

      return
      end
c-----------------------------------------------------------------------

      subroutine splfit_tm(r,l,ic,vpot)
c get spline_fit at r of TM potential for center ic and angular momentum l
c stored on shifted exponential grid

      use pseudo_mod, only: MPS_GRID
      use pseudo_tm, only: arg, d2pot, r0, vpseudo

      implicit real*8(a-h,o-z)



      dlogag=dlog(arg(ic))
      xr=(dlog((r+r0(ic))/r0(ic)))/dlogag+1.d0
      jx=int(xr)
      if(jx.lt.1) call fatal_error('SPLFIT_TM: index < 1')
      if(jx.gt.MPS_GRID) call fatal_error('SPLFIT_TM: index > MPS_GRID')

      ref0=r0(ic)*arg(ic)**(jx-1)-r0(ic)
      ref1=(ref0+r0(ic))*arg(ic)-r0(ic)
      delh=ref1-ref0

c cubic spline interpolation

      bb=(r-ref0)/delh
      aa=(ref1-r)/delh
      cc=aa*(aa**2-1.d0)*delh**2/6.d0
      dd=bb*(bb**2-1.d0)*delh**2/6.d0
      vpot=aa*vpseudo(jx,ic,l)+bb*vpseudo(jx+1,ic,l)+
     &   cc*d2pot(jx,ic,l)+dd*d2pot(jx+1,ic,l)

      return
      end
