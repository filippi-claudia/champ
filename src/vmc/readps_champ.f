      subroutine readps_champ
c Written by Cyrus Umrigar, modified by Claudia Filippi
c for the moment always assume that highest component is the local

c Read CHAMP-formatted Troullier-Martins pseudopotentials.
c Reads in v in Hartrees and subtracts out local part from all except
c the lpot component.
c Also initializes quadrature pts.
c rmax_coul is the point at which the psp. becomes -Z/r to within eps.
c rmax_nloc is the point at which the psp. becomes local to within eps.
c For TM psps. rmax_nloc is considerably smaller than rmax_coul and so
c considerable computer time can be saved by evaluating the nonlocal
c components only for r < rmax_nloc rather than r < rmax_coul.

c Can use 3 different grids:
c igrid_ps=1, linear,              r(i)=r0_ps+(i-1)*h_ps
c         =2, exponential,         r(i)=r0_ps*exp((i-1)*h_ps)
c         =3, shifted exponential, r(i)=r0_ps*(exp((i-1)*h_ps)-1)
c The prefered grid is 3.

      use atom, only: znuc, cent, pecent, iwctype, nctype, ncent
      use const, only: pi, hb, etrial, delta, deltai, fbias, nelec, imetro, ipr
      use pseudo_champ, only: igrid_ps, rmax_coul, rmax_nloc
      use pseudo_tm, only: arg_ps, d2pot, nr_ps, r0_ps, rmax_ps, vpseudo

      use pseudo, only: lpot, nloc, vps, vpso

      implicit real*8(a-h,o-z)




      include 'vmc.h'
      include 'pseudo.h'
      include 'force.h'

      character*20 atomtyp,atomsymbol
      character*256 filename,pooldir,pp_id
      character*80 title




      common /qua/ xq0(MPS_QUAD),yq0(MPS_QUAD),zq0(MPS_QUAD)
     &,xq(MPS_QUAD),yq(MPS_QUAD),zq(MPS_QUAD),wq(MPS_QUAD),nquad

      dimension r(MPS_GRID),work(MPS_GRID),lpot_max(MPS_L)

c pool directory for pseudopotentials
      call p2gtad('general:pool',pooldir,'.',1)
      call stripquotes(pooldir)
      call p2gtad('general:pseudopot',pp_id,'none',1)
      call stripquotes(pp_id)
CVARDOC String to identify pseudopotential. If set, fancy names for 
CVARDOC the pseudopotential files will be used.  

      do 200 ict=1,nctype

        if(ict.lt.10) then
          write(atomtyp,'(i1)') ict
         elseif(ict.lt.100) then
          write(atomtyp,'(i2)') ict
         else
          call fatal_error('READPS_CHAMP: nctype>100')
        endif

        if(pp_id.eq.'none') then
c old naming convention
          filename=pooldir(1:index(pooldir,' ')-1)//'/'//
     &               'pseudopot_champ'//atomtyp(1:index(atomtyp,' ')-1)
         else
c new naming convention
          call p2gtad('atom_types:'//atomtyp(1:index(atomtyp,' ')-1)
     &         ,atomsymbol,'X',1)
          filename=pooldir(1:index(pooldir,' ')-1)//
     &             '/'//
     &             pp_id(1:index(pp_id,' ')-1)//
     &             '.pseudopot_champ.'//
     &             atomsymbol(1:index(atomsymbol,' ')-1)
        endif

        open(1,file=filename,status='old',form='formatted',err=999)
        write(45,'(''Reading CHAMP format pseudopotential file '',a20)') filename

c position file to skip comments
        title(1:1)='#'
        do while(title(1:1).eq.'#')
          read(1,'(a80)') title
        enddo

c The TM psp. format has npotd and npotu for down and up, but we just use one of them (lpot_max)
c They are the number of different l components of the psp.

        write(6,'(''Reading psp in champ format'')')
        read(title,*) lpot_max(ict),zion,r_asymp
        write(6,'(''ict,lpot_max(ict),zion,r_asymp'',2i2,f4.0,f8.3)') ict,lpot_max(ict),zion,r_asymp
        if(lpot_max(ict).le.0 .or. lpot_max(ict).gt.MPS_L) call fatal('READTM_CHAMP: lpot_max must be > 0 and <= MPS_L')

c Presently, only allow highest component to be the local
c       if(lpot(ict).gt.lpot_max(ict)) then
c         write(6,'(''lpot(ict),lpot_max(ict)='',2i3)') lpot(ict),lpot_max(ict)
c         call fatal('READTM_CHAMP: Cannot choose local psp. to be > number of l components, lpot(ict) > lpot_max(ict)')
c       endif
c       if(lpot_max(ict).gt.MPS_L) call fatal('READTM_CHAMP: lpot_max(ict).gt.MPS_L')

c If the local pseudopot component is not set in input, set it here
c       if(lpot(ict).le.0) then
          lpot(ict)=lpot_max(ict)
          write(6,'(''Center type'',i4,'' local pseudopot component reset to'',i3)') ict,lpot(ict)
c       endif

        write(6,'(''Center type'',i2,'' has'',i2,'' pseudopotential L components, and component''
     &  ,i2,'' is chosen to be local'')') ict,lpot_max(ict),lpot(ict)

        if(znuc(ict).ne.zion) then
          write(6,'(''znuc(ict) != zion in readps_tm'',2f6.1)') znuc(ict),zion
          call fatal('READTM_CHAMP: znuc(ict) != zion in readps_tm')
        endif

        read(1,*) igrid_ps(ict),nr_ps(ict),r0_ps(ict),h_ps
        nr=nr_ps(ict)
        write(6,'(''igrid_ps(ict),nr_ps(ict),r0_ps(ict),h_ps='',i2,i5,1pd22.15,0pf8.5)')
     &  igrid_ps(ict),nr_ps(ict),r0_ps(ict),h_ps
        arg_ps(ict)=exp(h_ps)

        if(igrid_ps(ict).lt.1 .or. igrid_ps(ict).gt.3) stop 'igrid_ps(ict) must be 1 or 2 or 3'
        if(igrid_ps(ict).lt.1 .and. r0_ps(ict).ne.0.d0) stop 'if igrid_ps(ict)=1 r0_ps(ict) must be 0'

        if(nr.lt.100) then
          write(6,'(''nr in psp grid too small'',2i6)') nr
          stop 'nr in psp grid too small'
        endif
        if(nr.gt.MPS_GRID .or. igrid_ps(ict).eq.2.and.nr.gt.MPS_GRID-1) then
          write(6,'(''nr > MPS_GRID'',2i6)') nr,MPS_GRID
          stop 'nr > MPS_GRID'
        endif

        if(igrid_ps(ict).eq.1 .or. igrid_ps(ict).eq.3) then
          nr=nr_ps(ict)
          do 10 ir=1,nr_ps(ict)
   10       read(1,*) r(ir),(vpseudo(ir,ict,i),i=1,lpot_max(ict))
          if(r(1).ne.0.d0) stop 'if igrid_ps is 1 or 3, r(1) must be 0'
         else
          nr_ps(ict)=nr_ps(ict)+1
          nr=nr_ps(ict)
          nrm1=nr-1
          do 20 ir=2,nr_ps(ict)
   20       read(1,*) r(ir),(vpseudo(ir,ict,i),i=1,lpot_max(ict))
          r(1)=0
          if(r0_ps(ict).le.0.d0 .or. h_ps.le.0.d0) then
            r0_ps(ict)=r(2)
            arg_ps(ict)=r(3)/r(2)
            h_ps=dlog(arg_ps(ict))
            write(6,'(''Grid parameters deduced from grid values are, r0_ps(ict),h_ps,arg_ps(ict)='',9f10.5)')
     &      r0_ps(ict),h_ps,arg_ps(ict)
          endif
          do 30 i=1,lpot_max(ict)
            call intpol(r(2),vpseudo(2,ict,i),nrm1,r(1),vpseudo(1,ict,i),1,3)
   30       write(6,'(''Interpolated psp'',9f16.12)') (vpseudo(ir,ict,i),ir=1,5)
        endif

        if(r0_ps(ict).lt.0.d0 .or. r0_ps(ict).gt.1.d-2) stop 'r0_ps in psp grid is not reasonable'
        if(h_ps.le.0.d0 .or. h_ps.gt.1.d-1) stop 'h_ps in psp grid is not reasonable'

        close(1)

c Find the point beyond which r*v differs from zion by no more than .5*d-6.
c irmax_coul is used for the endpoint of the spline where it is assumed that
c the derivative of the local component is zion/r(irmax_coul)**2 and that of 
c the local component is 0.  Also, rmax_coul is used in splfit_champ when it is
c called in a calculation of a periodic system.
        rmax_coul(ict)=0.d0
        irmax_coul=0
        do 50 ir=nr,1,-1
          do 50 i=1,lpot_max(ict)
            if(dabs(r(ir)*vpseudo(ir,ict,i)+zion).gt..5d-6) then
              rmax_coul(ict)=max(rmax_coul(ict),r(ir))
              irmax_coul=ir
              goto 60
            endif
   50   continue
   60   irmax_coul=min(irmax_coul+5,nr)

c Find the point beyond which the various v components differ from each other by no more than .5*d-6
        rmax_nloc(ict)=0.d0
        irmax_nloc=0
        do 70 ir=nr,1,-1
          do 70 i=2,lpot_max(ict)
            do 70 j=1,i-1
              if(dabs(vpseudo(ir,ict,i)-vpseudo(ir,ict,j)).gt..5d-6) then
                rmax_nloc(ict)=max(rmax_nloc(ict),r(ir))
                irmax_nloc=ir
                goto 80
              endif
   70   continue
   80   irmax_nloc=irmax_nloc+1
        rmax_nloc(ict)=r(irmax_nloc)

c rmax_nloc is used in getvps_champ to decide whether to calculate calculate nonloc part of psp
c or to set it to zero.  irmax_coul is used to decide how far out to spline the psp. components
c so irmax_coul must be >= irmax_nloc.
        irmax_coul=max(irmax_coul,irmax_nloc)
        rmax_coul(ict)=r(irmax_coul)

        write(6,'(''center '',i3,'' pseudopot rmax_coul,irmax_coul,rmax_nloc,irmax_nloc= '',2(f6.2,i5))')
     &  ict,rmax_coul(ict),irmax_coul,rmax_nloc(ict),irmax_nloc

        if(ipr.ge.1) then
          write(38,'(''r(j)  (vpseudo(j,ict,i),i=1,lpot_max(ict))  -znuc(ict)/r(j)'')')
          do 104 j=2,nr
            if(r(j).gt.0.d0) then
              write(38,'(1pd12.6,9d14.6)') r(j),(vpseudo(j,ict,i),i=1,lpot_max(ict)),-znuc(ict)/r(j)
             else
              write(38,'(1pd12.6,9d14.6)') r(j),(vpseudo(j,ict,i),i=1,lpot_max(ict))
            endif
  104     continue
        endif

        do 110 i=1,lpot_max(ict)
          if(i.ne.lpot(ict)) then
            do 105 j=1,nr
  105         vpseudo(j,ict,i)=vpseudo(j,ict,i)-vpseudo(j,ict,lpot(ict))
          endif

  110   continue

        if(rmax_coul(ict).eq.0.d0) goto 200

        do 190 i=1,lpot_max(ict)

c Construct the spline

c Warning: the next line is correct only for pseudopotentials that have zero derivative at the origin.
c At present this routine is only used for such pseudopotentials, since for the Dolg pseudopotentials
c we call readps_gauss.f
          dpot1=0.d0

c Set derivative at end point equal to 0 for nonlocal components and Z/r^2 for local
          dpotn=0.d0
          if(i.eq.lpot(ict)) dpotn=zion/r(irmax_coul)**2
c         if(i.eq.lpot(ict)) call deriv_intpol(r,vpseudo(1,ict,i),nr,r(irmax_coul),dpotn,irmax_coul,3)

          write(6,'(''dpot1,dpotn'',1p2e15.5)') dpot1,dpotn

c get second derivative for spline fit
          call spline2(r,vpseudo(1,ict,i),irmax_coul,dpot1,dpotn,d2pot(1,ict,i),work)

          do 190 j=1,nr
            if(ipr.ge.2) then
              if(i.eq.1) write(35,'(1p5d14.6)') r(j),vpseudo(j,ict,i),d2pot(j,ict,i),-znuc(ict)/r(j),-2*znuc(ict)/r(j)**3
              if(i.eq.2) write(36,'(1p5d14.6)') r(j),vpseudo(j,ict,i),d2pot(j,ict,i),-znuc(ict)/r(j),-2*znuc(ict)/r(j)**3
              if(i.eq.3) write(37,'(1p5d14.6)') r(j),vpseudo(j,ict,i),d2pot(j,ict,i),-znuc(ict)/r(j),-2*znuc(ict)/r(j)**3
            endif
  190   continue

  200 continue

      call gesqua(nquad,xq0,yq0,zq0,wq)
c     call gesqua(nquad,xq,yq,zq,wq)

      write(6,'(''quadrature points'')')
      do 210 i=1,nquad
  210   write(6,'(''xyz,w'',4f10.5)') xq0(i),yq0(i),zq0(i),wq(i)

      return

  999 write(6,'(''Error: Pseudopot. file '',a20,'' is missing'')') filename
      stop 'Pseudopot. file is missing'

      end

c-----------------------------------------------------------------------
      subroutine getvps_champ(r_en,iel)
c compute pseudopotential for electron iel

      use atom, only: znuc, cent, pecent, iwctype, nctype, ncent
      use pseudo_champ, only: igrid_ps, rmax_coul, rmax_nloc
      use pseudo_tm, only: arg_ps, d2pot, nr_ps, r0_ps, rmax_ps, vpseudo

      use pseudo, only: lpot, nloc, vps, vpso

      implicit real*8(a-h,o-z)



      include 'vmc.h'
      include 'pseudo.h'
      include 'force.h'




      dimension r_en(MELEC,MCENT)

      do 10 ic=1,ncent
        ict=iwctype(ic)

        r=r_en(iel,ic)
c local potential
        if(r.lt.rmax_coul(ict)) then
          call splfit_champ(r,lpot(ict),ict,vpot)
          vps(iel,ic,lpot(ict))=vpot
         else
          vps(iel,ic,lpot(ict))=-znuc(ict)/r
        endif
c       write(6,'(''ic,iel,r,vpot='',2i3,f6.3,f9.5)') ic,iel,r,vps(iel,ic,lpot(ict))
c non-local pseudopotential
        do 10 l=1,lpot(ict)
          if(l.ne.lpot(ict)) then
            if(r.lt.rmax_nloc(ict)) then
              call splfit_champ(r,l,ict,vpot)
              vps(iel,ic,l)=vpot
             else
              vps(iel,ic,l)=0.d0
            endif
          endif
   10 continue

      return
      end
c-----------------------------------------------------------------------

      subroutine splfit_champ(r,l,ict,vpot)
c get spline fit at r=r(iel,ic) of pseudopotential for center-type ict and
c angular momentum l stored on shifted exponential grid
c Note: I check if r < rmax_coul(ict) because this routine is called from
c ewald without going through getvps_tm.
c We assume that rmax_nloc(ict) <= rmax_coul(ict).

      use atom, only: znuc, cent, pecent, iwctype, nctype, ncent
      use pseudo_champ, only: igrid_ps, rmax_coul, rmax_nloc
      use pseudo_tm, only: arg_ps, d2pot, nr_ps, r0_ps, rmax_ps, vpseudo

      use pseudo, only: lpot, nloc, vps, vpso

      implicit real*8(a-h,o-z)



      include 'vmc.h'
      include 'pseudo.h'
      include 'force.h'



      if(r.lt.rmax_coul(ict)) then

        if(igrid_ps(ict).eq.1)then
c Warning: this needs fixing to h_ps
          stop 'May need fixing in splfit_ps'
          xr=(r-r0_ps(ict))/arg_ps(ict)+1
          jx=int(xr)
          ref0=r0_ps(ict)+arg_ps(ict)*(jx-1)
          ref1=ref0+arg_ps(ict)
          delh=arg_ps(ict)
         elseif(igrid_ps(ict).eq.2)then
          h_ps=dlog(arg_ps(ict))
          xr=dlog(r/r0_ps(ict))/h_ps+1
          jx=int(xr)
          ref0=r0_ps(ict)*arg_ps(ict)**(jx-1)
          ref1=ref0*arg_ps(ict)
          delh=ref1-ref0
         elseif(igrid_ps(ict).eq.3)then
          h_ps=dlog(arg_ps(ict))
          xr=(dlog((r+r0_ps(ict))/r0_ps(ict)))/h_ps+1
          jx=int(xr)
          ref0=r0_ps(ict)*arg_ps(ict)**(jx-1)-r0_ps(ict)
          ref1=(ref0+r0_ps(ict))*arg_ps(ict)-r0_ps(ict)
          delh=ref1-ref0
        endif

        if(jx.lt.1) then
          write(6,'(''ict,jx,xr,r,r0_ps(ict),arg_ps(ict),h_ps='',2i3,5d12.4)')
     &    ict,jx,xr,r,r0_ps(ict),arg_ps(ict),h_ps
          write(6,'(''Warning: index < 1 in splfit_champ, r,xr='',2d12.4)') r,xr
          jx=1
        endif

        if(jx.gt.MPS_GRID) then
          write(6,'(''ict,jx,xr,r,r0_ps(ict),arg_ps(ict),h_ps='',2i3,5d12.4)')
     &    ict,jx,xr,r,r0_ps(ict),arg_ps(ict),h_ps
          stop 'index > MPS_GRID in splfit_champ'
        endif

c cubic spline interpolation

        bb=(r-ref0)/delh
        aa=(ref1-r)/delh
        cc=aa*(aa**2-1.d0)*delh**2/6.d0
        dd=bb*(bb**2-1.d0)*delh**2/6.d0
        vpot=aa*vpseudo(jx,ict,l)+bb*vpseudo(jx+1,ict,l)+
     &  cc*d2pot(jx,ict,l)+dd*d2pot(jx+1,ict,l)
       else
        if(l.eq.lpot(ict)) then
          vpot=-znuc(ict)/r
         else
          vpot=0
        endif
      endif

      return
      end
