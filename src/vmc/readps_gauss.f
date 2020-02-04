      subroutine readps_gauss
c read 'Quantum-chemist' gauss pseudopotentials
c file format: one text file with basename gauss_ecp.dat
c              for each atom type 
c first line : arbitrary label (written to log-file)
c second line: number of projectors + 1 (i.e. total number of components)
c remaining lines: components in the order (local,L=0,L=1 ...)
c     repeated for each component
c        number terms 
c        repeated for each term in this component 
c          coefficient power exponent 
c
c NOTE: as usual power n means r**(n-2)
c
      implicit real*8(a-h,o-z)
      include 'vmc.h'
      include 'pseudo.h'
      include 'force.h'

      character*20 atomtyp,atomsymbol
      character*80 label
      character*256 filename,pooldir,pp_id

      parameter (ncoef=5)

      common /const/ pi,hb,etrial,delta,deltai,fbias,nelec,imetro,ipr
      common /atom/ znuc(MCTYPE),cent(3,MCENT),pecent
     &,iwctype(MCENT),nctype,ncent

      common /pseudo/ vps(MELEC,MCENT,MPS_L),vpso(MELEC,MCENT,MPS_L,MFORCE)
     &,lpot(MCTYPE),nloc

      common /gauss_ecp/ necp_term(MPS_L,MCTYPE),necp_power(MGAUSS
     &     ,MPS_L,MCTYPE),ecp_coef(MGAUSS,MPS_L,MCTYPE)
     &     ,ecp_exponent(MGAUSS,MPS_L,MCTYPE)

      common /qua/ xq0(MPS_QUAD),yq0(MPS_QUAD),zq0(MPS_QUAD)
     &,xq(MPS_QUAD),yq(MPS_QUAD),zq(MPS_QUAD),wq(MPS_QUAD),nquad

c pool directory for pseudopotentials
      call p2gtad('general:pool',pooldir,'.',1)
      call stripquotes(pooldir)
      call p2gtad('general:pseudopot',pp_id,'none',1)
      call stripquotes(pp_id)
CVARDOC String to identify pseudopotential. If set, fancy names for 
CVARDOC the pseudopotential files will be used.  

      do 300 ic=1,nctype

        if(ic.lt.10) then
          write(atomtyp,'(i1)') ic
         elseif(ic.lt.100) then
          write(atomtyp,'(i2)') ic
         else
          call fatal_error('READPS_GAUSS: nctype>100')
        endif


        if(pp_id.eq.'none') then
          call fatal_error('READPS_GAUSS: ECP name missing')
         else
          call p2gtad('atom_types:'//atomtyp(1:index(atomtyp,' ')-1)
     &         ,atomsymbol,'X',1)
          filename=pooldir(1:index(pooldir,' ')-1)//
     &         '/'//
     &         pp_id(1:index(pp_id,' ')-1)//
     &         '.gauss_ecp.dat.'//
     &         atomsymbol(1:index(atomsymbol,' ')-1)
       endif

        open(1,file=filename(1:index(filename,' ')),status='old')
 
        write(45,'(''Reading pseudopotential file '',a)')
     &      filename(1:index(filename,' '))

c label 
        read(1,100,err=999,end=1000)   label
        write(45,101) ic,label

c max projector
        read(1,*,err=999,end=1000) lpot(ic)
        write(45,111) lpot(ic)

        if(lpot(ic).gt.MPS_L)
     &  call fatal_error('READPS_GAUSS: increase MPS_L')

c read terms of local part and all non-local parts
c local part first in file, but stored at index lpot
c non-local l=0 at index 1 etc, up to lpot-1 

        do 200 l=1,lpot(ic)
          if(l.eq.1)then
           idx=lpot(ic)
           else
            idx=l-1
          endif
          read(1,*,err=999,end=1000) necp_term(idx,ic)

          if(necp_term(idx,ic).gt.MGAUSS) 
     &     call fatal_error('READPS_GAUSS: increase MGAUSS')
          write(45,112) l,necp_term(idx,ic)
          do 200 i=1,necp_term(idx,ic)
            read(1,*,err=999,end=1000) ecp_coef(i,idx,ic),
     &        necp_power(i,idx,ic),ecp_exponent(i,idx,ic)
            write(45,113)  ecp_coef(i,idx,ic),necp_power(i,idx,ic)
     &        ,ecp_exponent(i,idx,ic)
 200  continue

      close(1)
 300  continue

      call gesqua(nquad,xq0,yq0,zq0,wq)

 100  format(a80)
 101  format('ECP for atom type ',i4,' label= ',a80)
 111  format('                         lpot = ',i3)
 112  format('    component, #terms ',2i6)
 113  format('    coef, power, expo ',f16.8,i2,f16.8)


      return

 999  continue
      call fatal_error('READPS_GAUSS: error while reading potential')

 1000 continue
      call fatal_error('READPS_GAUSS: end of file while reading potential')

      return
      end
c-----------------------------------------------------------------------

c compute gauss-pseudopotential for electron iel
      subroutine getvps_gauss(rvec_en,r_en,iel)

      implicit real*8(a-h,o-z)
      include 'vmc.h'
      include 'pseudo.h'
      include 'force.h'

      common /atom/ znuc(MCTYPE),cent(3,MCENT),pecent
     &,iwctype(MCENT),nctype,ncent

      common /pseudo/ vps(MELEC,MCENT,MPS_L),vpso(MELEC,MCENT,MPS_L,MFORCE)
     &,lpot(MCTYPE),nloc

      common /da_pseudo/ da_pecent(3,MCENT),da_vps(3,MELEC,MCENT,MPS_L),
     & da_nonloc(3,MCENT)

      dimension r_en(MELEC,MCENT),rvec_en(3,MELEC,MCENT)

      do 10 ic=1,ncent
        ict=iwctype(ic)

        r=max(1.0d-10,r_en(iel,ic))
        ri=1.d0/r
        ri2=ri*ri
c local potential
        vpot=0.d0
        dvpot=0.d0
        call gauss_pot(r,lpot(ict),ict,vpot,dvpot)
        vpot=vpot-znuc(ict)*ri
        vps(iel,ic,lpot(ict))=vpot
        do 5 k=1,3
          reni=rvec_en(k,iel,ic)*ri
    5     da_vps(k,iel,ic,lpot(ict))=-(dvpot+znuc(ict)*ri2)*reni
c non-local pseudopotential
        do 10 l=1,lpot(ict)-1
         vpot=0.d0
         dvpot=0.d0
         call gauss_pot(r,l,ict,vpot,dvpot)
         vps(iel,ic,l)=vpot
         do 10 k=1,3
          reni=rvec_en(k,iel,ic)*ri
          da_vps(k,iel,ic,l)=-dvpot*reni
   10 continue

c     do ic=1,ncent
c       write(6,*) 'HELLO_GAUSS',da_vps(1,iel,ic,lpot(iwctype(ic)))
c     enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine gauss_pot(r,l,ict,vpot,dvpot)
      implicit real*8(a-h,o-z)
      include 'vmc.h'
      include 'pseudo.h'
      common /gauss_ecp/ necp_term(MPS_L,MCTYPE),necp_power(MGAUSS
     &     ,MPS_L,MCTYPE),ecp_coef(MGAUSS,MPS_L,MCTYPE)
     &     ,ecp_exponent(MGAUSS,MPS_L,MCTYPE)

      v = 0.d0
      dv = 0.d0
      rsq=r**2

      do i=1, necp_term(l,ict)
       if(necp_power(i,l,ict).ne.2)then
        p = r**(necp_power(i,l,ict)-2)
        dp = (necp_power(i,l,ict)-2)*p/r
       else
        p = 1.d0
        dp= 0.d0
       endif
       e = ecp_coef(i,l,ict)*exp(-ecp_exponent(i,l,ict)*rsq)
       v = v + p*e
       dv = dv + (dp -2*p*ecp_exponent(i,l,ict)*r)*e
      enddo

      vpot=v
      dvpot=dv
      return
      end
