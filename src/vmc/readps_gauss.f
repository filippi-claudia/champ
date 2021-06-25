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
      use pseudo_mod, only: MPS_L, MGAUSS
      use atom, only: nctype
      use gauss_ecp, only: ecp_coef, ecp_exponent, necp_power, necp_term
      use pseudo, only: lpot
      use qua, only: nquad, wq, xq0, yq0, zq0
      use general, only: filename, filenames_ps_gauss

      use precision_kinds, only: dp
      implicit none

      integer :: i, ic, identify, idx, index
      integer :: l, names

      character*80 label

CVARDOC String to identify pseudopotential. If set, fancy names for
CVARDOC the pseudopotential files will be used.

      do 300 ic=1,nctype
        if (nctype.gt.100) call fatal_error('READPS_GAUSS: nctype>100')

        filename = filenames_ps_gauss(ic)
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

      use vmc_mod, only: MELEC, MCENT
      use atom, only: znuc, iwctype, ncent, ncent_tot
      use const, only: nelec
      use pseudo, only: lpot, vps

      use da_pseudo, only: da_vps

      use precision_kinds, only: dp
      implicit none

      integer :: ic, ict, iel, k, l
      real(dp) :: dvpot, r, reni, ri, ri2
      real(dp) :: vpot
      real(dp), dimension(nelec,ncent_tot) :: r_en
      real(dp), dimension(3,nelec,ncent_tot) :: rvec_en






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
      use gauss_ecp, only: ecp_coef, ecp_exponent, necp_power, necp_term

      use precision_kinds, only: dp
      implicit none

      integer :: i, ict, l
      real(dp) :: dp_val, dv, dvpot, e, p
      real(dp) :: r, rsq, v, vpot


      v = 0.d0
      dv = 0.d0
      rsq=r**2

      do i=1, necp_term(l,ict)
       if(necp_power(i,l,ict).ne.2)then
        p = r**(necp_power(i,l,ict)-2)
        dp_val = (necp_power(i,l,ict)-2)*p/r
       else
        p = 1.d0
        dp_val= 0.d0
       endif
       if (ecp_exponent(i,l,ict)*rsq.gt.7.0D2) then
          e = 0.0d0
       else
          e = ecp_coef(i,l,ict)*exp(-ecp_exponent(i,l,ict)*rsq)
       endif
       v = v + p*e
       dv = dv + (dp_val -2*p*ecp_exponent(i,l,ict)*r)*e
      enddo

      vpot=v
      dvpot=dv
      return
      end
