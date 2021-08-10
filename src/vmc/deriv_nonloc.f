      subroutine deriv_nonlocj(iel,x,rshift,rvec_en,r_en,rr_en,rr_en2,dd1,value,gn,vjn,da_ratio_jn)

c Written by Claudia Filippi, modified by Cyrus Umrigar
      use atom, only: iwctype, nctype, ncent, ncent_tot

      use jaspar, only: nspin2, sspinn, is
      use const, only: nelec
      use da_jastrow4val, only: da_j
      use derivjas, only: go
      use elec, only: nup
      use jaso, only: fso
      use jaspointer, only: npoint, npointa
      use optwf_nparmj, only: nparma, nparmb, nparmc
      use optwf_parms, only: nparmj
      use bparm, only: nocuspb, nspin2b
      use contr2, only: ijas
      use contr2, only: isc
      use contrl_per, only: iperiodic
      use force_analy, only: iforce_analy

      use precision_kinds, only: dp
      implicit none

      interface
      function dpsibnl(u,isb,ipar)
        use precision_kinds, only: dp
        real(dp), intent(in) :: u
        integer, intent(in) :: isb
        integer, intent(in) :: ipar
        real(dp) :: dpsibnl
      endfunction dpsibnl

      function dpsianl(rr,it)
        use precision_kinds, only: dp
        real(dp), intent(in) :: rr
        integer, intent(in) :: it
        real(dp) :: dpsianl
      endfunction dpsianl

      end interface

      integer :: i, ic, iel, ipar, ipara
      integer :: iparm, iparm0, isb, it
      integer :: j, jj, jparm, k
      real(dp) :: dd1u, deriv_psianl, deriv_psibnl, deriv_psinl, dum
      real(dp) :: dumk, fsumn, rij, u
      real(dp) :: value
      real(dp), dimension(3,*) :: x
      real(dp), dimension(3,nelec,ncent_tot) :: rshift
      real(dp), dimension(3,nelec,*) :: rvec_en
      real(dp), dimension(nelec,ncent_tot) :: r_en
      real(dp), dimension(nelec,ncent_tot) :: rr_en
      real(dp), dimension(nelec,ncent_tot) :: rr_en2
      real(dp), dimension(*) :: gn
      real(dp), dimension(nelec,nelec) :: fsn
      real(dp), dimension(3) :: dx
      real(dp), dimension(nelec,ncent_tot) :: dd1
      real(dp), dimension(3) :: vjn
      real(dp), dimension(3,ncent_tot) :: da_ratio_jn
      real(dp), parameter :: half = .5d0




      fsumn=0

       do 1 k=1,3
   1     vjn(k)=0.d0

C TMP
c     do 5 iparm=1,nparmj
c   5   gn(iparm)=gvalue(iparm)
      do 5 iparm=1,nparmj
    5   gn(iparm)=0

      if (nelec.lt.2) goto 47

      ipara=nparma(1)
      if(ijas.ge.4.and.ijas.le.6) then
        do 7 it=2,nctype
    7     ipara=ipara+nparma(it)
      endif

c     do 15 i=1,nelec
c       if(i.ne.iel) then
c         rvec_ee(1)=x(1,i)-x(1,iel)
c         rvec_ee(2)=x(2,i)-x(2,iel)
c         rvec_ee(3)=x(3,i)-x(3,iel)
c         rij=rvec_ee(1)**2+rvec_ee(2)**2+rvec_ee(3)**2
c         rij=dsqrt(rij)
c         call scale_dist(rij,u(i),1)
c       endif
c  15 continue

      do 45 jj=1,nelec

        if(jj.eq.iel) goto 45
        if(jj.lt.iel) then
          i=iel
          j=jj
         else
          i=jj
          j=iel
        endif

        sspinn=1
        ipar=0
        if(i.le.nup .or. j.gt.nup) then
          if(nspin2.ge.2) then
            is=2
            isb=is
            if(nspin2.eq.3 .and. j.gt.nup) then
             is=3
             isb=is
            endif
           else
            is=1
            isb=is
            if(nspin2b.eq.2) then
              isb=2
             elseif(nocuspb.eq.0) then
              sspinn=half
            endif
          endif
          ipar=1
         else
          is=1
          isb=1
        endif

        do 10 k=1,3
   10     dx(k)=x(k,jj)-x(k,iel)
        if(iperiodic.eq.0) then
          rij=0
          do 20 k=1,3
   20       rij=rij+dx(k)**2
          rij=dsqrt(rij)
         else
          call find_image3(dx,rij)
        endif

c e-e terms
        if(iforce_analy.eq.0) then
          call scale_dist(rij,u,1)
         else
          call scale_dist1(rij,u,dd1u,1)
          dum=dpsibnl(u,isb,ipar)*dd1u/rij
          do 30 k=1,3
            dumk=-dum*dx(k)
   30       vjn(k)=vjn(k)+dumk
        endif

        iparm0=ipara
        if(isb.eq.2) iparm0=iparm0+nparmb(1)
        fsn(i,j)=deriv_psibnl(u,gn(iparm0+1),isb,ipar)
c       fsn(i,j)=fsn(i,j) + deriv_psibnl(u(jj),gn(iparm0+1),isb,ipar)

        do 25 jparm=1,nparmb(isb)
          iparm=iparm0+jparm
   25     gn(iparm)=gn(iparm)-go(i,j,iparm)

c e-e-n terms
c The scaling is switched in deriv_psinl, so do not do it here.
      if(isc.ge.12) call scale_dist(rij,u,3)

        do 40 ic=1,ncent
          it=iwctype(ic)
c         ri=r_en(i,ic)
c         rj=r_en(j,ic)
          iparm0=npoint(it)
   40     fsn(i,j)=fsn(i,j) +
     &    deriv_psinl(u,rshift(1,i,ic),rshift(1,j,ic),rr_en2(i,ic),rr_en2(j,ic),gn(iparm0+1),it)
c    &    deriv_psinl(u(jj),rr_en2(i,ic),rr_en2(j,ic),gn(iparm0+1),it)

        do 42 it=1,nctype
          iparm0=npoint(it)
          do 42 jparm=1,nparmc(it)
            iparm=iparm0+jparm
   42       gn(iparm)=gn(iparm)-go(i,j,iparm)

        fsumn=fsumn+fsn(i,j)-fso(i,j)
   45 continue

c e-n terms
   47 fsn(iel,iel)=0

      if(ijas.ge.4.and.ijas.le.6) then
        do 55 ic=1,ncent
          it=iwctype(ic)
          iparm0=npointa(it)
   55     fsn(iel,iel)=fsn(iel,iel)+
     &                 deriv_psianl(rr_en(iel,ic),gn(iparm0+1),it)
        do 56 it=1,nctype
          iparm0=npointa(it)
          do 56 jparm=1,nparma(it)
            iparm=iparm0+jparm
   56       gn(iparm)=gn(iparm)-go(iel,iel,iparm)
      endif

      fsumn=fsumn+fsn(iel,iel)-fso(iel,iel)
      value=fsumn

      if(iforce_analy.eq.0) return

       do 70 ic=1,ncent
        it=iwctype(ic)
        dum=dpsianl(rr_en(iel,ic),it)*dd1(iel,ic)/r_en(iel,ic)
        do 70 k=1,3
          dumk=dum*rvec_en(k,iel,ic)
          vjn(k)=vjn(k)+dumk
   70     da_ratio_jn(k,ic)=-dumk-da_j(k,iel,ic)


      return
      end
