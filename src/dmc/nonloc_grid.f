      subroutine nonloc_grid(iel,iw,x,psid,imove)

      use const, only: nelec
      use contrl_per, only: iperiodic
      use contrldmc, only: tau
      use atom, only: cent, ncent
      use config, only: xold_dmc
      use optjas, only: MPARMJ
      use qua, only: nquad, xq, yq, zq
      use casula, only: icasula, t_vpsp
      use distance_mod, only: r_en, rshift, rvec_en
      use wfsec, only: iwf, iwftype
      use optwf_contrl, only: ioptci, ioptjas, ioptorb

      use precision_kinds, only: dp
      implicit none

      interface
         function rannyu(idum)
          use precision_kinds, only: dp
         implicit none
         integer,intent(in) :: idum
         real(dp) :: rannyu
         end function rannyu
      end interface

      integer :: i, i1, i2, ic, ic_good
      integer :: iel, iel_good, ii, imove
      integer :: ioptci_sav, ioptjas_sav, ioptorb_sav, iq
      integer :: iq_good, iw
      real(dp) :: costh, one, p, pe, psid
      real(dp) :: psidi, ri, t_cum
      real(dp) :: t_norm, t_normi, tauprim
      real(dp), dimension(2) :: vpsp_det
      real(dp), dimension(MPARMJ) :: dvpsp_dj
      real(dp), dimension(*) :: x

c here vpsp_det and dvpsp_det are dummy

      iwf=iwftype(1)

      tauprim=tau
      if(icasula.gt.0)then
        call distances(iel,xold_dmc(1,1,iw,1))

        ioptjas_sav=ioptjas
        ioptorb_sav=ioptorb
        ioptci_sav=ioptci
        ioptjas=0
        ioptorb=0
        ioptci=0

        call nonloc_pot(xold_dmc(1,1,iw,1),rshift,rvec_en,r_en,pe,vpsp_det,dvpsp_dj,t_vpsp,iel,1)

        call multideterminant_tmove(psid,iel)

        ioptjas=ioptjas_sav
        ioptorb=ioptorb_sav
        ioptci=ioptci_sav

        i1=iel
        i2=iel
       else

        i1=1
        i2=nelec
      endif
      imove=0

c     do i=i1,i2
c     write(6,'(''t_vpsp before = '',100e10.2)') ((t_vpsp(ic,iq,i),ic=1,ncent),iq=1,nquad)
c     enddo

      t_norm=0.d0
      psidi=1.d0/psid
      do 10 ii=i1,i2
        i=ii
        if(i.gt.nelec) i=i-nelec
        if(i.lt.1) i=i+nelec
        do 10 iq=1,nquad
          do 10 ic=1,ncent
            t_vpsp(ic,iq,i)=t_vpsp(ic,iq,i)
            if(t_vpsp(ic,iq,i).gt.0.d0) t_vpsp(ic,iq,i)=0.d0
 10         t_norm=t_norm-t_vpsp(ic,iq,i)
      t_norm=1.d0+t_norm*tauprim
      t_normi=1.d0/t_norm
c     write(6,*) 'tnormi=',t_normi
c     do i=i1,i2
c     write(6,'(''t_vpsp after = '',100f14.6)') ((-tauprim*t_vpsp(ic,iq,i)*t_normi,ic=1,ncent),iq=1,nquad)
c     enddo

      if(t_norm.eq.1.d0) return

      t_cum=0.d0
      p=rannyu(0)
      do 20 ii=i1,i2
        i=ii
        if(i.gt.nelec) i=i-nelec
        if(i.lt.1) i=i+nelec
        do 20 iq=1,nquad
          do 20 ic=1,ncent
            t_cum=t_cum-tauprim*t_vpsp(ic,iq,i)*t_normi
            if(t_cum.gt.p)then
              ic_good=ic
              iq_good=iq
              iel_good=i
              imove=1
              go to 30
            endif
 20   continue

 30   if(imove.eq.1)then
        iq=iq_good
        ic=ic_good
        iel=iel_good
        if(icasula.lt.0) call distances(iel,xold_dmc(1,1,iw,1))
        ri=one/r_en(iel,ic)
        costh=rvec_en(1,iel,ic)*xq(iq)
     &       +rvec_en(2,iel,ic)*yq(iq)
     &       +rvec_en(3,iel,ic)*zq(iq)
        costh=costh*ri

        if(iperiodic.eq.0) then
          x(1)=r_en(iel,ic)*xq(iq)+cent(1,ic)
          x(2)=r_en(iel,ic)*yq(iq)+cent(2,ic)
          x(3)=r_en(iel,ic)*zq(iq)+cent(3,ic)
         else
          x(1)=r_en(iel,ic)*xq(iq)+cent(1,ic)+rshift(1,iel,ic)
          x(2)=r_en(iel,ic)*yq(iq)+cent(2,ic)+rshift(2,iel,ic)
          x(3)=r_en(iel,ic)*zq(iq)+cent(3,ic)+rshift(3,iel,ic)
        endif
c       write(6,*) 'moved B',iw,iel,(xold_dmc(kk,iel,iw,1),kk=1,3)
c       write(6,*) 'moved A',iw,iel,(x(kk),kk=1,3)
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine t_vpsp_sav


      use const, only: nelec
      use atom, only: ncent. ncent_tot
      use qua, only: nquad
      use pseudo_mod, only: MPS_QUAD

      use casula, only: t_vpsp
      use precision_kinds, only: dp
      implicit none

      integer :: i, ic, iq

      real(dp), dimension(ncent_tot, MPS_QUAD, nelec) :: t_vpsp_save




      save t_vpsp_save

      do 10 i=1,nelec
        do 10 iq=1,nquad
          do 10 ic=1,ncent
   10       t_vpsp_save(ic,iq,i)=t_vpsp(ic,iq,i)

      return

      entry t_vpsp_get

      do 20 i=1,nelec
        do 20 iq=1,nquad
          do 20 ic=1,ncent
   20       t_vpsp(ic,iq,i)=t_vpsp_save(ic,iq,i)

      return
      end
