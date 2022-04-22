      module nonlpsi
      use error, only: fatal_error
      contains
      function psinl(u,rshifti,rshiftj,rri,rrj,it)
c Written by Claudia Filippi, modified by Cyrus Umrigar

      use vmc_mod, only: nordj
      use jaspar3, only: c
      use jaspar4, only: nordc
      use jaspar6, only: asymp_r
      use wfsec, only: iwf
      use contr2, only: ijas
      use precision_kinds, only: dp
      use scale_dist_mod, only: switch_scale

      implicit none

      integer :: it, jp, k, l, l_hi
      integer :: ll, m, n
      real(dp) :: rri, rrj, rrri, rrrj, u
      real(dp) :: uuu
      real(dp) :: psinl
      real(dp), dimension(0:nordj) :: uu
      real(dp), dimension(0:nordj) :: ss
      real(dp), dimension(0:nordj) :: tt
      real(dp), dimension(3) :: rshifti
      real(dp), dimension(3) :: rshiftj
      real(dp), parameter :: one = 1.d0
      real(dp), parameter :: two = 2.d0
      real(dp), parameter :: half = 0.5d0
      real(dp), parameter :: eps = 1.d-12

c Not updated for ijas=5,6 because we will probably stay with ijas=4
c If we want to use ijas=5,6 update this routine similarly to psi.f
      if(ijas.ge.5) call fatal_error('PSINL: ijas >= 5 not implemented')

      psinl=0.d0
      if(nordc.le.1) return

      if(rri.eq.asymp_r .or. rrj.eq.asymp_r) return
      do k=1,3
        if(abs(rshifti(k)-rshiftj(k)).gt.eps) return
      enddo

      uuu=u
      rrri=rri
      rrrj=rrj
      call switch_scale(uuu)
      call switch_scale(rrri)
      call switch_scale(rrrj)

      uu(0)=one
      ss(0)=two
      tt(0)=one
      do jp=1,nordc
        uu(jp)=uuu**jp
        ss(jp)=rrri**jp+rrrj**jp
        tt(jp)=(rrri*rrrj)**jp
      enddo

      ll=0
      do n=2,nordc
        do k=n-1,0,-1
          if(k.eq.0) then
            l_hi=n-k-2
           else
            l_hi=n-k
          endif
          do l=l_hi,0,-1
            m=(n-k-l)/2
            if(2*m.eq.n-k-l) then
              ll=ll+1
              psinl=psinl+c(ll,it,iwf)*uu(k)*ss(l)*tt(m)
            endif
          enddo
        enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      function psianl(rri,it)


      use jaspar4, only: a4, norda
      use jaspar6, only: asymp_jasa, asymp_r
      use wfsec, only: iwf
      use contr2, only: ijas
      use precision_kinds, only: dp
      implicit none

      integer :: i, it
      real(dp) :: rri
      real(dp) :: psianl


c Not updated for ijas=5,6 because we will probably stay with ijas=4
c If we want to use ijas=5,6 update this routine similarly to psi.f
      if(ijas.ge.5) call fatal_error('PSINL: ijas >= 5 not implemented')

      psianl=0.d0
      if(rri.eq.asymp_r) return

      psianl=a4(1,it,iwf)*rri/(1.d0+a4(2,it,iwf)*rri)-asymp_jasa(it)
      do i=2,norda
        psianl=psianl+a4(i+1,it,iwf)*rri**i
      enddo

      return
      end
c-----------------------------------------------------------------------

      function psibnl(u,isb,ipar)

      use jaspar, only: sspinn
      use jaspar3, only: b
      use jaspar4, only: nordb
      use jaspar6, only: asymp_jasb, asymp_r
      use wfsec, only: iwf
      use contr2, only: ijas
      use precision_kinds, only: dp

      implicit none

      integer :: i, ipar, isb
      real(dp) :: fee, psibnl, u

c Not updated for ijas=5,6 because we will probably stay with ijas=4
c If we want to use ijas=5,6 update this routine similarly to psi.f
      if(ijas.ge.5) call fatal_error('PSINL: ijas >= 5 not implemented')

      psibnl=0.d0
      if(ijas.le.2) return
      if(u.eq.asymp_r) return

      fee=b(1,isb,iwf)*u/(1+b(2,isb,iwf)*u)
      psibnl=sspinn*fee

      psibnl=psibnl-asymp_jasb(ipar+1)
      do i=2,nordb
        psibnl=psibnl+b(i+1,isb,iwf)*u**i
      enddo

      return
      end
c-----------------------------------------------------------------------
      function dpsianl(rri,it)

      use jaspar4, only: a4, norda
      use jaspar6, only: asymp_r
      use wfsec, only: iwf
      use contr2, only: ijas
      use precision_kinds, only: dp

      implicit none

      integer :: i, it
      real(dp) :: rri
      real(dp) :: dpsianl

c Not updated for ijas=5,6 because we will probably stay with ijas=4
c If we want to use ijas=5,6 update this routine similarly to psi.f
      if(ijas.ge.5) call fatal_error('PSINL: ijas >= 5 not implemented')

      dpsianl=0.d0
      if(rri.eq.asymp_r) return

      do i=2,norda
        dpsianl=dpsianl+i*a4(i+1,it,iwf)*rri**(i-1)
      enddo

      return
      end

c-----------------------------------------------------------------------
      function dpsibnl(u,isb,ipar)

      use jaspar, only: sspinn
      use jaspar3, only: b
      use jaspar4, only: nordb
      use jaspar6, only: asymp_r
      use wfsec, only: iwf
      use contr2, only: ijas
      use precision_kinds, only: dp

      implicit none

      integer :: i, isb, ipar
      real(dp) :: bot, boti, dbot, dfee, dtop
      real(dp) :: top, u
      real(dp) :: dpsibnl

c Not updated for ijas=5,6 because we will probably stay with ijas=4
c If we want to use ijas=5,6 update this routine similarly to psi.f
      if(ijas.ge.5) call fatal_error('PSINL: ijas >= 5 not implemented')

      dpsibnl=0.d0
      if(u.eq.asymp_r) return

      top=b(1,isb,iwf)*u
      dtop=b(1,isb,iwf)
      bot=1+b(2,isb,iwf)*u
      dbot=b(2,isb,iwf)
      boti=1.d0/bot

      dfee=dtop*boti-top*boti*boti*dbot
      dpsibnl=sspinn*dfee

      do i=2,nordb
        dpsibnl=dpsibnl+i*b(i+1,isb,iwf)*u**(i-1)
      enddo

      return
      end
      end module 
