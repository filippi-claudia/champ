      module deriv_nonlpsi
      contains
      function deriv_psinl(u,rshifti,rshiftj,rri,rrj,gn,it,iwfjas)
c Written by Claudia Filippi, modified by Cyrus Umrigar

      use vmc_mod, only: nwftypejas
      use jastrow, only: c, nordc, ijas, nordj
      use jaspar6, only: asymp_r
      use jastrow, only: c,ijas,nordj
      use multiple_geo, only: iwf
      use optwf_wjas, only: iwjasc
      use vardep, only: cdep, iwdepend, nvdepend
      use cuspmat4, only: d, iwc4
      use precision_kinds, only: dp
      use scale_dist_mod, only: switch_scale
      use vardep,  only: cdep,iwdepend,nvdepend
      implicit none

      integer :: id, ideriv, iparm, it, jj
      integer :: jp, jparm, k, l, iwfjas
      integer :: l_hi, ll, m, n
      real(dp) :: deriv_psinl, p, rri, rrj, rrri
      real(dp) :: rrrj, u
      real(dp), dimension(3) :: rshifti
      real(dp), dimension(3) :: rshiftj
      real(dp), dimension(*) :: gn
      real(dp), dimension(0:nordj) :: uu
      real(dp), dimension(0:nordj) :: ss
      real(dp), dimension(0:nordj) :: tt
      real(dp), parameter :: zero = 0.d0
      real(dp), parameter :: one = 1.d0
      real(dp), parameter :: two = 2.d0
      real(dp), parameter :: half = 0.5d0
      real(dp), parameter :: eps = 1.d-12

      deriv_psinl = 0.0




      if(ijas.ge.4.and.ijas.le.6) then

        deriv_psinl=0
        if(nordc.eq.0) return

        if(rri.eq.asymp_r .or. rrj.eq.asymp_r) return
        do k=1,3
          if(abs(rshifti(k)-rshiftj(k)).gt.eps) return
        enddo

        if(nwftypejas.gt.1) iwf=iwfjas

        uu(1)=u
        rrri=rri
        rrrj=rrj
        call switch_scale(uu(1))
        call switch_scale(rrri)
        call switch_scale(rrrj)

        uu(0)=1
        ss(0)=2
        tt(0)=1
        do jp=1,nordc
          uu(jp)=uu(1)*uu(jp-1)
          ss(jp)=rrri**jp+rrrj**jp
          tt(jp)=(rrri*rrrj)**jp
        enddo

        ll=0
        jj=1
        jparm=1
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
                p=uu(k)*ss(l)*tt(m)
                deriv_psinl=deriv_psinl+c(ll,it,iwf)*p

                ideriv=0
                if(ll.eq.iwjasc(jparm,it)) then
                  ideriv=2
                 else
                  do id=1,2*(nordc-1)
                    if(ll.eq.iwc4(id)) then
                      jj=id
                      if(nvdepend(jj,it).gt.0) ideriv=1
                    endif
                  enddo
                endif

                if(ideriv.eq.1) then
                  do id=1,nvdepend(jj,it)
                    iparm=iwdepend(jj,id,it)
                    gn(iparm)=gn(iparm)+cdep(jj,id,it)*p
                  enddo
c                 jj=jj+1
                 elseif(ideriv.eq.2) then
                  gn(jparm)=gn(jparm)+p
                  jparm=jparm+1
                endif
              endif
            enddo
          enddo
        enddo

      endif

      return
      end

c-----------------------------------------------------------------------
      function deriv_psianl(rri,gn,it,iwfjas)


      use jastrow, only: a4, norda, asymp_jasa, ijas
      use jaspar6, only: asymp_r
      use optwf_nparmj, only: nparma
      use optwf_wjas, only: iwjasa
      use multiple_geo, only: iwf
      use precision_kinds, only: dp
      use vmc_mod, only: nwftypejas

      implicit none

      integer :: i, iord, it, jparm, iwfjas
      real(dp) :: bot, deriv_psianl, gen, rri, top
      real(dp), dimension(*) :: gn
      real(dp), parameter :: one = 1.d0

      deriv_psianl = 0.0






c Note: This routine is only called with iwf=1, but parts of it are
c written for general iwf, whereas others (asymp_r) assume iwf=1.

      if(rri.eq.asymp_r) then
        deriv_psianl=0
        return
      endif

      if(nwftypejas.gt.1) iwf=iwfjas

      if(ijas.ge.4.and.ijas.le.6) then

        deriv_psianl=a4(1,it,iwf)*rri/(one+a4(2,it,iwf)*rri)-asymp_jasa(it,iwf)
        do i=2,norda
          deriv_psianl=deriv_psianl+a4(i+1,it,iwf)*rri**i
        enddo
        do jparm=1,nparma(it)
            if(iwjasa(jparm,it).eq.1) then
              top=rri
              bot=one+a4(2,it,iwf)*rri
              gen=top/bot-asymp_r/(1+a4(2,it,iwf)*asymp_r)
             elseif(iwjasa(jparm,it).eq.2) then
              top=-a4(1,it,iwf)*rri*rri
              bot=one+a4(2,it,iwf)*rri
              bot=bot*bot
              gen=top/bot+a4(1,it,iwf)*asymp_r**2/(1+a4(2,it,iwf)*asymp_r)**2
             else
              iord=iwjasa(jparm,it)-1
              gen=rri**iord-asymp_r**iord
            endif
            gn(jparm)=gn(jparm)+gen
        enddo
      endif

      return
      end

c-----------------------------------------------------------------------
      function deriv_psibnl(u,gn,isb,ipar,iwfjas)


      use jastrow, only: nordb
      use jaspar6, only: asymp_r
      use jastrow, only: asymp_jasb,b,ijas,sspinn
      use multiple_geo, only: iwf
      use optwf_nparmj, only: nparmb
      use optwf_wjas, only: iwjasb
      use precision_kinds, only: dp
      use vmc_mod, only: nwftypejas

      implicit none

      integer :: i, iord, ipar, isb, jparm, iwfjas
      real(dp) :: bot, deriv_psibnl, fee, gee, top
      real(dp) :: u
      real(dp), dimension(*) :: gn
      real(dp), parameter :: one = 1.d0















c Note: This routine is only called with iwf=1, but parts of it are
c written for general iwf, whereas others (asymp_r) assume iwf=1.

c Not updated for ijas=5,6 because we will probably stay with ijas=4
c If we want to use ijas=5,6 update this routine similarly to psi.f
      if(ijas.ge.5) stop 'ijas >= 5 not implemented in psibnl'

      if(u.eq.asymp_r) then
        deriv_psibnl=0
        return
      endif

      if(nwftypejas.gt.1) iwf=iwfjas

      fee=b(1,isb,iwf)*u/(one+b(2,isb,iwf)*u)

      deriv_psibnl=sspinn*fee-asymp_jasb(ipar+1,iwf)
      if(ijas.ge.4.and.ijas.le.6) then
        do i=2,nordb
          deriv_psibnl=deriv_psibnl+b(i+1,isb,iwf)*u**i
        enddo
      endif

      do jparm=1,nparmb(isb)
        if(iwjasb(jparm,isb).eq.1) then
          top=u
          bot=one+b(2,isb,iwf)*u
          gee=sspinn*(top/bot-asymp_r/(1+b(2,isb,iwf)*asymp_r))
         elseif(iwjasb(jparm,isb).eq.2) then
          top=-b(1,isb,iwf)*u*u
          bot=one+b(2,isb,iwf)*u
          bot=bot*bot
          gee=sspinn*(top/bot+b(1,isb,iwf)*asymp_r**2/(1+b(2,isb,iwf)*asymp_r)**2)
         else
          iord=iwjasb(jparm,isb)-1
          gee=u**iord-asymp_r**iord
        endif
        gn(jparm)=gn(jparm)+gee
      enddo

      return
      end
      end module 
