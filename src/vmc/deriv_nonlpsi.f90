      module deriv_nonlpsi
      contains
      function deriv_psinl(u,rri,rrj,gn,it,iwfjas)
! Written by Claudia Filippi, modified by Cyrus Umrigar

      use vmc_mod, only: nwftypejas
      use jastrow, only: c, nordc, ijas, nordj
      use jaspar6, only: asymp_r
      use jastrow, only: cutjas_en,cutjas_eni
      use jastrow, only: a4,c,ijas,nordj
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
      real(dp) :: rrrj, u, term, xi,xj
      real(dp), dimension(*) :: gn
      real(dp), dimension(0:nordj) :: uu
      real(dp), dimension(0:nordj) :: ri
      real(dp), dimension(0:nordj) :: rj
      real(dp), dimension(0:nordj) :: ss
      real(dp), dimension(0:nordj) :: tt
      real(dp), parameter :: eps = 1.d-12

      deriv_psinl = 0.0
      if(nordc.eq.0) return



      if(nwftypejas.gt.1) iwf=iwfjas

      uu(1)=u
      ri(1)=rri
      rj(1)=rrj

      if(ijas.eq.4) then
         if(rri.eq.asymp_r .or. rrj.eq.asymp_r) return
         call switch_scale(uu(1))
         call switch_scale(ri(1))
         call switch_scale(rj(1))
         term=1.d0
      elseif(ijas.eq.1) then
         if(rri.gt.cutjas_en(it,iwf).or. rrj.gt.cutjas_en(it,iwf)) return
         xi=rri*cutjas_eni(it,iwf)
         xj=rrj*cutjas_eni(it,iwf)
         term=((1.d0-xi)*(1.d0-xj))**3
      endif

      uu(0)=1
      ri(0)=1
      rj(0)=1
      ss(0)=2
      tt(0)=1
      do jp=1,nordc
        uu(jp)=uu(1)*uu(jp-1)
        ri(jp)=ri(1)*ri(jp-1)
        rj(jp)=rj(1)*rj(jp-1)
        ss(jp)=ri(jp)+rj(jp)
        tt(jp)=ri(jp)*rj(jp)
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
                  p=uu(k)*ss(l)*tt(m)*term
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
!     jj=jj+1
                  elseif(ideriv.eq.2) then
                     gn(jparm)=gn(jparm)+p
                     jparm=jparm+1
                  endif
               endif
            enddo
         enddo
      enddo



      return
      end

!-----------------------------------------------------------------------
      function deriv_psianl(rri,gn,it,iwfjas)


      use jastrow, only: a4, norda, asymp_jasa, ijas
      use jaspar6, only: asymp_r
      use jastrow, only: cutjas_en,cutjas_eni
      use optwf_nparmj, only: nparma
      use optwf_wjas, only: iwjasa
      use multiple_geo, only: iwf
      use precision_kinds, only: dp
      use vmc_mod, only: nwftypejas

      implicit none

      integer :: i, iord, it, jparm, iwfjas
      real(dp) :: a1_cusp, bot, deriv_psianl, gen, rri, top, term, xi
      real(dp) :: da1_cusp
      real(dp), dimension(*) :: gn
      real(dp), dimension(norda) :: ri
      real(dp), parameter :: one = 1.d0

      deriv_psianl = 0.0d0
      if(ijas.eq.4.and.rri.eq.asymp_r) return
! Note: This routine is only called with iwf=1, but parts of it are
! written for general iwf, whereas others (asymp_r) assume iwf=1.
      if(nwftypejas.gt.1) iwf=iwfjas
      if(ijas.eq.1.and.rri.gt.cutjas_en(it,iwf))  return


      ri(1)=rri
      if(ijas.ge.4) then
        deriv_psianl=a4(1,it,iwf)*rri/(one+a4(2,it,iwf)*rri)-asymp_jasa(it,iwf)
        do i=2,norda
           ri(i)=ri(1)*ri(i-1)
           deriv_psianl=deriv_psianl+a4(i+1,it,iwf)*ri(i)
        enddo
        do jparm=1,nparma(it)
            if(iwjasa(jparm,it).eq.1) then
              top=rri
              bot=one+a4(2,it,iwf)*rri
              gen=top/bot-asymp_r/(1+a4(2,it,iwf)*asymp_r)
             elseif(iwjasa(jparm,it).eq.2) then
              top=-a4(1,it,iwf)*ri(2)
              bot=one+a4(2,it,iwf)*rri
              bot=bot*bot
              gen=top/bot+a4(1,it,iwf)*asymp_r**2/(1+a4(2,it,iwf)*asymp_r)**2
             else
              iord=iwjasa(jparm,it)-1
              gen=ri(iord)-asymp_r**iord
            endif
            gn(jparm)=gn(jparm)+gen
         enddo
      elseif(ijas.eq.1) then
         xi=ri(1)*cutjas_eni(it,iwf)
         term=(1.d0-xi)**3
         a1_cusp=3.d0*a4(1,it,iwf)*cutjas_eni(it,iwf)
         deriv_psianl=a1_cusp*ri(1)+a4(1,it,iwf)
         do i=2,norda
            ri(i)=ri(1)*ri(i-1)
            deriv_psianl=deriv_psianl+a4(i,it,iwf)*ri(i)
         enddo
         deriv_psianl=deriv_psianl*term
         do jparm=1,nparma(it)
            if(iwjasa(jparm,it).eq.1) then
               da1_cusp=3.d0*cutjas_eni(it,iwf)
               gen=da1_cusp*ri(1)+1.d0
            else
               iord=iwjasa(jparm,it)
               gen=ri(iord)
            endif
            gn(jparm)=gn(jparm)+gen*term
         enddo
      endif

      return
      end

!-----------------------------------------------------------------------
      function deriv_psibnl(u,gn,isb,ipar,iwfjas)


      use jastrow, only: nordb
      use jaspar6, only: asymp_r
      use jastrow, only: cutjas_ee,cutjas_eei
      use jastrow, only: asymp_jasb,b,ijas,sspinn
      use multiple_geo, only: iwf
      use optwf_nparmj, only: nparmb
      use optwf_wjas, only: iwjasb
      use precision_kinds, only: dp
      use vmc_mod, only: nwftypejas

      implicit none

      integer :: i, iord, ipar, isb, jparm, iwfjas
      real(dp) :: a1_cusp, b1_cusp, bot, deriv_psibnl, fee, gee, term, top
      real(dp) :: u, xij
      real(dp), dimension(*) :: gn
      real(dp), dimension(nordb) :: rij
      real(dp), parameter :: one = 1.d0


      deriv_psibnl=0
      if(ijas.eq.4.and.u.eq.asymp_r) return

!     Note: This routine is only called with iwf=1, but parts of it are
!     written for general iwf, whereas others (asymp_r) assume iwf=1.
      if(nwftypejas.gt.1) iwf=iwfjas
      if(ijas.eq.1.and.u.gt.cutjas_ee(isb,iwf))  return

      rij(1)=u

      if(ijas.eq.4) then

         fee=b(1,isb,iwf)*u/(one+b(2,isb,iwf)*u)

         deriv_psibnl=sspinn*fee-asymp_jasb(ipar+1,iwf)
         do i=2,nordb
            rij(i)=rij(1)*rij(i-1)
            deriv_psibnl=deriv_psibnl+b(i+1,isb,iwf)*rij(i)
         enddo


         do jparm=1,nparmb(isb)
            if(iwjasb(jparm,isb).eq.1) then
               top=u
               bot=one+b(2,isb,iwf)*u
               gee=sspinn*(top/bot-asymp_r/(1+b(2,isb,iwf)*asymp_r))
            elseif(iwjasb(jparm,isb).eq.2) then
               top=-b(1,isb,iwf)*rij(2)
               bot=one+b(2,isb,iwf)*u
               bot=bot*bot
               gee=sspinn*(top/bot+b(1,isb,iwf)*asymp_r**2/(1+b(2,isb,iwf)*asymp_r)**2)
            else
               iord=iwjasb(jparm,isb)-1
               gee=rij(iord)-asymp_r**iord
            endif
            gn(jparm)=gn(jparm)+gee
         enddo

      elseif(ijas.eq.1) then

         xij=rij(1)*cutjas_eei(isb,iwf)
         term=(1.d0-xij)**3

         b1_cusp=sspinn*0.5+3.d0*b(1,isb,iwf)*cutjas_eei(isb,iwf)
         deriv_psibnl=b1_cusp*rij(1)+b(1,isb,iwf)
         do i=2,nordb
            rij(i)=rij(1)*rij(i-1)
            deriv_psibnl=deriv_psibnl+b(i,isb,iwf)*rij(i)
         enddo
         deriv_psibnl=deriv_psibnl*term

         do jparm=1,nparmb(isb)
            if(iwjasb(jparm,isb).eq.1) then
               gee=3.d0*xij+1.d0
            else
               iord=iwjasb(jparm,isb)
               gee=rij(iord)
            endif
            gn(jparm)=gn(jparm)+gee*term
         enddo

      endif


      return
      end
end module
