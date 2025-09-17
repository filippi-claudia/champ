      module deriv_jastrow4_mod
      contains
      subroutine deriv_jastrow4(x,fjo,d2o,fsumo,fso,fijo,d2ijo,g,go,d2g,gvalue)
! Written by Cyrus Umrigar and Claudia Filippi
      use bparm,   only: nocuspb,nspin2b
      use contrl_file, only: ounit
      use cuspmat4, only: d,iwc4
      use da_jastrow, only: da_d2j,da_j,da_vj
      use distance_mod, only: r_ee,r_en,rvec_ee,rvec_en
      use ijasnonlin, only: d1d2a,d1d2b,d2d2a,d2d2b
      use jastrow, only: norda,nordb,nordc
      use jastrow, only: asymp_r
      use jaspointer, only: npoint,npointa
      use jastrow, only: a4,asymp_jasa,asymp_jasb,b,c,nordj
      use jastrow, only: sspinn
      use jastrow4_mod, only: da_jastrow4_een, da_jastrow4_en
      use m_force_analytic, only: iforce_analy
      use multiple_geo, only: iwf
      use optwf_control, only: ioptjas
      use optwf_nparmj, only: nparma,nparmb
      use optwf_parms, only: nparmj
      use optwf_wjas, only: iwjasa,iwjasb,iwjasc
      use precision_kinds, only: dp
      use scale_dist_mod, only: scale_dist2,switch_scale2, scale_dist3, switch_scale3
      use system,  only: iwctype,ncent,nctype,nelec,nup
      use vardep,  only: cdep,iwdepend,nvdepend
      implicit none

      integer :: i, ic, id, ideriv, ij
      integer :: im1, iord, ipar, iparm
      integer :: iparm0, iparma, isb, it
      integer :: j, jj, jparm, k
      integer :: l, l_hi, ll, m, n 
      real(dp) :: bot, bot0, bot2, boti, botii
      real(dp) :: botu, botuu, cd, d2o
      real(dp) :: dd1, dd2, dd3, dd7, dd8, dd9, dd10, dd11, dd12
      real(dp) :: fc, fee, feeu, feeuu
      real(dp) :: fen, feni, fenii, feniii
      real(dp) :: fi, fii, fiii, fiij, fij, fj, fjj, fjji, fjjj, fu, fui, fuii, fuij, fuj, fujj, fuu, fuui, fuuj
      real(dp) :: gee=0, geeu=0, geeuu=0, gen
      real(dp) :: geni, genii, gi, gii
      real(dp) :: gj, gjj, gp, gu
      real(dp) :: gui, guj, guu, pc
      real(dp) :: pii, pj, pjj, ppi
      real(dp) :: pu, pui, puj, puu
      real(dp) :: ri, rij, rj, s
      real(dp) :: t, top, topi, topii
      real(dp) :: topu, topuu, u2mst, u2pst
      real(dp) :: value
      real(dp), dimension(3, *) :: x
      real(dp), dimension(-3:nordj) :: uu, ss, tt, rri, rrj
      real(dp) :: fsumo
      real(dp), dimension(3, *) :: fjo
      real(dp), dimension(nelec, *) :: fso
      real(dp), dimension(3, nelec, *) :: fijo
      real(dp), dimension(nelec, *) :: d2ijo
      real(dp), dimension(*) :: gvalue
      real(dp), dimension(*) :: d2g
      real(dp), dimension(3, nelec, *) :: g
      real(dp), dimension(nelec, nelec, *) :: go

      real(dp), parameter :: zero = 0.d0
      real(dp), parameter :: one = 1.d0
      real(dp), parameter :: two = 2.d0
      real(dp), parameter :: half = .5d0
      real(dp), parameter :: eps = 1.d-12

      iparma=nparma(1)
      do it=2,nctype
       iparma=iparma+nparma(it)
      enddo

      fsumo=0.d0
      d2o=0.d0
      do i=1,nelec
        fjo(1,i)=0.d0
        fjo(2,i)=0.d0
        fjo(3,i)=0.d0
      enddo

      if(iforce_analy.gt.0) then
        da_j=0.d0
        da_d2j=0.d0
        da_vj=0.d0
      endif

      do i=-3,-1
        uu(i)=0
        ss(i)=0
        tt(i)=0
        rri(i)=0
        rrj(i)=0
      enddo
      uu(0)=1
      ss(0)=2
      tt(0)=1
      rri(0)=1
      rrj(0)=1

      do iparm=1,nparmj
        gvalue(iparm)=0
        d2g(iparm)=0
        do i=1,nelec
          g(1,i,iparm)=0
          g(2,i,iparm)=0
          g(3,i,iparm)=0
        enddo
      enddo
      do isb=1,nspin2b
        d1d2b(isb)=0
        d2d2b(isb)=0
      enddo
      do it=1,nctype
        d1d2a(it)=0
        d2d2a(it)=0
      enddo

      if (nelec.lt.2) goto 65

! e-e and e-e-n terms
      ij=0
      do i=2,nelec
      im1=i-1
      do j=1,im1
      ij=ij+1

      fso(i,j)=0
      d2ijo(i,j)=0
      do k=1,3
        fijo(k,i,j)=0
        fijo(k,j,i)=0
      enddo
      do iparm=1,nparmj
        go(i,j,iparm)=0
      enddo

      sspinn=1
      isb=1
      ipar=0
      if(i.le.nup .or. j.gt.nup) then
        if(nspin2b.eq.2) then
          isb=2
         elseif(nocuspb.eq.0) then
          sspinn=half
        endif
        ipar=1
      endif

      rij=r_ee(ij)

      call scale_dist2(rij,uu(1),dd1,dd2)

      top=sspinn*b(1,isb,iwf)*uu(1)
      topu=sspinn*b(1,isb,iwf)
      topuu=0

      bot=1+b(2,isb,iwf)*uu(1)
      botu=b(2,isb,iwf)
      botuu=0
      bot2=bot*bot

      fee=top/bot-asymp_jasb(ipar+1,iwf)
      feeu=topu/bot-botu*top/bot2
      feeuu=topuu-(botuu*top+2*botu*topu)/bot+2*botu**2*top/bot2
      feeuu=feeuu/bot

      do iord=2,nordb
        uu(iord)=uu(1)*uu(iord-1)
        fee=fee+b(iord+1,isb,iwf)*uu(iord)
        feeu=feeu+b(iord+1,isb,iwf)*iord*uu(iord-1)
        feeuu=feeuu+b(iord+1,isb,iwf)*iord*(iord-1)*uu(iord-2)
      enddo

      feeuu=feeuu*dd1*dd1+feeu*dd2
      feeu=feeu*dd1/rij

      fso(i,j)=fso(i,j)+fee
      do k=1,3
        fijo(k,i,j)= fijo(k,i,j) + feeu*rvec_ee(k,ij)
        fijo(k,j,i)= fijo(k,j,i) - feeu*rvec_ee(k,ij)
      enddo
      d2ijo(i,j)=d2ijo(i,j)+2*(feeuu+2*feeu)

! derivatives of wave function wrt b(1),b(2) and rest of b(i)
      iparm0=iparma
      if(isb.eq.2) iparm0=iparm0+nparmb(1)
      do jparm=1,nparmb(isb)
        iparm=iparm0+jparm

        if(iwjasb(jparm,isb).eq.1) then
          top=sspinn*uu(1)
          topu=sspinn
          topuu=zero

          bot=one+b(2,isb,iwf)*uu(1)
          botu=b(2,isb,iwf)
          botuu=zero
          bot2=bot*bot

          gee=top/bot-sspinn*asymp_r/(1+b(2,isb,iwf)*asymp_r)
          geeu=topu/bot-botu*top/bot2
          geeuu=topuu-(botuu*top+2*botu*topu)/bot+2*botu**2*top/bot2
          geeuu=geeuu/bot

         elseif(iwjasb(jparm,isb).eq.2) then
          top=-sspinn*b(1,isb,iwf)*uu(1)*uu(1)
          topu=-sspinn*2*b(1,isb,iwf)*uu(1)
          topuu=-sspinn*2*b(1,isb,iwf)

          bot0=one+b(2,isb,iwf)*uu(1)
          bot=bot0*bot0
          botu=2*bot0*b(2,isb,iwf)
          botuu=2*b(2,isb,iwf)*b(2,isb,iwf)
          bot2=bot*bot

          gee=top/bot+sspinn*b(1,isb,iwf)*asymp_r**2/(1+b(2,isb,iwf)*asymp_r)**2
          geeu=topu/bot-botu*top/bot2
          geeuu=topuu-(botuu*top+2*botu*topu)/bot+2*botu**2*top/bot2
          geeuu=geeuu/bot

          if(ioptjas.gt.0) then
            d2d2b(isb)=d2d2b(isb)-2*top*uu(1)/bot/bot0-2*sspinn*b(1,isb,iwf)*asymp_r**3/(1+b(2,isb,iwf)*asymp_r)**3
            d1d2b(isb)=d1d2b(isb)-sspinn*uu(1)*uu(1)/bot+sspinn*asymp_r**2/(1+b(2,isb,iwf)*asymp_r)**2
          endif
         else
          iord=iwjasb(jparm,isb)-1
          gee=uu(iord)-asymp_r**iord
          geeu=iord*uu(iord-1)
          geeuu=iord*(iord-1)*uu(iord-2)

        endif

        geeuu=geeuu*dd1*dd1+geeu*dd2
        geeu=geeu*dd1/rij

        go(i,j,iparm)=go(i,j,iparm)+gee
        gvalue(iparm)=gvalue(iparm)+gee

        g(1,i,iparm)=g(1,i,iparm) + geeu*rvec_ee(1,ij)
        g(2,i,iparm)=g(2,i,iparm) + geeu*rvec_ee(2,ij)
        g(3,i,iparm)=g(3,i,iparm) + geeu*rvec_ee(3,ij)
        g(1,j,iparm)=g(1,j,iparm) - geeu*rvec_ee(1,ij)
        g(2,j,iparm)=g(2,j,iparm) - geeu*rvec_ee(2,ij)
        g(3,j,iparm)=g(3,j,iparm) - geeu*rvec_ee(3,ij)

        d2g(iparm)=d2g(iparm) + two*(geeuu+two*geeu)

      enddo

! There are no C terms to order 1.
   30 if(nordc.le.1) goto 58

      call switch_scale2(uu(1),dd1,dd2)
      do iord=2,nordc
        uu(iord)=uu(1)*uu(iord-1)
      enddo

      do ic=1,ncent
        it=iwctype(ic)

        ri=r_en(i,ic)
        rj=r_en(j,ic)

        if(iforce_analy.eq.0) then
          call scale_dist2(ri,rri(1),dd7,dd9)
          call scale_dist2(rj,rrj(1),dd8,dd10)
          if(rri(1).eq.asymp_r .and. rrj(1).eq.asymp_r) cycle
          call switch_scale2(rri(1),dd7,dd9)
          call switch_scale2(rrj(1),dd8,dd10)
        else
          call scale_dist3(ri,rri(1),dd7,dd9,dd11)
          call scale_dist3(rj,rrj(1),dd8,dd10,dd12)
          if(rri(1).eq.asymp_r .and. rrj(1).eq.asymp_r) cycle
          call switch_scale3(rri(1),dd7,dd9,dd11)
          call switch_scale3(rrj(1),dd8,dd10,dd12)
       endif

        s=ri+rj
        t=ri-rj
!       u2mt2=rij*rij-t*t
        u2pst=rij*rij+s*t
        u2mst=rij*rij-s*t
!       s2mu2=s*s-rij*rij
!       s2mt2=s*s-t*t

        do iord=1,nordc
          rri(iord)=rri(1)*rri(iord-1)
          rrj(iord)=rrj(1)*rrj(iord-1)
          ss(iord)=rri(iord)+rrj(iord)
          tt(iord)=rri(iord)*rrj(iord)
        enddo

        fc=0
        fu=0
        fui=0
        fuj=0
        fuu=0
        fuui=0
        fuuj=0
        fi=0
        fii=0
        fij=0
        fiii=0
        fiij=0
        fj=0
        fjj=0
        fjjj=0
        fjji=0
        fui=0
        fuii=0
        fuij=0
        fuj=0
        fujj=0
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
                pc=uu(k)*ss(l)*tt(m)
                pu=k*uu(k-1)*ss(l)*tt(m)
                puu=k*(k-1)*uu(k-2)*ss(l)*tt(m)
                ppi=uu(k) &
     &          *((l+m)*rri(l+m-1)*rrj(m)+m*rri(m-1)*rrj(l+m))
                pii=uu(k) &
     &          *((l+m)*(l+m-1)*rri(l+m-2)*rrj(m) &
     &          +m*(m-1)*rri(m-2)*rrj(l+m))
                pj=uu(k) &
     &          *((l+m)*rrj(l+m-1)*rri(m)+m*rrj(m-1)*rri(l+m))
                pjj=uu(k) &
     &          *((l+m)*(l+m-1)*rrj(l+m-2)*rri(m) &
     &          +m*(m-1)*rrj(m-2)*rri(l+m))
                pui=k*uu(k-1) &
     &          *((l+m)*rri(l+m-1)*rrj(m)+m*rri(m-1)*rrj(l+m))
                puj=k*uu(k-1) &
     &          *((l+m)*rrj(l+m-1)*rri(m)+m*rrj(m-1)*rri(l+m))

                fc=fc+c(ll,it,iwf)*pc
                fu=fu+c(ll,it,iwf)*pu
                fuu=fuu+c(ll,it,iwf)*puu
                fi=fi+c(ll,it,iwf)*ppi
                fii=fii+c(ll,it,iwf)*pii
                fj=fj+c(ll,it,iwf)*pj
                fjj=fjj+c(ll,it,iwf)*pjj
                fui=fui+c(ll,it,iwf)*pui
                fuj=fuj+c(ll,it,iwf)*puj
                if(iforce_analy.gt.0) then
                  fuii=fuii+c(ll,it,iwf)*k*uu(k-1) &
                       *((l+m)*(l+m-1)*rri(l+m-2)*rrj(m)+m*(m-1)*rri(m-2)*rrj(l+m))
                  fuij=fuij+c(ll,it,iwf)*k*uu(k-1) &
                       *((l+m)*m*rri(l+m-1)*rrj(m-1)+m*(l+m)*rri(m-1)*rrj(l+m-1))
                  fujj=fujj+c(ll,it,iwf)*k*uu(k-1) &
                       *((l+m)*(l+m-1)*rrj(l+m-2)*rri(m)+m*(m-1)*rrj(m-2)*rri(l+m))
                  fuui=fuui+c(ll,it,iwf)*k*(k-1)*uu(k-2) &
                       *((l+m)*rri(l+m-1)*rrj(m)+m*rri(m-1)*rrj(l+m))
                  fuuj=fuuj+c(ll,it,iwf)*k*(k-1)*uu(k-2) &
                       *((l+m)*rrj(l+m-1)*rri(m)+m*rrj(m-1)*rri(l+m))
                  fij=fij+c(ll,it,iwf)*uu(k) &
                       *((l+m)*m*rri(l+m-1)*rrj(m-1)+m*(l+m)*rri(m-1)*rrj(l+m-1))
                  fiii=fiii+c(ll,it,iwf)*uu(k) &
                       *((l+m)*(l+m-1)*(l+m-2)*rri(l+m-3)*rrj(m) &
                       +m*(m-1)*(m-2)*rri(m-3)*rrj(l+m))
                  fiij=fiij+c(ll,it,iwf)*uu(k) &
                       *((l+m)*(l+m-1)*m*rri(l+m-2)*rrj(m-1) &
                       +m*(m-1)*(l+m)*rri(m-2)*rrj(l+m-1))
                  fjjj=fjjj+c(ll,it,iwf)*uu(k) &
                       *((l+m)*(l+m-1)*(l+m-2)*rrj(l+m-3)*rri(m) &
                       +m*(m-1)*(m-2)*rrj(m-3)*rri(l+m))
                  fjji=fjji+c(ll,it,iwf)*uu(k) &
                       *((l+m)*(l+m-1)*m*rrj(l+m-2)*rri(m-1) &
                       +m*(m-1)*(l+m)*rrj(m-2)*rri(l+m-1))
                endif

! derivatives of wave function wrt c-parameters
!               ideriv=0
!               if(ll.eq.iwc4(jj)) then
!                 if(nvdepend(jj,it).gt.0) then
!                   ideriv=1
!                  else
!                   jj=jj+1
!                 endif
!                elseif(ll.eq.iwjasc(jparm,it)) then
!                 ideriv=2
!               endif

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

                if(ideriv.gt.0) then
                  gp=pc
                  gu=pu
                  guu=puu
                  gi=ppi
                  gii=pii
                  gj=pj
                  gjj=pjj
                  gui=pui
                  guj=puj

                  guu=guu*dd1*dd1+gu*dd2
                  gu=gu*dd1/rij

                  gui=gui*dd1*dd7
                  guj=guj*dd1*dd8

                  gii=gii*dd7*dd7+gi*dd9
                  gjj=gjj*dd8*dd8+gj*dd10
                  gi=gi*dd7/ri
                  gj=gj*dd8/rj

                  if(ideriv.eq.1) then

                  do id=1,nvdepend(jj,it)

                  iparm=npoint(it)+iwdepend(jj,id,it)
                  cd=cdep(jj,id,it)

                  go(i,j,iparm)=go(i,j,iparm)+cd*gp
                  gvalue(iparm)=gvalue(iparm)+cd*gp

                  g(1,i,iparm)=g(1,i,iparm)+cd*(gi*rvec_en(1,i,ic)+gu*rvec_ee(1,ij))
                  g(2,i,iparm)=g(2,i,iparm)+cd*(gi*rvec_en(2,i,ic)+gu*rvec_ee(2,ij))
                  g(3,i,iparm)=g(3,i,iparm)+cd*(gi*rvec_en(3,i,ic)+gu*rvec_ee(3,ij))
                  g(1,j,iparm)=g(1,j,iparm)+cd*(gj*rvec_en(1,j,ic)-gu*rvec_ee(1,ij))
                  g(2,j,iparm)=g(2,j,iparm)+cd*(gj*rvec_en(2,j,ic)-gu*rvec_ee(2,ij))
                  g(3,j,iparm)=g(3,j,iparm)+cd*(gj*rvec_en(3,j,ic)-gu*rvec_ee(3,ij))

                  d2g(iparm)=d2g(iparm) + cd*(2*(guu + 2*gu) &
     &            + gui*u2pst/(ri*rij) + guj*u2mst/(rj*rij) &
     &            + gii + 2*gi + gjj + 2*gj)
                  enddo

!                 jj=jj+1

                  elseif(ideriv.eq.2) then

                  iparm=npoint(it)+jparm

                  go(i,j,iparm)=go(i,j,iparm)+gp
                  gvalue(iparm)=gvalue(iparm)+gp

                  g(1,i,iparm)=g(1,i,iparm)+gi*rvec_en(1,i,ic)+gu*rvec_ee(1,ij)
                  g(2,i,iparm)=g(2,i,iparm)+gi*rvec_en(2,i,ic)+gu*rvec_ee(2,ij)
                  g(3,i,iparm)=g(3,i,iparm)+gi*rvec_en(3,i,ic)+gu*rvec_ee(3,ij)
                  g(1,j,iparm)=g(1,j,iparm)+gj*rvec_en(1,j,ic)-gu*rvec_ee(1,ij)
                  g(2,j,iparm)=g(2,j,iparm)+gj*rvec_en(2,j,ic)-gu*rvec_ee(2,ij)
                  g(3,j,iparm)=g(3,j,iparm)+gj*rvec_en(3,j,ic)-gu*rvec_ee(3,ij)

                  d2g(iparm)=d2g(iparm) + 2*(guu + 2*gu) &
     &            + gui*u2pst/(ri*rij) + guj*u2mst/(rj*rij) &
     &            + gii + 2*gi + gjj + 2*gj

                  jparm=jparm+1
                  endif

                endif
              endif
            enddo
          enddo
        enddo

        if(iforce_analy.eq.1) call da_jastrow4_een(i,j,ic,rvec_en(1,1,ic),rvec_ee(1,ij),ri,rj,rij,u2pst,u2mst, &
          fi,fii,fiii,fij,fiij,fj,fjj,fjji,fjjj,fui,fuii,fuij,fuj,fujj,fuui,fuuj,dd1,dd2,dd7,dd8,dd9,dd10,dd11,dd12)

        fuu=fuu*dd1*dd1+fu*dd2
        fu=fu*dd1/rij

        fui=fui*dd1*dd7
        fuj=fuj*dd1*dd8

        fii=fii*dd7*dd7+fi*dd9
        fjj=fjj*dd8*dd8+fj*dd10
        fi=fi*dd7/ri
        fj=fj*dd8/rj

        fso(i,j)=fso(i,j) + fc

        fijo(1,i,j)=fijo(1,i,j) + fi*rvec_en(1,i,ic)+fu*rvec_ee(1,ij)
        fijo(2,i,j)=fijo(2,i,j) + fi*rvec_en(2,i,ic)+fu*rvec_ee(2,ij)
        fijo(3,i,j)=fijo(3,i,j) + fi*rvec_en(3,i,ic)+fu*rvec_ee(3,ij)
        fijo(1,j,i)=fijo(1,j,i) + fj*rvec_en(1,j,ic)-fu*rvec_ee(1,ij)
        fijo(2,j,i)=fijo(2,j,i) + fj*rvec_en(2,j,ic)-fu*rvec_ee(2,ij)
        fijo(3,j,i)=fijo(3,j,i) + fj*rvec_en(3,j,ic)-fu*rvec_ee(3,ij)
!       write(ounit,'(''i,j,fijo2='',2i5,9d12.4)') i,j,(fijo(k,i,j),k=1,3)

        d2ijo(i,j)=d2ijo(i,j) + 2*(fuu + 2*fu) + fui*u2pst/(ri*rij) &
     &  + fuj*u2mst/(rj*rij) + fii + 2*fi + fjj + 2*fj

   57 continue
      enddo

   58 fsumo=fsumo+fso(i,j)
      fjo(1,i)=fjo(1,i)+fijo(1,i,j)
      fjo(2,i)=fjo(2,i)+fijo(2,i,j)
      fjo(3,i)=fjo(3,i)+fijo(3,i,j)
      fjo(1,j)=fjo(1,j)+fijo(1,j,i)
      fjo(2,j)=fjo(2,j)+fijo(2,j,i)
      fjo(3,j)=fjo(3,j)+fijo(3,j,i)
!     div_vj(i)=div_vj(i)+d2ijo(i,j)/2
!     div_vj(j)=div_vj(j)+d2ijo(i,j)/2
      d2o=d2o+d2ijo(i,j)
      enddo
      enddo

! e-n terms
   65 do i=1,nelec

        fso(i,i)=0
        fijo(1,i,i)=0
        fijo(2,i,i)=0
        fijo(3,i,i)=0
        d2ijo(i,i)=0

        do iparm=1,nparmj
          go(i,i,iparm)=zero
        enddo

        do ic=1,ncent
          it=iwctype(ic)

          ri=r_en(i,ic)

          if(iforce_analy.eq.0) then
            call scale_dist2(ri,rri(1),dd7,dd9)
          else
            call scale_dist3(ri,rri(1),dd7,dd9,dd11)
          endif

          top=a4(1,it,iwf)*rri(1)
          topi=a4(1,it,iwf)
          topii=0

          bot=1+a4(2,it,iwf)*rri(1)
          boti=a4(2,it,iwf)
          botii=0
          bot2=bot*bot

          fen=top/bot-asymp_jasa(it,iwf)
          feni=topi/bot-boti*top/bot2
          fenii=topii-(botii*top+2*boti*topi)/bot+2*boti**2*top/bot2
          fenii=fenii/bot

          do iord=2,norda
            rri(iord)=rri(1)*rri(iord-1)
            fen=fen+a4(iord+1,it,iwf)*rri(iord)
            feni=feni+a4(iord+1,it,iwf)*iord*rri(iord-1)
            fenii=fenii+a4(iord+1,it,iwf)*iord*(iord-1)*rri(iord-2)
          enddo

          if(iforce_analy.eq.1) then
! temporarely assuming no Pade' term in A Jastrow
            feniii=0.d0
            do iord=3,norda
              feniii=feniii+a4(iord+1,it,iwf)*iord*(iord-1)*(iord-2)*rri(iord-3)
            enddo
            call da_jastrow4_en(i,ic,rvec_en(1,i,ic),ri,feni,fenii,feniii,dd7,dd9,dd11)
          endif

          fenii=fenii*dd7*dd7+feni*dd9
          feni=feni*dd7/ri

          fso(i,i)=fso(i,i)+fen

          fijo(1,i,i)=fijo(1,i,i) + feni*rvec_en(1,i,ic)
          fijo(2,i,i)=fijo(2,i,i) + feni*rvec_en(2,i,ic)
          fijo(3,i,i)=fijo(3,i,i) + feni*rvec_en(3,i,ic)
!         write(ounit,'(''fijo='',9d12.4)') (fijo(k,i,i),k=1,3),feni,rvec_en(1,i,ic)

          d2ijo(i,i) = d2ijo(i,i) + fenii + 2*feni

          do jparm=1,nparma(it)
            iparm=npointa(it)+jparm

            if(iwjasa(jparm,it).eq.1) then
              top=rri(1)
              topi=one
              topii=zero

              bot=one+a4(2,it,iwf)*rri(1)
              boti=a4(2,it,iwf)
              botii=zero
              bot2=bot*bot

              gen=top/bot-asymp_r/(1+a4(2,it,iwf)*asymp_r)
              geni=topi/bot-boti*top/bot2
              genii=topii-(botii*top+2*boti*topi)/bot+2*boti**2*top/bot2
              genii=genii/bot

             elseif(iwjasa(jparm,it).eq.2) then
              top=-a4(1,it,iwf)*rri(1)*rri(1)
              topi=-2*a4(1,it,iwf)*rri(1)
              topii=-2*a4(1,it,iwf)

              bot0=one+a4(2,it,iwf)*rri(1)
              bot=bot0*bot0
              boti=2*bot0*a4(2,it,iwf)
              botii=2*a4(2,it,iwf)*a4(2,it,iwf)
              bot2=bot*bot

              gen=top/bot+a4(1,it,iwf)*asymp_r**2/(1+a4(2,it,iwf)*asymp_r)**2
              geni=topi/bot-boti*top/bot2
              genii=topii-(botii*top+2*boti*topi)/bot+2*boti**2*top/bot2
              genii=genii/bot

              if(ioptjas.gt.0) then
                d2d2a(it)=d2d2a(it)-2*top*rri(1)/bot/bot0-2*a4(1,it,iwf)*asymp_r**3/(1+a4(2,it,iwf)*asymp_r)**3
                d1d2a(it)=d1d2a(it)-rri(1)*rri(1)/bot+asymp_r**2/(1+a4(2,it,iwf)*asymp_r)**2
              endif
             else
              iord=iwjasa(jparm,it)-1
              gen=rri(iord)-asymp_r**iord
              geni=iord*rri(iord-1)
              genii=iord*(iord-1)*rri(iord-2)

            endif
!           write(ounit,*) 'CIAO',iord,rri(iord),genii,d77,geni,dd9

            genii=genii*dd7*dd7+geni*dd9
            geni=geni*dd7/r_en(i,ic)

            go(i,i,iparm)=go(i,i,iparm)+gen
            gvalue(iparm)=gvalue(iparm)+gen

            g(1,i,iparm)=g(1,i,iparm)+geni*rvec_en(1,i,ic)
            g(2,i,iparm)=g(2,i,iparm)+geni*rvec_en(2,i,ic)
            g(3,i,iparm)=g(3,i,iparm)+geni*rvec_en(3,i,ic)

            d2g(iparm)=d2g(iparm)+genii+two*geni
          enddo
   80     continue
        enddo

        fsumo=fsumo+fso(i,i)
        fjo(1,i)=fjo(1,i)+fijo(1,i,i)
        fjo(2,i)=fjo(2,i)+fijo(2,i,i)
        fjo(3,i)=fjo(3,i)+fijo(3,i,i)
!       write(ounit,'(''v='',9d12.4)') (v(k,i),k=1,3)
!       div_vj(i)=div_vj(i)+d2ijo(i,i)
        d2o=d2o+d2ijo(i,i)
      enddo

!     write(ounit,*) (d2g(iparm),iparm=1,nparmj)

!     write(ounit,*) 'd2d2',nspin2b,d2d2b(1),d2d2b(2),asymp_r
!     write(ounit,*) 'asym',asymp_r,asymp_jasa(1,1),asymp_jasb(1,1)

      return
      end
      end module

