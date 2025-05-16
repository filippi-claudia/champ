      module jastrow4_mod
      contains
      subroutine jastrow_factor4(x,fjo,d2o,fsumo,fso,fijo,d2ijo)
! Written by Cyrus Umrigar, modified by C. Filippi
      use fragments, only: eloc_i, elocfrag, ifragelec, nfrag
      use constants, only: hb
      use contrldmc, only: icut_e
      use da_jastrow, only: da_d2j, da_j, da_vj
      use system, only: iwctype, ncent, nelec, nup
      use jastrow, only: sspinn, b, c, scalek, a4, norda, nordb, nordc, asymp_jasa, asymp_jasb, nordj
      use multiple_geo, only: iwf
      use bparm, only: nocuspb, nspin2b
      use scale_dist_mod, only: scale_dist2, switch_scale2, scale_dist3, switch_scale3
      use m_force_analytic, only: iforce_analy
      use distance_mod, only: r_en, rvec_en, r_ee, rvec_ee
      use precision_kinds, only: dp
      use contrl_file,    only: ounit
      use vmc_mod, only: nwftypejas
      implicit none

      integer :: i, ic, ij, im1, iord
      integer :: ipar, isb, it, j
      integer :: k, l, l_hi, ll
      integer :: m, n
      real(dp) :: bot, bot2, boti, botii, botu, botuu
      real(dp) :: dd1, dd2, dd3, dd7, dd8, dd9, dd10, dd11, dd12
      real(dp) :: fc, fee, feeu, feeuu
      real(dp) :: fen, feni, fenii, feniii
      real(dp) :: fi, fii, fiii, fiij, fij, fj, fjj, fjji, fjjj, fu, fui, fuii, fuij, fuj, fujj, fuu, fuui, fuuj
      real(dp) :: ri, rij, rj, s, t, term
      real(dp) :: top, topi, topii, topu
      real(dp) :: topuu, u2mst, u2pst
      real(dp), dimension(3, *) :: x
      real(dp), dimension(-3:nordj) :: uu, ss, tt, rri, rrj
      real(dp) :: fsumo, d2o
      real(dp), dimension(3, *) :: fjo
      real(dp), dimension(nelec, *) :: fso
      real(dp), dimension(3, nelec, *) :: fijo
      real(dp), dimension(nelec, *) :: d2ijo
      real(dp), parameter :: half = .5d0
      real(dp), parameter :: eps = 1.d-12
      real(dp), dimension(3,ncent) :: tempjas

      fsumo=0.d0
      d2o=0.d0
      do j=1,nelec
        do k=1,3
          fjo(k,j)=0
        enddo
      enddo

      if(iforce_analy.gt.0) then
        da_j=0.d0
        da_vj=0.d0
        da_d2j=0.d0
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
!      write(ounit,'(''rij,u in ee'',2f9.5)') rij,uu(1)

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

! There are no C terms to order 1.
      22 if(nordc.le.1) goto 55

      call switch_scale2(uu(1),dd1,dd2)
      do iord=2,nordc
        uu(iord)=uu(1)*uu(iord-1)
      enddo
!      write(ounit,'(''rij,u in een'',2f12.9)') rij,uu(1)

      do ic=1,ncent
        it=iwctype(ic)

        ri=r_en(i,ic)
        rj=r_en(j,ic)

        if(iforce_analy.eq.0) then
          call scale_dist2(ri,rri(1),dd7,dd9)
          call scale_dist2(rj,rrj(1),dd8,dd10)

          call switch_scale2(rri(1),dd7,dd9)
          call switch_scale2(rrj(1),dd8,dd10)
        else
          call scale_dist3(ri,rri(1),dd7,dd9,dd11)
          call scale_dist3(rj,rrj(1),dd8,dd10,dd12)

          call switch_scale3(rri(1),dd7,dd9,dd11)
          call switch_scale3(rrj(1),dd8,dd10,dd12)
       endif
!        write(ounit,'(''ri,rri in een'',2f12.9)') ri,rri(1)

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
!               write(ounit,*) 'order u s t',k,l,m
                fc=fc+c(ll,it,iwf)*uu(k)*ss(l)*tt(m)
                fu=fu+c(ll,it,iwf)*k*uu(k-1)*ss(l)*tt(m)
                fuu=fuu+c(ll,it,iwf)*k*(k-1)*uu(k-2)*ss(l)*tt(m)
                fi=fi+c(ll,it,iwf)*uu(k) &
                *((l+m)*rri(l+m-1)*rrj(m)+m*rri(m-1)*rrj(l+m))
                fii=fii+c(ll,it,iwf)*uu(k) &
                *((l+m)*(l+m-1)*rri(l+m-2)*rrj(m) &
                +m*(m-1)*rri(m-2)*rrj(l+m))
                fj=fj+c(ll,it,iwf)*uu(k) &
                *((l+m)*rrj(l+m-1)*rri(m)+m*rrj(m-1)*rri(l+m))
                fjj=fjj+c(ll,it,iwf)*uu(k) &
                *((l+m)*(l+m-1)*rrj(l+m-2)*rri(m) &
                +m*(m-1)*rrj(m-2)*rri(l+m))
                fui=fui+c(ll,it,iwf)*k*uu(k-1) &
                *((l+m)*rri(l+m-1)*rrj(m)+m*rri(m-1)*rrj(l+m))
                fuj=fuj+c(ll,it,iwf)*k*uu(k-1) &
                *((l+m)*rrj(l+m-1)*rri(m)+m*rrj(m-1)*rri(l+m))
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
              endif
!        write(ounit,'(''rij,ri,rj'',9f10.5)') rij,ri,rj,uu(1),rri(1),rrj(1)
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
!        write(ounit,'(''i,j,fijo2='',2i5,9d12.4)') i,j,(fijo(k,i,j),k=1,3)

        d2ijo(i,j)=d2ijo(i,j) + 2*(fuu + 2*fu) + fui*u2pst/(ri*rij) &
        + fuj*u2mst/(rj*rij) + fii + 2*fi + fjj + 2*fj

      50 continue
      enddo

      55 fsumo=fsumo+fso(i,j)
      fjo(1,i)=fjo(1,i)+fijo(1,i,j)
      fjo(2,i)=fjo(2,i)+fijo(2,i,j)
      fjo(3,i)=fjo(3,i)+fijo(3,i,j)
      fjo(1,j)=fjo(1,j)+fijo(1,j,i)
      fjo(2,j)=fjo(2,j)+fijo(2,j,i)
      fjo(3,j)=fjo(3,j)+fijo(3,j,i)
      d2o=d2o+d2ijo(i,j)
      
      if (icut_e.lt.0) then
        eloc_i(i)=eloc_i(i)-hb*0.5d0*d2ijo(i,j)
        eloc_i(j)=eloc_i(j)-hb*0.5d0*d2ijo(i,j)
      endif
      if (nfrag.gt.1) then
        elocfrag(ifragelec(i)) = elocfrag(ifragelec(i)) - hb*0.5d0*d2ijo(i,j)
        elocfrag(ifragelec(j)) = elocfrag(ifragelec(j)) - hb*0.5d0*d2ijo(i,j)
      endif  
      enddo
      enddo

! e-n terms
      65 do i=1,nelec

        fso(i,i)=0
        fijo(1,i,i)=0
        fijo(2,i,i)=0
        fijo(3,i,i)=0
        d2ijo(i,i)=0

        do ic=1,ncent
          it=iwctype(ic)

          ri=r_en(i,ic)

          if(iforce_analy.eq.0) then
            call scale_dist2(ri,rri(1),dd7,dd9)
           else
            call scale_dist3(ri,rri(1),dd7,dd9,dd11)
          endif
!          write(ounit,'(''ri,rri in en'',2f9.5)') ri,rri(1)

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
            rri(iord)=rri(1)**iord
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
!          write(ounit,'(''fijo='',9d12.4)') (fijo(k,i,i),k=1,3),feni,rvec_en(1,i,ic)

          d2ijo(i,i) = d2ijo(i,i) + fenii + 2*feni

      80   continue
        enddo

        fsumo=fsumo+fso(i,i)
        fjo(1,i)=fjo(1,i)+fijo(1,i,i)
        fjo(2,i)=fjo(2,i)+fijo(2,i,i)
        fjo(3,i)=fjo(3,i)+fijo(3,i,i)
!        write(ounit,'(''v='',9d12.4)') (fjo(k,i),k=1,3)
        d2o=d2o+d2ijo(i,i)
        if (icut_e.lt.0) then
          eloc_i(i) = eloc_i(i) - hb * d2ijo(i,i)
        endif 
        if (nfrag.gt.1) then
          elocfrag(ifragelec(i)) = elocfrag(ifragelec(i)) - hb * d2ijo(i,i)
        endif
      enddo

  !     tempjas = 0.d0
  !     do ic=1,ncent
  !     do k = 1,3
  !     do i = 1, nelec
  !       do j =1,i
  !         tempjas(k,ic) = tempjas(k,ic) + da_j(k,i,j,ic)
  !       enddo
  !     enddo
  !   enddo
  ! enddo

      ! do ic=1,ncent
      !   write(ounit,'(''da_j,ic='',i4,1000f15.11)') ic, (tempjas(k,ic),k=1,3)
      !   do i=1,1
      !    write(ounit,'(''da_vj,ic,i='',2i4,1000f15.11)') ic, i, ((da_vj(k,l,i,ic),k=1,3),l=1,3)
      !   enddo
      !   write(ounit,'(''da_d2j,ic='',i4,1000f15.11)') ic, (da_d2j(k,ic),k=1,3)
      ! enddo

      return
      end

!-----------------------------------------------------------------------
      function nterms4(nord)
! Written by Cyrus Umrigar

      implicit none

      integer :: i, k, l, l_hi, m
      integer :: n, nord, nterms4

      i=0
      do n=2,nord
        do k=n-1,0,-1
          if(k.eq.0) then
            l_hi=n-k-2
           else
            l_hi=n-k
          endif
          do l=l_hi,0,-1
            m=(n-k-l)/2
            if(2*m.eq.n-k-l) then
              i=i+1
            endif
          enddo
        enddo
      enddo
      nterms4=i
!     write(ounit,'(''nterms4='',i5)') nterms4
      return
      end
!-----------------------------------------------------------------------
      subroutine da_jastrow4_en(i,ic,rvec_en,r,feni,fenii,feniii,dd1,dd2,dd3)

      use da_jastrow, only: da_d2j, da_j, da_vj
      use contrl_file, only: ounit
      use precision_kinds, only: dp

      implicit none

      integer :: i, ic, k, l
      real(dp) :: dd1, dd2, dd3, feni, fenii, feniii
      real(dp) :: r, ri, ri2
      real(dp), dimension(3) :: rvec_en

      ri=1.d0/r
      ri2=ri*ri

! can use iwf and loop over nwftypejas, since iwf is passed, this
! will only be called with 1 jas I believe, using iwf's original
! purpose

      do k=1,3
        da_j(k,i,i,ic)=da_j(k,i,i,ic)-rvec_en(k)*ri*feni*dd1
        da_d2j(k,ic)=da_d2j(k,ic)-rvec_en(k)*ri*(feniii*dd1*dd1*dd1+fenii*dd1*(3*dd2+2*dd1*ri)+feni*(dd3+2*dd2*ri-2*dd1*ri2))
        do l=1,3
          da_vj(k,l,i,ic)=da_vj(k,l,i,ic)-rvec_en(k)*rvec_en(l)*ri2*(fenii*dd1*dd1+feni*dd2-feni*dd1*ri)
        enddo
        da_vj(k,k,i,ic)=da_vj(k,k,i,ic)-feni*dd1*ri
!       write(ounit,*) 'CIAO',k,ic,da_d2j(k,ic)
      enddo

      return
      end
!-----------------------------------------------------------------------
      subroutine da_jastrow4_een(i,j,ic,rvec_en,rvec_ee,ri,rj,rij,u2pst,u2mst, &
          fi,fii,fiii,fij,fiij,fj,fjj,fjji,fjjj,fui,fuii,fuij,fuj,fujj,fuui,fuuj,dd1,dd2,dd7,dd8,dd9,dd10,dd11,dd12)
               
      use da_jastrow, only: da_d2j, da_j, da_vj
      use contrl_file,    only: ounit
      use contrl_file,    only: ounit
      use precision_kinds, only: dp
      use system, only: nelec

      implicit none

      integer :: i, ic, j, k, l
      real(dp) :: dd1,dd2,dd7,dd8,dd9,dd10,dd11,dd12
      real(dp) :: fi,fii,fiii,fij,fiij,fj,fjj,fjji,fjjj,fui,fuii,fuij,fuj,fujj,fuui,fuuj
      real(dp) :: ri, ri_i, ri_i2, rj, rj_i, rj_i2, rij, rij_i,u2pst,u2mst
      real(dp) :: dum,dum1i,dum1j,dum1ij,dum2i,dum2j
      real(dp), dimension(3,nelec) :: rvec_en
      real(dp), dimension(3) :: rvec_ee

      ri_i=1.d0/ri
      rj_i=1.d0/rj
      rij_i=1.d0/rij
      ri_i2=ri_i*ri_i
      rj_i2=rj_i*rj_i

      !fi*dd7*ri_i*rvec_en(1,i,ic)+fu*dd1*rij_i*rvec_ee(1,ij)
      !fj*dd8*rj_i*rvec_en(1,j,ic)-fu*dd1*rij_i*rvec_ee(1,ij)

      !d2ijo(i,j)=d2ijo(i,j) + 2*(fuu*dd1*dd1 + fu*dd2 + 2*fu*dd1/rij) + fui*dd1*dd7*(rij2+(ri2-rj2))/(ri*rij) &
      !+ fuj*dd1*dd8*(rij2-(ri2-rj2))/(rj*rij) + fii*dd7*dd7 + fi*dd9 + 2*fi*dd7/ri + fjj*dd8*dd8 + fj*dd10 + 2*fj*dd8/rj

      dum1i=ri_i2*(fii*dd7*dd7+fi*dd9-fi*dd7*ri_i)
      dum1j=rj_i2*(fjj*dd8*dd8+fj*dd10-fj*dd8*rj_i)
      dum1ij=rj_i*ri_i*fij*dd7*dd8
      dum2i=ri_i*(fiii*dd7*dd7*dd7+fii*dd7*(3*dd9+2*dd7*ri_i) +fi*(dd11+2*dd9*ri_i-2*dd7*ri_i2) &
                                                   +dd7*(fjji*dd8*dd8+fij*(dd10+2*dd8*rj_i)) &
                                                   +2*dd7*(fuui*dd1*dd1 + fui*dd2 + 2*fui*dd1*rij_i) &
                                                   +dd1*(fuii*dd7*dd7*u2pst*ri_i &
                                                   +fui*(u2pst*(dd9*ri_i-dd7*ri_i2) + 2*dd7) &
                                                   +dd8*(fuij*dd7*u2mst - 2*fuj*ri)*rj_i)*rij_i )
      dum2j=rj_i*(fjjj*dd8*dd8*dd8+fjj*dd8*(3*dd10+2*dd8*rj_i)+fj*(dd12+2*dd10*rj_i-2*dd8*rj_i2) &
                                                   +dd8*(fiij*dd7*dd7+fij*(dd9+2*dd7*ri_i)) &
                                                   +2*dd8*(fuuj*dd1*dd1 + fuj*dd2 + 2*fuj*dd1*rij_i) &
                                                   +dd1*(fujj*dd8*dd8*u2mst*rj_i &
                                                   +fuj*(u2mst*(dd10*rj_i-dd8*rj_i2)+ 2*dd8) &
                                                   +dd7*(fuij*dd8*u2pst -2*fui*rj)*ri_i)*rij_i)
      do k=1,3
        da_j(k,i,j,ic)=da_j(k,i,j,ic)-rvec_en(k,i)*ri_i*fi*dd7-rvec_en(k,j)*rj_i*fj*dd8
        da_d2j(k,ic)=da_d2j(k,ic)-rvec_en(k,i)*dum2i -rvec_en(k,j)*dum2j

        dum=dd1*rij_i*(rvec_en(k,i)*fui*dd7*ri_i+rvec_en(k,j)*fuj*dd8*rj_i)
        do l=1,3
          da_vj(k,l,i,ic)=da_vj(k,l,i,ic)-rvec_en(l,i)*(rvec_en(k,i)*dum1i+rvec_en(k,j)*dum1ij)&
                                         -rvec_ee(l)*dum

          da_vj(k,l,j,ic)=da_vj(k,l,j,ic)-rvec_en(l,j)*(rvec_en(k,j)*dum1j+rvec_en(k,i)*dum1ij) &
                                         +rvec_ee(l)*dum
        enddo
        da_vj(k,k,i,ic)=da_vj(k,k,i,ic)-fi*dd7*ri_i
        da_vj(k,k,j,ic)=da_vj(k,k,j,ic)-fj*dd8*rj_i

!        da_d2j(k,ic)=da_d2j(k,ic)-rvec_en(k,i)*ri_i*(fiii*dd7*dd7*dd7+fii*dd7*(3*dd9+2*dd7*ri_i) +fi*(dd11+2*dd9*ri_i-2*dd7*ri_i2) &
!                                                   +dd7*(fjji*dd8*dd8+fij*(dd10+2*dd8*rj_i)) &
!                                                   +2*dd7*(fuui*dd1*dd1 + fui*dd2 + 2*fui*dd1*rij_i) &
!                                                   +fuii*dd1*dd7*dd7*(u2pst)*ri_i*rij_i &
!                                                   +fui*dd1*dd9*(u2pst)*ri_i*rij_i &
!                                                   +fui*dd1*dd7*2*rij_i &
!                                                   -fui*dd1*dd7*(u2pst)*ri_i2*rij_i &
!                                                   +fuij*dd1*dd7*dd8*(u2mst)*rj_i*rij_i &
!                                                   -fuj*dd1*dd8*2*ri*rj_i*rij_i) &
!                                 -rvec_en(k,j)*rj_i*(fjjj*dd8*dd8*dd8+fjj*dd8*(3*dd10+2*dd8*rj_i)+fj*(dd12+2*dd10*rj_i-2*dd8*rj_i2) &
!                                                   +dd8*(fiij*dd7*dd7+fij*(dd9+2*dd7*ri_i)) &
!                                                   +2*dd8*(fuuj*dd1*dd1 + fuj*dd2 + 2*fuj*dd1*rij_i) &
!                                                   +fujj*dd1*dd8*dd8*(u2mst)*rj_i*rij_i &
!                                                   +fuj*dd1*dd10*(u2mst)*rj_i*rij_i &
!                                                   +fuj*dd1*dd8*2*rij_i &
!                                                   -fuj*dd1*dd8*(u2mst)*rj_i2*rij_i &
!                                                   +fuij*dd1*dd7*dd8*(u2pst)*ri_i*rij_i &
!                                                   -fui*dd1*dd7*2*rj*ri_i*rij_i)
!       do l=1,3
!         da_vj(k,l,i,ic)=da_vj(k,l,i,ic)-rvec_en(l,i)*rvec_en(k,i)*ri_i2*(fii*dd7*dd7+fi*dd9-fi*dd7*ri_i) &
!                                        -rvec_en(l,i)*rvec_en(k,j)*ri_i*rj_i*fij*dd7*dd8 &
!                                        -rvec_ee(l)*dd1*rij_i*(rvec_en(k,i)*fui*dd7*ri_i+rvec_en(k,j)*fuj*dd8*rj_i)

!         da_vj(k,l,j,ic)=da_vj(k,l,j,ic)-rvec_en(l,j)*rvec_en(k,j)*rj_i2*(fjj*dd8*dd8+fj*dd10-fj*dd8*rj_i) &
!                                        -rvec_en(l,j)*rvec_en(k,i)*rj_i*ri_i*fij*dd7*dd8 &
!                                        +rvec_ee(l)*dd1*rij_i*(rvec_en(k,i)*fui*dd7*ri_i+rvec_en(k,j)*fuj*dd8*rj_i)
!       enddo
!       da_vj(k,k,i,ic)=da_vj(k,k,i,ic)-fi*dd7*ri_i
!       da_vj(k,k,j,ic)=da_vj(k,k,j,ic)-fj*dd8*rj_i
!
!       write(ounit,*) k,i,ic,da_d2j(k,ic),(da_vj(k,l,i,ic),l=1,3)

      enddo

      return
      end

end module
