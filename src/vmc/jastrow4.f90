      module jastrow4_mod
      contains
      subroutine jastrow_factor4(x,fjo,d2o,fsumo,fso,fijo,d2ijo)
! Written by Cyrus Umrigar, modified by C. Filippi
      use da_jastrow, only: da_d2j, da_vj
      use system, only: iwctype, ncent, nelec, nup
      use jastrow, only: sspinn, b, c, scalek, a4, norda, nordb, nordc, asymp_jasa, asymp_jasb, nordj
      use multiple_geo, only: iwf
      use bparm, only: nocuspb, nspin2b
      use scale_dist_mod, only: scale_dist2, switch_scale2
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
      real(dp) :: bot, bot2, boti, botii, botu
      real(dp) :: botuu, dd1, dd10
      real(dp) :: dd2, dd7, dd8, dd9
      real(dp) :: fc, fee, feeu, feeuu
      real(dp) :: fen, feni, feni_save, fenii
      real(dp) :: fenii_save, fi, fii, fj
      real(dp) :: fjj, fu, fui
      real(dp) :: fuj, fuu, ri, rij
      real(dp) :: rj, s, t, term
      real(dp) :: top, topi, topii, topu
      real(dp) :: topuu, u2mst, u2pst
      real(dp), dimension(3, *) :: x
      real(dp), dimension(-2:nordj) :: uu, ss, tt, rri, rrj
      !real(dp), dimension(-2:nordj) :: ss
      !real(dp), dimension(-2:nordj) :: tt
      !real(dp), dimension(-2:nordj) :: rri
      !real(dp), dimension(-2:nordj) :: rrj
      real(dp), parameter :: half = .5d0
      real(dp), parameter :: eps = 1.d-12
! replace global variables of fsumo, fjo, fso, fijo, and d2ijo
! with locals of the same name, with one less dimension
      real(dp) :: fsumo, d2o
      real(dp), dimension(3, *) :: fjo
      real(dp), dimension(nelec, *) :: fso
      real(dp), dimension(3, nelec, *) :: fijo
      real(dp), dimension(nelec, *) :: d2ijo


      fsumo=0.d0
      d2o=0.d0
      do j=1,nelec
        do k=1,3
          fjo(k,j)=0
        enddo
      enddo
      if(iforce_analy.gt.0) then
        da_d2j=0.d0
        da_vj=0.d0
      endif

      do i=-2,-1
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

        call scale_dist2(ri,rri(1),dd7,dd9)
        call scale_dist2(rj,rrj(1),dd8,dd10)

        call switch_scale2(rri(1),dd7,dd9)
        call switch_scale2(rrj(1),dd8,dd10)
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
        fuu=0
        fi=0
        fii=0
        fj=0
        fjj=0
        fui=0
        fuj=0
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
              endif
!        write(ounit,'(''rij,ri,rj'',9f10.5)') rij,ri,rj,uu(1),rri(1),rrj(1)
            enddo
          enddo
        enddo

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

          call scale_dist2(ri,rri(1),dd7,dd9)
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

          feni_save=feni
          fenii_save=fenii
          fenii=fenii*dd7*dd7+feni*dd9
          feni=feni*dd7/ri

          fso(i,i)=fso(i,i)+fen

          fijo(1,i,i)=fijo(1,i,i) + feni*rvec_en(1,i,ic)
          fijo(2,i,i)=fijo(2,i,i) + feni*rvec_en(2,i,ic)
          fijo(3,i,i)=fijo(3,i,i) + feni*rvec_en(3,i,ic)
!          write(ounit,'(''fijo='',9d12.4)') (fijo(k,i,i),k=1,3),feni,rvec_en(1,i,ic)

          d2ijo(i,i) = d2ijo(i,i) + fenii + 2*feni

          if(iforce_analy.eq.1) call da_jastrow4(iwf,i,ic,it,rvec_en(1,i,ic),ri,rri,feni_save,fenii_save,dd7,dd9)
      80   continue
        enddo

        fsumo=fsumo+fso(i,i)
        fjo(1,i)=fjo(1,i)+fijo(1,i,i)
        fjo(2,i)=fjo(2,i)+fijo(2,i,i)
        fjo(3,i)=fjo(3,i)+fijo(3,i,i)
!        write(ounit,'(''v='',9d12.4)') (fjo(k,i),k=1,3)
        d2o=d2o+d2ijo(i,i)
      enddo

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
      subroutine da_jastrow4(iwf,i,ic,it,rvec_en,r,rr,feni,fenii,dd1,dd2)

      use da_jastrow, only: da_d2j, da_j, da_vj
      use jastrow, only: a4, norda, nordj
      use scale_more, only: dd3
      use precision_kinds, only: dp
      use scale_more, only: dd3

      implicit none

      integer :: i, ic, iord, it, iwf
      integer :: k, l
      real(dp) :: dd1, dd2, feni, fenii, feniii
      real(dp) :: r, ri, ri2
      real(dp), dimension(3) :: rvec_en
      real(dp), dimension(-2:nordj) :: rr

      feniii=0.d0
      do iord=3,norda
        feniii=feniii+a4(iord+1,it,iwf)*iord*(iord-1)*(iord-2)*rr(iord-3)
      enddo

      ri=1.d0/r
      ri2=ri*ri

! can use iwf and loop over nwftypejas, since iwf is passed, this
! will only be called with 1 jas I believe, using iwf's original
! purpose

      do k=1,3
        da_j(k,i,ic)=-rvec_en(k)*ri*feni*dd1
        da_d2j(k,ic)=da_d2j(k,ic)-rvec_en(k)*ri*(feniii*dd1*dd1*dd1+fenii*dd1*(3*dd2+2*dd1*ri)+feni*(dd3+2*dd2*ri-2*dd1*ri2))
        do l=1,3
          da_vj(k,l,i,ic)=da_vj(k,l,i,ic)-rvec_en(k)*rvec_en(l)*ri2*(fenii*dd1*dd1+feni*dd2-feni*dd1*ri)
        enddo
        da_vj(k,k,i,ic)=da_vj(k,k,i,ic)-feni*dd1*ri
      enddo

      return
      end
end module
