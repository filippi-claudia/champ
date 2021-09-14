      subroutine jastrow4(x,v,d2,div_vj,value)
c Written by Cyrus Umrigar, modified by C. Filippi
c Jastrow 4,5 must be used with one of isc=2,4,6,7,12,14,16,17
c Jastrow 6   must be used with one of isc=6,7

      use vmc_mod, only: nordj
      use atom, only: iwctype, ncent
      use jaspar, only: sspinn
      use const, only: nelec
      use elec, only: nup
      use jaso, only: d2ijo, d2o, fijo, fjo, fso, fsumo
      use jaspar3, only: b, c, scalek
      use jaspar4, only: a4, norda, nordb, nordc
      use jaspar6, only: asymp_jasa, asymp_jasb, c1_jas6
      use jaspar6, only: cutjas
      use wfsec, only: iwf
      use bparm, only: nocuspb, nspin2b
      use contr2, only: ijas
      use contr2, only: isc
      use force_analy, only: iforce_analy
      use distance_mod, only: rshift, r_en, rvec_en, r_ee, rvec_ee
      use precision_kinds, only: dp
      use contrl_file,    only: ounit
      implicit none

      integer :: i, ic, ij, im1, iord
      integer :: ipar, isb, it, j
      integer :: k, l, l_hi, ll
      integer :: m, n
      real(dp) :: bot, bot2, boti, botii, botu
      real(dp) :: botuu, d2, dd1, dd10
      real(dp) :: dd2, dd7, dd8, dd9
      real(dp) :: fc, fee, feeu, feeuu
      real(dp) :: fen, feni, feni_save, fenii
      real(dp) :: fenii_save, fi, fii, fj
      real(dp) :: fjj, fsum, fu, fui
      real(dp) :: fuj, fuu, ri, rij
      real(dp) :: rj, s, t, term
      real(dp) :: top, topi, topii, topu
      real(dp) :: topuu, u2mst, u2pst, value
      real(dp), dimension(3, *) :: x
      real(dp), dimension(3, *) :: v
      real(dp), dimension(*) :: div_vj
      real(dp), dimension(-2:nordj) :: uu
      real(dp), dimension(-2:nordj) :: ss
      real(dp), dimension(-2:nordj) :: tt
      real(dp), dimension(-2:nordj) :: rri
      real(dp), dimension(-2:nordj) :: rrj
      real(dp), parameter :: half = .5d0
      real(dp), parameter :: eps = 1.d-12

      fsum=0
      do 5 i=-2,-1
        uu(i)=0
        ss(i)=0
        tt(i)=0
        rri(i)=0
    5   rrj(i)=0
      uu(0)=1
      ss(0)=2
      tt(0)=1
      rri(0)=1
      rrj(0)=1

      if (nelec.lt.2) goto 65

c e-e and e-e-n terms
      ij=0
      do 60 i=2,nelec
      im1=i-1
      do 60 j=1,im1
      ij=ij+1

      fso(i,j)=0
      d2ijo(i,j)=0
      do 10 k=1,3
        fijo(k,i,j)=0
   10   fijo(k,j,i)=0

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

      call scale_dist2(rij,uu(1),dd1,dd2,1)
c     write(ounit,'(''rij,u in ee'',2f9.5)') rij,uu(1)

c Check rij after scaling because uu(1) used in e-e-n terms too
      if(rij.gt.cutjas) goto 22

      top=sspinn*b(1,isb,iwf)*uu(1)
      topu=sspinn*b(1,isb,iwf)
      topuu=0

      if(ijas.eq.4.or.ijas.eq.5) then
        bot=1+b(2,isb,iwf)*uu(1)
        botu=b(2,isb,iwf)
       elseif(ijas.eq.6) then
        bot=1+b(2,isb,iwf)*(1-uu(1))
        botu=-b(2,isb,iwf)
      endif
      botuu=0
      bot2=bot*bot

      fee=top/bot-asymp_jasb(ipar+1)
      feeu=topu/bot-botu*top/bot2
      feeuu=topuu-(botuu*top+2*botu*topu)/bot+2*botu**2*top/bot2
      feeuu=feeuu/bot

      do 20 iord=2,nordb
        uu(iord)=uu(1)*uu(iord-1)
        if(ijas.eq.4) then
          fee=fee+b(iord+1,isb,iwf)*uu(iord)
          feeu=feeu+b(iord+1,isb,iwf)*iord*uu(iord-1)
          feeuu=feeuu+b(iord+1,isb,iwf)*iord*(iord-1)*uu(iord-2)
         elseif(ijas.eq.5.or.ijas.eq.6) then
          fee=fee+sspinn*b(iord+1,isb,iwf)*uu(iord)
          feeu=feeu+sspinn*b(iord+1,isb,iwf)*iord*uu(iord-1)
          feeuu=feeuu+sspinn*b(iord+1,isb,iwf)*iord*(iord-1)*uu(iord-2)
        endif
   20 continue

      feeuu=feeuu*dd1*dd1+feeu*dd2
      feeu=feeu*dd1/rij

      fso(i,j)=fso(i,j)+fee
      do 21 k=1,3
        fijo(k,i,j)= fijo(k,i,j) + feeu*rvec_ee(k,ij)
   21   fijo(k,j,i)= fijo(k,j,i) - feeu*rvec_ee(k,ij)
      d2ijo(i,j)=d2ijo(i,j)+2*(feeuu+2*feeu)

c There are no C terms to order 1.
   22 if(nordc.le.1) goto 55

      if(isc.ge.12) call scale_dist2(rij,uu(1),dd1,dd2,3)
      if(ijas.eq.4.or.ijas.eq.5) then
        call switch_scale2(uu(1),dd1,dd2)
        do 25 iord=2,nordc
   25     uu(iord)=uu(1)*uu(iord-1)
      endif
c     write(ounit,'(''rij,u in een'',2f12.9)') rij,uu(1)

      do 50 ic=1,ncent
        it=iwctype(ic)

        ri=r_en(i,ic)
        rj=r_en(j,ic)

        if(ri.gt.cutjas .or. rj.gt.cutjas) goto 50
        do 27 k=1,3
   27     if(abs(rshift(k,i,ic)-rshift(k,j,ic)).gt.eps) goto 50

        call scale_dist2(ri,rri(1),dd7,dd9,2)
        call scale_dist2(rj,rrj(1),dd8,dd10,2)

        if(ijas.eq.4.or.ijas.eq.5) then
          call switch_scale2(rri(1),dd7,dd9)
          call switch_scale2(rrj(1),dd8,dd10)
        endif
c     write(ounit,'(''ri,rri in een'',2f12.9)') ri,rri(1)

        s=ri+rj
        t=ri-rj
c       u2mt2=rij*rij-t*t
        u2pst=rij*rij+s*t
        u2mst=rij*rij-s*t
c       s2mu2=s*s-rij*rij
c       s2mt2=s*s-t*t

        do 30 iord=1,nordc
          rri(iord)=rri(1)*rri(iord-1)
          rrj(iord)=rrj(1)*rrj(iord-1)
          ss(iord)=rri(iord)+rrj(iord)
   30     tt(iord)=rri(iord)*rrj(iord)

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
        do 40 n=2,nordc
          do 40 k=n-1,0,-1
            if(k.eq.0) then
              l_hi=n-k-2
             else
              l_hi=n-k
            endif
            do 40 l=l_hi,0,-1
              m=(n-k-l)/2
              if(2*m.eq.n-k-l) then
                ll=ll+1
                fc=fc+c(ll,it,iwf)*uu(k)*ss(l)*tt(m)
                fu=fu+c(ll,it,iwf)*k*uu(k-1)*ss(l)*tt(m)
                fuu=fuu+c(ll,it,iwf)*k*(k-1)*uu(k-2)*ss(l)*tt(m)
                fi=fi+c(ll,it,iwf)*uu(k)
     &          *((l+m)*rri(l+m-1)*rrj(m)+m*rri(m-1)*rrj(l+m))
                fii=fii+c(ll,it,iwf)*uu(k)
     &          *((l+m)*(l+m-1)*rri(l+m-2)*rrj(m)
     &          +m*(m-1)*rri(m-2)*rrj(l+m))
                fj=fj+c(ll,it,iwf)*uu(k)
     &          *((l+m)*rrj(l+m-1)*rri(m)+m*rrj(m-1)*rri(l+m))
                fjj=fjj+c(ll,it,iwf)*uu(k)
     &          *((l+m)*(l+m-1)*rrj(l+m-2)*rri(m)
     &          +m*(m-1)*rrj(m-2)*rri(l+m))
                fui=fui+c(ll,it,iwf)*k*uu(k-1)
     &          *((l+m)*rri(l+m-1)*rrj(m)+m*rri(m-1)*rrj(l+m))
                fuj=fuj+c(ll,it,iwf)*k*uu(k-1)
     &          *((l+m)*rrj(l+m-1)*rri(m)+m*rrj(m-1)*rri(l+m))
              endif
c     write(ounit,'(''rij,ri,rj'',9f10.5)') rij,ri,rj,uu(1),rri(1),rrj(1)
   40   continue

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
c       write(ounit,'(''i,j,fijo2='',2i5,9d12.4)') i,j,(fijo(k,i,j),k=1,3)

        d2ijo(i,j)=d2ijo(i,j) + 2*(fuu + 2*fu) + fui*u2pst/(ri*rij)
     &  + fuj*u2mst/(rj*rij) + fii + 2*fi + fjj + 2*fj

  50  continue

  55  fsum=fsum+fso(i,j)
      v(1,i)=v(1,i)+fijo(1,i,j)
      v(2,i)=v(2,i)+fijo(2,i,j)
      v(3,i)=v(3,i)+fijo(3,i,j)
      v(1,j)=v(1,j)+fijo(1,j,i)
      v(2,j)=v(2,j)+fijo(2,j,i)
      v(3,j)=v(3,j)+fijo(3,j,i)
      div_vj(i)=div_vj(i)+d2ijo(i,j)/2
      div_vj(j)=div_vj(j)+d2ijo(i,j)/2
  60  d2=d2+d2ijo(i,j)

c e-n terms
  65  do 90 i=1,nelec

        fso(i,i)=0
        fijo(1,i,i)=0
        fijo(2,i,i)=0
        fijo(3,i,i)=0
        d2ijo(i,i)=0

        do 80 ic=1,ncent
          it=iwctype(ic)

          ri=r_en(i,ic)
          if(ri.gt.cutjas) goto 80

          call scale_dist2(ri,rri(1),dd7,dd9,1)
c     write(ounit,'(''ri,rri in en'',2f9.5)') ri,rri(1)

          top=a4(1,it,iwf)*rri(1)
          topi=a4(1,it,iwf)
          topii=0

          if(ijas.eq.4.or.ijas.eq.5) then
            bot=1+a4(2,it,iwf)*rri(1)
            boti=a4(2,it,iwf)
           elseif(ijas.eq.6) then
            bot=1+a4(2,it,iwf)*(1-rri(1))
            boti=-a4(2,it,iwf)
          endif
          botii=0
          bot2=bot*bot

          fen=top/bot-asymp_jasa(it)
          feni=topi/bot-boti*top/bot2
          fenii=topii-(botii*top+2*boti*topi)/bot+2*boti**2*top/bot2
          fenii=fenii/bot

        do 70 iord=2,norda
          rri(iord)=rri(1)**iord
          fen=fen+a4(iord+1,it,iwf)*rri(iord)
          feni=feni+a4(iord+1,it,iwf)*iord*rri(iord-1)
   70     fenii=fenii+a4(iord+1,it,iwf)*iord*(iord-1)*rri(iord-2)

          feni_save=feni
          fenii_save=fenii
          fenii=fenii*dd7*dd7+feni*dd9
          feni=feni*dd7/ri

          fso(i,i)=fso(i,i)+fen

          fijo(1,i,i)=fijo(1,i,i) + feni*rvec_en(1,i,ic)
          fijo(2,i,i)=fijo(2,i,i) + feni*rvec_en(2,i,ic)
          fijo(3,i,i)=fijo(3,i,i) + feni*rvec_en(3,i,ic)
c         write(ounit,'(''fijo='',9d12.4)') (fijo(k,i,i),k=1,3),feni,rvec_en(1,i,ic)

          d2ijo(i,i) = d2ijo(i,i) + fenii + 2*feni

          if(iforce_analy.eq.1) call da_jastrow4(iwf,i,ic,it,rvec_en(1,i,ic),ri,rri,feni_save,fenii_save,dd7,dd9)
   80   continue

        fsum=fsum+fso(i,i)
        v(1,i)=v(1,i)+fijo(1,i,i)
        v(2,i)=v(2,i)+fijo(2,i,i)
        v(3,i)=v(3,i)+fijo(3,i,i)
c       write(ounit,'(''v='',9d12.4)') (v(k,i),k=1,3)
        div_vj(i)=div_vj(i)+d2ijo(i,i)
   90   d2=d2+d2ijo(i,i)

      if(ijas.eq.6) then
        term=1/(c1_jas6*scalek(iwf))
        fsum=term*fsum
        d2=term*d2
        do 100 i=1,nelec
          div_vj(i)=term*div_vj(i)
          do 95 k=1,3
   95       v(k,i)=term*v(k,i)
          do 100 j=1,nelec
            d2ijo(i,j)=term*d2ijo(i,j)
            do 100 k=1,3
  100         fijo(k,i,j)=term*fijo(k,i,j)
      endif

      fsumo=fsum
      d2o=d2
      do 110 i=1,nelec
        fjo(1,i)=v(1,i)
        fjo(2,i)=v(2,i)
  110   fjo(3,i)=v(3,i)

      value=fsum

      return
      end

c-----------------------------------------------------------------------
      function nterms4(nord)
c Written by Cyrus Umrigar

      implicit none

      integer :: i, k, l, l_hi, m
      integer :: n, nord, nterms4

      i=0
      do 20 n=2,nord
        do 20 k=n-1,0,-1
          if(k.eq.0) then
            l_hi=n-k-2
           else
            l_hi=n-k
          endif
          do 20 l=l_hi,0,-1
            m=(n-k-l)/2
            if(2*m.eq.n-k-l) then
              i=i+1
            endif
   20 continue
      nterms4=i
c     write(ounit,'(''nterms4='',i5)') nterms4
      return
      end
c-----------------------------------------------------------------------
      subroutine da_jastrow4(iwf,i,ic,it,rvec_en,r,rr,feni,fenii,dd1,dd2)

      use vmc_mod, only: nordj
      use da_jastrow4val, only: da_d2j, da_j, da_vj
      use jaspar4, only: a4, norda
      use scale_more, only: dd3
      use precision_kinds, only: dp

      implicit none

      integer :: i, ic, iord, it, iwf
      integer :: k, l
      real(dp) :: dd1, dd2, feni, fenii, feniii
      real(dp) :: r, ri, ri2
      real(dp), dimension(3) :: rvec_en
      real(dp), dimension(-2:nordj) :: rr

      feniii=0.d0
      do 10 iord=3,norda
   10   feniii=feniii+a4(iord+1,it,iwf)*iord*(iord-1)*(iord-2)*rr(iord-3)

      ri=1.d0/r
      ri2=ri*ri

      do 30 k=1,3
        da_j(k,i,ic)=-rvec_en(k)*ri*feni*dd1
        da_d2j(k,i,ic)=-rvec_en(k)*ri*(feniii*dd1*dd1*dd1+fenii*dd1*(3*dd2+2*dd1*ri)+feni*(dd3+2*dd2*ri-2*dd1*ri2))
        do 20 l=1,3
   20     da_vj(k,l,i,ic)=-rvec_en(k)*rvec_en(l)*ri2*(fenii*dd1*dd1+feni*dd2-feni*dd1*ri)
   30   da_vj(k,k,i,ic)=da_vj(k,k,i,ic)-feni*dd1*ri

      return
      end
