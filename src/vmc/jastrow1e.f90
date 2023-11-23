      module jastrow1e_mod
      contains
      subroutine jastrow1e(iel,x,fjn,d2n,fsumn,fsn,fijn,d2ijn, &
                              fjo,d2o,fsumo,fso,fijo,d2ijo,iflag)
! Written by Cyrus Umrigar and Claudia Filippi
! Jastrow 4,5 must be used with one of isc=2,4,6,7,12,14,16,17
! Jastrow 6   must be used with one of isc=6,7
      use system, only: iwctype, ncent, nelec, nup
      use jastrow, only: sspinn, b, c, a4, norda, nordb, nordc, asymp_jasa, asymp_jasb, ijas, isc, nordj
      use jaspar6, only: cutjas
      use jastrow, only: cutjas_ee, cutjas_en, cutjas_eei, cutjas_eni
      use multiple_geo, only: iwf
      use bparm, only: nocuspb, nspin2b
      use distance_mod, only: r_en, rvec_en, r_ee, rvec_ee
      use precision_kinds, only: dp
      use scale_dist_mod, only: scale_dist1,scale_dist2,switch_scale1
      use scale_dist_mod, only: switch_scale2
      use system,  only: iwctype,ncent,nelec,nup


      implicit none

      integer :: i, ic, iel, iflag, ij
      integer :: iord, ipar, isb, it
      integer :: j, jj, k, l
      integer :: l_hi, ll, m, n
      real(dp) :: a1_cusp, b1_cusp, d2, dd1
      real(dp) :: fc, fee, feeu, feeu_save, feeuu
      real(dp) :: fen, feni, feni_save, fenii, fi, fi_save
      real(dp) :: fii, fj, fj_save, fjj, fu
      real(dp) :: fui, fuj, fuu
      real(dp) :: s, t, term, termi, termii, termj, termjj
      real(dp) :: termu, termuu, term1, term2
      real(dp) :: u2mst, u2pst, value, xi, xij, xj
      real(dp), dimension(3, *) :: x
      real(dp), dimension(-2:nordj) :: rij
      real(dp), dimension(-2:nordj) :: ss
      real(dp), dimension(-2:nordj) :: tt
      real(dp), dimension(-2:nordj) :: ri
      real(dp), dimension(-2:nordj) :: rj
      real(dp), parameter :: half = .5d0
      real(dp), parameter :: eps = 1.d-12
      real(dp) :: fsumo, d2o
      real(dp), dimension(3, *) :: fjo
      real(dp), dimension(nelec, *) :: fso
      real(dp), dimension(3, nelec, *) :: fijo
      real(dp), dimension(nelec, *) :: d2ijo
      real(dp) :: fsumn, d2n
      real(dp), dimension(3, *) :: fjn
      real(dp), dimension(nelec, *) :: fsn
      real(dp), dimension(3, nelec, *) :: fijn
      real(dp), dimension(nelec, *) :: d2ijn

      do i=-2,-1
        rij(i)=0
        ss(i)=0
        tt(i)=0
        ri(i)=0
        rj(i)=0
      enddo
      rij(0)=1
      ss(0)=2
      tt(0)=1
      ri(0)=1
      rj(0)=1
      if (nelec.lt.2) goto 65

      do i=1,nelec
        do k=1,3
        fjn(k,i)=fjo(k,i)
        enddo
      enddo
      fsumn=fsumo
      d2n=d2o

      do jj=1,nelec

      if(jj.eq.iel) goto 60
      if(jj.lt.iel) then
        i=iel
        j=jj
       else
        i=jj
        j=iel
      endif
      ij=((i-1)*(i-2))/2+j

      fijn(1,i,j)=0
      fijn(2,i,j)=0
      fijn(3,i,j)=0
      fijn(1,j,i)=0
      fijn(2,j,i)=0
      fijn(3,j,i)=0
      fsn(i,j)=0
      d2ijn(i,j)=0

      sspinn=1
      ipar=0
      isb=1
      if(i.le.nup .or. j.gt.nup) then
         sspinn=half
        if(nspin2b.eq.2) then
          isb=2
         elseif(nocuspb.gt.0) then
          sspinn=1
        endif
        ipar=1
      endif

      rij(1)=r_ee(ij)
      if(rij(1).gt.cutjas_ee(isb,iwf)) goto 22

      b1_cusp=sspinn*0.5+3.d0*b(1,isb,iwf)*cutjas_eei(isb,iwf)
      fee=b1_cusp*rij(1)+b(1,isb,iwf)
      feeu=b1_cusp

      do iord=2,nordb
        rij(iord)=rij(1)*rij(iord-1)
        fee=fee+b(iord,isb,iwf)*rij(iord)
        feeu=feeu+b(iord,isb,iwf)*iord*rij(iord-1)
      enddo

      if(iflag.gt.0) then
         feeuu=0
         do iord=2,nordb
            feeuu=feeuu+b(iord,isb,iwf)*iord*(iord-1)*rij(iord-2)
         enddo
      endif

      xij=rij(1)*cutjas_eei(isb,iwf)
      term=(1.d0-xij)**3
      termu=-3*(1.d0-xij)**2*cutjas_eei(isb,iwf)
      termuu=6*(1.d0-xij)*cutjas_eei(isb,iwf)*cutjas_eei(isb,iwf)

      feeu_save=feeu*term+fee*termu
      feeuu=feeuu*term+2*feeu*termu+fee*termuu

      feeu=feeu_save/rij(1)

      fsn(i,j)=fsn(i,j) + fee*term


      fijn(1,i,j)=fijn(1,i,j) + feeu*rvec_ee(1,ij)
      fijn(2,i,j)=fijn(2,i,j) + feeu*rvec_ee(2,ij)
      fijn(3,i,j)=fijn(3,i,j) + feeu*rvec_ee(3,ij)
      fijn(1,j,i)=fijn(1,j,i) - feeu*rvec_ee(1,ij)
      fijn(2,j,i)=fijn(2,j,i) - feeu*rvec_ee(2,ij)
      fijn(3,j,i)=fijn(3,j,i) - feeu*rvec_ee(3,ij)

      if(iflag.gt.0) d2ijn(i,j)=d2ijn(i,j) + 2*(feeuu+2*feeu)

! There are no C terms to order 1.
   22 if(nordc.le.1) goto 55


      do iord=2,nordc
         rij(iord)=rij(1)*rij(iord-1)
      enddo


      do ic=1,ncent
        it=iwctype(ic)

        ri(1)=r_en(i,ic)
        rj(1)=r_en(j,ic)

        if(ri(1).gt.cutjas_en(it,iwf) .or. rj(1).gt.cutjas_en(it,iwf)) goto 50

        do iord=1,nordc
          ri(iord)=ri(1)*ri(iord-1)
          rj(iord)=rj(1)*rj(iord-1)
          ss(iord)=ri(iord)+rj(iord)
          tt(iord)=ri(iord)*rj(iord)
        enddo

        xi=ri(1)*cutjas_eni(it,iwf)
        xj=rj(1)*cutjas_eni(it,iwf)
        term1=(1.d0-xi)**3
        term2=(1.d0-xj)**3
        term=term1*term2
        termi=-3*(1.d0-xi)**2*cutjas_eni(it,iwf)*term2
        termii=6*(1.d0-xi)*cutjas_eni(it,iwf)*cutjas_eni(it,iwf)*term2
        termj=-3*(1.d0-xj)**2*cutjas_eni(it,iwf)*term1
        termjj=6*(1.d0-xj)*cutjas_eni(it,iwf)*cutjas_eni(it,iwf)*term1


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
                fc=fc+c(ll,it,iwf)*rij(k)*ss(l)*tt(m)
                fu=fu+c(ll,it,iwf)*k*rij(k-1)*ss(l)*tt(m)
                fi=fi+c(ll,it,iwf)*rij(k) &
                  *((l+m)*ri(l+m-1)*rj(m)+m*ri(m-1)*rj(l+m))
                fj=fj+c(ll,it,iwf)*rij(k) &
                  *((l+m)*rj(l+m-1)*ri(m)+m*rj(m-1)*ri(l+m))

                if(iflag.gt.0) then
                  fuu=fuu+c(ll,it,iwf)*k*(k-1)*rij(k-2)*ss(l)*tt(m)
                  fii=fii+c(ll,it,iwf)*rij(k) &
                    *((l+m)*(l+m-1)*ri(l+m-2)*rj(m) &
                    +m*(m-1)*ri(m-2)*rj(l+m))
                  fjj=fjj+c(ll,it,iwf)*rij(k) &
                    *((l+m)*(l+m-1)*rj(l+m-2)*ri(m) &
                    +m*(m-1)*rj(m-2)*ri(l+m))
                  fui=fui+c(ll,it,iwf)*k*rij(k-1) &
                    *((l+m)*ri(l+m-1)*rj(m)+m*ri(m-1)*rj(l+m))
                  fuj=fuj+c(ll,it,iwf)*k*rij(k-1) &
                    *((l+m)*rj(l+m-1)*ri(m)+m*rj(m-1)*ri(l+m))
                endif

              endif
            enddo
          enddo
        enddo

        if(iflag.gt.0) then
          s=ri(1)+rj(1)
          t=ri(1)-rj(1)
          u2pst=rij(1)*rij(1)+s*t
          u2mst=rij(1)*rij(1)-s*t

        endif


        fu=fu*term/rij(1)

        fi_save=fi*term+fc*termi
        fj_save=fj*term+fc*termj
        fii=fii*term+2*fi*termi+fc*termii
        fjj=fjj*term+2*fj*termj+fc*termjj

        fi=fi_save/ri(1)
        fj=fj_save/rj(1)


        fsn(i,j)=fsn(i,j) + fc

        fijn(1,i,j)=fijn(1,i,j) + fi*rvec_en(1,i,ic)+fu*rvec_ee(1,ij)
        fijn(2,i,j)=fijn(2,i,j) + fi*rvec_en(2,i,ic)+fu*rvec_ee(2,ij)
        fijn(3,i,j)=fijn(3,i,j) + fi*rvec_en(3,i,ic)+fu*rvec_ee(3,ij)
        fijn(1,j,i)=fijn(1,j,i) + fj*rvec_en(1,j,ic)-fu*rvec_ee(1,ij)
        fijn(2,j,i)=fijn(2,j,i) + fj*rvec_en(2,j,ic)-fu*rvec_ee(2,ij)
        fijn(3,j,i)=fijn(3,j,i) + fj*rvec_en(3,j,ic)-fu*rvec_ee(3,ij)

        if(iflag.gt.0) &
          d2ijn(i,j)=d2ijn(i,j) + 2*(fuu + 2*fu) + fui*u2pst/(ri(1)*rij(1)) &
           + fuj*u2mst/(rj(1)*rij(1)) + fii + 2*fi + fjj + 2*fj
   50 continue
      enddo

   55 fsumn=fsumn+fsn(i,j)-fso(i,j)
      fjn(1,i)=fjn(1,i)+fijn(1,i,j)-fijo(1,i,j)
      fjn(2,i)=fjn(2,i)+fijn(2,i,j)-fijo(2,i,j)
      fjn(3,i)=fjn(3,i)+fijn(3,i,j)-fijo(3,i,j)
      fjn(1,j)=fjn(1,j)+fijn(1,j,i)-fijo(1,j,i)
      fjn(2,j)=fjn(2,j)+fijn(2,j,i)-fijo(2,j,i)
      fjn(3,j)=fjn(3,j)+fijn(3,j,i)-fijo(3,j,i)
      d2n=d2n+d2ijn(i,j)-d2ijo(i,j)
   60 continue
      enddo

! e-n terms

   65 fijn(1,iel,iel)=0
      fijn(2,iel,iel)=0
      fijn(3,iel,iel)=0
      fsn(iel,iel)=0
      d2ijn(iel,iel)=0

      do ic=1,ncent
        it=iwctype(ic)

        ri(1)=r_en(iel,ic)
        if(ri(1).gt.cutjas_en(it,iwf)) goto 80

        a1_cusp=3.d0*a4(1,it,iwf)*cutjas_eni(it,iwf)
        fen=a1_cusp*ri(1)+a4(1,it,iwf)
        feni=a1_cusp


        do iord=2,norda
          ri(iord)=ri(1)*ri(iord-1)
          fen=fen+a4(iord,it,iwf)*ri(iord)
          feni=feni+a4(iord,it,iwf)*iord*ri(iord-1)
        enddo

        if(iflag.gt.0) then
           fenii=0
           do iord=2,norda
              fenii=fenii+a4(iord+1,it,iwf)*iord*(iord-1)*ri(iord-2)
           enddo
        endif

        xi=ri(1)*cutjas_eni(it,iwf)
        term=(1.d0-xi)**3
        termi=-3*(1.d0-xi)**2*cutjas_eni(it,iwf)
        termii=6*(1.d0-xi)*cutjas_eni(it,iwf)*cutjas_eni(it,iwf)

        feni_save=feni*term+fen*termi
        fenii=fenii*term+2*feni*termi+fen*termii

        feni=feni_save/ri(1)

        fsn(iel,iel)=fsn(iel,iel)+fen*term

        fijn(1,iel,iel)=fijn(1,iel,iel) + feni*rvec_en(1,iel,ic)
        fijn(2,iel,iel)=fijn(2,iel,iel) + feni*rvec_en(2,iel,ic)
        fijn(3,iel,iel)=fijn(3,iel,iel) + feni*rvec_en(3,iel,ic)

        if(iflag.gt.0) d2ijn(iel,iel) = d2ijn(iel,iel) + fenii + 2*feni
   80 continue
      enddo

      fsumn=fsumn+fsn(iel,iel)-fso(iel,iel)
      fjn(1,iel)=fjn(1,iel)+fijn(1,iel,iel)-fijo(1,iel,iel)
      fjn(2,iel)=fjn(2,iel)+fijn(2,iel,iel)-fijo(2,iel,iel)
      fjn(3,iel)=fjn(3,iel)+fijn(3,iel,iel)-fijo(3,iel,iel)
      d2n=d2n+d2ijn(iel,iel)-d2ijo(iel,iel)


      return
      end
      end module
