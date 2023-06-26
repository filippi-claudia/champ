      module deriv_jastrow1_mod
      contains
      subroutine deriv_jastrow1(x,v,d2,value)

      use bparm, only: nocuspb,nspin2b
      use contrl_file, only: ounit
      use cuspmat4, only: d,iwc4
      use derivjas, only: d2g,g,go,gvalue
      use distance_mod, only: r_ee,r_en,rshift,rvec_ee,rvec_en
      use jastrow, only: cutjas_ee, cutjas_eei, cutjas_en, cutjas_eni
      use jastrow, only: norda,nordb,nordc
      use jastrow, only: a4,b,c,ijas,isc,nordj
      use jastrow, only: sspinn
      use jaspointer, only: npoint,npointa
      use jastrow1_mod, only: da_jastrow1
      use jastrow4_mod, only: da_jastrow4
      use jastrow_update, only: d2ijo,d2o,fijo,fjo,fso,fsumo
      use m_force_analytic, only: iforce_analy
      use multiple_geo, only: iwf
      use optwf_control, only: ioptjas
      use optwf_nparmj, only: nparma,nparmb
      use optwf_parms, only: nparmj
      use optwf_wjas, only: iwjasa,iwjasb,iwjasc
      use precision_kinds, only: dp
      use system,  only: iwctype,ncent,nctype,nelec,nup
      use vardep,  only: cdep,iwdepend,nvdepend
      implicit none

      integer :: i, ic, id, ideriv, ij
      integer :: im1, iord, ipar, iparm
      integer :: iparm0, iparma, iparmb, isb, it
      integer :: j, jj, jparm, k
      integer :: l, l_hi, ll, m, n
      real(dp) :: a1_cusp, da1_cusp, bot, bot0, bot2, boti, botii
      real(dp) :: botu, botuu, b1_cusp, db1_cusp, cd, d2
      real(dp) :: fc, fee, feeu, feeu_save, feeuu, fen, feni
      real(dp) :: feni_save, fenii, fenii_save, fi, fi_save
      real(dp) :: fii, fj, fj_save, fjj, fsum
      real(dp) :: fu, fui, fuj, fuu
      real(dp) :: gee, geeu, geeu_save, geeuu, gen
      real(dp) :: geni, geni_save, genii, gi, gi_save, gii
      real(dp) :: gj, gj_save, gjj, gp, gu
      real(dp) :: gui, guj, guu, pc
      real(dp) :: pii, pj, pjj, ppi
      real(dp) :: pu, pui, puj, puu
      real(dp) :: s, t, term1, term2
      real(dp) :: term, termi, termii, termj, termjj, termu, termuu
      real(dp) :: u2mst, u2pst
      real(dp) :: value, xi, xij, xj
      real(dp), dimension(3, *) :: x
      real(dp), dimension(3, *) :: v
      real(dp), dimension(-2:nordj) :: ri
      real(dp), dimension(-2:nordj) :: rj
      real(dp), dimension(-2:nordj) :: rij
      real(dp), dimension(-2:nordj) :: ss
      real(dp), dimension(-2:nordj) :: tt
      real(dp), parameter :: zero = 0.d0
      real(dp), parameter :: one = 1.d0
      real(dp), parameter :: two = 2.d0
      real(dp), parameter :: half = .5d0
      real(dp), parameter :: eps = 1.d-12

      iparma=nparma(1)
      do it=2,nctype
       iparma=iparma+nparma(it)
      enddo

      fsum=0
      do i=-2,-1
        rij(i)=0
        ri(i)=0
        rj(i)=0
        ss(i)=0
        tt(i)=0
      enddo
      rij(0)=1
      ri(0)=1
      rj(0)=1
      ss(0)=2
      tt(0)=1

      do iparm=1,nparmj
        gvalue(iparm)=0
        d2g(iparm)=0
        do i=1,nelec
          g(1,i,iparm)=0
          g(2,i,iparm)=0
          g(3,i,iparm)=0
        enddo
      enddo

      if (nelec.lt.2) goto 65

c e-e and e-e-n terms
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
      feeuu=0

      do iord=2,nordb
        rij(iord)=rij(1)*rij(iord-1)

        fee=fee+b(iord,isb,iwf)*rij(iord)
        feeu=feeu+b(iord,isb,iwf)*iord*rij(iord-1)
        feeuu=feeuu+b(iord,isb,iwf)*iord*(iord-1)*rij(iord-2)
      enddo

      xij=rij(1)*cutjas_eei(isb,iwf)
      term=(1.d0-xij)**3
      termu=-3*(1.d0-xij)**2*cutjas_eei(isb,iwf)
      termuu=6*(1.d0-xij)*cutjas_eei(isb,iwf)*cutjas_eei(isb,iwf)

      feeu_save=feeu*term+fee*termu
      feeuu=feeuu*term+2*feeu*termu+fee*termuu

      feeu=feeu_save/rij(1)

      fso(i,j)=fso(i,j)+fee*term
      do k=1,3
        fijo(k,i,j)= fijo(k,i,j) + feeu*rvec_ee(k,ij)
        fijo(k,j,i)= fijo(k,j,i) - feeu*rvec_ee(k,ij)
      enddo
      d2ijo(i,j)=d2ijo(i,j)+2*(feeuu+2*feeu)

c derivatives of wave function wrt b(i)
      iparm0=iparma
      if(isb.eq.2) iparm0=iparm0+nparmb(1)
      do jparm=1,nparmb(isb)
        iparm=iparm0+jparm

        if(iwjasb(jparm,isb).eq.1) then

          db1_cusp=3.d0*cutjas_eei(isb,iwf)
          gee=db1_cusp*rij(1)+1.d0
          geeu=db1_cusp
          geeuu=0.d0

        else

          iord=iwjasb(jparm,isb)
          gee=rij(iord)
          geeu=iord*rij(iord-1)
          geeuu=iord*(iord-1)*rij(iord-2)

        endif

        geeu_save=geeu*term+gee*termu
        geeuu=geeuu*term+2*geeu*termu+gee*termuu

        geeu=geeu_save/rij(1)

        go(i,j,iparm)=go(i,j,iparm)+gee*term
        gvalue(iparm)=gvalue(iparm)+gee*term

        g(1,i,iparm)=g(1,i,iparm) + geeu*rvec_ee(1,ij)
        g(2,i,iparm)=g(2,i,iparm) + geeu*rvec_ee(2,ij)
        g(3,i,iparm)=g(3,i,iparm) + geeu*rvec_ee(3,ij)
        g(1,j,iparm)=g(1,j,iparm) - geeu*rvec_ee(1,ij)
        g(2,j,iparm)=g(2,j,iparm) - geeu*rvec_ee(2,ij)
        g(3,j,iparm)=g(3,j,iparm) - geeu*rvec_ee(3,ij)

        d2g(iparm)=d2g(iparm) + two*(geeuu+two*geeu)
      enddo

c MISSING -> derivative wrt iparmb parm = cutjas_ee

c There are no C terms to order 1.
   22 if(nordc.le.1) goto 55

      do ic=1,ncent
        it=iwctype(ic)

        ri(1)=r_en(i,ic)
        rj(1)=r_en(j,ic)

        if(ri(1).gt.cutjas_en(it,iwf) .or. rj(1).gt.cutjas_en(it,iwf)) goto 50

        s=ri(1)+rj(1)
        t=ri(1)-rj(1)
        u2pst=rij(1)*rij(1)+s*t
        u2mst=rij(1)*rij(1)-s*t

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
        termi=-3*(1.d0-xi)**2*cutjas_eei(it,iwf)*term2
        termii=6*(1.d0-xi)*cutjas_eei(it,iwf)*cutjas_eei(it,iwf)*term2
        termj=-3*(1.d0-xj)**2*cutjas_eei(it,iwf)*term1
        termjj=6*(1.d0-xj)*cutjas_eei(it,iwf)*cutjas_eei(it,iwf)*term1

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
                pc=rij(k)*ss(l)*tt(m)
                pu=k*rij(k-1)*ss(l)*tt(m)
                puu=k*(k-1)*rij(k-2)*ss(l)*tt(m)
                ppi=rij(k)
     &          *((l+m)*ri(l+m-1)*rj(m)+m*ri(m-1)*rj(l+m))
                pii=rij(k)
     &          *((l+m)*(l+m-1)*ri(l+m-2)*rj(m)
     &          +m*(m-1)*ri(m-2)*rj(l+m))
                pj=rij(k)
     &          *((l+m)*rj(l+m-1)*ri(m)+m*rj(m-1)*ri(l+m))
                pjj=rij(k)
     &          *((l+m)*(l+m-1)*rj(l+m-2)*ri(m)
     &          +m*(m-1)*rj(m-2)*ri(l+m))
                pui=k*rij(k-1)
     &          *((l+m)*ri(l+m-1)*rj(m)+m*ri(m-1)*rj(l+m))
                puj=k*rij(k-1)
     &          *((l+m)*rj(l+m-1)*ri(m)+m*rj(m-1)*ri(l+m))

                fc=fc+c(ll,it,iwf)*pc
                fu=fu+c(ll,it,iwf)*pu
                fuu=fuu+c(ll,it,iwf)*puu
                fi=fi+c(ll,it,iwf)*ppi
                fii=fii+c(ll,it,iwf)*pii
                fj=fj+c(ll,it,iwf)*pj
                fjj=fjj+c(ll,it,iwf)*pjj
                fui=fui+c(ll,it,iwf)*pui
                fuj=fuj+c(ll,it,iwf)*puj

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

                  gui=gu*termi
                  guj=gu*termj
                  gu=gu*term/rij(1)
                  guu=guu*term

                  gi_save=gi*term+gp*termi
                  gj_save=gj*term+gp*termj
                  gii=gii*term+2*gi*termi+gp*termii
                  gjj=gjj*term+2*gj*termj+gp*termjj

                  gi=gi_save/ri(1)
                  gj=gj_save/rj(1)

                  if(ideriv.eq.1) then

                  do id=1,nvdepend(jj,it)

                  iparm=npoint(it)+iwdepend(jj,id,it)
                  cd=cdep(jj,id,it)

                  go(i,j,iparm)=go(i,j,iparm)+cd*gp*term
                  gvalue(iparm)=gvalue(iparm)+cd*gp*term

                  g(1,i,iparm)=g(1,i,iparm)+cd*(gi*rvec_en(1,i,ic)+gu*rvec_ee(1,ij))
                  g(2,i,iparm)=g(2,i,iparm)+cd*(gi*rvec_en(2,i,ic)+gu*rvec_ee(2,ij))
                  g(3,i,iparm)=g(3,i,iparm)+cd*(gi*rvec_en(3,i,ic)+gu*rvec_ee(3,ij))
                  g(1,j,iparm)=g(1,j,iparm)+cd*(gj*rvec_en(1,j,ic)-gu*rvec_ee(1,ij))
                  g(2,j,iparm)=g(2,j,iparm)+cd*(gj*rvec_en(2,j,ic)-gu*rvec_ee(2,ij))
                  g(3,j,iparm)=g(3,j,iparm)+cd*(gj*rvec_en(3,j,ic)-gu*rvec_ee(3,ij))

                  d2g(iparm)=d2g(iparm) + cd*(2*(guu + 2*gu)
     &            + gui*u2pst/(ri(1)*rij(1)) + guj*u2mst/(rj(1)*rij(1))
     &            + gii + 2*gi + gjj + 2*gj)
                  enddo

                  elseif(ideriv.eq.2) then

                  iparm=npoint(it)+jparm

                  go(i,j,iparm)=go(i,j,iparm)+gp*term
                  gvalue(iparm)=gvalue(iparm)+gp*term

                  g(1,i,iparm)=g(1,i,iparm)+gi*rvec_en(1,i,ic)+gu*rvec_ee(1,ij)
                  g(2,i,iparm)=g(2,i,iparm)+gi*rvec_en(2,i,ic)+gu*rvec_ee(2,ij)
                  g(3,i,iparm)=g(3,i,iparm)+gi*rvec_en(3,i,ic)+gu*rvec_ee(3,ij)
                  g(1,j,iparm)=g(1,j,iparm)+gj*rvec_en(1,j,ic)-gu*rvec_ee(1,ij)
                  g(2,j,iparm)=g(2,j,iparm)+gj*rvec_en(2,j,ic)-gu*rvec_ee(2,ij)
                  g(3,j,iparm)=g(3,j,iparm)+gj*rvec_en(3,j,ic)-gu*rvec_ee(3,ij)

                  d2g(iparm)=d2g(iparm) + 2*(guu + 2*gu)
     &            + gui*u2pst/(ri(1)*rij(1)) + guj*u2mst/(rj(1)*rij(1))
     &            + gii + 2*gi + gjj + 2*gj

                  jparm=jparm+1
                  endif

                endif
              endif
            enddo
          enddo
        enddo

        fui=fu*termi
        fuj=fu*termj
        fu=fu*term/rij(1)
        fuu=fuu*term

        fi_save=fi*term+fc*termi
        fj_save=fj*term+fc*termj
        fii=fii*term+2*fi*termi+fc*termii
        fjj=fjj*term+2*fj*termj+fc*termjj

        fi=fi_save/ri(1)
        fj=fj_save/rj(1)

        fso(i,j)=fso(i,j) + fc*term

        fijo(1,i,j)=fijo(1,i,j) + fi*rvec_en(1,i,ic)+fu*rvec_ee(1,ij)
        fijo(2,i,j)=fijo(2,i,j) + fi*rvec_en(2,i,ic)+fu*rvec_ee(2,ij)
        fijo(3,i,j)=fijo(3,i,j) + fi*rvec_en(3,i,ic)+fu*rvec_ee(3,ij)
        fijo(1,j,i)=fijo(1,j,i) + fj*rvec_en(1,j,ic)-fu*rvec_ee(1,ij)
        fijo(2,j,i)=fijo(2,j,i) + fj*rvec_en(2,j,ic)-fu*rvec_ee(2,ij)
        fijo(3,j,i)=fijo(3,j,i) + fj*rvec_en(3,j,ic)-fu*rvec_ee(3,ij)

        d2ijo(i,j)=d2ijo(i,j) + 2*(fuu + 2*fu) + fui*u2pst/(ri(1)*rij(1))
     &  + fuj*u2mst/(rj(1)*rij(1)) + fii + 2*fi + fjj + 2*fj

  50  continue
      enddo

  55  fsum=fsum+fso(i,j)
      v(1,i)=v(1,i)+fijo(1,i,j)
      v(2,i)=v(2,i)+fijo(2,i,j)
      v(3,i)=v(3,i)+fijo(3,i,j)
      v(1,j)=v(1,j)+fijo(1,j,i)
      v(2,j)=v(2,j)+fijo(2,j,i)
      v(3,j)=v(3,j)+fijo(3,j,i)
      d2=d2+d2ijo(i,j)

      enddo
      enddo

c e-n terms
  65  do i=1,nelec

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

          ri(1)=r_en(i,ic)
          if(ri(1).gt.cutjas_en(it,iwf)) goto 80

          a1_cusp=3.d0*a4(1,it,iwf)*cutjas_eni(it,iwf)
          fen=a1_cusp*ri(1)+a4(1,it,iwf)
          feni=a1_cusp
          fenii=0

          do iord=2,norda
            ri(iord)=ri(1)*ri(iord-1)
            fen=fen+a4(iord,it,iwf)*ri(iord)
            feni=feni+a4(iord,it,iwf)*iord*ri(iord-1)
            fenii=fenii+a4(iord,it,iwf)*iord*(iord-1)*ri(iord-2)
          enddo

          xi=ri(1)*cutjas_eni(it,iwf)
          term=(1.d0-xi)**3
          termi=-3*(1.d0-xi)**2*cutjas_eni(it,iwf)
          termii=6*(1.d0-xi)*cutjas_eni(it,iwf)*cutjas_eni(it,iwf)

          feni_save=feni*term+fen*termi
          fenii=fenii*term+2*feni*termi+fen*termii

          feni=feni_save/ri(1)

          fso(i,i)=fso(i,i)+fen*term

          fijo(1,i,i)=fijo(1,i,i) + feni*rvec_en(1,i,ic)
          fijo(2,i,i)=fijo(2,i,i) + feni*rvec_en(2,i,ic)
          fijo(3,i,i)=fijo(3,i,i) + feni*rvec_en(3,i,ic)

          d2ijo(i,i) = d2ijo(i,i) + fenii + 2*feni

          do jparm=1,nparma(it)
            iparm=npointa(it)+jparm

            if(iwjasa(jparm,it).eq.1) then
              da1_cusp=3.d0*cutjas_eni(it,iwf)
              fen=a1_cusp*ri(1)+1.d0
              feni=a1_cusp
              fenii=0
             else

              iord=iwjasa(jparm,it)
              gen=ri(iord)
              geni=iord*ri(iord-1)
              genii=iord*(iord-1)*ri(iord-2)
            endif

            geni_save=geni*term+gen*termi
            genii=genii*term+2*geni*termi+gen*termii

            geni=geni_save/ri(1)

            go(i,i,iparm)=go(i,i,iparm)+gen*term
            gvalue(iparm)=gvalue(iparm)+gen*term

            g(1,i,iparm)=g(1,i,iparm)+geni*rvec_en(1,i,ic)
            g(2,i,iparm)=g(2,i,iparm)+geni*rvec_en(2,i,ic)
            g(3,i,iparm)=g(3,i,iparm)+geni*rvec_en(3,i,ic)

            d2g(iparm)=d2g(iparm)+genii+two*geni
          enddo

c MISSING a4(norda+1,it) as paramters

          if(iforce_analy.eq.1) call da_jastrow1(iwf,i,ic,it,rvec_en(1,i,ic),ri,feni_save,fenii)
   80   continue
        enddo


        fsum=fsum+fso(i,i)
        v(1,i)=v(1,i)+fijo(1,i,i)
        v(2,i)=v(2,i)+fijo(2,i,i)
        v(3,i)=v(3,i)+fijo(3,i,i)
        d2=d2+d2ijo(i,i)
      enddo

      fsumo=fsum
      d2o=d2
      do i=1,nelec
        fjo(1,i)=v(1,i)
        fjo(2,i)=v(2,i)
        fjo(3,i)=v(3,i)
      enddo

      value=fsum

      return
      end

      end module deriv_jastrow1_mod
