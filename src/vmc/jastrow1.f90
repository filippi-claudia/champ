      module jastrow1_mod
      contains
      subroutine jastrow_factor1(x,fjo,d2o,fsumo,fso,fijo,d2ijo)
! Written by Cyrus Umrigar, modified by C. Filippi

      use da_jastrow, only: da_d2j, da_vj
      use system, only: iwctype, ncent, nelec, nup
      use jastrow, only: sspinn, b, c, scalek, a4, norda, nordb, nordc, nordj
      use multiple_geo, only: iwf
      use bparm, only: nocuspb, nspin2b
      use jastrow, only: cutjas_ee, cutjas_en, cutjas_eei, cutjas_eni
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
      real(dp) :: a1_cusp, b1_cusp, d2, dd2
      real(dp) :: fc, fee, feeu, feeu_save,feeuu
      real(dp) :: fen, feni, feni_save, fenii
      real(dp) :: fenii_save, fi, fi_save, fii, fj, fj_save
      real(dp) :: fjj, fsum, fu, fui, fuj, fuu
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

      if (nelec.lt.2) goto 65

!     e-e and e-e-n terms
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

!     There are no C terms to order 1.
 22         if(nordc.le.1) goto 55

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
                           fc=fc+c(ll,it,iwf)*rij(k)*ss(l)*tt(m)
                           fu=fu+c(ll,it,iwf)*k*rij(k-1)*ss(l)*tt(m)
                           fuu=fuu+c(ll,it,iwf)*k*(k-1)*rij(k-2)*ss(l)*tt(m)
                           fi=fi+c(ll,it,iwf)*rij(k) &
                              *((l+m)*ri(l+m-1)*rj(m)+m*ri(m-1)*rj(l+m))
                           fii=fii+c(ll,it,iwf)*rij(k) &
                              *((l+m)*(l+m-1)*ri(l+m-2)*rj(m) &
                              +m*(m-1)*ri(m-2)*rj(l+m))
                           fj=fj+c(ll,it,iwf)*rij(k) &
                                *((l+m)*rj(l+m-1)*ri(m)+m*rj(m-1)*ri(l+m))
                           fjj=fjj+c(ll,it,iwf)*rij(k) &
                                *((l+m)*(l+m-1)*rj(l+m-2)*ri(m) &
                                +m*(m-1)*rj(m-2)*ri(l+m))
                           fui=fui+c(ll,it,iwf)*k*rij(k-1) &
                                *((l+m)*ri(l+m-1)*rj(m)+m*ri(m-1)*rj(l+m))
                           fuj=fuj+c(ll,it,iwf)*k*rij(k-1) &
                                *((l+m)*rj(l+m-1)*ri(m)+m*rj(m-1)*ri(l+m))
                        endif
!     write(ounit,'(''rij,ri,rj'',9f10.5)') rij,ri,rj,rij(1),ri(1),rj(1)
                     enddo
                  enddo
               enddo

               fi_save=fi*term+fc*termi
               fj_save=fj*term+fc*termj
               fii=fii*term+2*fi*termi+fc*termii
               fjj=fjj*term+2*fj*termj+fc*termjj

               fui=fu*termi
               fuj=fu*termj
               fu=fu*term/rij(1)
               fuu=fuu*term

               fi=fi_save/ri(1)
               fj=fj_save/rj(1)

               fso(i,j)=fso(i,j) + fc*term


               fijo(1,i,j)=fijo(1,i,j) + fi*rvec_en(1,i,ic)+fu*rvec_ee(1,ij)
               fijo(2,i,j)=fijo(2,i,j) + fi*rvec_en(2,i,ic)+fu*rvec_ee(2,ij)
               fijo(3,i,j)=fijo(3,i,j) + fi*rvec_en(3,i,ic)+fu*rvec_ee(3,ij)
               fijo(1,j,i)=fijo(1,j,i) + fj*rvec_en(1,j,ic)-fu*rvec_ee(1,ij)
               fijo(2,j,i)=fijo(2,j,i) + fj*rvec_en(2,j,ic)-fu*rvec_ee(2,ij)
               fijo(3,j,i)=fijo(3,j,i) + fj*rvec_en(3,j,ic)-fu*rvec_ee(3,ij)
!     write(ounit,'(''i,j,fijo2='',2i5,9d12.4)') i,j,(fijo(k,i,j),k=1,3)

               d2ijo(i,j)=d2ijo(i,j) + 2*(fuu + 2*fu) + fui*u2pst/(ri(1)*rij(1)) &
                    + fuj*u2mst/(rj(1)*rij(1)) + fii + 2*fi + fjj + 2*fj

 50            continue
            enddo

 55         fsumo=fsumo+fso(i,j)
            fjo(1,i)=fjo(1,i)+fijo(1,i,j)
            fjo(2,i)=fjo(2,i)+fijo(2,i,j)
            fjo(3,i)=fjo(3,i)+fijo(3,i,j)
            fjo(1,j)=fjo(1,j)+fijo(1,j,i)
            fjo(2,j)=fjo(2,j)+fijo(2,j,i)
            fjo(3,j)=fjo(3,j)+fijo(3,j,i)
            d2o=d2o+d2ijo(i,j)
         enddo
      enddo

!     e-n terms
 65   do i=1,nelec

         fso(i,i)=0
         fijo(1,i,i)=0
         fijo(2,i,i)=0
         fijo(3,i,i)=0
         d2ijo(i,i)=0

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
!     write(ounit,'(''fijo='',9d12.4)') (fijo(k,i,i),k=1,3),feni,rvec_en(1,i,ic)

            d2ijo(i,i) = d2ijo(i,i) + fenii + 2*feni

            if(iforce_analy.eq.1) call da_jastrow1(iwf,i,ic,it,rvec_en(1,i,ic),ri,feni_save,fenii)
 80         continue
         enddo

         fsumo=fsumo+fso(i,i)
         fjo(1,i)=fjo(1,i)+fijo(1,i,i)
         fjo(2,i)=fjo(2,i)+fijo(2,i,i)
         fjo(3,i)=fjo(3,i)+fijo(3,i,i)
!     write(ounit,'(''v='',9d12.4)') (fjo(k,i),k=1,3)
         d2o=d2o+d2ijo(i,i)
      enddo


      return
      end

!-----------------------------------------------------------------------
      subroutine da_jastrow1(iwf,i,ic,it,rvec_en,r,feni,fenii)

      use da_jastrow, only: da_d2j, da_j, da_vj
      use error, only: fatal_error
      use jastrow, only: a4, norda, nordj
      use scale_more, only: dd3
      use precision_kinds, only: dp
      use scale_more, only: dd3

      implicit none

      integer :: i, ic, iord, it, iwf
      integer :: k, l
      real(dp) :: feni, fenii, feniii
      real(dp) :: ri, ri2
      real(dp), dimension(3) :: rvec_en
      real(dp), dimension(-2:nordj) :: r

      call fatal_error('DA_JAS1: da_jastrow1 to fix')

      feniii=0.d0
      do iord=3,norda
        feniii=feniii+a4(iord,it,iwf)*iord*(iord-1)*(iord-2)*r(iord-3)
      enddo

      ri=1.d0/r(1)
      ri2=ri*ri

      ! can use iwf and loop over nwftypejas, since iwf is passed, this
      ! will only be called with 1 jas I believe, using iwf's original
      ! purpose

      do k=1,3
        da_j(k,i,ic)=-rvec_en(k)*ri*feni
        da_d2j(k,ic)=da_d2j(k,ic)-rvec_en(k)*ri*(feniii+fenii*(2*ri)+feni*(-2*ri2))
        do l=1,3
          da_vj(k,l,i,ic)=da_vj(k,l,i,ic)-rvec_en(k)*rvec_en(l)*ri2*(fenii-feni**ri)
        enddo
        da_vj(k,k,i,ic)=da_vj(k,k,i,ic)-feni*ri
      enddo

      return
      end
      end module
