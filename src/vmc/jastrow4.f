      subroutine jastrow4(x,v,d2,div_vj,value,istate)
c     Written by Cyrus Umrigar, modified by C. Filippi
c     Jastrow 4,5 must be used with one of isc=2,4,6,7,12,14,16,17
c     Jastrow 6   must be used with one of isc=6,7
      use force_mod, only: MFORCE, MFORCE_WT_PRD, MWF
      use vmc_mod, only: MELEC, MORB, MBASIS, MDET, MCENT, MCTYPE, MCTYP3X
      use vmc_mod, only: NSPLIN, nrad, MORDJ, MORDJ1, MMAT_DIM, MMAT_DIM2, MMAT_DIM20
      use vmc_mod, only: radmax, delri
      use vmc_mod, only: NEQSX, MTERMS
      use vmc_mod, only: MCENT3, NCOEF, MEXCIT
      use atom, only: iwctype, ncent
      use jaspar, only: sspinn
      use const, only: nelec
      use elec, only: nup
      use jaso, only: d2ijo, d2o, fijo, fjo, fso, fsumo
      use jaspar3, only: a, b, c, scalek
      use jaspar4, only: a4, norda, nordb, nordc
      use jaspar6, only: asymp_jasa, asymp_jasb, c1_jas6
      use jaspar6, only: cutjas
      use wfsec, only: iwf
      use bparm, only: nocuspb, nspin2b
      use contr2, only: ijas
      use contr2, only: isc
      use force_analy, only: iforce_analy
      use distance_mod, only: rshift, r_en, rvec_en, r_ee, rvec_ee

      implicit real*8(a-h,o-z)
      
      parameter (half=.5d0,eps=1.d-12)
      dimension x(3,*),v(3,*),div_vj(*)
      dimension uu(-2:MORDJ),ss(-2:MORDJ),tt(-2:MORDJ),rri(-2:MORDJ)
     &,rrj(-2:MORDJ)

      fsum=0.0d0
      do i=-2,-1
         uu(i)=0.0d0
         ss(i)=0.0d0
         tt(i)=0.0d0
         rri(i)=0.0d0
         rrj(i)=0.0d0
      enddo
      uu(0)=1.0d0
      ss(0)=2.0d0
      tt(0)=1.0d0
      rri(0)=1.0d0
      rrj(0)=1.0d0

      if (nelec.lt.2) goto 65

c e-e and e-e-n terms
      ij=0
      do 60 i=2,nelec
      im1=i-1
      do 60 j=1,im1
      ij=ij+1

      fso(i,j,istate)=0.0d0
      d2ijo(i,j,istate)=0.0d0
      do k=1,3
        fijo(k,i,j,istate)=0.0d0
        fijo(k,j,i,istate)=0.0d0
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

      call scale_dist2(rij,uu(1),dd1,dd2,1)
c     write(6,'(''rij,u in ee'',2f9.5)') rij,uu(1)

c Check rij after scaling because uu(1) used in e-e-n terms too
      if(rij.gt.cutjas) goto 22

      top=sspinn*b(1,isb,istate,iwf)*uu(1)
      topu=sspinn*b(1,isb,istate,iwf)
      topuu=0.0d0

      if(ijas.eq.4.or.ijas.eq.5) then
         bot=1+b(2,isb,istate,iwf)*uu(1)
         botu=b(2,isb,istate,iwf)
       elseif(ijas.eq.6) then
         bot=1+b(2,isb,istate,iwf)*(1-uu(1))
         botu=-b(2,isb,istate,iwf)
      endif
      botuu=0.0d0
      bot2=bot*bot

      fee=top/bot-asymp_jasb(ipar+1,istate)
      feeu=topu/bot-botu*top/bot2
      feeuu=topuu-(botuu*top+2*botu*topu)/bot+2*botu**2*top/bot2
      feeuu=feeuu/bot

      do iord=2,nordb
         uu(iord)=uu(1)*uu(iord-1)
         if(ijas.eq.4) then
            fee=fee+b(iord+1,isb,istate,iwf)*uu(iord)
            feeu=feeu+b(iord+1,isb,istate,iwf)*iord*uu(iord-1)
            feeuu=feeuu+b(iord+1,isb,istate,iwf)*iord*(iord-1)*uu(iord-2)
         elseif(ijas.eq.5.or.ijas.eq.6) then
            fee=fee+sspinn*b(iord+1,isb,istate,iwf)*uu(iord)
            feeu=feeu+sspinn*b(iord+1,isb,istate,iwf)*iord*uu(iord-1)
            feeuu=feeuu+sspinn*b(iord+1,isb,istate,iwf)*iord*(iord-1)*uu(iord-2)
         endif
      enddo

      feeuu=feeuu*dd1*dd1+feeu*dd2
      feeu=feeu*dd1/rij

      fso(i,j,istate)=fso(i,j,istate)+fee
      do k=1,3
         fijo(k,i,j,istate)= fijo(k,i,j,istate) + feeu*rvec_ee(k,ij)
         fijo(k,j,i,istate)= fijo(k,j,i,istate) - feeu*rvec_ee(k,ij)
      enddo
      d2ijo(i,j,istate)=d2ijo(i,j,istate)+2*(feeuu+2*feeu)

c There are no C terms to order 1.
   22 if(nordc.le.1) goto 55

      if(isc.ge.12) call scale_dist2(rij,uu(1),dd1,dd2,3)
      if(ijas.eq.4.or.ijas.eq.5) then
         call switch_scale2(uu(1),dd1,dd2)
         do iord=2,nordc
            uu(iord)=uu(1)*uu(iord-1)
         enddo
      endif

      do 50 ic=1,ncent
        it=iwctype(ic)
        ri=r_en(i,ic)
        rj=r_en(j,ic)

        if(ri.gt.cutjas .or. rj.gt.cutjas) goto 50
        do k=1,3
           if(abs(rshift(k,i,ic)-rshift(k,j,ic)).gt.eps) goto 50
        enddo

        call scale_dist2(ri,rri(1),dd7,dd9,2)
        call scale_dist2(rj,rrj(1),dd8,dd10,2)

        if(ijas.eq.4.or.ijas.eq.5) then
           call switch_scale2(rri(1),dd7,dd9)
           call switch_scale2(rrj(1),dd8,dd10)
        endif

        s=ri+rj
        t=ri-rj
        u2pst=rij*rij+s*t
        u2mst=rij*rij-s*t

        do iord=1,nordc
           rri(iord)=rri(1)*rri(iord-1)
           rrj(iord)=rrj(1)*rrj(iord-1)
           ss(iord)=rri(iord)+rrj(iord)
           tt(iord)=rri(iord)*rrj(iord)
        enddo

        fc=0.0d0
        fu=0.0d0
        fuu=0.0d0
        fi=0.0d0
        fii=0.0d0
        fj=0.0d0
        fjj=0.0d0
        fui=0.0d0
        fuj=0.0d0

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
                    fc=fc+c(ll,it,istate,iwf)*uu(k)*ss(l)*tt(m)
                    fu=fu+c(ll,it,istate,iwf)*k*uu(k-1)*ss(l)*tt(m)
                    fuu=fuu+c(ll,it,istate,iwf)*k*(k-1)*uu(k-2)*ss(l)*tt(m)
                    fi=fi+c(ll,it,istate,iwf)*uu(k)
     &              *((l+m)*rri(l+m-1)*rrj(m)+m*rri(m-1)*rrj(l+m))
                    fii=fii+c(ll,it,istate,iwf)*uu(k)
     &              *((l+m)*(l+m-1)*rri(l+m-2)*rrj(m)
     &              +m*(m-1)*rri(m-2)*rrj(l+m))
                    fj=fj+c(ll,it,istate,iwf)*uu(k)
     &              *((l+m)*rrj(l+m-1)*rri(m)+m*rrj(m-1)*rri(l+m))
                    fjj=fjj+c(ll,it,istate,iwf)*uu(k)
     &              *((l+m)*(l+m-1)*rrj(l+m-2)*rri(m)
     &              +m*(m-1)*rrj(m-2)*rri(l+m))
                    fui=fui+c(ll,it,istate,iwf)*k*uu(k-1)
     &              *((l+m)*rri(l+m-1)*rrj(m)+m*rri(m-1)*rrj(l+m))
                    fuj=fuj+c(ll,it,istate,iwf)*k*uu(k-1)
     &              *((l+m)*rrj(l+m-1)*rri(m)+m*rrj(m-1)*rri(l+m))
                 endif
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

        fso(i,j,istate)=fso(i,j,istate) + fc

        fijo(1,i,j,istate)=fijo(1,i,j,istate) + fi*rvec_en(1,i,ic)+fu*rvec_ee(1,ij)
        fijo(2,i,j,istate)=fijo(2,i,j,istate) + fi*rvec_en(2,i,ic)+fu*rvec_ee(2,ij)
        fijo(3,i,j,istate)=fijo(3,i,j,istate) + fi*rvec_en(3,i,ic)+fu*rvec_ee(3,ij)
        fijo(1,j,i,istate)=fijo(1,j,i,istate) + fj*rvec_en(1,j,ic)-fu*rvec_ee(1,ij)
        fijo(2,j,i,istate)=fijo(2,j,i,istate) + fj*rvec_en(2,j,ic)-fu*rvec_ee(2,ij)
        fijo(3,j,i,istate)=fijo(3,j,i,istate) + fj*rvec_en(3,j,ic)-fu*rvec_ee(3,ij)
c       write(6,'(''i,j,fijo2='',2i5,9d12.4)') i,j,(fijo(k,i,j),k=1,3)

        d2ijo(i,j,istate)=d2ijo(i,j,istate) + 2*(fuu + 2*fu) + fui*u2pst/(ri*rij)
     &  + fuj*u2mst/(rj*rij) + fii + 2*fi + fjj + 2*fj

  50  continue

  55  fsum=fsum+fso(i,j,istate)
      v(1,i)=v(1,i)+fijo(1,i,j,istate)
      v(2,i)=v(2,i)+fijo(2,i,j,istate)
      v(3,i)=v(3,i)+fijo(3,i,j,istate)
      v(1,j)=v(1,j)+fijo(1,j,i,istate)
      v(2,j)=v(2,j)+fijo(2,j,i,istate)
      v(3,j)=v(3,j)+fijo(3,j,i,istate)
      div_vj(i)=div_vj(i)+d2ijo(i,j,istate)/2
      div_vj(j)=div_vj(j)+d2ijo(i,j,istate)/2
  60  d2=d2+d2ijo(i,j,istate)

c e-n terms
  65  do 90 i=1,nelec

        fso(i,i,istate)=0.0d0
        fijo(1,i,i,istate)=0.0d0
        fijo(2,i,i,istate)=0.0d0
        fijo(3,i,i,istate)=0.0d0
        d2ijo(i,i,istate)=0.0d0

        do ic=1,ncent
           it=iwctype(ic)
           ri=r_en(i,ic)

           if(ri.gt.cutjas) cycle

           call scale_dist2(ri,rri(1),dd7,dd9,1)

           top=a4(1,it,istate,iwf)*rri(1)
           topi=a4(1,it,istate,iwf)
           topii=0.0d0

           if(ijas.eq.4.or.ijas.eq.5) then
              bot=1+a4(2,it,istate,iwf)*rri(1)
              boti=a4(2,it,istate,iwf)
            elseif(ijas.eq.6) then
              bot=1+a4(2,it,istate,iwf)*(1-rri(1))
              boti=-a4(2,it,istate,iwf)
           endif

           botii=0.0d0
           bot2=bot*bot

           fen=top/bot-asymp_jasa(it,istate)
           feni=topi/bot-boti*top/bot2
           fenii=topii-(botii*top+2*boti*topi)/bot+2*boti**2*top/bot2
           fenii=fenii/bot

           do iord=2,norda
              rri(iord)=rri(1)**iord
              fen=fen+a4(iord+1,it,istate,iwf)*rri(iord)
              feni=feni+a4(iord+1,it,istate,iwf)*iord*rri(iord-1)
              fenii=fenii+a4(iord+1,it,istate,iwf)*iord*(iord-1)*rri(iord-2)
           enddo

           feni_save=feni
           fenii_save=fenii
           fenii=fenii*dd7*dd7+feni*dd9
           feni=feni*dd7/ri

           fso(i,i,istate)=fso(i,i,istate)+fen

           fijo(1,i,i,istate)=fijo(1,i,i,istate) + feni*rvec_en(1,i,ic)
           fijo(2,i,i,istate)=fijo(2,i,i,istate) + feni*rvec_en(2,i,ic)
           fijo(3,i,i,istate)=fijo(3,i,i,istate) + feni*rvec_en(3,i,ic)

           d2ijo(i,i,istate) = d2ijo(i,i,istate) + fenii + 2*feni

           if(iforce_analy.eq.1) then
              call da_jastrow4(iwf,i,ic,it,rvec_en(1,i,ic),
     &                        ri,rri,feni_save,fenii_save,dd7,dd9,istate)
           endif
        enddo

        fsum=fsum+fso(i,i,istate)
        v(1,i)=v(1,i)+fijo(1,i,i,istate)
        v(2,i)=v(2,i)+fijo(2,i,i,istate)
        v(3,i)=v(3,i)+fijo(3,i,i,istate)
        div_vj(i)=div_vj(i)+d2ijo(i,i,istate)
   90   d2=d2+d2ijo(i,i,istate)

      if(ijas.eq.6) then
         term=1/(c1_jas6*scalek(iwf))
         fsum=term*fsum
         d2=term*d2
         do i=1,nelec
            div_vj(i)=term*div_vj(i)
            do k=1,3
               v(k,i)=term*v(k,i)
            enddo
            do j=1,nelec
               d2ijo(i,j,istate)=term*d2ijo(i,j,istate)
               do k=1,3
                  fijo(k,i,j,istate)=term*fijo(k,i,j,istate)
               enddo
            enddo
         enddo
      endif

      fsumo(istate)=fsum
      d2o(istate)=d2
      do i=1,nelec
         fjo(1,i,istate)=v(1,i)
         fjo(2,i,istate)=v(2,i)
         fjo(3,i,istate)=v(3,i)
      enddo

      value=fsum

      end subroutine

c-----------------------------------------------------------------------

      function nterms4(nord)
c     Written by Cyrus Umrigar
      implicit real*8(a-h,o-z)

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

      return

      end function

c-----------------------------------------------------------------------

      subroutine da_jastrow4(iwf,i,ic,it,rvec_en,r,rr,feni,fenii,dd1,dd2,istate)
      use force_mod, only: MFORCE, MFORCE_WT_PRD, MWF
      use vmc_mod, only: MELEC, MORB, MBASIS, MDET, MCENT, MCTYPE, MCTYP3X
      use vmc_mod, only: NSPLIN, nrad, MORDJ, MORDJ1, MMAT_DIM, MMAT_DIM2, MMAT_DIM20
      use vmc_mod, only: radmax, delri
      use vmc_mod, only: NEQSX, MTERMS
      use vmc_mod, only: MCENT3, NCOEF, MEXCIT
      use da_jastrow4val, only: da_d2j, da_j, da_vj
      use jaspar4, only: a4, norda
      use scale_more, only: dd3

      implicit real*8(a-h,o-z)

      dimension rvec_en(3),rr(-2:MORDJ)

      feniii=0.0d0
      do iord=3,norda
         feniii=feniii+a4(iord+1,it,istate,iwf)*iord*(iord-1)*(iord-2)*rr(iord-3)
      enddo

      ri=1.0d0/r
      ri2=ri*ri

      do k=1,3
         da_j(k,i,ic,istate)=-rvec_en(k)*ri*feni*dd1
         da_d2j(k,i,ic,istate)=-rvec_en(k)*ri*(feniii*dd1*dd1*dd1+
     &        fenii*dd1*(3*dd2+2*dd1*ri)+feni*(dd3+2*dd2*ri-2*dd1*ri2))
         do l=1,3
            da_vj(k,l,i,ic,istate)=-rvec_en(k)*rvec_en(l)*ri2*(fenii*dd1*dd1+feni*dd2-feni*dd1*ri)
         enddo
         da_vj(k,k,i,ic,istate)=da_vj(k,k,i,ic,istate)-feni*dd1*ri
      enddo

      end subroutine
