      subroutine jastrow4e(iel,x,v,d2,value,iflag,istate) 
c     Written by Cyrus Umrigar and Claudia Filippi
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
      use jaspar3, only: a, b, c
      use jaspar4, only: a4, norda, nordb, nordc
      use jaspar6, only: asymp_jasa, asymp_jasb
      use jaspar6, only: cutjas
      use wfsec, only: iwf
      use bparm, only: nocuspb, nspin2b
      use contr2, only: ijas
      use contr2, only: isc
      use jasn, only: d2ijn, d2n, fijn, fjn, fsn, fsumn
      use distance_mod, only: rshift, r_en, rvec_en, r_ee, rvec_ee

      implicit real*8(a-h,o-z)

      parameter (half=.5d0,eps=1.d-12)
      dimension x(3,*),v(3,*)
      dimension uu(-2:MORDJ),ss(-2:MORDJ),tt(-2:MORDJ),rri(-2:MORDJ)
     &,rrj(-2:MORDJ)

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

      do i=1,nelec
         do k=1,3
            fjn(k,i,istate)=fjo(k,i,istate)
         enddo
      enddo

      fsumn(istate)=fsumo(istate)
      d2n(istate)=d2o(istate)

      do 60 jj=1,nelec

      if(jj.eq.iel) goto 60
      if(jj.lt.iel) then
        i=iel
        j=jj
       else
        i=jj
        j=iel
      endif
      ij=((i-1)*(i-2))/2+j

      fijn(1,i,j,istate)=0.0d0
      fijn(2,i,j,istate)=0.0d0
      fijn(3,i,j,istate)=0.0d0
      fijn(1,j,i,istate)=0.0d0
      fijn(2,j,i,istate)=0.0d0
      fijn(3,j,i,istate)=0.0d0
      fsn(i,j,istate)=0.0d0
      d2ijn(i,j,istate)=0.0d0

      sspinn=1.0d0
      ipar=0
      isb=1
      if(i.le.nup .or. j.gt.nup) then
         if(nspin2b.eq.2) then
            isb=2
          elseif(nocuspb.eq.0) then
            sspinn=half
         endif
         ipar=1
      endif

      rij=r_ee(ij)

      if(iflag.eq.0) then
        call scale_dist1(rij,uu(1),dd1,1)
       else
        call scale_dist2(rij,uu(1),dd1,dd2,1)
      endif

      if(rij.gt.cutjas) goto 22

      top=sspinn*b(1,isb,istate,iwf)*uu(1)
      topu=sspinn*b(1,isb,istate,iwf)

      bot=1+b(2,isb,istate,iwf)*uu(1)
      botu=b(2,isb,istate,iwf)
      bot2=bot*bot

      fee=top/bot-asymp_jasb(ipar+1,istate)
      feeu=topu/bot-botu*top/bot2

      do iord=2,nordb
         uu(iord)=uu(1)*uu(iord-1)
         fee=fee+b(iord+1,isb,istate,iwf)*uu(iord)
         feeu=feeu+b(iord+1,isb,istate,iwf)*iord*uu(iord-1)
      enddo

      if(iflag.gt.0) then
         topuu=0.0d0
         botuu=0.0d0
         feeuu=topuu-(botuu*top+2*botu*topu)/bot+2*botu**2*top/bot2
         feeuu=feeuu/bot
         do iord=2,nordb
            feeuu=feeuu+b(iord+1,isb,istate,iwf)*iord*(iord-1)*uu(iord-2)
         enddo
         feeuu=feeuu*dd1*dd1+feeu*dd2
      endif

      feeu=feeu*dd1/rij
      fsn(i,j,istate)=fsn(i,j,istate) + fee
      fijn(1,i,j,istate)=fijn(1,i,j,istate)+feeu*rvec_ee(1,ij)
      fijn(2,i,j,istate)=fijn(2,i,j,istate)+feeu*rvec_ee(2,ij)
      fijn(3,i,j,istate)=fijn(3,i,j,istate)+feeu*rvec_ee(3,ij)
      fijn(1,j,i,istate)=fijn(1,j,i,istate)-feeu*rvec_ee(1,ij)
      fijn(2,j,i,istate)=fijn(2,j,i,istate)-feeu*rvec_ee(2,ij)
      fijn(3,j,i,istate)=fijn(3,j,i,istate)-feeu*rvec_ee(3,ij)

      if(iflag.gt.0) d2ijn(i,j,istate)=d2ijn(i,j,istate) + 2*(feeuu+2*feeu)

c There are no C terms to order 1.
   22 if(nordc.le.1) goto 55

      if(iflag.eq.0) then
         if(isc.ge.12) call scale_dist1(rij,uu(1),dd1,3)
         if(ijas.eq.4.or.ijas.eq.5) call switch_scale1(uu(1),dd1)
      else
         if(isc.ge.12) call scale_dist2(rij,uu(1),dd1,dd2,3)
         if(ijas.eq.4.or.ijas.eq.5) call switch_scale2(uu(1),dd1,dd2)
      endif

      if(ijas.eq.4.or.ijas.eq.5) then
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

        if(iflag.eq.0) then
           call scale_dist1(ri,rri(1),dd7,2)
           call scale_dist1(rj,rrj(1),dd8,2)
           if(ijas.eq.4.or.ijas.eq.5) then
              call switch_scale1(rri(1),dd7)
              call switch_scale1(rrj(1),dd8)
           endif
        else
           call scale_dist2(ri,rri(1),dd7,dd9,2)
           call scale_dist2(rj,rrj(1),dd8,dd10,2)
           if(ijas.eq.4.or.ijas.eq.5) then
              call switch_scale2(rri(1),dd7,dd9)
              call switch_scale2(rrj(1),dd8,dd10)
           endif
        endif

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
                    fi=fi+c(ll,it,istate,iwf)*uu(k)
     &              *((l+m)*rri(l+m-1)*rrj(m)+m*rri(m-1)*rrj(l+m))
                    fj=fj+c(ll,it,istate,iwf)*uu(k)
     &              *((l+m)*rrj(l+m-1)*rri(m)+m*rrj(m-1)*rri(l+m))

                    if(iflag.gt.0) then
                       fuu=fuu+c(ll,it,istate,iwf)*k*(k-1)*uu(k-2)*ss(l)*tt(m)
                       fii=fii+c(ll,it,istate,iwf)*uu(k)
     &                 *((l+m)*(l+m-1)*rri(l+m-2)*rrj(m)
     &                 +m*(m-1)*rri(m-2)*rrj(l+m))
                       fjj=fjj+c(ll,it,istate,iwf)*uu(k)
     &                 *((l+m)*(l+m-1)*rrj(l+m-2)*rri(m)
     &                 +m*(m-1)*rrj(m-2)*rri(l+m))
                       fui=fui+c(ll,it,istate,iwf)*k*uu(k-1)
     &                 *((l+m)*rri(l+m-1)*rrj(m)+m*rri(m-1)*rrj(l+m))
                       fuj=fuj+c(ll,it,istate,iwf)*k*uu(k-1)
     &                 *((l+m)*rrj(l+m-1)*rri(m)+m*rrj(m-1)*rri(l+m))
                    endif
                 endif
              enddo
           enddo
        enddo

        if(iflag.gt.0) then
           s=ri+rj
           t=ri-rj
           u2pst=rij*rij+s*t
           u2mst=rij*rij-s*t

           fuu=fuu*dd1*dd1+fu*dd2
           fui=fui*dd1*dd7
           fuj=fuj*dd1*dd8
           fii=fii*dd7*dd7+fi*dd9
           fjj=fjj*dd8*dd8+fj*dd10
        endif

        fu=fu*dd1/rij
        fi=fi*dd7/ri
        fj=fj*dd8/rj

        fsn(i,j,istate)=fsn(i,j,istate) + fc

        fijn(1,i,j,istate)=fijn(1,i,j,istate)+fi*rvec_en(1,i,ic)+fu*rvec_ee(1,ij)
        fijn(2,i,j,istate)=fijn(2,i,j,istate)+fi*rvec_en(2,i,ic)+fu*rvec_ee(2,ij)
        fijn(3,i,j,istate)=fijn(3,i,j,istate)+fi*rvec_en(3,i,ic)+fu*rvec_ee(3,ij)
        fijn(1,j,i,istate)=fijn(1,j,i,istate)+fj*rvec_en(1,j,ic)-fu*rvec_ee(1,ij)
        fijn(2,j,i,istate)=fijn(2,j,i,istate)+fj*rvec_en(2,j,ic)-fu*rvec_ee(2,ij)
        fijn(3,j,i,istate)=fijn(3,j,i,istate)+fj*rvec_en(3,j,ic)-fu*rvec_ee(3,ij)

        if(iflag.gt.0) 
     &    d2ijn(i,j,istate)=d2ijn(i,j,istate) + 2*(fuu + 2*fu) + fui*u2pst/(ri*rij)
     &    + fuj*u2mst/(rj*rij) + fii + 2*fi + fjj + 2*fj
   50 continue

   55 fsumn(istate)=fsumn(istate)+fsn(i,j,istate)-fso(i,j,istate)
      fjn(1,i,istate)=fjn(1,i,istate)+fijn(1,i,j,istate)-fijo(1,i,j,istate)
      fjn(2,i,istate)=fjn(2,i,istate)+fijn(2,i,j,istate)-fijo(2,i,j,istate)
      fjn(3,i,istate)=fjn(3,i,istate)+fijn(3,i,j,istate)-fijo(3,i,j,istate)
      fjn(1,j,istate)=fjn(1,j,istate)+fijn(1,j,i,istate)-fijo(1,j,i,istate)
      fjn(2,j,istate)=fjn(2,j,istate)+fijn(2,j,i,istate)-fijo(2,j,i,istate)
      fjn(3,j,istate)=fjn(3,j,istate)+fijn(3,j,i,istate)-fijo(3,j,i,istate)
      d2n(istate)=d2n(istate)+d2ijn(i,j,istate)-d2ijo(i,j,istate)
   60 continue

c e-n terms

   65 fijn(1,iel,iel,istate)=0.0d0
      fijn(2,iel,iel,istate)=0.0d0
      fijn(3,iel,iel,istate)=0.0d0
      fsn(iel,iel,istate)=0.0d0
      d2ijn(iel,iel,istate)=0.0d0

      do ic=1,ncent
         it=iwctype(ic)
         ri=r_en(iel,ic)

         if(ri.gt.cutjas) cycle

         if(iflag.eq.0) then
           call scale_dist1(ri,rri(1),dd7,1)
          else
           call scale_dist2(ri,rri(1),dd7,dd9,1)
         endif

         top=a4(1,it,istate,iwf)*rri(1)
         topi=a4(1,it,istate,iwf)

         bot=a4(2,it,istate,iwf)*rri(1)
         boti=a4(2,it,istate,iwf)

         bot=1+bot
         bot2=bot*bot
         fen=top/bot-asymp_jasa(it,istate)
         feni=topi/bot-boti*top/bot2

         do iord=2,norda
            rri(iord)=rri(1)**iord
            fen=fen+a4(iord+1,it,istate,iwf)*rri(iord)
            feni=feni+a4(iord+1,it,istate,iwf)*iord*rri(iord-1)
         enddo

         if(iflag.gt.0) then
            topii=0.0d0
            botii=0.0d0
            fenii=topii-(botii*top+2*boti*topi)/bot+2*boti**2*top/bot2
            fenii=fenii/bot
            do iord=2,norda
               fenii=fenii+a4(iord+1,it,istate,iwf)*iord*(iord-1)*rri(iord-2)
            enddo
            fenii=fenii*dd7*dd7+feni*dd9
         endif

         feni=feni*dd7/ri

         fsn(iel,iel,istate)=fsn(iel,iel,istate)+fen
         fijn(1,iel,iel,istate)=fijn(1,iel,iel,istate)+feni*rvec_en(1,iel,ic)
         fijn(2,iel,iel,istate)=fijn(2,iel,iel,istate)+feni*rvec_en(2,iel,ic)
         fijn(3,iel,iel,istate)=fijn(3,iel,iel,istate)+feni*rvec_en(3,iel,ic)

         if(iflag.gt.0) d2ijn(iel,iel,istate) = d2ijn(iel,iel,istate) + fenii + 2*feni
      enddo

      fsumn(istate)=fsumn(istate)+fsn(iel,iel,istate)-fso(iel,iel,istate)
      fjn(1,iel,istate)=fjn(1,iel,istate)+fijn(1,iel,iel,istate)-fijo(1,iel,iel,istate)
      fjn(2,iel,istate)=fjn(2,iel,istate)+fijn(2,iel,iel,istate)-fijo(2,iel,iel,istate)
      fjn(3,iel,istate)=fjn(3,iel,istate)+fijn(3,iel,iel,istate)-fijo(3,iel,iel,istate)
      d2n(istate)=d2n(istate)+d2ijn(iel,iel,istate)-d2ijo(iel,iel,istate)

      do i=1,nelec
         v(1,i)=fjn(1,i,istate)
         v(2,i)=fjn(2,i,istate)
         v(3,i)=fjn(3,i,istate)
      enddo
      value=fsumn(istate)
      d2=d2n(istate)

      end subroutine
