      function psi(rij,ri,rj,it,istate)
c     Written by Cyrus Umrigar, modified by Claudia Filippi
c     **Warning** This routine needs to be upgraded to check rshifts
c     if we add in the capability to use numerical Laplacian for
c     periodic systems.
      use force_mod, only: MFORCE, MFORCE_WT_PRD, MWF
      use vmc_mod, only: MELEC, MORB, MBASIS, MDET, MCENT, MCTYPE, MCTYP3X
      use vmc_mod, only: NSPLIN, nrad, MORDJ, MORDJ1, MMAT_DIM, MMAT_DIM2, MMAT_DIM20
      use vmc_mod, only: radmax, delri
      use vmc_mod, only: NEQSX, MTERMS
      use vmc_mod, only: MCENT3, NCOEF, MEXCIT
      use jaspar3, only: a, c
      use jaspar4, only: nordc
      use jaspar6, only: cutjas
      use wfsec, only: iwf
      use contr2, only: ijas

      implicit real*8(a-h,o-z)

      parameter (zero=0.d0,one=1.d0,two=2.d0)
      dimension uu(0:MORDJ),ss(0:MORDJ),tt(0:MORDJ)

      psi=0.0d0

      u=rij
      s=ri+rj
      t=dabs(ri-rj)

      call scale_dist(rij,u,3)
      call scale_dist(ri,rri,2)
      call scale_dist(rj,rrj,2)

      if(nordc.le.1) return

      if(ri.gt.cutjas .or. rj.gt.cutjas) return

      if(ijas.eq.4.or.ijas.eq.5) then
         call switch_scale(u)
         call switch_scale(rri)
         call switch_scale(rrj)
      endif

      uu(0)=one
      ss(0)=2.0d0
      tt(0)=one
      do jp=1,nordc
         uu(jp)=u**jp
         ss(jp)=rri**jp+rrj**jp
         tt(jp)=(rri*rrj)**jp
      enddo

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
                  psi=psi+c(ll,it,istate,iwf)*uu(k)*ss(l)*tt(m)
               endif
            enddo
         enddo
      enddo

      return

      end function

c-----------------------------------------------------------------------

      function psia(ri,it,istate)
      use force_mod, only: MFORCE, MFORCE_WT_PRD, MWF
      use vmc_mod, only: MELEC, MORB, MBASIS, MDET, MCENT, MCTYPE, MCTYP3X
      use vmc_mod, only: NSPLIN, nrad, MORDJ, MORDJ1, MMAT_DIM, MMAT_DIM2, MMAT_DIM20
      use vmc_mod, only: radmax, delri
      use vmc_mod, only: NEQSX, MTERMS
      use vmc_mod, only: MCENT3, NCOEF, MEXCIT
      use jaspar3, only: a
      use jaspar4, only: a4, norda
      use jaspar6, only: asymp_jasa
      use jaspar6, only: cutjas
      use wfsec, only: iwf
      use contr2, only: ijas

      implicit real*8(a-h,o-z)

      parameter(zero=0.d0,one=1.d0)

      psia=zero

      if(ijas.lt.4.or.ijas.gt.6) return

      if(ri.gt.cutjas) return

      call scale_dist(ri,rri,1)

      if(ijas.eq.4.or.ijas.eq.5) then
         psia=a4(1,it,istate,iwf)*rri/(one+a4(2,it,istate,iwf)*rri)-asymp_jasa(it,istate)
         do i=2,norda
            psia=psia+a4(i+1,it,istate,iwf)*rri**i
         enddo
      elseif(ijas.eq.6) then
         psia=a4(1,it,istate,iwf)*rri/(one+a4(2,it,istate,iwf)*(1-rri))
         do i=2,norda
            psia=psia+a4(i+1,it,istate,iwf)*rri**i
         enddo
      endif

      return

      end function

c-----------------------------------------------------------------------

      function psib(rij,isb,ipar)
      use force_mod, only: MFORCE, MFORCE_WT_PRD, MWF
      use vmc_mod, only: MELEC, MORB, MBASIS, MDET, MCENT, MCTYPE, MCTYP3X
      use vmc_mod, only: NSPLIN, nrad, MORDJ, MORDJ1, MMAT_DIM, MMAT_DIM2, MMAT_DIM20
      use vmc_mod, only: radmax, delri
      use vmc_mod, only: NEQSX, MTERMS
      use vmc_mod, only: MCENT3, NCOEF, MEXCIT
      use jaspar, only: sspinn
      use jaspar3, only: a, b
      use jaspar4, only: nordb
      use jaspar6, only: asymp_jasb
      use jaspar6, only: cutjas
      use wfsec, only: iwf
      use contr2, only: ijas

      implicit real*8(a-h,o-z)

      parameter(zero=0.d0,one=1.d0)

      psib=zero
      if(ijas.lt.4.or.ijas.gt.6) return

      if(rij.gt.cutjas) return

      u=rij
      call scale_dist(rij,u,1)

      if(ijas.eq.4) then
         psib=sspinn*b(1,isb,istate,iwf)*u/(one+b(2,isb,istate,iwf)*u)-asymp_jasb(ipar+1,istate)
         do i=2,nordb
            psib=psib+b(i+1,isb,istate,iwf)*u**i
         enddo
      elseif(ijas.eq.5) then
         psib=b(1,isb,istate,iwf)*u/(one+b(2,isb,istate,iwf)*u)-asymp_jasb(ipar+1,istate)
         do i=2,nordb
            psib=psib+b(i+1,isb,istate,iwf)*u**i
         enddo
         psib=sspinn*psib
      elseif(ijas.eq.6) then
         psib=b(1,isb,istate,iwf)*u/(one+b(2,isb,istate,iwf)*(1-u))
         do i=2,nordb
            psib=psib+b(i+1,isb,istate,iwf)*u**i
         enddo
         psib=sspinn*psib
      endif

      return

      end function
