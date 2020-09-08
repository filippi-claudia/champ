      function psi(rij,ri,rj,it)
c Written by Cyrus Umrigar, modified by Claudia Filippi
c **Warning** This routine needs to be upgraded to check rshifts
c if we add in the capability to use numerical Laplacian for
c periodic systems.

      use force_mod, only: MFORCE, MFORCE_WT_PRD, MWF
      use vmc, only: MELEC, MORB, MBASIS, MDET, MCENT, MCTYPE, MCTYP3X
      use vmc, only: NSPLIN, nrad, MORDJ, MORDJ1, MMAT_DIM, MMAT_DIM2, MMAT_DIM20
      use vmc, only: radmax, delri
      use vmc, only: NEQSX, MTERMS
      use vmc, only: MCENT3, NCOEF, MEXCIT
      use jaspar3, only: a, c

      use jaspar4, only: nordc
      use jaspar6, only: cutjas
      use wfsec, only: iwf
      use contr2, only: ijas
      implicit real*8(a-h,o-z)












      parameter (zero=0.d0,one=1.d0,two=2.d0)





      dimension uu(0:MORDJ),ss(0:MORDJ),tt(0:MORDJ)

      psi=0

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
      ss(0)=2
      tt(0)=one
      do 40 jp=1,nordc
        uu(jp)=u**jp
        ss(jp)=rri**jp+rrj**jp
   40   tt(jp)=(rri*rrj)**jp

      ll=0
      do 50 n=2,nordc
        do 50 k=n-1,0,-1
          if(k.eq.0) then
            l_hi=n-k-2
           else
            l_hi=n-k
          endif
          do 50 l=l_hi,0,-1
            m=(n-k-l)/2
            if(2*m.eq.n-k-l) then
              ll=ll+1
              psi=psi+c(ll,it,iwf)*uu(k)*ss(l)*tt(m)
            endif
   50 continue

      return
      end

c-----------------------------------------------------------------------
      function psia(ri,it)

      use force_mod, only: MFORCE, MFORCE_WT_PRD, MWF
      use vmc, only: MELEC, MORB, MBASIS, MDET, MCENT, MCTYPE, MCTYP3X
      use vmc, only: NSPLIN, nrad, MORDJ, MORDJ1, MMAT_DIM, MMAT_DIM2, MMAT_DIM20
      use vmc, only: radmax, delri
      use vmc, only: NEQSX, MTERMS
      use vmc, only: MCENT3, NCOEF, MEXCIT
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
        psia=a4(1,it,iwf)*rri/(one+a4(2,it,iwf)*rri)-asymp_jasa(it)
        do 10 i=2,norda
   10     psia=psia+a4(i+1,it,iwf)*rri**i
       elseif(ijas.eq.6) then
        psia=a4(1,it,iwf)*rri/(one+a4(2,it,iwf)*(1-rri))
        do 20 i=2,norda
   20     psia=psia+a4(i+1,it,iwf)*rri**i
      endif

      return
      end
c-----------------------------------------------------------------------
      function psib(rij,isb,ipar)

      use force_mod, only: MFORCE, MFORCE_WT_PRD, MWF
      use vmc, only: MELEC, MORB, MBASIS, MDET, MCENT, MCTYPE, MCTYP3X
      use vmc, only: NSPLIN, nrad, MORDJ, MORDJ1, MMAT_DIM, MMAT_DIM2, MMAT_DIM20
      use vmc, only: radmax, delri
      use vmc, only: NEQSX, MTERMS
      use vmc, only: MCENT3, NCOEF, MEXCIT
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
        psib=sspinn*b(1,isb,iwf)*u/(one+b(2,isb,iwf)*u)-asymp_jasb(ipar+1)
        do 10 i=2,nordb
   10     psib=psib+b(i+1,isb,iwf)*u**i
       elseif(ijas.eq.5) then
        psib=b(1,isb,iwf)*u/(one+b(2,isb,iwf)*u)-asymp_jasb(ipar+1)
        do 20 i=2,nordb
   20     psib=psib+b(i+1,isb,iwf)*u**i
        psib=sspinn*psib
       elseif(ijas.eq.6) then
        psib=b(1,isb,iwf)*u/(one+b(2,isb,iwf)*(1-u))
        do 30 i=2,nordb
   30     psib=psib+b(i+1,isb,iwf)*u**i
        psib=sspinn*psib
      endif

      return
      end
