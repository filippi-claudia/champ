      function psinl(u,rshifti,rshiftj,rri,rrj,it)
c Written by Claudia Filippi, modified by Cyrus Umrigar

      use jaspar, only: nspin1, nspin2, sspin, sspinn, is
      use jaspar1, only: cjas1, cjas2
      use elec, only: ndn, nup
      use jaspar3, only: a, b, c, fck, nord, scalek

      use jaspar4, only: a4, norda, nordb, nordc
      implicit real*8(a-h,o-z)





      include 'vmc.h'
      include 'force.h'

      parameter (one=1.d0,two=2.d0,half=0.5d0,eps=1.d-12)

      common /pars/ a00,a20,a21,eps_fock,c0000,c1110,c2000,
     &   xm1,xm2,xm12,xms,xma,Z
      common /contr2/ ijas,icusp,icusp2,isc,ianalyt_lap
     &,ifock,i3body,irewgt,iaver,istrch


      common /jaspar6/ cutjas,cutjasi,c1_jas6i,c1_jas6,c2_jas6,
     &asymp_r,asymp_jasa(MCTYPE),asymp_jasb(2)
      common /wfsec/ iwftype(MFORCE),iwf,nwftype
      common /chck/ bot

      dimension uu(0:MORDJ),ss(0:MORDJ),tt(0:MORDJ),rshifti(3),rshiftj(3)

c Not updated for ijas=5,6 because we will probably stay with ijas=4
c If we want to use ijas=5,6 update this routine similarly to psi.f
      if(ijas.ge.5) call fatal_error('PSINL: ijas >= 5 not implemented')

      psinl=0.d0
      if(nordc.le.1) return

      if(rri.eq.asymp_r .or. rrj.eq.asymp_r) return
      do 37 k=1,3
   37   if(abs(rshifti(k)-rshiftj(k)).gt.eps) return

      uuu=u
      rrri=rri
      rrrj=rrj
      call switch_scale(uuu)
      call switch_scale(rrri)
      call switch_scale(rrrj)

      uu(0)=one
      ss(0)=two
      tt(0)=one
      do 40 jp=1,nordc
        uu(jp)=uuu**jp
        ss(jp)=rrri**jp+rrrj**jp
   40   tt(jp)=(rrri*rrrj)**jp

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
              psinl=psinl+c(ll,it,iwf)*uu(k)*ss(l)*tt(m)
            endif
   50 continue

      return
      end
c-----------------------------------------------------------------------
      function psianl(rri,it)

      use jaspar3, only: a, b, c, fck, nord, scalek

      use jaspar4, only: a4, norda, nordb, nordc
      implicit real*8(a-h,o-z)


      include 'vmc.h'
      include 'force.h'
      common /contr2/ ijas,icusp,icusp2,isc,ianalyt_lap
     &,ifock,i3body,irewgt,iaver,istrch

      common /jaspar6/ cutjas,cutjasi,c1_jas6i,c1_jas6,c2_jas6,
     &asymp_r,asymp_jasa(MCTYPE),asymp_jasb(2)
      common /wfsec/ iwftype(MFORCE),iwf,nwftype

c Not updated for ijas=5,6 because we will probably stay with ijas=4
c If we want to use ijas=5,6 update this routine similarly to psi.f
      if(ijas.ge.5) call fatal_error('PSINL: ijas >= 5 not implemented')

      psianl=0.d0
      if(rri.eq.asymp_r) return

      psianl=a4(1,it,iwf)*rri/(1.d0+a4(2,it,iwf)*rri)-asymp_jasa(it)
      do 10 i=2,norda
   10   psianl=psianl+a4(i+1,it,iwf)*rri**i

      return
      end
c-----------------------------------------------------------------------
      function psibnl(u,isb,ipar)

      use jaspar, only: nspin1, nspin2, sspin, sspinn, is
      use jaspar3, only: a, b, c, fck, nord, scalek

      use jaspar4, only: a4, norda, nordb, nordc
      implicit real*8(a-h,o-z)



      include 'vmc.h'
      include 'force.h'
      common /contr2/ ijas,icusp,icusp2,isc,ianalyt_lap
     &,ifock,i3body,irewgt,iaver,istrch

      common /jaspar6/ cutjas,cutjasi,c1_jas6i,c1_jas6,c2_jas6,
     &asymp_r,asymp_jasa(MCTYPE),asymp_jasb(2)
      common /wfsec/ iwftype(MFORCE),iwf,nwftype

c Not updated for ijas=5,6 because we will probably stay with ijas=4
c If we want to use ijas=5,6 update this routine similarly to psi.f
      if(ijas.ge.5) call fatal_error('PSINL: ijas >= 5 not implemented')

      psibnl=0.d0
      if(ijas.le.2) return
      if(u.eq.asymp_r) return

      fee=b(1,isb,iwf)*u/(1+b(2,isb,iwf)*u)
      psibnl=sspinn*fee

      psibnl=psibnl-asymp_jasb(ipar+1)
      do 10 i=2,nordb
   10   psibnl=psibnl+b(i+1,isb,iwf)*u**i

      return
      end
c-----------------------------------------------------------------------
      function dpsianl(rri,it)

      use jaspar3, only: a, b, c, fck, nord, scalek

      use jaspar4, only: a4, norda, nordb, nordc
      implicit real*8(a-h,o-z)


      include 'vmc.h'
      include 'force.h'
      common /contr2/ ijas,icusp,icusp2,isc,ianalyt_lap
     &,ifock,i3body,irewgt,iaver,istrch

      common /jaspar6/ cutjas,cutjasi,c1_jas6i,c1_jas6,c2_jas6,
     &asymp_r,asymp_jasa(MCTYPE),asymp_jasb(2)
      common /wfsec/ iwftype(MFORCE),iwf,nwftype

c Not updated for ijas=5,6 because we will probably stay with ijas=4
c If we want to use ijas=5,6 update this routine similarly to psi.f
      if(ijas.ge.5) call fatal_error('PSINL: ijas >= 5 not implemented')

      dpsianl=0.d0
      if(rri.eq.asymp_r) return

      do 10 i=2,norda
   10   dpsianl=dpsianl+i*a4(i+1,it,iwf)*rri**(i-1)

      return
      end

c-----------------------------------------------------------------------
      function dpsibnl(u,isb,ipar)

      use jaspar, only: nspin1, nspin2, sspin, sspinn, is
      use jaspar3, only: a, b, c, fck, nord, scalek

      use jaspar4, only: a4, norda, nordb, nordc
      implicit real*8(a-h,o-z)



      include 'vmc.h'
      include 'force.h'
      common /contr2/ ijas,icusp,icusp2,isc,ianalyt_lap
     &,ifock,i3body,irewgt,iaver,istrch

      common /jaspar6/ cutjas,cutjasi,c1_jas6i,c1_jas6,c2_jas6,
     &asymp_r,asymp_jasa(MCTYPE),asymp_jasb(2)
      common /wfsec/ iwftype(MFORCE),iwf,nwftype

c Not updated for ijas=5,6 because we will probably stay with ijas=4
c If we want to use ijas=5,6 update this routine similarly to psi.f
      if(ijas.ge.5) call fatal_error('PSINL: ijas >= 5 not implemented')

      dpsibnl=0.d0
      if(u.eq.asymp_r) return

      top=b(1,isb,iwf)*u
      dtop=b(1,isb,iwf)
      bot=1+b(2,isb,iwf)*u
      dbot=b(2,isb,iwf)
      boti=1.d0/bot

      dfee=dtop*boti-top*boti*boti*dbot
      dpsibnl=sspinn*dfee

      do 10 i=2,nordb
   10   dpsibnl=dpsibnl+i*b(i+1,isb,iwf)*u**(i-1)

      return
      end
