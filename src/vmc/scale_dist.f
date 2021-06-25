      subroutine set_scale_dist(ipr)
c Written by Cyrus Umrigar
      use force_mod, only: MFORCE, MFORCE_WT_PRD, MWF
      use vmc_mod, only: MELEC, MORB, MBASIS, MDET, MCENT, MCTYPE, MCTYP3X
      use vmc_mod, only: NSPLIN, nrad, MORDJ, MORDJ1, MMAT_DIM, MMAT_DIM2, MMAT_DIM20
      use vmc_mod, only: radmax, delri
      use vmc_mod, only: NEQSX, MTERMS
      use vmc_mod, only: MCENT3, NCOEF, MEXCIT
      use atom, only: nctype

      use jaspar, only: sspinn, is
      use jaspar3, only: a, b, c, scalek

      use jaspar4, only: a4, norda, nordb
      use jaspar6, only: asymp_jasa, asymp_jasb, asymp_r, c1_jas6, c1_jas6i, c2_jas6
      use jaspar6, only: cutjas
      use bparm, only: nocuspb, nspin2b
      use contr2, only: ijas
      use contr2, only: isc
      use precision_kinds, only: dp
      implicit none

      integer :: i, iord, ipr, isp, it
      real(dp) :: val_cutjas
      real(dp), parameter :: third = 1.d0/3.d0










c isc = 2,3 are exponential scalings
c isc = 4,5 are inverse power scalings
c isc = 3,5 have zero 2nd and 3rd derivatives at 0
c isc = 6,7 are exponential and power law scalings resp. that have
c       zero, value and 1st 2 derivatives at cut_jas.
c isc = 12,14, are  are similar 2,4 resp. but are quadratic at 0 so as
c       to not contribute to the cusp without having to impose cusp conditions.
c isc = 16,17, iflag=1 are the same as 6,7 resp.
c       This is used to scale r_en in J_en and r_ee in J_ee for solids.
c isc = 16,17, iflag=2 are similar 6,7 resp. but are quadratic at 0 so as
c       to not contribute to the cusp without having to impose cusp conditions.
c       This is used to scale r_en in J_een for solids.
c isc = 16,17, iflag=3 have infinite range and are quadratic at 0 so as
c       to not contribute to the cusp without having to impose cusp conditions.
c       This is used to scale r_ee in J_een for solids.

c Evaluate constants that need to be reset if scalek is being varied.
c Warning: At present we are assuming that the same scalek is used
c for primary and secondary wavefns.  Otherwise c1_jas6i, c1_jas6, c2_jas6
c should be dimensioned to MWF
c Note val_cutjas is set isc=6,7 when isc=16,17 because isc=6,7 scalings are
c used for J_en and J_ee when isc=16,17.
      if(isc.eq.6 .or. isc.eq.16) then
        val_cutjas=exp(-third*scalek(1)*cutjas)
       elseif(isc.eq.7 .or. isc.eq.17) then
        val_cutjas=1/(1+third*scalek(1)*cutjas)
       else
        val_cutjas=0
        cutjas=1.d99
      endif
      c1_jas6i=1-val_cutjas
      c1_jas6=1/c1_jas6i
      c2_jas6=val_cutjas*c1_jas6
      if(ipr.gt.1) write(6,'(''cutjas,c1_jas6,c2_jas6='',9d12.5)')
     &cutjas,c1_jas6,c2_jas6

c Calculate asymptotic value of A and B terms
      asymp_r=c1_jas6i/scalek(1)
      do 10 it=1,nctype
        asymp_jasa(it)=a4(1,it,1)*asymp_r/(1+a4(2,it,1)*asymp_r)
        do 10 iord=2,norda
   10     asymp_jasa(it)=asymp_jasa(it)+a4(iord+1,it,1)*asymp_r**iord

      if(ijas.eq.4) then
        do 20 i=1,2
          if(i.eq.1) then
            sspinn=1
            isp=1
           else
            if(nspin2b.eq.1.and.nocuspb.eq.0) then
              sspinn=0.5d0
             else
              sspinn=1
            endif
            isp=nspin2b
          endif
          asymp_jasb(i)=sspinn*b(1,isp,1)*asymp_r/(1+b(2,isp,1)*asymp_r)
          do 20 iord=2,nordb
   20       asymp_jasb(i)=asymp_jasb(i)+b(iord+1,isp,1)*asymp_r**iord
       elseif(ijas.eq.5) then
        do 35 i=1,2
          if(i.eq.1) then
            sspinn=1
            isp=1
           else
            if(nspin2b.eq.1.and.nocuspb.eq.0) then
              sspinn=0.5d0
             else
              sspinn=1
            endif
            isp=nspin2b
          endif
          asymp_jasb(i)=b(1,isp,1)*asymp_r/(1+b(2,isp,1)*asymp_r)
          do 30 iord=2,nordb
   30       asymp_jasb(i)=asymp_jasb(i)+b(iord+1,isp,1)*asymp_r**iord
   35     asymp_jasb(i)=sspinn*asymp_jasb(i)
      endif
      if((ijas.eq.4.or.ijas.eq.5).and.ipr.gt.1) then
        write(6,'(''asymp_r='',f10.6)') asymp_r
        write(6,'(''asympa='',10f10.6)') (asymp_jasa(it),it=1,nctype)
        write(6,'(''asympb='',10f10.6)') (asymp_jasb(i),i=1,2)
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine scale_dist(r,rr,iflag)
c Written by Cyrus Umrigar
c Scale interparticle distances.

      use force_mod, only: MFORCE, MFORCE_WT_PRD, MWF
      use vmc_mod, only: MELEC, MORB, MBASIS, MDET, MCENT, MCTYPE, MCTYP3X
      use vmc_mod, only: NSPLIN, nrad, MORDJ, MORDJ1, MMAT_DIM, MMAT_DIM2, MMAT_DIM20
      use vmc_mod, only: radmax, delri
      use vmc_mod, only: NEQSX, MTERMS
      use vmc_mod, only: MCENT3, NCOEF, MEXCIT
      use jaspar3, only: a, c, scalek

      use jaspar6, only: asymp_r, c1_jas6, c2_jas6
      use jaspar6, only: cutjas, cutjasi
      use wfsec, only: iwf
      use contr2, only: ijas
      use contr2, only: isc
      use precision_kinds, only: dp
      implicit none

      integer :: iflag
      real(dp) :: deni, exprij, r, r_by_cut, r_by_cut2
      real(dp) :: rr, rsc, rsc2, rsc3
      real(dp) :: term
      real(dp), parameter :: zero = 0.d0
      real(dp), parameter :: one = 1.d0
      real(dp), parameter :: two = 2.d0
      real(dp), parameter :: half = 0.5d0
      real(dp), parameter :: third = 1.d0/3.d0
      real(dp), parameter :: d4b3 = 4.d0/3.d0








c isc = 2,3 are exponential scalings
c isc = 4,5 are inverse power scalings
c isc = 3,5 have zero 2nd and 3rd derivatives at 0
c isc = 6,7 are exponential and power law scalings resp. that have
c       zero, value and 1st 2 derivatives at cut_jas.
c isc = 12,14, are  are similar 2,4 resp. but are quadratic at 0 so as
c       to not contribute to the cusp without having to impose cusp conditions.
c isc = 16,17, iflag=1 are the same as 6,7 resp.
c       This is used to scale r_en in J_en and r_ee in J_ee for solids.
c isc = 16,17, iflag=2 are similar 6,7 resp. but are quadratic at 0 so as
c       to not contribute to the cusp without having to impose cusp conditions.
c       This is used to scale r_en in J_een for solids.
c isc = 16,17, iflag=3 have infinite range and are quadratic at 0 so as
c       to not contribute to the cusp without having to impose cusp conditions.
c       This is used to scale r_ee in J_een for solids.

      if(scalek(iwf).ne.zero) then
        if(isc.eq.2 .or. (isc.eq.12.and.iflag.eq.1)) then
          rr=(one-dexp(-scalek(iwf)*r))/scalek(iwf)
c     write(6,'(''r,rr='',9d12.4)') r,rr,scalek(iwf)
	 elseif(isc.eq.3) then
	  rsc=scalek(iwf)*r
          rsc2=rsc*rsc
          rsc3=rsc*rsc2
	  exprij=dexp(-rsc-half*rsc2-third*rsc3)
	  rr=(one-exprij)/scalek(iwf)
	 elseif(isc.eq.4 .or. (isc.eq.14.and.iflag.eq.1)) then
	  deni=one/(one+scalek(iwf)*r)
	  rr=r*deni
 	 elseif(isc.eq.5) then
	  deni=one/(one+(scalek(iwf)*r)**3)**third
	  rr=r*deni
         elseif(isc.eq.6 .or. (isc.eq.16.and.iflag.eq.1)) then
          if(r.gt.cutjas) then
            if(ijas.eq.4.or.ijas.eq.5) then
              rr=asymp_r
             elseif(ijas.eq.6) then
              rr=0
            endif
           else
            r_by_cut=r*cutjasi
            r_by_cut2=r_by_cut**2
            term=(1-r_by_cut+third*r_by_cut2)
            if(ijas.eq.4.or.ijas.eq.5) then
              rr=(1-dexp(-scalek(iwf)*r*term))/scalek(iwf)
             elseif(ijas.eq.6) then
              rr=c1_jas6*dexp(-scalek(iwf)*r*term)-c2_jas6
            endif
          endif
         elseif(isc.eq.7 .or. (isc.eq.17.and.iflag.eq.1)) then
          if(r.gt.cutjas) then
            if(ijas.eq.4.or.ijas.eq.5) then
              rr=asymp_r
             elseif(ijas.eq.6) then
              rr=0
            endif
           else
            r_by_cut=r*cutjasi
            r_by_cut2=r_by_cut**2
            term=(1-r_by_cut+third*r_by_cut2)
            deni=1/(1+scalek(iwf)*r*term)
            if(ijas.eq.4.or.ijas.eq.5) then
              rr=r*term*deni
             elseif(ijas.eq.6) then
              rr=c1_jas6*deni-c2_jas6
	    endif
	  endif
         elseif(isc.eq.16 .and. iflag.eq.2) then
          if(r.gt.cutjas) then
            if(ijas.eq.4.or.ijas.eq.5) then
              rr=asymp_r
             elseif(ijas.eq.6) then
              rr=0
            endif
           else
            r_by_cut=r*cutjasi
            r_by_cut2=r_by_cut**2
            term=(1-d4b3*r_by_cut+half*r_by_cut2)
            if(ijas.eq.4.or.ijas.eq.5) then
              rr=(1-dexp(-scalek(iwf)*r**2*term))/scalek(iwf)
             elseif(ijas.eq.6) then
              rr=c1_jas6*dexp(-scalek(iwf)*r**2*term)-c2_jas6
            endif
          endif
         elseif((isc.eq.12 .and. iflag.ge.2) .or. (isc.eq.16 .and. iflag.eq.3)) then
          if(ijas.eq.4.or.ijas.eq.5) then
            rr=(1-dexp(-scalek(iwf)*r**2))/scalek(iwf)
           elseif(ijas.eq.6) then
            rr=c1_jas6*dexp(-scalek(iwf)*r**2)-c2_jas6
          endif
         elseif(isc.eq.17 .and. iflag.eq.2) then
          if(r.gt.cutjas) then
            if(ijas.eq.4.or.ijas.eq.5) then
              rr=asymp_r
             elseif(ijas.eq.6) then
              rr=0
            endif
           else
            r_by_cut=r*cutjasi
            r_by_cut2=r_by_cut**2
            term=(1-d4b3*r_by_cut+half*r_by_cut2)
            deni=1/(1+scalek(iwf)*r**2*term)
            if(ijas.eq.4.or.ijas.eq.5) then
              rr=r**2*term*deni
             elseif(ijas.eq.6) then
              rr=c1_jas6*deni-c2_jas6
	    endif
	  endif
         elseif((isc.eq.14 .and. iflag.ge.2) .or. (isc.eq.17 .and. iflag.eq.3)) then
          deni=1/(1+scalek(iwf)*r**2)
          if(ijas.eq.4.or.ijas.eq.5) then
            rr=r**2*deni
           elseif(ijas.eq.6) then
            rr=c1_jas6*deni-c2_jas6
	  endif
	endif
       else
        rr=r
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine scale_dist1(r,rr,dd1,iflag)
c Written by Cyrus Umrigar
c Scale interparticle distances and calculate the 1st derivative
c of the scaled distances wrt the unscaled ones for calculating the
c gradient and laplacian.

      use force_mod, only: MFORCE, MFORCE_WT_PRD, MWF
      use vmc_mod, only: MELEC, MORB, MBASIS, MDET, MCENT, MCTYPE, MCTYP3X
      use vmc_mod, only: NSPLIN, nrad, MORDJ, MORDJ1, MMAT_DIM, MMAT_DIM2, MMAT_DIM20
      use vmc_mod, only: radmax, delri
      use vmc_mod, only: NEQSX, MTERMS
      use vmc_mod, only: MCENT3, NCOEF, MEXCIT
      use jaspar3, only: a, c, scalek

      use jaspar6, only: asymp_r, c1_jas6, c2_jas6
      use jaspar6, only: cutjas, cutjasi
      use wfsec, only: iwf
      use contr2, only: ijas
      use contr2, only: isc
      use precision_kinds, only: dp
      implicit none

      integer :: iflag
      real(dp) :: dd1, deni, exprij, r, r_by_cut
      real(dp) :: r_by_cut2, rr, rsc, rsc2
      real(dp) :: rsc3, term, term2
      real(dp), parameter :: zero = 0.d0
      real(dp), parameter :: one = 1.d0
      real(dp), parameter :: two = 2.d0
      real(dp), parameter :: half = 0.5d0
      real(dp), parameter :: third = 1.d0/3.d0
      real(dp), parameter :: d4b3 = 4.d0/3.d0








c isc = 2,3 are exponential scalings
c isc = 4,5 are inverse power scalings
c isc = 3,5 have zero 2nd and 3rd derivatives at 0
c isc = 6,7 are exponential and power law scalings resp. that have
c       zero, value and 1st 2 derivatives at cut_jas.
c isc = 12,14, are  are similar 2,4 resp. but are quadratic at 0 so as
c       to not contribute to the cusp without having to impose cusp conditions.
c isc = 16,17, iflag=1 are the same as 6,7 resp.
c       This is used to scale r_en in J_en and r_ee in J_ee for solids.
c isc = 16,17, iflag=2 are similar 6,7 resp. but are quadratic at 0 so as
c       to not contribute to the cusp without having to impose cusp conditions.
c       This is used to scale r_en in J_een for solids.
c isc = 16,17, iflag=3 have infinite range and are quadratic at 0 so as
c       to not contribute to the cusp without having to impose cusp conditions.
c       This is used to scale r_ee in J_een for solids.

      if(scalek(iwf).ne.zero) then
        if(isc.eq.2 .or. (isc.eq.12.and.iflag.eq.1)) then
          rr=(one-dexp(-scalek(iwf)*r))/scalek(iwf)
          dd1=one-scalek(iwf)*rr
c     write(6,'(''r,rr='',9d12.4)') r,rr,dd1,scalek(iwf)
	 elseif(isc.eq.3) then
	  rsc=scalek(iwf)*r
          rsc2=rsc*rsc
          rsc3=rsc*rsc2
	  exprij=dexp(-rsc-half*rsc2-third*rsc3)
	  rr=(one-exprij)/scalek(iwf)
	  dd1=(one+rsc+rsc2)*exprij
	 elseif(isc.eq.4 .or. (isc.eq.14.and.iflag.eq.1)) then
	  deni=one/(one+scalek(iwf)*r)
	  rr=r*deni
	  dd1=deni*deni
 	 elseif(isc.eq.5) then
	  deni=one/(one+(scalek(iwf)*r)**3)**third
	  rr=r*deni
	  dd1=deni**4
         elseif(isc.eq.6 .or. (isc.eq.16.and.iflag.eq.1)) then
          if(r.gt.cutjas) then
            if(ijas.eq.4.or.ijas.eq.5) then
              rr=asymp_r
             elseif(ijas.eq.6) then
              rr=0
            endif
            dd1=0
           else
            r_by_cut=r*cutjasi
            r_by_cut2=r_by_cut**2
            term=(1-r_by_cut+third*r_by_cut2)
            term2=(1-2*r_by_cut+r_by_cut2)
            if(ijas.eq.4.or.ijas.eq.5) then
              rr=(1-dexp(-scalek(iwf)*r*term))/scalek(iwf)
c             rr_plus_c2=rr+c2_jas6
              dd1=(one-scalek(iwf)*rr)*term2
             elseif(ijas.eq.6) then
              rr=c1_jas6*dexp(-scalek(iwf)*r*term)-c2_jas6
c             rr_plus_c2=rr+c2_jas6
c             dd1=-scalek(iwf)*rr_plus_c2*term2
            endif
          endif
         elseif(isc.eq.7 .or. (isc.eq.17.and.iflag.eq.1)) then
          if(r.gt.cutjas) then
            if(ijas.eq.4.or.ijas.eq.5) then
              rr=asymp_r
             elseif(ijas.eq.6) then
              rr=0
            endif
            dd1=0
           else
            r_by_cut=r*cutjasi
            r_by_cut2=r_by_cut**2
            term=(1-r_by_cut+third*r_by_cut2)
            term2=(1-2*r_by_cut+r_by_cut2)
            deni=1/(1+scalek(iwf)*r*term)
            if(ijas.eq.4.or.ijas.eq.5) then
              rr=r*term*deni
              dd1=term2*deni*deni
             elseif(ijas.eq.6) then
              rr=c1_jas6*deni-c2_jas6
c             rr_plus_c2=rr+c2_jas6
c             dd1=-scalek(iwf)*rr_plus_c2*deni*term2
	    endif
	  endif
         elseif(isc.eq.16 .and. iflag.eq.2) then
          if(r.gt.cutjas) then
            if(ijas.eq.4.or.ijas.eq.5) then
              rr=asymp_r
             elseif(ijas.eq.6) then
              rr=0
            endif
            dd1=0
           else
            r_by_cut=r*cutjasi
            r_by_cut2=r_by_cut**2
            term=(1-d4b3*r_by_cut+half*r_by_cut2)
            term2=(1-2*r_by_cut+r_by_cut2)
            if(ijas.eq.4.or.ijas.eq.5) then
              rr=(1-dexp(-scalek(iwf)*r**2*term))/scalek(iwf)
              dd1=2*r*(one-scalek(iwf)*rr)*term2
             elseif(ijas.eq.6) then
              rr=c1_jas6*dexp(-scalek(iwf)*r**2*term)-c2_jas6
c             rr_plus_c2=rr+c2_jas6
c             dd1=-scalek(iwf)*rr_plus_c2*term2
            endif
          endif
         elseif((isc.eq.12 .and. iflag.ge.2) .or. (isc.eq.16 .and. iflag.eq.3)) then
          if(ijas.eq.4.or.ijas.eq.5) then
            rr=(1-dexp(-scalek(iwf)*r**2))/scalek(iwf)
            dd1=2*r*(one-scalek(iwf)*rr)
           elseif(ijas.eq.6) then
            rr=c1_jas6*dexp(-scalek(iwf)*r**2)-c2_jas6
c           rr_plus_c2=rr+c2_jas6
c           dd1=-scalek(iwf)*rr_plus_c2*term2
          endif
         elseif(isc.eq.17 .and. iflag.eq.2) then
          if(r.gt.cutjas) then
            if(ijas.eq.4.or.ijas.eq.5) then
              rr=asymp_r
             elseif(ijas.eq.6) then
              rr=0
            endif
            dd1=0
           else
            r_by_cut=r*cutjasi
            r_by_cut2=r_by_cut**2
            term=(1-d4b3*r_by_cut+half*r_by_cut2)
            term2=(1-2*r_by_cut+r_by_cut2)
            deni=1/(1+scalek(iwf)*r**2*term)
            if(ijas.eq.4.or.ijas.eq.5) then
              rr=r**2*term*deni
              dd1=2*r*term2*deni*deni
             elseif(ijas.eq.6) then
              rr=c1_jas6*deni-c2_jas6
c             rr_plus_c2=rr+c2_jas6
c             dd1=-scalek(iwf)*rr_plus_c2*deni*term2
	    endif
	  endif
         elseif((isc.eq.14 .and. iflag.ge.2) .or. (isc.eq.17 .and. iflag.eq.3)) then
          deni=1/(1+scalek(iwf)*r**2)
          if(ijas.eq.4.or.ijas.eq.5) then
            rr=r**2*deni
            dd1=2*r*deni*deni
           elseif(ijas.eq.6) then
            rr=c1_jas6*deni-c2_jas6
c           rr_plus_c2=rr+c2_jas6
c           dd1=-scalek(iwf)*rr_plus_c2*deni*term2
	  endif
	endif
       else
        rr=r
        dd1=one
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine scale_dist2(r,rr,dd1,dd2,iflag)
c Written by Cyrus Umrigar
c Scale interparticle distances and calculate the 1st and 2nd derivs
c of the scaled distances wrt the unscaled ones for calculating the
c gradient and laplacian.

      use force_mod, only: MFORCE, MFORCE_WT_PRD, MWF
      use vmc_mod, only: MELEC, MORB, MBASIS, MDET, MCENT, MCTYPE, MCTYP3X
      use vmc_mod, only: NSPLIN, nrad, MORDJ, MORDJ1, MMAT_DIM, MMAT_DIM2, MMAT_DIM20
      use vmc_mod, only: radmax, delri
      use vmc_mod, only: NEQSX, MTERMS
      use vmc_mod, only: MCENT3, NCOEF, MEXCIT
      use jaspar3, only: a, c, scalek

      use jaspar6, only: asymp_r, c1_jas6, c2_jas6
      use jaspar6, only: cutjas, cutjasi
      use scale_more, only: dd3
      use wfsec, only: iwf
      use contr2, only: ijas
      use contr2, only: isc
      use precision_kinds, only: dp
      implicit none

      integer :: iflag
      real(dp) :: dd1, dd2, deni, exprij, r
      real(dp) :: r_by_cut, r_by_cut2, rr, rsc
      real(dp) :: rsc2, rsc3, term, term2
      real(dp), parameter :: zero = 0.d0
      real(dp), parameter :: one = 1.d0
      real(dp), parameter :: two = 2.d0
      real(dp), parameter :: half = 0.5d0
      real(dp), parameter :: third = 1.d0/3.d0
      real(dp), parameter :: d4b3 = 4.d0/3.d0










c isc = 2,3 are exponential scalings
c isc = 4,5 are inverse power scalings
c isc = 3,5 have zero 2nd and 3rd derivatives at 0
c isc = 6,7 are exponential and power law scalings resp. that have
c       zero, value and 1st 2 derivatives at cut_jas.
c isc = 12,14, are  are similar 2,4 resp. but are quadratic at 0 so as
c       to not contribute to the cusp without having to impose cusp conditions.
c isc = 16,17, iflag=1 are the same as 6,7 resp.
c       This is used to scale r_en in J_en and r_ee in J_ee for solids.
c isc = 16,17, iflag=2 are similar 6,7 resp. but are quadratic at 0 so as
c       to not contribute to the cusp without having to impose cusp conditions.
c       This is used to scale r_en in J_een for solids.
c isc = 16,17, iflag=3 have infinite range and are quadratic at 0 so as
c       to not contribute to the cusp without having to impose cusp conditions.
c       This is used to scale r_ee in J_een for solids.

      if(scalek(iwf).ne.zero) then
        if(isc.eq.2 .or. (isc.eq.12.and.iflag.eq.1)) then
          rr=(one-dexp(-scalek(iwf)*r))/scalek(iwf)
          dd1=one-scalek(iwf)*rr
          dd2=-scalek(iwf)*dd1
          dd3=-scalek(iwf)*dd2
c     write(6,'(''r,rr='',9d12.4)') r,rr,dd1,dd2,scalek(iwf)
	 elseif(isc.eq.3) then
	  rsc=scalek(iwf)*r
          rsc2=rsc*rsc
          rsc3=rsc*rsc2
	  exprij=dexp(-rsc-half*rsc2-third*rsc3)
	  rr=(one-exprij)/scalek(iwf)
	  dd1=(one+rsc+rsc2)*exprij
	  dd2=-scalek(iwf)*rsc2*(3+2*rsc+rsc2)*exprij
	 elseif(isc.eq.4 .or. (isc.eq.14.and.iflag.eq.1)) then
	  deni=one/(one+scalek(iwf)*r)
	  rr=r*deni
	  dd1=deni*deni
	  dd2=-two*scalek(iwf)*deni*dd1
 	 elseif(isc.eq.5) then
	  deni=one/(one+(scalek(iwf)*r)**3)**third
	  rr=r*deni
	  dd1=deni**4
	  dd2=-4*(scalek(iwf)*r)**2*scalek(iwf)*dd1*dd1/deni
         elseif(isc.eq.6 .or. (isc.eq.16.and.iflag.eq.1)) then
          if(r.gt.cutjas) then
            if(ijas.eq.4.or.ijas.eq.5) then
              rr=asymp_r
             elseif(ijas.eq.6) then
              rr=0
            endif
            dd1=0
            dd2=0
           else
            r_by_cut=r*cutjasi
            r_by_cut2=r_by_cut**2
            term=(1-r_by_cut+third*r_by_cut2)
            term2=(1-2*r_by_cut+r_by_cut2)
            if(ijas.eq.4.or.ijas.eq.5) then
              rr=(1-dexp(-scalek(iwf)*r*term))/scalek(iwf)
c             rr_plus_c2=rr+c2_jas6
              dd1=(one-scalek(iwf)*rr)*term2
              dd2=-scalek(iwf)*dd1*term2+2*(one-scalek(iwf)*rr)*cutjasi*(r_by_cut-1)
             elseif(ijas.eq.6) then
              rr=c1_jas6*dexp(-scalek(iwf)*r*term)-c2_jas6
c             rr_plus_c2=rr+c2_jas6
c             dd1=-scalek(iwf)*rr_plus_c2*term2
c             dd2=-scalek(iwf)*(dd1*term2+rr_plus_c2*2*cutjasi*(r_by_cut-1))
            endif
          endif
         elseif(isc.eq.7 .or. (isc.eq.17.and.iflag.eq.1)) then
          if(r.gt.cutjas) then
            if(ijas.eq.4.or.ijas.eq.5) then
              rr=asymp_r
             elseif(ijas.eq.6) then
              rr=0
            endif
            dd1=0
            dd2=0
           else
            r_by_cut=r*cutjasi
            r_by_cut2=r_by_cut**2
            term=(1-r_by_cut+third*r_by_cut2)
            term2=(1-2*r_by_cut+r_by_cut2)
            deni=1/(1+scalek(iwf)*r*term)
            if(ijas.eq.4.or.ijas.eq.5) then
              rr=r*term*deni
              dd1=term2*deni*deni
              dd2=-2*deni**2*(scalek(iwf)*deni*term2**2-cutjasi*(r_by_cut-1))
             elseif(ijas.eq.6) then
              rr=c1_jas6*deni-c2_jas6
c             rr_plus_c2=rr+c2_jas6
c             dd1=-scalek(iwf)*rr_plus_c2*deni*term2
c             dd2=-2*scalek(iwf)*deni*(dd1*term2+rr_plus_c2*cutjasi*(r_by_cut-1))
	    endif
	  endif
         elseif(isc.eq.16 .and. iflag.eq.2) then
          if(r.gt.cutjas) then
            if(ijas.eq.4.or.ijas.eq.5) then
              rr=asymp_r
             elseif(ijas.eq.6) then
              rr=0
            endif
            dd1=0
            dd2=0
           else
            r_by_cut=r*cutjasi
            r_by_cut2=r_by_cut**2
            term=(1-d4b3*r_by_cut+half*r_by_cut2)
            term2=(1-2*r_by_cut+r_by_cut2)
            if(ijas.eq.4.or.ijas.eq.5) then
              rr=(1-dexp(-scalek(iwf)*r**2*term))/scalek(iwf)
              dd1=2*r*(one-scalek(iwf)*rr)*term2
c             dd2=-scalek(iwf)*dd1*term2+2*(one-scalek(iwf)*rr)*cutjasi*(r_by_cut-1)
              dd2=-2*(scalek(iwf)*r*dd1*term2
     &        -(one-scalek(iwf)*rr)*(1-4*r_by_cut+3*r_by_cut2))
             elseif(ijas.eq.6) then
              rr=c1_jas6*dexp(-scalek(iwf)*r**2*term)-c2_jas6
c             rr_plus_c2=rr+c2_jas6
c             dd1=-scalek(iwf)*rr_plus_c2*term2
c             dd2=-scalek(iwf)*(dd1*term2+rr_plus_c2*2*cutjasi*(r_by_cut-1))
            endif
          endif
         elseif((isc.eq.12 .and. iflag.ge.2) .or. (isc.eq.16 .and. iflag.eq.3)) then
          if(ijas.eq.4.or.ijas.eq.5) then
            rr=(1-dexp(-scalek(iwf)*r**2))/scalek(iwf)
            dd1=2*r*(one-scalek(iwf)*rr)
            dd2=-2*(scalek(iwf)*r*dd1-(one-scalek(iwf)*rr))
           elseif(ijas.eq.6) then
            rr=c1_jas6*dexp(-scalek(iwf)*r**2)-c2_jas6
c           rr_plus_c2=rr+c2_jas6
c           dd1=-scalek(iwf)*rr_plus_c2*term2
c           dd2=-scalek(iwf)*(dd1*term2+rr_plus_c2*2*cutjasi*(r_by_cut-1))
          endif
         elseif(isc.eq.17 .and. iflag.eq.2) then
          if(r.gt.cutjas) then
            if(ijas.eq.4.or.ijas.eq.5) then
              rr=asymp_r
             elseif(ijas.eq.6) then
              rr=0
            endif
            dd1=0
            dd2=0
           else
            r_by_cut=r*cutjasi
            r_by_cut2=r_by_cut**2
            term=(1-d4b3*r_by_cut+half*r_by_cut2)
            term2=(1-2*r_by_cut+r_by_cut2)
            deni=1/(1+scalek(iwf)*r**2*term)
            if(ijas.eq.4.or.ijas.eq.5) then
              rr=r**2*term*deni
              dd1=2*r*term2*deni*deni
c             dd2=-2*deni**2*(scalek(iwf)*deni*term2**2-cutjasi*(r_by_cut-1))
              dd2=-2*deni**2*(4*scalek(iwf)*r**2*deni*term2**2
     &        -(1-4*r_by_cut+3*r_by_cut2))
             elseif(ijas.eq.6) then
              rr=c1_jas6*deni-c2_jas6
c             rr_plus_c2=rr+c2_jas6
c             dd1=-scalek(iwf)*rr_plus_c2*deni*term2
c             dd2=-2*scalek(iwf)*deni*(dd1*term2+rr_plus_c2*cutjasi*(r_by_cut-1))
	    endif
	  endif
         elseif((isc.eq.14 .and. iflag.ge.2) .or. (isc.eq.17 .and. iflag.eq.3)) then
          deni=1/(1+scalek(iwf)*r**2)
          if(ijas.eq.4.or.ijas.eq.5) then
            rr=r**2*deni
            dd1=2*r*deni*deni
            dd2=-2*deni**2*(4*scalek(iwf)*r**2*deni-1)
           elseif(ijas.eq.6) then
            rr=c1_jas6*deni-c2_jas6
c           rr_plus_c2=rr+c2_jas6
c           dd1=-scalek(iwf)*rr_plus_c2*deni*term2
c           dd2=-2*scalek(iwf)*deni*(dd1*term2+rr_plus_c2*cutjasi*(r_by_cut-1))
	  endif
	endif
       else
        rr=r
        dd1=one
        dd2=0
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine switch_scale(rr)
c Written by Cyrus Umrigar
c Switch scaling for ijas=4,5 from that appropriate for A,B terms to
c that appropriate for C terms, for dist.

      use force_mod, only: MFORCE, MFORCE_WT_PRD, MWF
      use vmc_mod, only: MELEC, MORB, MBASIS, MDET, MCENT, MCTYPE, MCTYP3X
      use vmc_mod, only: NSPLIN, nrad, MORDJ, MORDJ1, MMAT_DIM, MMAT_DIM2, MMAT_DIM20
      use vmc_mod, only: radmax, delri
      use vmc_mod, only: NEQSX, MTERMS
      use vmc_mod, only: MCENT3, NCOEF, MEXCIT
      use jaspar3, only: a, c, scalek

      use jaspar6, only: c1_jas6
      use wfsec, only: iwf
      use precision_kinds, only: dp
      implicit none


      real(dp) :: rr






      rr=1-c1_jas6*scalek(iwf)*rr

      return
      end
c-----------------------------------------------------------------------
      subroutine switch_scale1(rr,dd1)
c Written by Cyrus Umrigar
c Switch scaling for ijas=4,5 from that appropriate for A,B terms to
c that appropriate for C terms, for dist and 1st deriv.

      use force_mod, only: MFORCE, MFORCE_WT_PRD, MWF
      use vmc_mod, only: MELEC, MORB, MBASIS, MDET, MCENT, MCTYPE, MCTYP3X
      use vmc_mod, only: NSPLIN, nrad, MORDJ, MORDJ1, MMAT_DIM, MMAT_DIM2, MMAT_DIM20
      use vmc_mod, only: radmax, delri
      use vmc_mod, only: NEQSX, MTERMS
      use vmc_mod, only: MCENT3, NCOEF, MEXCIT
      use jaspar3, only: a, c, scalek

      use jaspar6, only: c1_jas6
      use wfsec, only: iwf
      use precision_kinds, only: dp
      implicit none


      real(dp) :: dd1, rr






      rr=1-c1_jas6*scalek(iwf)*rr
      dd1=-c1_jas6*scalek(iwf)*dd1

      return
      end
c-----------------------------------------------------------------------
      subroutine switch_scale2(rr,dd1,dd2)
c Written by Cyrus Umrigar
c Switch scaling for ijas=4,5 from that appropriate for A,B terms to
c that appropriate for C terms, for dist and 1st two derivs.

      use force_mod, only: MFORCE, MFORCE_WT_PRD, MWF
      use vmc_mod, only: MELEC, MORB, MBASIS, MDET, MCENT, MCTYPE, MCTYP3X
      use vmc_mod, only: NSPLIN, nrad, MORDJ, MORDJ1, MMAT_DIM, MMAT_DIM2, MMAT_DIM20
      use vmc_mod, only: radmax, delri
      use vmc_mod, only: NEQSX, MTERMS
      use vmc_mod, only: MCENT3, NCOEF, MEXCIT
      use jaspar3, only: a, c, scalek

      use jaspar6, only: c1_jas6
      use wfsec, only: iwf
      use precision_kinds, only: dp
      implicit none


      real(dp) :: dd1, dd2, rr






      rr=1-c1_jas6*scalek(iwf)*rr
      dd1=-c1_jas6*scalek(iwf)*dd1
      dd2=-c1_jas6*scalek(iwf)*dd2

      return
      end
c-----------------------------------------------------------------------
