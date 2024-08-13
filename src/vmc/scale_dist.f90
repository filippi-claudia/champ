module scale_dist_mod
contains
      subroutine set_scale_dist(iadiag,ipr)
! Written by Cyrus Umrigar
      use bparm,   only: nocuspb,nspin2b
      use contrl_file, only: ounit
      use jastrow, only: a4,asymp_jasa,asymp_jasb,asymp_r,b,ijas,isc,norda,nordb
      use jastrow, only: scalek,sspinn
      use multiple_geo, only: nforce,nwftype
      use precision_kinds, only: dp
      use system,  only: nctype
      use vmc_mod, only: nwftypejas

      implicit none

      integer :: i, iadiag, i0, i1, j, iord, ipr, isp, it, nwf

! isc = 2,3 are exponential scalings
! isc = 4,5 are inverse power scalings
! isc = 3,5 have zero 2nd and 3rd derivatives at 0

! Evaluate constants that need to be reset if scalek is being varied

! Warning: we are assuming that same scalek is used for primary and secondary wavefns

! Calculate asymptotic value of A and B terms
      asymp_r=1.d0/scalek(1)

      j=iadiag
      do it=1,nctype
        asymp_jasa(it,j)=a4(1,it,j)*asymp_r/(1+a4(2,it,j)*asymp_r)
        do iord=2,norda
          asymp_jasa(it,j)=asymp_jasa(it,j)+a4(iord+1,it,j)*asymp_r**iord
        enddo
      enddo

      do i=1,2
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
        asymp_jasb(i,j)=sspinn*b(1,isp,j)*asymp_r/(1+b(2,isp,j)*asymp_r)
        do iord=2,nordb
          asymp_jasb(i,j)=asymp_jasb(i,j)+b(iord+1,isp,j)*asymp_r**iord
        enddo
      enddo

      if(ipr.gt.1) then
        write(ounit,'(''Jastrow type='',i4)') j
        write(ounit,'(''asymp_r='',f10.6)') asymp_r
        write(ounit,'(''asympa='',10f10.6)') (asymp_jasa(it,1),it=1,nctype)
        write(ounit,'(''asympb='',10f10.6)') (asymp_jasb(i,1),i=1,2)
      endif

      return
      end
!-----------------------------------------------------------------------
      subroutine scale_dist(r,rr)
! Written by Cyrus Umrigar
! Scale interparticle distances.

      use jastrow, only: ijas,isc,scalek
      use multiple_geo, only: iwf
      use precision_kinds, only: dp
      implicit none

      real(dp) :: deni, exprij, r
      real(dp) :: rr, rsc, rsc2, rsc3
      real(dp), parameter :: zero = 0.d0
      real(dp), parameter :: one = 1.d0
      real(dp), parameter :: two = 2.d0
      real(dp), parameter :: half = 0.5d0
      real(dp), parameter :: third = 1.d0/3.d0

! isc = 2,3 are exponential scalings
! isc = 4,5 are inverse power scalings
! isc = 3,5 have zero 2nd and 3rd derivatives at 0

      if(scalek(iwf).ne.zero) then
        if(isc.eq.2) then
          rr=(one-dexp(-scalek(iwf)*r))/scalek(iwf)
         elseif(isc.eq.3) then
          rsc=scalek(iwf)*r
          rsc2=rsc*rsc
          rsc3=rsc*rsc2
          exprij=dexp(-rsc-half*rsc2-third*rsc3)
          rr=(one-exprij)/scalek(iwf)
         elseif(isc.eq.4) then
          deni=one/(one+scalek(iwf)*r)
          rr=r*deni
         elseif(isc.eq.5) then
          deni=one/(one+(scalek(iwf)*r)**3)**third
          rr=r*deni
        endif
       else
        rr=r
      endif
!     write(ounit,'(''r,rr='',9d12.4)') r,rr,scalek(iwf)

      return
      end
!-----------------------------------------------------------------------
      subroutine scale_dist1(r,rr,dd1)
! Written by Cyrus Umrigar
! Scale interparticle distances and calculate the 1st derivative
! of the scaled distances wrt the unscaled ones for calculating the
! gradient and laplacian.

      use jastrow, only: asymp_r,ijas,isc,scalek
      use multiple_geo, only: iwf
      use precision_kinds, only: dp
      implicit none

      real(dp) :: dd1, deni, exprij, r
      real(dp) :: rr, rsc, rsc2,  rsc3
      real(dp), parameter :: zero = 0.d0
      real(dp), parameter :: one = 1.d0
      real(dp), parameter :: two = 2.d0
      real(dp), parameter :: half = 0.5d0
      real(dp), parameter :: third = 1.d0/3.d0

! isc = 2,3 are exponential scalings
! isc = 4,5 are inverse power scalings
! isc = 3,5 have zero 2nd and 3rd derivatives at 0

      if(scalek(iwf).ne.zero) then
        if(isc.eq.2) then
          rr=(one-dexp(-scalek(iwf)*r))/scalek(iwf)
          dd1=one-scalek(iwf)*rr
         elseif(isc.eq.3) then
          rsc=scalek(iwf)*r
          rsc2=rsc*rsc
          rsc3=rsc*rsc2
          exprij=dexp(-rsc-half*rsc2-third*rsc3)
          rr=(one-exprij)/scalek(iwf)
          dd1=(one+rsc+rsc2)*exprij
         elseif(isc.eq.4) then
          deni=one/(one+scalek(iwf)*r)
          rr=r*deni
          dd1=deni*deni
         elseif(isc.eq.5) then
          deni=one/(one+(scalek(iwf)*r)**3)**third
          rr=r*deni
          dd1=deni**4
        endif
       else
        rr=r
        dd1=one
      endif
!     write(ounit,'(''r,rr='',9d12.4)') r,rr,dd1,scalek(iwf)

      return
      end
!-----------------------------------------------------------------------
      subroutine scale_dist2(r,rr,dd1,dd2)
! Written by Cyrus Umrigar
! Scale interparticle distances and calculate the 1st and 2nd derivs
! of the scaled distances wrt the unscaled ones for calculating the
! gradient and laplacian.

      use jastrow, only: asymp_r,ijas,isc,scalek
      use multiple_geo, only: iwf
      use precision_kinds, only: dp
      use scale_more, only: dd3
      implicit none

      real(dp) :: dd1, dd2, deni, exprij, r
      real(dp) :: rr, rsc, rsc2, rsc3, term, term2
      real(dp), parameter :: zero = 0.d0
      real(dp), parameter :: one = 1.d0
      real(dp), parameter :: two = 2.d0
      real(dp), parameter :: half = 0.5d0
      real(dp), parameter :: third = 1.d0/3.d0

! isc = 2,3 are exponential scalings
! isc = 4,5 are inverse power scalings
! isc = 3,5 have zero 2nd and 3rd derivatives at 0

      if(scalek(iwf).ne.zero) then
        if(isc.eq.2) then
          rr=(one-dexp(-scalek(iwf)*r))/scalek(iwf)
          dd1=one-scalek(iwf)*rr
          dd2=-scalek(iwf)*dd1
          dd3=-scalek(iwf)*dd2
         elseif(isc.eq.3) then
          rsc=scalek(iwf)*r
          rsc2=rsc*rsc
          rsc3=rsc*rsc2
          exprij=dexp(-rsc-half*rsc2-third*rsc3)
          rr=(one-exprij)/scalek(iwf)
          dd1=(one+rsc+rsc2)*exprij
          dd2=-scalek(iwf)*rsc2*(3+2*rsc+rsc2)*exprij
         elseif(isc.eq.4) then
          deni=one/(one+scalek(iwf)*r)
          rr=r*deni
          dd1=deni*deni
          dd2=-two*scalek(iwf)*deni*dd1
         elseif(isc.eq.5) then
          deni=one/(one+(scalek(iwf)*r)**3)**third
          rr=r*deni
          dd1=deni**4
          dd2=-4*(scalek(iwf)*r)**2*scalek(iwf)*dd1*dd1/deni
        endif
       else
        rr=r
        dd1=one
        dd2=0
      endif
!     write(ounit,'(''r,rr='',9d12.4)') r,rr,dd1,dd2,scalek(iwf)

      return
      end
!-----------------------------------------------------------------------
      subroutine switch_scale(rr)
! Written by Cyrus Umrigar
! Switch scaling for ijas=4 from that appropriate for A,B terms to
! that appropriate for C terms, for dist.

      use jastrow, only: scalek
      use multiple_geo, only: iwf
      use precision_kinds, only: dp
      implicit none

      real(dp) :: rr

      rr=1-scalek(iwf)*rr

      return
      end
!-----------------------------------------------------------------------
      subroutine switch_scale1(rr,dd1)
! Written by Cyrus Umrigar
! Switch scaling for ijas=4 from that appropriate for A,B terms to
! that appropriate for C terms, for dist and 1st deriv.

      use jastrow, only: scalek
      use multiple_geo, only: iwf
      use precision_kinds, only: dp
      implicit none

      real(dp) :: dd1, rr

      rr=1-scalek(iwf)*rr
      dd1=-scalek(iwf)*dd1

      return
      end
!-----------------------------------------------------------------------
      subroutine switch_scale2(rr,dd1,dd2)
! Written by Cyrus Umrigar
! Switch scaling for ijas=4 from that appropriate for A,B terms to
! that appropriate for C terms, for dist and 1st two derivs.

      use jastrow, only: scalek
      use multiple_geo, only: iwf
      use precision_kinds, only: dp
      implicit none

      real(dp) :: dd1, dd2, rr

      rr=1-scalek(iwf)*rr
      dd1=-scalek(iwf)*dd1
      dd2=-scalek(iwf)*dd2

      return
      end
!-----------------------------------------------------------------------
end module 
