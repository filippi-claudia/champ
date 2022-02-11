      subroutine splfit(r,irb,ic,iwf,wfv,ider)
c Written by Claudia Filippi
c get spline_fit at r of basis fn irb of center ic and force iwf
c 1st and 2nd derivs also calculated if ider=1.

      use vmc_mod, only: NCOEF
      use numbas, only: arg, d2rwf, igrid, nr, r0, rwf

      use numexp, only: ae, ce
      use precision_kinds, only: dp
      implicit none

      integer :: ic, ider, ii, irb, iwf
      integer :: jx
      real(dp) :: aa, bb, cc, dd, delh
      real(dp) :: delhi, dlogag, dlrad0, r
      real(dp) :: ref0, ref1, xr
      real(dp), dimension(4) :: wfv
      real(dp), parameter :: sixth = 1.d0/6.d0

      if(igrid(ic).eq.1)then
        xr=(r-r0(ic))/arg(ic)+1
        jx=int(xr)
        ref0=r0(ic)+arg(ic)*(jx-1)
        ref1=ref0+arg(ic)
        delh=arg(ic)
       elseif(igrid(ic).eq.2)then
        dlrad0=dlog(r0(ic))
        dlogag=dlog(arg(ic))
        xr=(dlog(r)-dlrad0)/dlogag+1
        jx=int(xr)
        ref0=r0(ic)*arg(ic)**(jx-1)
        ref1=ref0*arg(ic)
        delh=ref1-ref0
       elseif(igrid(ic).eq.3)then
        dlogag=dlog(arg(ic))
        xr=(dlog((r+r0(ic))/r0(ic)))/dlogag+1.d0
        jx=int(xr)
        ref0=r0(ic)*arg(ic)**(jx-1)-r0(ic)
        ref1=(ref0+r0(ic))*arg(ic)-r0(ic)
        delh=ref1-ref0
      endif

      if(jx.lt.1) then
c use small radius polynomial

        wfv(1)=ce(1,irb,ic,iwf)
        do ii=2,NCOEF
          wfv(1)=wfv(1)+ce(ii,irb,ic,iwf)*r**(ii-1)
        enddo
        if(ider.ge.1) then
          wfv(2)=0.d0
          wfv(3)=0.d0
          wfv(4)=0.d0
          do ii=2,NCOEF
            wfv(2)=wfv(2)+(ii-1)*ce(ii,irb,ic,iwf)*r**(ii-2)
            wfv(3)=wfv(3)+(ii-1)*(ii-2)*ce(ii,irb,ic,iwf)*r**(ii-3)
            wfv(4)=wfv(4)+(ii-1)*(ii-2)*(ii-3)*ce(ii,irb,ic,iwf)*r**(ii-4)
          enddo
        endif

       elseif (jx.ge.nr(ic)) then
c use large radius exponential

        wfv(1)=ae(1,irb,ic,iwf)*dexp(-ae(2,irb,ic,iwf)*r)
        if(ider.ge.1) then
          wfv(2)=-ae(2,irb,ic,iwf)*wfv(1)
          wfv(3)=-ae(2,irb,ic,iwf)*wfv(2)
          wfv(4)=-ae(2,irb,ic,iwf)*wfv(3)
        endif

       else
c cubic spline interpolation

        delhi=1.d0/delh
        bb=(r-ref0)*delhi
        aa=(ref1-r)*delhi
        cc=aa*(aa**2-1.d0)*delh**2*sixth
        dd=bb*(bb**2-1.d0)*delh**2*sixth
        wfv(1)=aa*rwf(jx,irb,ic,iwf)+bb*rwf(jx+1,irb,ic,iwf)+
     &         cc*d2rwf(jx,irb,ic,iwf)+dd*d2rwf(jx+1,irb,ic,iwf)

        if(ider.ge.1) then
          wfv(2)=(rwf(jx+1,irb,ic,iwf)-rwf(jx,irb,ic,iwf))*delhi+
     &    (-(3.d0*aa**2-1.d0)*d2rwf(jx,irb,ic,iwf)
     &     +(3.d0*bb**2-1.d0)*d2rwf(jx+1,irb,ic,iwf))*delh*sixth
          wfv(3)=aa*d2rwf(jx,irb,ic,iwf)+bb*d2rwf(jx+1,irb,ic,iwf)
          wfv(4)=(-d2rwf(jx,irb,ic,iwf)+d2rwf(jx+1,irb,ic,iwf))*delhi
        endif
      endif

      return
      end
