      function grater(fn,xmin,xmax,abserr,relerr,nsing,xsing,error
     1,numcal,mregions)
      implicit real*8(a-h,o-z)
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c This is Cyrus Umrigar's rewrite of Mike Teter's integration routine.  
c This is a general purpose integration routine usable whenever the
c function can be evaluated at arbitrary points.
c If the non-analyticites of the function are known, then it helps to
c specify their locations.
c It is not necessary to list non-analyticities at end-points.
c It uses a 3-pt Gauss-Legendre formula and a 9-pt formula.
c The 9-pt formula integrates polynomials upto order 9 correctly while
c The 3-pt formula integrates polynomials upto order 5 correctly.
c For comparison Simpson is good to order 3 and Bode to order 5.
c Is it possible to do better, while still using the same points for
c the 9-pt and 3-pt formulae?  We could integrate polynomials upto
c order 18 if we used the 9-pt Gauss-Legendre formula.  But then we could
c not reuse points and the 3-pt formula would be only good to order 3.

c fn       = function to be integrated (declare external in MAIN)
c xmin     = lower limit
c xmax     = upper limit
c abserr   = absolute tolerable error
c relerr   = relative tolerable error
c nsing    = number of singularities
c xsing    = array of locations of non-analyticities
c error    = estimated error (conservative)
c numcal   = the number of times fn was called
c mregions = the maximum number of regions occuring simultaneously
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      external fn
      parameter (MBOUND=1001,MREGION=MBOUND-1,zero=0.d0,eps=1.d-15)
      dimension x(MBOUND),fval(3,MREGION),dx(3),wt3(3),wt9(9)
      dimension xsing(*)
c     save dx,wt3,wt9
      data dx/0.1127016653792583d0,0.5d0,0.8872983346207417d0/
      data wt3/0.2777777777777777778d0,0.4444444444444444444d0,
     1         0.2777777777777777778d0/
      data wt9/0.0616938806304841571d0,0.108384229110206161d0,
     1         0.0398463603260281088d0,0.175209035316976464d0,
     2         0.229732989232610220d0, 0.175209035316976464d0,
     3         0.0398463603260281088d0,0.108384229110206161d0,
     4         0.0616938806304841571d0/

      if(abserr.le.zero .and. relerr.lt.eps) stop 'abserr and relerr =0'

      safety=.1d0
      isafety=0
      numcal_sav=999999999
      error_sav=9.d99
    5 error=0.d0
      grater=0.d0

c The array x stores the boundary points of the regions.
c Place initial boundary points at end-pts of region and at those
c singular pts that are within the region.
c nregion is the number of different intervals into which the 
c integration region is currently divided. The number of regions can
c grow if more accuracy is needed by dividing the right-most region
c into three regions. The number of regions shrinks when the integral
c over the right-most region is accurate enough, in which case that
c integral is added to the total (stored in grater) and the region
c is removed from consideration (and a new region is the right-most).
      x(1)=xmin
      nsing_in=0
      do 10 j=1,nsing
        if((xsing(j)-xmin)*(xmax-xsing(j)).gt.zero) then
          nsing_in=nsing_in+1
          x(nsing_in+1)=xsing(j)
        endif
   10 continue
      x(nsing_in+2)=xmax
      nregion=nsing_in+1
      mregions=nregion

c For each region, calculate the function, store at three selected points
c and calculate initial estimate of integral.
      estim=0.d0
      do 20 i=1,nregion
        del3=x(i+1)-x(i)
        do 20 j=1,3
          fval(j,i)=fn(x(i)+del3*dx(j))
   20     estim=estim+fval(j,i)*wt3(j)*del3
      absestim=abs(estim)
      numcal=3*nregion

   30 continue
        if(nregion+3.ge.MBOUND) then
          write(*,*) 'TOO MANY REGIONS IN GRATER'
          write(6,'(9d15.8)') frac,abserr,relerr
          write(6,'(9d15.8)') value3,value9,value3-value9
          do 35 i=1,nregion
   35       write(6,'(f15.10,9d15.8)') x(i),(fval(j,i),j=1,3)
          stop
        endif

c Divide the rightmost region into three subregions.  
        del9=x(nregion+1)-x(nregion)
        x(nregion+3)=x(nregion+1)
        x(nregion+1)=x(nregion)+del9*dx(1)*2.
        x(nregion+2)=x(nregion+3)-del9*dx(1)*2.

c The three data points already found for the region become the 
c middle data points (number 2 in first index of fval) for each region.
        fval(2,nregion+2)=fval(3,nregion)
        fval(2,nregion+1)=fval(2,nregion)
        fval(2,nregion)=fval(1,nregion)

c Now do the integral over the right-most region in two different ways-
c a 3-point integral (value3) over each of the 3 subregions 
c and a more accurate 9-point integral (value9) over the whole region.
c value3 is used only for the error estimate.
        icount=0
        value9=0.d0
        value3=0.d0
        do 40 j=nregion,nregion+2
          del3=x(j+1)-x(j)
          fval(1,j)=fn(x(j)+dx(1)*del3)
          fval(3,j)=fn(x(j)+dx(3)*del3)
          do 40 k=1,3
            icount=icount+1
            value9=value9+wt9(icount)*fval(k,j)*del9
   40       value3=value3+fval(k,j)*wt3(k)*del3
        numcal=numcal+6
        dif=abs(value9-value3)

c If one of the following condition is true, add in this integral to the
c total, and reduce the number of regions under consideration.
c The 3rd condition can be dangerous if the integrand oscillates and
c yields a small integral, but without it sometimes takes for ever
c calculating over tiny intervals.  Because of this dangerous condition
c we check at the end and loop back to the beginning if accuracy criteria
c are not satisfied.
        frac=del9/(xmax-xmin)
        if(dif.le.abserr*frac .or. dif.le.relerr*abs(value9)*safety
     1  .or.  dif.le.relerr*absestim*safety) then
c       if(dif.le.abserr*frac .or. dif.le.relerr*abs(value9) .or.
c    1  (frac .le. 1.d-14)) then
c The following commented out lines are Steve White's error criterion.
c       if(dif .le. abserr*frac .or. dif.le.relerr*abs(value9) .or. 
c    1   (frac .le. 1.d-15) .or.
c    2   (frac .le. 1.d-8 .and. dif .le. abserr*0.1d0)) then
          grater=grater+value9
          error=error+abs(dif)
          nregion=nregion-1
c If no more regions, we are done.
          if(nregion.le.0) then
            if(error.le.5*relerr*dabs(grater).or. error.le.5*abserr)
     &      then
              if(isafety.ge.1)
     &        write(6,'(''num,val,error'',i8,2d21.14,d10.2)'
     &        ) numcal,grater,absestim,error
              return
             else
              write(6,'(''num,val,error'',i8,2d21.14,d10.2)'
     &        ) numcal,grater,absestim,error
              write(6,'(''relerr,abserr='',2d12.6)') relerr,abserr
              safety=.5*safety
              if(safety.lt.1.d-3) stop 'safety2 < 1.d-3'
              if(isafety.ge.1 .and. numcal*error.gt.numcal_sav*error_sav
     &          ) then
                write(6,*) 'Warning: round-off error in grater'
                return
              endif
              isafety=isafety+1
              numcal_sav=numcal
              error_sav=error
              goto 5
            endif
          endif
         else
c If the integration is insufficiently accurate, make each of the 
c three subregions of the right-most region into regions.
c On next pass the right-most of these is the new current region.
          nregion=nregion+2
          mregions = max(mregions,nregion)
        endif
      goto 30
      end
