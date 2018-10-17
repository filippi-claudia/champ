c  bcspline -- dmc 30 May 1996

c  set up coefficients for bicubic spline with following BC's:
c  FULL BC CONTROL at all bdys

c  algorithm note -- handling of y & z explicit BC's while maintaining
c  full C2 differentiability is delicate.  Basic method:  use a fully
c  C2 method based on the "not-a-knot" BC, and then, correct to meet
c  each user BC by calculating a C2 spline that is zero at all grid
c  points but satisfies a BC which is the difference btw the user spec
c  and the not-a-knot result; add the coeffs of this into the original.

c  for this more workspace is needed: nwk .ge. 4*inx*inth +5*max(inx,inth)


      subroutine r8bcspline(x,inx,th,inth,fspl,inf3,
     >                    ibcxmin,bcxmin,ibcxmax,bcxmax,
     >                    ibcthmin,bcthmin,ibcthmax,bcthmax,
     >                    wk,nwk,ilinx,ilinth,ier)

!============
! idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
c     INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER inth,inf3,ibcxmin,ibcxmax,ibcthmin,ibcthmax,nwk,ilinx
      INTEGER ilinth,ier,inx,iflg2,ix,itest,ierx,ierth,inxo,ith
      INTEGER intho,ic,ibcthmina,ibcthmaxa,iasc,iinc,iawk,jx,jth,ii
      INTEGER iadr,ia5w,iaspl
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 xo2,xo6,zhxn,zhth,zcur,zdiff1,zdiff2
!============
      REAL*8 x(inx),th(inth),fspl(4,4,inf3,inth),wk(nwk)
      REAL*8 bcxmin(inth),bcxmax(inth)
      REAL*8 bcthmin(inx),bcthmax(inx)

c  input:
c    x(1...inx) -- abscissae, first dimension of data
c   th(1...inth) -- abscissae, second (periodic) dimension of data
c   fspl(1,1,1..inx,1..inth) -- function values
c   inf3 -- fspl dimensioning, inf3.ge.inx required.

c  boundary conditions input:
c   ibcxmin -- indicator for boundary condition at x(1):
c    bcxmin(...) -- boundary condition data
c     =-1 -- periodic boundary condition
c     =0 -- use "not a knot", bcxmin(...) ignored
c     =1 -- match slope, specified at x(1),th(ith) by bcxmin(ith)
c     =2 -- match 2nd derivative, specified at x(1),th(ith) by bcxmin(ith)
c     =3 -- boundary condition is slope=0 (df/dx=0) at x(1), all th(j)
c     =4 -- boundary condition is d2f/dx2=0 at x(1), all th(j)
c     =5 -- match 1st derivative to 1st divided difference
c     =6 -- match 2nd derivative to 2nd divided difference
c     =7 -- match 3rd derivative to 3rd divided difference
c           (for more detailed definition of BCs 5-7, see the
c           comments of subroutine mkspline)
c   NOTE bcxmin(...) referenced ONLY if ibcxmin=1 or ibcxmin=2

c   ibcxmax -- indicator for boundary condition at x(nx):
c    bcxmax(...) -- boundary condition data
c     (interpretation as with ibcxmin, bcxmin)
c   NOTE:  if ibcxmin=-1, ibcxmax is ignored! ...and the BC is periodic.

c   ibcthmin -- indicator for boundary condition at th(1):
c    bcthmin(...) -- boundary condition data
c     (interpretation as with ibcxmin, bcxmin)
c   ibcthmax -- indicator for boundary condition at th(inth):
c    bcthmax(...) -- boundary condition data
c     (interpretation as with ibcxmin, bcxmin)
c   NOTE:  if ibcthmin=-1, ibcthmax is ignored! ...and the BC is periodic.

c   NOTE the bcxmin,bcxmax,bcthmin,bcthmax arrays are only used if the
c     corresponding boundary condition flags are set to 1 or 2.
c     Carefully note the dimensioning of these arrays!

c  output:
c   fspl(*,*,1..inx,1..inth) -- bicubic spline coeffs (4x4)
c   ...fspl(1,1,*,*) is not replaced.

c   ilinx -- =1 on output if x(inx) pts are nearly evenly spaced (tol=1e-3)
c   ilinth-- =1 on output if th(inth) evenly spaced (tol=1e-3)

c   ier -- completion code, 0 for normal

c  workspace:
c   wk -- must be at least 5*max(inx,inth) large
c                          5*max(inx,inth) + 4*inx*inth large
c                          if explicit non-zero d/dth or d2/dth2 BC's
c                          are supplied.
c  nwk -- size of workspace of workspace provided

c---------------------------------
c  in what follows, "f" is an abbreviation for "fspl"...

c  compute bicubic spline of 2d function, given values at the grid
c  grid crossing points, f(1,1,i,j)=f(x(i),th(j)).

c  on evaluation:  for point x btw x(i) and x(i+1), dx=x-x(i)
c                       and th btw th(j) and th(j+1), dt=th-th(j),

c      spline =
c        f(1,1,i,j) + dx*f(2,1,i,j) + dx**2*f(3,1,i,j) + dx**3*f(4,1,i,j)
c   +dt*(f(1,2,i,j) + dx*f(2,2,i,j) + dx**2*f(3,2,i,j) + dx**3*f(4,2,i,j))
c   +d2*(f(1,3,i,j) + dx*f(2,3,i,j) + dx**2*f(3,3,i,j) + dx**3*f(4,3,i,j))
c   +d3*(f(1,4,i,j) + dx*f(2,4,i,j) + dx**2*f(3,4,i,j) + dx**3*f(4,4,i,j))

c      where d2=dt**2 and d3=dt**3.

      integer iselect1(10)
      integer iselect2(10)

c---------------------------------

c  see if 2nd pass is needed due to "non-linear" d/dth bdy cond.

      iflg2=0
      if(ibcthmin.ne.-1) then
         if((ibcthmin.eq.1).or.(ibcthmin.eq.2)) then
            do ix=1,inx
               if (bcthmin(ix).ne.0.0d0) iflg2=1
            enddo
         endif
         if((ibcthmax.eq.1).or.(ibcthmax.eq.2)) then
            do ix=1,inx
               if (bcthmax(ix).ne.0.0d0) iflg2=1
            enddo
         endif
      endif

      ier=0
      itest=5*max(inx,inth)
      if(iflg2.eq.1) then
         itest=itest +4*inx*inth
      endif

      if(nwk.lt.itest) then
         write(6,9901) nwk,itest
 9901    format(' ?bcspline:  workspace too small:'/
     >          '  user supplied:  nwk=',i6,'; need at least:  ',i6/
     >          '  nwk=4*nx*ny +5*max(nx,ny) will work for any user'/
     >          '  choice of bdy conditions.')
         ier=1
      endif
      if(inx.lt.4) then
         write(6,'('' ?bcspline:  at least 4 x points required.'')')
         ier=1
      endif
      if(inth.lt.4) then
         write(6,'('' ?bcspline:  need at least 4 theta points.'')')
         ier=1
      endif

      call ibc_ck(ibcxmin,'bcspline','xmin',-1,7,ier)
      if(ibcxmin.ge.0) call ibc_ck(ibcxmax,'bcspline','xmax',0,7,ier)
      call ibc_ck(ibcthmin,'bcspline','thmin',-1,7,ier)
      if(ibcthmin.ge.0) call ibc_ck(ibcthmax,'bcspline','thmax',0,7,ier)

c  check ilinx & x vector

      call r8splinck(x,inx,ilinx,1.0D-3,ierx)
      if(ierx.ne.0) ier=2

      if(ier.eq.2) then
         write(6,'('' ?bcspline:  x axis not strict ascending'')')
      endif

c  check ilinth & th vector

      call r8splinck(th,inth,ilinth,1.0D-3,ierth)
      if(ierth.ne.0) ier=3

      if(ier.eq.3) then
         write(6,'('' ?bcspline:  th axis not strict ascending'')')
      endif

      if(ier.ne.0) return

c------------------------------------

      xo2=0.5d0
      xo6=1.0d0/6.0d0

c  spline in x direction

      inxo=4*(inx-1)
      do ith=1,inth

c  copy the function in

         do ix=1,inx
            wk(4*(ix-1)+1)=fspl(1,1,ix,ith)
         enddo

         if(ibcxmin.eq.1) then
            wk(2)=bcxmin(ith)
         else if(ibcxmin.eq.2) then
            wk(3)=bcxmin(ith)
         endif

         if(ibcxmax.eq.1) then
            wk(inxo+2)=bcxmax(ith)
         else if(ibcxmax.eq.2) then
            wk(inxo+3)=bcxmax(ith)
         endif

c  use Wayne's routine

         call r8v_spline(ibcxmin,ibcxmax,inx,x,wk,wk(4*inx+1))

c  copy the coefficients out

         do ix=1,inx-1
            fspl(2,1,ix,ith)=wk(4*(ix-1)+2)
            fspl(3,1,ix,ith)=wk(4*(ix-1)+3)*xo2
            fspl(4,1,ix,ith)=wk(4*(ix-1)+4)*xo6
         enddo

      enddo

c-----------------------------------

c  spline in theta direction

      intho=4*(inth-1)
      do ix=1,inx-1

c  spline each x coeff

         do ic=1,4

c  copy ordinates in

            do ith=1,inth
               wk(4*(ith-1)+1)=fspl(ic,1,ix,ith)
            enddo

c  first pass:  use a linear BC -- if flag indicates BC correction
c  will be needed, it will be done later

            wk(2)=0.0d0
            wk(3)=0.0d0
            wk(intho+2)=0.0d0
            wk(intho+3)=0.0d0

            ibcthmina=ibcthmin
            ibcthmaxa=ibcthmax
            if(iflg2.eq.1) then
               if((ibcthmin.eq.1).or.(ibcthmin.eq.2)) ibcthmina=0
               if((ibcthmax.eq.1).or.(ibcthmax.eq.2)) ibcthmaxa=0
            endif

            call r8v_spline(ibcthmina,ibcthmaxa,inth,th,wk,wk(4*inth+1))

c  copy coeffs out

            do ith=1,inth-1
               fspl(ic,2,ix,ith)=wk(4*(ith-1)+2)
               fspl(ic,3,ix,ith)=wk(4*(ith-1)+3)*xo2
               fspl(ic,4,ix,ith)=wk(4*(ith-1)+4)*xo6
            enddo

         enddo

      enddo

c  now make correction for user BC's if needed

      if(iflg2.eq.1) then

         iasc=1                         ! wk addr for correction splines
         iinc=4*inth                    ! spacing btw correction splines
         iawk=iasc+4*inth*inx

c  last grid zone widths

         zhxn=x(inx)-x(inx-1)
         jx=inx-1
         zhth=th(inth)-th(inth-1)
         jth=inth-1

         do ii=1,10
            iselect1(ii)=0
            iselect2(ii)=0
         enddo
         if(ibcthmin.eq.1) iselect1(3)=1
         if(ibcthmin.eq.2) iselect1(5)=1
         if(ibcthmax.eq.1) iselect2(3)=1
         if(ibcthmax.eq.2) iselect2(5)=1

c  loop over BC's

         do ix=1,inx

c  (a) d/dth @ th(1) difference btw current BC and user request

            if(ibcthmin.eq.1) then
               if(ix.lt.inx) then
                  zcur=fspl(1,2,ix,1)   ! 1st deriv.
               else
                  zcur=fspl(1,2,jx,1)+zhxn*(fspl(2,2,jx,1)+zhxn*
     >               (fspl(3,2,jx,1)+zhxn*fspl(4,2,jx,1)))
               endif
               zdiff1=bcthmin(ix)-zcur
            else if(ibcthmin.eq.2) then
               if(ix.lt.inx) then
                  zcur=2.0d0*fspl(1,3,ix,1) ! 2nd deriv.
               else
                  zcur=2.0d0*(fspl(1,3,jx,1)+zhxn*(fspl(2,3,jx,1)+zhxn*
     >               (fspl(3,3,jx,1)+zhxn*fspl(4,3,jx,1))))
               endif
               zdiff1=bcthmin(ix)-zcur
            else
               zdiff1=0.0d0
            endif

c  (b) d/dth @ th(inth) difference btw current BC and user request

            if(ibcthmax.eq.1) then
               if(ix.lt.inx) then
c  1st deriv.
                  zcur=fspl(1,2,ix,jth)+zhth*(2.0d0*fspl(1,3,ix,jth)+
     >               zhth*3.0d0*fspl(1,4,ix,jth))
               else
                  call r8bcspeval(x(inx),th(inth),iselect2,  zcur,
     >               x,inx,th,inth,ilinx,ilinth,fspl,inf3,ier)
                  if(ier.ne.0) return
               endif
               zdiff2=bcthmax(ix)-zcur
            else if(ibcthmax.eq.2) then
               if(ix.lt.inx) then
c  2nd deriv.
                  zcur=2.0d0*fspl(1,3,ix,jth)+
     >               6.0d0*zhth*fspl(1,4,ix,jth)
               else
                  call r8bcspeval(x(inx),th(inth),iselect2,  zcur,
     >               x,inx,th,inth,ilinx,ilinth,fspl,inf3,ier)
                  if(ier.ne.0) return
               endif
               zdiff2=bcthmax(ix)-zcur
            else
               zdiff2=0.0d0
            endif

c  ok compute the theta spline with BC's to span the difference(s)
c  these theta "correction splines" are zero at all the grid points
c  but have at least one non-zero 1st or 2nd derivative BC

            iadr=iasc+(ix-1)*iinc
            do ith=1,inth
               wk(iadr+4*(ith-1))=0.0d0
            enddo

            wk(iadr+1)=0.0d0
            wk(iadr+2)=0.0d0
            wk(iadr+intho+1)=0.0d0
            wk(iadr+intho+2)=0.0d0

            if(ibcthmin.eq.1) then
               wk(iadr+1)=zdiff1
            else if(ibcthmin.eq.2) then
               wk(iadr+2)=zdiff1
            endif

            if(ibcthmax.eq.1) then
               wk(iadr+intho+1)=zdiff2
            else if(ibcthmax.eq.2) then
               wk(iadr+intho+2)=zdiff2
            endif

            call r8v_spline(ibcthmin,ibcthmax,inth,th,wk(iadr),wk(iawk))
         enddo

c  add in results to main array -- th spline coef corrections

         do ix=1,inx
            iadr=iasc+(ix-1)*iinc
            do ith=1,inth-1
               wk(iadr+4*(ith-1)+2)=wk(iadr+4*(ith-1)+2)*xo2
               wk(iadr+4*(ith-1)+3)=wk(iadr+4*(ith-1)+3)*xo6
               if(ix.lt.inx) then
                  fspl(1,2,ix,ith)=fspl(1,2,ix,ith)+wk(iadr+4*(ith-1)+1)
                  fspl(1,3,ix,ith)=fspl(1,3,ix,ith)+wk(iadr+4*(ith-1)+2)
                  fspl(1,4,ix,ith)=fspl(1,4,ix,ith)+wk(iadr+4*(ith-1)+3)
               endif
            enddo
         enddo

c  compute the x splines of the th spline correction coeffs

         ia5w=iawk+4*inx

         do ith=1,inth-1
            do ic=2,4
               do ix=1,inx
                  iaspl=iasc+iinc*(ix-1)
                  wk(iawk+4*(ix-1))=wk(iaspl+4*(ith-1)+(ic-1))
               enddo

c  use zero BCs for this correction spline

               wk(iawk+1)=0.0d0
               wk(iawk+2)=0.0d0
               wk(iawk+inxo+1)=0.0d0
               wk(iawk+inxo+2)=0.0d0

c  periodic spline of correction spline higher coeffs (1st coeffs are
c  all zero by defn of the correction spline

               call r8v_spline(ibcxmin,ibcxmax,inx,x,wk(iawk),wk(ia5w))

               do ix=1,inx-1
                  fspl(2,ic,ix,ith)=fspl(2,ic,ix,ith)+
     >               wk(iawk+4*(ix-1)+1)
                  fspl(3,ic,ix,ith)=fspl(3,ic,ix,ith)+
     >               wk(iawk+4*(ix-1)+2)*xo2
                  fspl(4,ic,ix,ith)=fspl(4,ic,ix,ith)+
     >               wk(iawk+4*(ix-1)+3)*xo6
               enddo

            enddo
         enddo                          ! ith

      endif                             ! BC correction needs test

      return
      end
