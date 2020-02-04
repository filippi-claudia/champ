c  tcspline -- dmc 20 Jan 1999

c  set up coefficients for bicubic spline with following BC's:
c  * LHS and RHS handled as in cubspl.for for 1st coordinate
c  * derivatives periodic in second coordinate (use pspline.for)

c workspace:
c  if phi bdy cond. is periodic, not-a-knot, df/dphi = 0 everywhere,
c  or d2f/dphi2 = 0 everywhere, then the phi boundary condition is
c  "linear" and a workspace of size at least:

c     nwk = 20*inx*inth + 10*max(inx,inth,inph)

c  will suffice.

c  if the phi bdy cond. involves specification of df/dphi .ne. 0 or
c  d2f/dphi .ne. 0 at any (x,theta) grid point, then, the phi boundary
c  condition is "non-linear", a correction step is needed, and a workspace
c  of size at least:

c     nwk = 16*inx*inth*inph

c  is required.

      subroutine r8tcspline(x,inx,th,inth,ph,inph,fspl,inf4,inf5,
     >                    ibcxmin,bcxmin,ibcxmax,bcxmax,inb1x,
     >                    ibcthmin,bcthmin,ibcthmax,bcthmax,inb1th,
     >                    ibcphmin,bcphmin,ibcphmax,bcphmax,inb1ph,
     >                    wk,nwk,ilinx,ilinth,ilinph,ier)

!============
! idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
c     INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER inth,inph,inf4,inf5,ibcxmin,ibcxmax,inb1x,ibcthmin
      INTEGER ibcthmax,inb1th,ibcphmin,ibcphmax,inb1ph,nwk,ilinx
      INTEGER ilinth,ilinph,ier,inx,iflg,ith,ix,itest,ierx,ierth
      INTEGER ierph,iaspl2,iabcx1,iabcx2,iabcth1,iabcth2,iawk,inwk
      INTEGER iph,inpho,ic1,ic2,ibcphmina,ibcphmaxa,iabcph1,iabcph2
      INTEGER iaccoef,i,iskip1,iskip2,ia1,ia2
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 xo2,xo6,zcur
!============
      REAL*8 x(inx),th(inth),ph(inph)
      REAL*8 fspl(4,4,4,inf4,inf5,inph),wk(nwk)
      REAL*8 bcxmin(inb1x,inph),bcxmax(inb1x,inph) ! inth x inph defined
      REAL*8 bcthmin(inb1th,inph),bcthmax(inb1th,inph) ! inx x inph defined
      REAL*8 bcphmin(inb1ph,inth),bcphmax(inb1ph,inth) ! inx x inth defined

c  input:
c    x(1...inx) -- abscissae, first dimension of data
c   th(1...inth) -- abscissae, second (periodic) dimension of data
c   ph(1...inph) -- abscissae, third (periodic) dimension of data
c   fspl(1,1,1,1..inx,1..inth,1..inph) -- function values
c   inf4 -- fspl dimensioning, inf4.ge.inx required.
c   inf5 -- fspl dimensioning, inf5.ge.inth required.

c  boundary conditions input:

c   bc data at xmin, xmax  vs.  theta,phi
c   bc data at thmin, thmax  vs.  x,phi
c   bc data at phmin, phmax  vs.  x,theta

c   ibcxmin -- indicator for boundary condition at x(1):
c    bcxmin(...) -- boundary condition data
c     =-1 -- use periodic boundary condition
c     =0 -- use "not a knot", bcxmin(...) ignored
c     =1 -- match slope, specified at x(1),th(ith),ph(iph) by bcxmin(ith,iph)
c     =2 -- match 2nd derivative, specified at x(1),th(ith),ph(iph)
c           by bcxmin(ith,iph
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
c     NOTE:  if ibcxmin=-1 then the periodic BC applies on both sides
c            and ibcxmax is ignored.
c   inb1x -- 1st dimension of bcxmin, bcxmax: if ibcxmin or ibcxmax .gt. 0
c            this must be .ge. inth:

c   interpretation of ibcthmin,bcthmin,ibcthmax,bcthmax,inb1th
c     is same as with ibcxmin,...

c   interpretation of ibcphmin,bcphmin,ibcphmax,bcphmax,inb1ph
c     is same as with ibcxmin,...

c   the explicit bdy condition arrays are referenced only if the
c     corresponding "ibc" flag values are set to 1 or 2.

c  output:
c   fspl(*,*,*,1..inx,1..inth,1..inph) -- bicubic spline coeffs (4x4)
c   ...fspl(1,1,1,*,*,*) is not replaced.

c   ilinx -- =1 on output if x(inx) pts are nearly evenly spaced (tol=1e-3)
c   ilinth-- =1 on output if th(inth) evenly spaced (tol=1e-3)
c   ilinph-- =1 on output if ph(inph) evenly spaced (tol=1e-3)

c   ier -- completion code, 0 for normal

c  workspace:
c   wk -- must be at least 5*max(inx,inth,inph) large -- or more, see
c         comments, above.
c  nwk -- size of workspace

c---------------------------------
c  ** in what follows, f is an abbreviation for fspl **

c  compute tricubic spline of 3d function, given values at the
c  grid crossing points, f(1,1,1,i,j,k)=f(x(i),th(j),ph(k)).

c  on evaluation:  for point x btw x(i) and x(i+1), dx=x-x(i)
c                       and th btw th(j) and th(j+1), dt=th-th(j),
c                       and ph btw ph(k) and ph(k+1), dp=ph-ph(k),

c      spline =
c        f(1,1,1,i,j,k)+dx*f(2,1,1,i,j,k)+dx2*f(3,1,1,i,j,k)+dx3*f(4,1,1,i,j,k)
c  +dt*(f(1,2,1,i,j,k)+dx*f(2,2,1,i,j,k)+dx2*f(3,2,1,i,j,k)+dx3*f(4,2,1,i,j,k))
c +dt2*(f(1,3,1,i,j,k)+dx*f(2,3,1,i,j,k)+dx2*f(3,3,1,i,j,k)+dx3*f(4,3,1,i,j,k))
c +dt3*(f(1,4,1,i,j,k)+dx*f(2,4,1,i,j,k)+dx2*f(3,4,1,i,j,k)+dx3*f(4,4,1,i,j,k))
c        +dp*(
c        f(1,1,2,i,j,k)+dx*f(2,1,2,i,j,k)+dx2*f(3,1,2,i,j,k)+dx3*f(4,1,2,i,j,k)
c  +dt*(f(1,2,2,i,j,k)+dx*f(2,2,2,i,j,k)+dx2*f(3,2,2,i,j,k)+dx3*f(4,2,2,i,j,k))
c +dt2*(f(1,3,2,i,j,k)+dx*f(2,3,2,i,j,k)+dx2*f(3,3,2,i,j,k)+dx3*f(4,3,2,i,j,k))
c +dt3*(f(1,4,2,i,j,k)+dx*f(2,4,2,i,j,k)+dx2*f(3,4,2,i,j,k)+dx3*f(4,4,2,i,j,k)))
c        +dp2*(
c        f(1,1,3,i,j,k)+dx*f(2,1,3,i,j,k)+dx2*f(3,1,3,i,j,k)+dx3*f(4,1,3,i,j,k)
c  +dt*(f(1,2,3,i,j,k)+dx*f(2,2,3,i,j,k)+dx2*f(3,2,3,i,j,k)+dx3*f(4,2,3,i,j,k))
c +dt2*(f(1,3,3,i,j,k)+dx*f(2,3,3,i,j,k)+dx2*f(3,3,3,i,j,k)+dx3*f(4,3,3,i,j,k))
c +dt3*(f(1,4,3,i,j,k)+dx*f(2,4,3,i,j,k)+dx2*f(3,4,3,i,j,k)+dx3*f(4,4,3,i,j,k)))
c        +dp3*(
c        f(1,1,4,i,j,k)+dx*f(2,1,4,i,j,k)+dx2*f(3,1,4,i,j,k)+dx3*f(4,1,4,i,j,k)
c  +dt*(f(1,2,4,i,j,k)+dx*f(2,2,4,i,j,k)+dx2*f(3,2,4,i,j,k)+dx3*f(4,2,4,i,j,k))
c +dt2*(f(1,3,4,i,j,k)+dx*f(2,3,4,i,j,k)+dx2*f(3,3,4,i,j,k)+dx3*f(4,3,4,i,j,k))
c +dt3*(f(1,4,4,i,j,k)+dx*f(2,4,4,i,j,k)+dx2*f(3,4,4,i,j,k)+dx3*f(4,4,4,i,j,k)))

c      where dx2=dx**2 and dx3=dx**3.
c      where dt2=dt**2 and dt3=dt**3.
c      where dp2=dp**2 and dp3=dp**3.

c---------------------------------
      integer iselect1(10)
      integer iselect2(10)

      REAL*8 z0,z1,ztol

      data z0/0.0D0/
      data z1/1.0D0/
      data ztol/1.0D-3/

c---------------------------------

      ier=0

      iflg=0

c  check phi bdy condition "linearity"

      if(ibcphmin.ne.-1) then
         if((ibcphmin.eq.1).or.(ibcphmin.eq.2)) then
            do ith=1,inth
               do ix=1,inx
                  if(bcphmin(ix,ith).ne.z0) iflg=1
               enddo
            enddo
         endif
         if((ibcphmax.eq.1).or.(ibcphmax.eq.2)) then
            do ith=1,inth
               do ix=1,inx
                  if(bcphmax(ix,ith).ne.z0) iflg=1
               enddo
            enddo
         endif
      endif

      itest=10*max(inx,inth,inph) + 20*inx*inth
      if(iflg.eq.1) then
         itest=16*inx*inth*inph
      endif

      if(nwk.lt.itest) then
         write(6,9901) nwk,itest
         ier=1
 9901    format(' ?tcspline:  workspace too small.'/
     >      '  user supplied nwk=',i7,'; need at least: ',i7/
     >      '  If no explicit df/dph boundary condition is set,'/
     >      '  nwk = 20*inx*inth + 10*max(inx,inth,inph) can be used.'/
     >      '  If an explicit df/dph or d2f/dph2 boundary condition'/
     >      '  is set, nwk=16*inx*inth*inph is required.')
      endif
      if(inx.lt.4) then
         write(6,'('' ?tcspline:  at least 4 x points required.'')')
         ier=1
      endif
      if(inth.lt.4) then
         write(6,'('' ?tcspline:  need at least 4 theta points.'')')
         ier=1
      endif
      if(inph.lt.4) then
         write(6,'('' ?tcspline:  need at least 4 phi points.'')')
         ier=1
      endif

      if((ibcxmin.eq.1).or.(ibcxmax.eq.1).or.(ibcxmin.eq.2).or.
     >   (ibcxmax.eq.2)) then
         if(inb1x.lt.inth) then
            ier=1
            write(6,
     >'('' ?tcspline:  1st dim of bcxmin/max arrays .lt. inth'')')
         endif
      endif

      if((ibcthmin.eq.1).or.(ibcthmax.eq.1).or.(ibcthmin.eq.2).or.
     >   (ibcthmax.eq.2)) then
         if(inb1th.lt.inx) then
            ier=1
            write(6,
     >'('' ?tcspline:  1st dim of bcthmin/max arrays .lt. inx'')')
         endif
      endif

      if((ibcphmin.eq.1).or.(ibcphmax.eq.1).or.(ibcphmin.eq.2).or.
     >   (ibcphmax.eq.2)) then
         if(inb1ph.lt.inx) then
            ier=1
            write(6,
     >'('' ?tcspline:  1st dim of bphmin/max arrays .lt. inx'')')
         endif
      endif

      call ibc_ck(ibcxmin,'tcspline','xmin',-1,7,ier)
      if(ibcxmin.ge.0) call ibc_ck(ibcxmax,'tcspline','xmax',0,7,ier)

      call ibc_ck(ibcthmin,'tcspline','thmin',-1,7,ier)
      if(ibcthmin.ge.0) call ibc_ck(ibcthmax,'tcspline','thmax',0,7,ier)

      call ibc_ck(ibcphmin,'tcspline','phmin',-1,7,ier)
      if(ibcphmax.ge.0) call ibc_ck(ibcphmax,'tcspline','phmax',0,7,ier)

c  check ilinx & x vector

      call r8splinck(x,inx,ilinx,ztol,ierx)
      if(ierx.ne.0) ier=2

      if(ier.eq.2) then
         write(6,'('' ?tcspline:  x axis not strict ascending'')')
      endif

c  check ilinth & th vector

      call r8splinck(th,inth,ilinth,ztol,ierth)
      if(ierth.ne.0) ier=3

      if(ier.eq.3) then
         write(6,'('' ?tcspline:  theta axis not strict ascending'')')
      endif

c  check ilinth & th vector

      call r8splinck(ph,inph,ilinph,ztol,ierph)
      if(ierph.ne.0) ier=4

      if(ier.eq.4) then
         write(6,'('' ?tcspline:  phi axis not strict ascending'')')
      endif

      if(ier.ne.0) return

c------------------------------------

c  part 1.  compute (x,theta) spline coeffs via an intermediate
c  routine that call bcspline

c  workspace addresses

      iaspl2=1
      iabcx1=iaspl2+16*inx*inth
      iabcx2=iabcx1+inth
      iabcth1=iabcx2+inth
      iabcth2=iabcth1+inx
      iawk=iabcth2+inx
      inwk=nwk-iawk+1

      do iph=1,inph

c  copy bc data

         do ix=1,inx
            wk(iabcth1+ix-1)=0.0d0
            wk(iabcth2+ix-1)=0.0d0
            if((ibcthmin.eq.1).or.(ibcthmin.eq.2)) then
               wk(iabcth1+ix-1)=bcthmin(ix,iph)
            endif
            if((ibcthmin.ne.-1).and.
     >         ((ibcthmax.eq.1).or.(ibcthmax.eq.2))) then
               wk(iabcth2+ix-1)=bcthmax(ix,iph)
            endif
         enddo
         do ith=1,inth
            wk(iabcx1+ith-1)=0.0d0
            wk(iabcx2+ith-1)=0.0d0
            if((ibcxmin.eq.1).or.(ibcxmin.eq.2)) then
               wk(iabcx1+ith-1)=bcxmin(ith,iph)
            endif
            if((ibcxmin.ne.-1).and.
     >         ((ibcxmax.eq.1).or.(ibcxmax.eq.2))) then
               wk(iabcx2+ith-1)=bcxmax(ith,iph)
            endif
         enddo

c  call 2d spline intermediary routine

         call r8tcsp23(x,inx,th,inth,fspl(1,1,1,1,1,iph),inf4,
     >      ibcxmin,wk(iabcx1),ibcxmax,wk(iabcx2),
     >      ibcthmin,wk(iabcth1),ibcthmax,wk(iabcth2),
     >      wk(iaspl2),wk(iawk),inwk,ilinx,ilinth,ier)

         if(ier.ne.0) then
            write(6,*) ' ?tcspline:  error in 2d spline, exiting.'
            return
         endif

      enddo

c  ok now fspl(*,*,1,*,*,*) have been evaluated and C2 in (x,theta)
c  now need to extend to coeffs in phi direction.

      xo2=0.5D0
      xo6=1.0D0/6.0D0

c  spline each (x,th) coeff in the phi direction

      inpho=4*(inph-1)
      do ith=1,inth-1
         do ix=1,inx-1

            do ic1=1,4
               do ic2=1,4

c  copy coeff. ordinates in

                  do iph=1,inph
                     wk(4*(iph-1)+1)=fspl(ic1,ic2,1,ix,ith,iph)
                  enddo

c  use linear BC on this first pass; will correct later if
c  necessary

                  wk(2)=0.0d0
                  wk(3)=0.0d0
                  wk(inpho+2)=0.0d0
                  wk(inpho+3)=0.0d0

                  ibcphmina=ibcphmin
                  ibcphmaxa=ibcphmax
                  if(iflg.eq.1) then
                     if((ibcphmin.eq.1).or.(ibcphmin.eq.2)) ibcphmina=0
                     if((ibcphmax.eq.1).or.(ibcphmax.eq.2)) ibcphmaxa=0
                  endif

                  call r8v_spline(ibcphmina,ibcphmaxa,inph,ph,wk,
     >               wk(4*inph+1))

c  copy coeffs out

                  do iph=1,inph-1
                     fspl(ic1,ic2,2,ix,ith,iph)=wk(4*(iph-1)+2)
                     fspl(ic1,ic2,3,ix,ith,iph)=wk(4*(iph-1)+3)*xo2
                     fspl(ic1,ic2,4,ix,ith,iph)=wk(4*(iph-1)+4)*xo6
                  enddo

               enddo                    ! ic2
            enddo                       ! ic1

         enddo                          ! ix
      enddo                             ! ith

c  if there are "non-linear" BCs requiring correction...

c  at each (x(ix),th(ith)) get the d/dph BC's right while preserving C2
c  everywhere...

      if(iflg.eq.1) then

c  first get BC correction numbers

         iabcph1=1
         iabcph2=iabcph1+inx*inth
         iaccoef=iabcph2+inx*inth
         iawk=iaccoef+12*inx*inth*inph
         inwk=nwk-iawk+1

         do i=1,10
            iselect1(i)=0
            iselect2(i)=0
         enddo

c  note because iflg=1, we know at least one of ibcphmin/max = 1 or 2

         iskip1=0
         if(ibcphmin.eq.1) then
            iselect1(4)=1               ! df/dph
         else if(ibcphmin.eq.2) then
            iselect1(7)=1               ! d2f/dph2
         else
            iskip1=1
         endif

         iskip2=0
         if(ibcphmax.eq.1) then
            iselect2(4)=1               ! df/dph
         else if(ibcphmax.eq.2) then
            iselect2(7)=1               ! d2f/dph2
         else
            iskip2=1
         endif

         ia1=iabcph1-1
         ia2=iabcph2-1
         do ith=1,inth
            do ix=1,inx
               ia1=ia1+1
               ia2=ia2+1

               if(iskip1.eq.0) then
                  call r8tcspeval(x(ix),th(ith),ph(1),iselect1, zcur,
     >               x,inx,th,inth,ph,inph,ilinx,ilinth,ilinph,
     >               fspl,inf4,inf5,ier)
                  if(ier.ne.0) then
                     write(6,*) ' ?? tcspline:  error in tcspeval call'
                     return
                  endif
                  wk(ia1)=bcphmin(ix,ith)-zcur ! correction needed
               else
                  wk(ia1)=z0
               endif

               if(iskip2.eq.0) then
                  call r8tcspeval(x(ix),th(ith),ph(inph),iselect2, zcur,
     >               x,inx,th,inth,ph,inph,ilinx,ilinth,ilinph,
     >               fspl,inf4,inf5,ier)
                  if(ier.ne.0) then
                     write(6,*) ' ?? tcspline:  error in tcspeval call'
                     return
                  endif
                  wk(ia2)=bcphmax(ix,ith)-zcur ! correction needed
               else
                  wk(ia2)=z0
               endif
            enddo
         enddo

         call r8tcspcorr(x,inx,th,inth,ph,inph,fspl,inf4,inf5,
     >      ibcxmin,ibcxmax,ibcthmin,ibcthmax,
     >      ibcphmin,wk(iabcph1),ibcphmax,wk(iabcph2),
     >      wk(iaccoef),wk(iawk),inwk)

      endif

      return
      end
c-----------------------------------------
      subroutine r8tcspcorr(x,inx,th,inth,ph,inph,fspl,inf4,inf5,
     >   ibcxmin,ibcxmax,ibcthmin,ibcthmax,
     >   ibcphmin,bcph1,ibcphmax,bcph2,ccorr,wk,nwk)

c  intermediary routine for tcspline:
c  do correction needed to get C2 3d spline with phi bdy conditions
c  matched.

c  all input unless noted:

!============
! idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
c     INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER inth,inph,inf4,inf5,nwk,inx,iawk2,inpho,ith,ix,iph
      INTEGER inxo,icph,intho,icx,icth
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 xo2,xo6,z0,zfac
!============
      REAL*8 x(inx)                       ! x axis
      REAL*8 th(inth)                     ! th axis
      REAL*8 ph(inph)                     ! ph axis

      REAL*8 fspl(4,4,4,inf4,inf5,inph)   ! spline coeffs -- adjusted

      integer ibcxmin,ibcxmax           ! x BC flags
      integer ibcthmin,ibcthmax         ! th BC flags
      integer ibcphmin,ibcphmax         ! ph BC flags

      REAL*8 bcph1(inx,inth)              ! ph BC correction array @ phi(1)
      REAL*8 bcph2(inx,inth)              ! ph BC correction array @ phi(inph)

c  workspaces:

      REAL*8 ccorr(3,4,inx,inth,inph)     ! correction coefficients (partial)

      REAL*8 wk(nwk)                      ! workspace

c---------------------------

      xo2=0.5D0
      xo6=1.0D0/6.0D0

c  1.  splines in phi -- fcns are zero everywhere but have non-zero BCs

      if(nwk.lt.10*max(inx,inth,inph)) then
         write(6,*) ' ?? programming error in tcspcorr (tcspline)'
         return
      endif

      z0=0.0d0
      iawk2=4*inph+1

      inpho=4*(inph-1)
      do ith=1,inth
         do ix=1,inx

            do iph=1,inph
               wk(4*(iph-1)+1)=z0
            enddo

c  set BC for this 1d spline correction

            if(ibcphmin.eq.1) then
               wk(2)=bcph1(ix,ith)
            else if(ibcphmin.eq.2) then
               wk(3)=bcph1(ix,ith)
            endif

            if(ibcphmax.eq.1) then
               wk(inpho+2)=bcph2(ix,ith)
            else if(ibcphmax.eq.2) then
               wk(inpho+3)=bcph2(ix,ith)
            endif

            call r8v_spline(ibcphmin,ibcphmax,inph,ph,wk,wk(iawk2))

c  copy non-zero coeffs out to ccorr

            do iph=1,inph-1
               ccorr(1,1,ix,ith,iph)=wk(4*(iph-1)+2)
               ccorr(2,1,ix,ith,iph)=wk(4*(iph-1)+3)*xo2
               ccorr(3,1,ix,ith,iph)=wk(4*(iph-1)+4)*xo6
            enddo

         enddo
      enddo

c  2. spline the coeffs in x -- use ibcx flags & zero for derivative
c  bc if necessary

      iawk2=4*inx+1

      inxo=4*(inx-1)
      do iph=1,inph-1
         do ith=1,inth

            do icph=1,3

               do ix=1,inx
                  wk(4*(ix-1)+1)=ccorr(icph,1,ix,ith,iph)
               enddo

c  zero BC:  correction spline

               wk(2)=0.0d0
               wk(3)=0.0d0
               wk(inxo+2)=0.0d0
               wk(inxo+3)=0.0d0

               call r8v_spline(ibcxmin,ibcxmax,inx,x,wk,wk(iawk2))

               do ix=1,inx-1
                  ccorr(icph,2,ix,ith,iph)=wk(4*(ix-1)+2)
                  ccorr(icph,3,ix,ith,iph)=wk(4*(ix-1)+3)*xo2
                  ccorr(icph,4,ix,ith,iph)=wk(4*(ix-1)+4)*xo6
               enddo

            enddo

         enddo
      enddo

c  3.  spline all the ccorr coefs in th -- use ibcth flags & zero for
c      derivative correction BC if necessary

c      add the results into fspl

      iawk2=4*inth+1

      intho=4*(inth-1)
      do iph=1,inph-1
         do ix=1,inx-1

            do icx=1,4
               do icph=1,3

                  do ith=1,inth
                     wk(4*(ith-1)+1)=ccorr(icph,icx,ix,ith,iph)
                  enddo

c  zero BC:  correction spline

                  wk(2)=0.0d0
                  wk(3)=0.0d0
                  wk(intho+2)=0.0d0
                  wk(intho+3)=0.0d0

                  call r8v_spline(ibcthmin,ibcthmax,inth,th,wk,
     >               wk(iawk2))

                  do ith=1,inth-1
                     do icth=1,4
                        zfac=1.0d0
                        if(icth.eq.3) zfac=xo2
                        if(icth.eq.4) zfac=xo6
                        fspl(icx,icth,icph+1,ix,ith,iph)=
     >                     fspl(icx,icth,icph+1,ix,ith,iph)+
     >                     wk(4*(ith-1)+icth)*zfac
                     enddo
                  enddo

               enddo
            enddo

         enddo
      enddo

      return
      end
c-----------------------------------------
      subroutine r8tcsp23(x,inx,th,inth,fspl,inf4,
     >   ibcxmin,bcxmin,ibcxmax,bcxmax,
     >   ibcthmin,bcthmin,ibcthmax,bcthmax,
     >   fspl2,wk,nwk,ilinx,ilinth,ier)

c  intermediary routines
c  call bcspline from tcspline loop
c  to set up 2d splines in each phi plane

!============
! idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
c     INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER inth,inf4,ibcxmin,ibcxmax,ibcthmin,ibcthmax,nwk,ilinx
      INTEGER ilinth,ier,inx,ith,ix,j,i
!============
      REAL*8 x(inx)                       ! x axis
      REAL*8 th(inth)                     ! th axis

      REAL*8 fspl(4,4,4,inf4,inth)        ! fspl array at one phi pt.
      REAL*8 fspl2(4,4,inx,inth)          ! temp fspl array for bcspline

      REAL*8 bcxmin(inth),bcxmax(inth)    ! d/dx BC's @ x(1),x(inx), th(*)
      REAL*8 bcthmin(inx),bcthmax(inx)    ! d/dth BC's @ th(1),th(inth), x(*)

      REAL*8 wk(nwk)

c--------------------

c  1.  copy spline data in

      do ith=1,inth
         do ix=1,inx
            fspl2(1,1,ix,ith)=fspl(1,1,1,ix,ith)
         enddo
      enddo

c  2.  compute the 2d spline

      call r8bcspline(x,inx,th,inth,fspl2,inx,
     >   ibcxmin,bcxmin,ibcxmax,bcxmax,
     >   ibcthmin,bcthmin,ibcthmax,bcthmax,
     >   wk,nwk,ilinx,ilinth,ier)
      if(ier.ne.0) return

c  3.  copy spline coeff results out

      do ith=1,inth-1
         do ix=1,inx-1
            do j=1,4
               do i=1,4
                  fspl(i,j,1,ix,ith)=fspl2(i,j,ix,ith)
               enddo
            enddo
         enddo
      enddo

      return
      end
