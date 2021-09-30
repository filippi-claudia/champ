
c  bcspeval -- eval bicubic spline function and/or derivatives

      subroutine r8bcspeval(xget,yget,iselect,fval,
     >                    x,nx,y,ny,ilinx,iliny,f,inf3,ier)

      IMPLICIT NONE
c     INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      integer iselect(6)
      integer ilinx,iliny,nx,ny,inf3,ier

      REAL*8 xget,yget
      REAL*8 fval(*)
      REAL*8 x(nx),y(ny),f(4,4,inf3,ny)

c  modification -- dmc 11 Jan 1999 -- remove SAVE stmts; break routine
C    into these parts:

C    bcspevxy -- find grid cell of target pt.
C    bcspevfn -- evaluate function using output of bcpsevxy

C    in cases where multiple functions are defined on the same grid,
C    time can be saved by using bcspevxy once and then bcspevfn
C    multiple times.

c  input:
c     (xget,yget)   location where interpolated value is desired
c                   x(1).le.xget.le.x(nx) expected
c                   y(1).le.yget.le.y(ny) expected

c     iselect       select desired output

c                     iselect(1)=1 -- want function value (f) itself
c                     iselect(2)=1 -- want  df/dx
c                     iselect(3)=1 -- want  df/dy
c                     iselect(4)=1 -- want  d2f/dx2
c                     iselect(5)=1 -- want  d2f/dy2
c                     iselect(6)=1 -- want  d2f/dxdy

c              example:  iselect(1)=iselect(2)=iselect(3)=1
c                            f, df/dx, and df/dy all evaluated
c                        iselect(4)=iselect(5)=iselect(6)=0
c                            2nd derivatives not evaluated.

c                   see fval (output) description.

c     x(1...nx)     independent coordinate x, strict ascending
c     y(1...ny)     independent coordinate y, strict ascending

c     ilinx  --  =1: flag that x is linearly spaced (avoid search for speed)
c     iliny  --  =1: flag that y is linearly spaced (avoid search for speed)

c  **CAUTION** actual even spacing of x, y is NOT CHECKED HERE!


c     f             the function values (at grid points) and spline coefs

c  evaluation formula:  for point x btw x(i) and x(i+1), dx=x-x(i)
c                             and y btw y(j) and y(j+1), dy=y-y(j),

c      spline value =
c        f(1,1,i,j) + dx*f(2,1,i,j) + dx**2*f(3,1,i,j) + dx**3*f(4,1,i,j)
c   +dy*(f(1,2,i,j) + dx*f(2,2,i,j) + dx**2*f(3,2,i,j) + dx**3*f(4,2,i,j))
c   +d2*(f(1,3,i,j) + dx*f(2,3,i,j) + dx**2*f(3,3,i,j) + dx**3*f(4,3,i,j))
c   +d3*(f(1,4,i,j) + dx*f(2,4,i,j) + dx**2*f(3,4,i,j) + dx**3*f(4,4,i,j))

c      where d2=dy**2 and d3=dy**3.

c  output:
c      up to 6 elements of fval, ordered as follows:
c        fval(1)=function value or lowest order derivative requested
c        fval(2)=next order derivative
c             etc
c        the ordering is a subset of the sequence given under the "iselect"
c        description.

c      ier = 0 -- successful completion; = 1 -- an error occurred.

c-------------------------------------------------------------------
c  local

      integer i(1),j(1)

      REAL*8 dx(1),dy(1)

c--------------------------
         i(1)=0
         j(1)=0

      call r8bcspevxy(xget,yget,x,nx,y,ny,ilinx,iliny,
     >   i(1),j(1),dx(1),dy(1),ier)
      if(ier.ne.0) return

      call r8bcspevfn(iselect,1,1,fval,i,j,dx,dy,f,inf3,ny)

      return
      end

c-------------------------------------------------------------------------
c  bcspevxy -- look up x-y zone

c  this is the "first part" of bcspeval, see comments, above.

      subroutine r8bcspevxy(xget,yget,x,nx,y,ny,ilinx,iliny,
     >   i,j,dx,dy,ier)

!============
! idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
c     INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER nxm,nym,ii,jj
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 zxget,zyget,zxtol,zytol
!============
      integer nx,ny                     ! array dimensions

      REAL*8 xget,yget                    ! target point
      REAL*8 x(nx),y(ny)                  ! indep. coords.

      integer ilinx                     ! =1:  assume x evenly spaced
      integer iliny                     ! =1:  assume y evenly spaced

c  output of bcspevxy

      integer i,j                       ! index to cell containing target pt
      REAL*8 dx,dy                        ! displacement of target pt w/in cell
                                        ! dx=x-x(i)  dy=y-y(j)

      integer ier                       ! return ier.ne.0 on error

c------------------------------------

      ier=0

c  range check

      zxget=xget
      zyget=yget

      if((xget.lt.x(1)).or.(xget.gt.x(nx))) then
         zxtol=4.0D-7*max(abs(x(1)),abs(x(nx)))
         if((xget.lt.x(1)-zxtol).or.(xget.gt.x(nx)+zxtol)) then
            ier=1
            write(6,1001) xget,x(1),x(nx)
 1001       format(' ?bcspeval:  xget=',1pe11.4,' out of range ',
     >         1pe11.4,' to ',1pe11.4)
         else
            if((xget.lt.x(1)-0.5d0*zxtol).or.
     >         (xget.gt.x(nx)+0.5d0*zxtol))
     >      write(6,1011) xget,x(1),x(nx)
 1011       format(' %bcspeval:  xget=',1pe15.8,' beyond range ',
     >         1pe15.8,' to ',1pe15.8,' (fixup applied)')
            if(xget.lt.x(1)) then
               zxget=x(1)
            else
               zxget=x(nx)
            endif
         endif
      endif
      if((yget.lt.y(1)).or.(yget.gt.y(ny))) then
         zytol=4.0D-7*max(abs(y(1)),abs(y(ny)))
         if((yget.lt.y(1)-zytol).or.(yget.gt.y(ny)+zytol)) then
            ier=1
            write(6,1002) yget,y(1),y(ny)
 1002       format(' ?bcspeval:  yget=',1pe11.4,' out of range ',
     >         1pe11.4,' to ',1pe11.4)
         else
         if((yget.lt.y(1)-0.5d0*zytol).or.(yget.gt.y(ny)+0.5d0*zytol))
     >      write(6,1012) yget,y(1),y(ny)
 1012       format(' %bcspeval:  yget=',1pe15.8,' beyond range ',
     >         1pe15.8,' to ',1pe15.8,' (fixup applied)')
            if(yget.lt.y(1)) then
               zyget=y(1)
            else
               zyget=y(ny)
            endif
         endif
      endif
      if(ier.ne.0) return

c  now find interval in which target point lies..

      nxm=nx-1
      nym=ny-1

      if(ilinx.eq.1) then
         ii=1+nxm*(zxget-x(1))/(x(nx)-x(1))
         i=min(nxm, ii)
         if(zxget.lt.x(i)) then
            i=i-1
         else if(zxget.gt.x(i+1)) then
            i=i+1
         endif
      else
         if((1.le.i).and.(i.lt.nxm)) then
            if((x(i).le.zxget).and.(zxget.le.x(i+1))) then
               continue  ! already have the zone
            else
               call r8zonfind(x,nx,zxget,i)
            endif
         else
            call r8zonfind(x,nx,zxget,i)
         endif
      endif

      if(iliny.eq.1) then
         jj=1+nym*(zyget-y(1))/(y(ny)-y(1))
         j=min(nym, jj)
         if(zyget.lt.y(j)) then
            j=j-1
         else if(zyget.gt.y(j+1)) then
            j=j+1
         endif
      else
         if((1.le.j).and.(j.lt.nym)) then
            if((y(j).le.zyget).and.(zyget.le.y(j+1))) then
               continue  ! already have the zone
            else
               call r8zonfind(y,ny,zyget,j)
            endif
         else
            call r8zonfind(y,ny,zyget,j)
         endif
      endif

      dx=zxget-x(i)
      dy=zyget-y(j)

      return
      end
c------------------------------------------------------------------------
c  bcspevfn -- OK now evaluate the bicubic spline

      subroutine r8bcspevfn(ict,ivec,ivd,fval,iv,jv,dxv,dyv,f,inf3,ny)

c  input:
!============
! idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
c     INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER ny,iaval,i,j
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 dx,dy
!============
      integer ict(6)                    ! selector:
c        ict(1)=1 for f      (don't evaluate f if ict(1)=0)
c        ict(2)=1 for df/dx   ""
c        ict(3)=1 for df/dy   ""
c        ict(4)=1 for d2f/dx2
c        ict(5)=1 for d2f/dy2
c        ict(6)=1 for d2f/dxdy

c    note:  if ict(1)=-1, evaluate f,d2f/dx2,d2f/dy2,d4f/dx2dy2

      integer ivec,ivd                  ! vector dimensioning

c    ivec-- number of vector pts (spline values to look up)
c    ivd -- 1st dimension of fval, .ge.ivec

c output:
      REAL*8 fval(ivd,6)                 ! output array

c    v = index to element in vector;
c  fval(v,1) = first item requested by ict(...),
c  fval(v,2) = 2nd item requested,  ...etc...

c  input:
      integer iv(ivec),jv(ivec)         ! grid cell indices -- vectors
      REAL*8 dxv(ivec),dyv(ivec)          ! displacements w/in cell -- vectors

      integer inf3                      ! 3rd dimension of f -- .ge. nx
      REAL*8 f(4,4,inf3,ny)               ! bicubic fcn spline coeffs array

c  usage example:

c  1.  for each element (xx(v),yy(v)) in a vector of (x,y) pairs,
c    find the x and y zone indices and displacements with respect
c    to the "lower left corner" of the zone; store these in vectors
c    iv,jv and dxv,dyv.

c  2.  set ict(1)=0, ict(2)=1, ict(3)=1, the rest zero -- get only
c      the 1st derivatives.

c  3.  ivec is the length of the vector; ivd is the 1st dimension
c      of the array fval to receive the output

c      real fval(ivd,6)
c      real xv(ivd),yv(ivd)
c      integer iv(ivd),jv(ivd)
c      real dxv(ivd),dyv(ivd)
c      integer ict(6)

c      real fspline(4,4,nx,ny)  ! spline coeffs
c      data ict/0,1,1,0,0,0/    ! this call:  want 1st derivatives
c                               ! only ... these will be output to
c                               ! fval(*,1) fval(*,2)
c      ...
c      do iv=1,ivec
c        ...                    ! find indices and displacements
c      enddo
c      call bcspevfn(ict,ivec,ivd,fval,iv,jv,dxv,dyv,fspline,nx,ny)

c-------------------------------------------------------------------
c  local:

      integer v                         ! vector element index

c  OK can now do evaluations

      iaval=0  ! fval addressing

      if((ict(1).gt.0).or.(ict(1).eq.-1)) then
c  evaluate f
         iaval=iaval+1
         do v=1,ivec
            i=iv(v)
            j=jv(v)
            dx=dxv(v)
            dy=dyv(v)
            fval(v,iaval)=
     >       f(1,1,i,j)+dy*(f(1,2,i,j)+dy*(f(1,3,i,j)+dy*f(1,4,i,j)))
     >  +dx*(f(2,1,i,j)+dy*(f(2,2,i,j)+dy*(f(2,3,i,j)+dy*f(2,4,i,j)))
     >  +dx*(f(3,1,i,j)+dy*(f(3,2,i,j)+dy*(f(3,3,i,j)+dy*f(3,4,i,j)))
     >  +dx*(f(4,1,i,j)+dy*(f(4,2,i,j)+dy*(f(4,3,i,j)+dy*f(4,4,i,j))))))
         enddo
      endif

      if((ict(2).gt.0).and.(ict(1).ne.-1)) then
c  evaluate df/dx
         iaval=iaval+1
         do v=1,ivec
            i=iv(v)
            j=jv(v)
            dx=dxv(v)
            dy=dyv(v)
            fval(v,iaval)=
     >         f(2,1,i,j)+dy*(f(2,2,i,j)+dy*(f(2,3,i,j)+dy*f(2,4,i,j)))
     >       +2.0d0*dx*(
     >         f(3,1,i,j)+dy*(f(3,2,i,j)+dy*(f(3,3,i,j)+dy*f(3,4,i,j)))
     >       +1.5d0*dx*(
     >         f(4,1,i,j)+dy*(f(4,2,i,j)+dy*(f(4,3,i,j)+dy*f(4,4,i,j)))
     >         ))
         enddo
      endif

      if((ict(3).gt.0).and.(ict(1).ne.-1)) then
c  evaluate df/dy
         iaval=iaval+1
         do v=1,ivec
            i=iv(v)
            j=jv(v)
            dx=dxv(v)
            dy=dyv(v)
            fval(v,iaval)=
     >         f(1,2,i,j)+dy*(2.0d0*f(1,3,i,j)+dy*3.0d0*f(1,4,i,j))
     >      +dx*(f(2,2,i,j)+dy*(2.0d0*f(2,3,i,j)+dy*3.0d0*f(2,4,i,j))
     >      +dx*(f(3,2,i,j)+dy*(2.0d0*f(3,3,i,j)+dy*3.0d0*f(3,4,i,j))
     >      +dx*(f(4,2,i,j)+dy*(2.0d0*f(4,3,i,j)+dy*3.0d0*f(4,4,i,j))
     >         )))
         enddo
      endif

      if((ict(4).gt.0).or.(ict(1).eq.-1)) then
c  evaluate d2f/dx2
         iaval=iaval+1
         do v=1,ivec
            i=iv(v)
            j=jv(v)
            dx=dxv(v)
            dy=dyv(v)
            fval(v,iaval)=
     >        2.0d0*(
     >         f(3,1,i,j)+dy*(f(3,2,i,j)+dy*(f(3,3,i,j)+dy*f(3,4,i,j))))
     >       +6.0d0*dx*(
     >         f(4,1,i,j)+dy*(f(4,2,i,j)+dy*(f(4,3,i,j)+dy*f(4,4,i,j))))
         enddo
      endif

      if((ict(5).gt.0).or.(ict(1).eq.-1)) then
c  evaluate d2f/dy2
         iaval=iaval+1
         do v=1,ivec
            i=iv(v)
            j=jv(v)
            dx=dxv(v)
            dy=dyv(v)
            fval(v,iaval)=
     >         2.0d0*f(1,3,i,j)+6.0d0*dy*f(1,4,i,j)
     >         +dx*(2.0d0*f(2,3,i,j)+6.0d0*dy*f(2,4,i,j)
     >         +dx*(2.0d0*f(3,3,i,j)+6.0d0*dy*f(3,4,i,j)
     >         +dx*(2.0d0*f(4,3,i,j)+6.0d0*dy*f(4,4,i,j))))
         enddo
      endif

      if((ict(6).gt.0).and.(ict(1).ne.-1)) then
c  evaluate d2f/dxdy
         iaval=iaval+1
         do v=1,ivec
            i=iv(v)
            j=jv(v)
            dx=dxv(v)
            dy=dyv(v)
            fval(v,iaval)=
     >         f(2,2,i,j)+dy*(2.0d0*f(2,3,i,j)+dy*3.0d0*f(2,4,i,j))
     > +2.d0*dx*(f(3,2,i,j)+dy*(2.0d0*f(3,3,i,j)+dy*3.0d0*f(3,4,i,j))
     >+1.5d0*dx*(f(4,2,i,j)+dy*(2.0d0*f(4,3,i,j)+dy*3.0d0*f(4,4,i,j))
     >         ))
         enddo
      endif

      if(ict(1).eq.-1) then
         iaval=iaval+1
         do v=1,ivec
            i=iv(v)
            j=jv(v)
            dx=dxv(v)
            dy=dyv(v)
            fval(v,iaval)=
     >         4.0d0*f(3,3,i,j)+12.0d0*dy*f(3,4,i,j)
     >         +dx*(12.0d0*f(4,3,i,j)+36.0d0*dy*f(4,4,i,j))
         enddo
      endif

      return
      end
c----------------------
