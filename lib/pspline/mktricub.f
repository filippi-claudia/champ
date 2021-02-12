      subroutine mktricub(x,nx,y,ny,z,nz,f,nf2,nf3,
     >                    ibcxmin,bcxmin,ibcxmax,bcxmax,inb1x,
     >                    ibcymin,bcymin,ibcymax,bcymax,inb1y,
     >                    ibczmin,bczmin,ibczmax,bczmax,inb1z,
     >                    ilinx,iliny,ilinz,ier)

c  setup a tricubic spline; store coefficients in compatct form
c  (as per suggestion of L. Zakharov, PPPL, Feb. 1999)
C  8 coeffs per (x,y,z) grid point:
C          f,fxx,fyy,fzz,fxxyy,fxxzz,fyyzz,fxxyyzz

C  input:
      integer nx                        ! length of x vector
      integer ny                        ! length of y vector
      integer nz                        ! length of z vector
      real x(nx)                        ! x vector, strict ascending
      real y(ny)                        ! y vector, strict ascending
      real z(nz)                        ! z vector, strict ascending

      integer nf2                       ! 2nd dim. of f array, nf2.ge.nx
      integer nf3                       ! 3rd dim. of f array, nf3.ge.ny

c  input/output:

      real f(8,nf2,nf3,nz)              ! data and spline coefficients

C  on input:  f(1,i,j,k) = f(x(i),y(j),z(k))
C  on output:  f(1,i,j,k) unchanged
C              f(2,i,j,k) = d2f/dx2(x(i),y(j),z(k))
C              f(3,i,j,k) = d2f/dy2(x(i),y(j),z(k))
C              f(4,i,j,k) = d2f/dz2(x(i),y(j),z(k))
C              f(5,i,j,k) = d4f/dx2dy2(x(i),y(j),z(k))
C              f(6,i,j,k) = d4f/dx2dz2(x(i),y(j),z(k))
C              f(7,i,j,k) = d4f/dy2dz2(x(i),y(j),z(k))
C              f(8,i,j,k) = d6f/dx2dy2dz2(x(i),y(j),z(k))

C  there is a rather Hermite like interpolation formula to go with
C  this-- see evtricub.for.  Also the bicubic formula is given in
C  mkbicubw.for; the tricubic formula is precisely analogous.

C  boundary condition data
C  inputs:
      integer inb1x                     ! 1st dim of xmin & xmax bc arrays
      integer inb1y                     ! 1st dim of ymin & ymax bc arrays
      integer inb1z                     ! 1st dim of zmin & zmax bc arrays

      integer ibcxmin,ibcxmax           ! BC type flag @xmin, xmax
      integer ibcymin,ibcymax           ! BC type flag @ymin, ymax
      integer ibczmin,ibczmax           ! BC type flag @zmin, zmax

      real bcxmin(inb1x,nz),bcxmax(inb1x,nz) ! xmin & xmax BC data, ny x nz
      real bcymin(inb1y,nz),bcymax(inb1y,nz) ! ymin & ymax BC data, nx x nz
      real bczmin(inb1z,ny),bczmax(inb1z,ny) ! zmin & zmax BC data, nx x ny

c  where BC data is not required, dummy scalars may be passed.
C  the ibc* flags determine whether BC data isneeded.

c  BC data:  bcxmin & bcxmax:  BC vs. y,z @xmin,xmax
C            bcymin & bcymax:  BC vs. x,z @ymin,ymax
C            bczmin & bczmax:  BC vs. x,y @zmin,zmax

c   ibcxmin -- indicator for boundary condition at xmin=x(1):
c    bcxmin(...) -- boundary condition data
c     =-1 -- use periodic boundary condition
c     =0 -- use "not a knot"
c     =1 -- match slope, specified at x(1),y(iy),z(iz) by bcxmin(iy,iz)
c     =2 -- match 2nd derivative, specified at x(1),y(iy),z(iz)
c           by bcxmin(iy,iz
c     =3 -- boundary condition is slope=0 (df/dx=0) at x(1), all y(j)
c     =4 -- boundary condition is d2f/dx2=0 at x(1), all y(j)
c   ***NOTE bcxmin(...) referenced ONLY if ibcxmin=1 or ibcxmin=2

c   ibcxmax -- indicator for boundary condition at x(nx):
c    bcxmax(...) -- boundary condition data
c     (interpretation as with ibcxmin, bcxmin)
c     NOTE:  if ibcxmin=-1 then the periodic BC applies on both sides
c            and ibcxmax, bcxmax are ignored.
c   inb1x -- 1st dimension of bcxmin, bcxmax: if ibcxmin or ibcxmax .gt. 0
c            this must be .ge. ny.

c   interpretation of ibcymin,bcymin,ibcymax,bcymax,inb1y
c     is same as with ibcxmin,...

c   interpretation of ibczmin,bczmin,ibczmax,bczmax,inb1z
c     is same as with ibcxmin,...

c   the explicit bdy condition arrays are referenced only if the
c     corresponding "ibc" flag values are set to 1 or 2.

c  output:
      integer ilinx                     ! x vector equal spacing flag
      integer iliny                     ! y vector equal spacing flag
      integer ilinz                     ! z vector equal spacing flag

c   ilinx -- =1 on output if x(nx) pts are nearly evenly spaced (tol=1e-3)
c   iliny -- =1 on output if y(ny) evenly spaced (tol=1e-3)
c   ilinz -- =1 on output if z(nz) evenly spaced (tol=1e-3)

      integer ier                       ! exit code
c   ier -- completion code, 0 for normal

C-----------------------------------------------------
c  workspace **dynamic allocation**
C  f90 dynamic array

      real zdummy
!!ftoken:skip
      integer, parameter :: rkind = kind(zdummy)
      real(rkind), dimension(:), allocatable :: wk
!!ftoken:resume
      integer istat

c  this version requires a very large workspace, nwk.ge.80*nx*ny*nz
c  so as to be able to use tcspline to calculate the spline coefficients.
c--- if dynamic allocation fails, use mktricubw with a statically allocated
c  workspace.

c-----------------------------------------------------
c  OK now just call mktricubw

      inwk=80*nx*ny*nz

!!ftoken:skip
      allocate(wk(inwk),stat=istat)
!!ftoken:resume
      if(istat.ne.0) then
         write(6,*) ' ??mktricub:  workspace memory allocation failure!'
         ier=99
         return
      endif

      call mktricubw(x,nx,y,ny,z,nz,f,nf2,nf3,
     >                    ibcxmin,bcxmin,ibcxmax,bcxmax,inb1x,
     >                    ibcymin,bcymin,ibcymax,bcymax,inb1y,
     >                    ibczmin,bczmin,ibczmax,bczmax,inb1z,
!!ftoken:skip
     >                    wk,inwk,ilinx,iliny,ilinz,ier)

      if(allocated(wk)) deallocate(wk)
!!ftoken:resume

c  that's all

      return
      end

