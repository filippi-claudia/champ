C---------------------------------------------------------------------
C  evaluate C1 cubic Hermite function interpolation -- 3d fcn
C   --vectorized-- dmc 10 Feb 1999

      subroutine r8fvtricub(ict,ivec,ivecd,
     >   fval,ii,jj,kk,xparam,yparam,zparam,
     >   hx,hxi,hy,hyi,hz,hzi,
     >   fin,inf2,inf3,nz)

C  use mktricub to set up spline coefficients...

!============
! idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
c     INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER inf3,nz,inf2,i,j,k,iadr
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 z36th,z216th,xp,xpi,xp2,xpi2,cx,cxi,hx2,cxd,cxdi,yp
      REAL*8 ypi,yp2,ypi2,cy,cyi,hy2,cyd,cydi,zp,zpi,zp2,zpi2,cz
      REAL*8 czi,hz2,czd,czdi
!============
      integer ict(10)                   ! requested output control
      integer ivec                      ! vector length
      integer ivecd                     ! vector dimension (1st dim of fval)

      integer ii(ivec),jj(ivec),kk(ivec) ! target cells (i,j,k)
      REAL*8 xparam(ivec),yparam(ivec),zparam(ivec)
                          ! normalized displacements from (i,j,k) corners

      REAL*8 hx(ivec),hy(ivec),hz(ivec)   ! grid spacing, and
      REAL*8 hxi(ivec),hyi(ivec),hzi(ivec) ! inverse grid spacing
           ! 1/(x(i+1)-x(i)) & 1/(y(j+1)-y(j)) & 1/(z(k+1)-z(i))

      REAL*8 fin(0:7,inf2,inf3,nz)        ! interpolant data (cf "evtricub")

      REAL*8 fval(ivecd,10)               ! output returned

C  for detailed description of fin, ict and fval see subroutine evtricub
C  comments.  Note ict is not vectorized; the same output
C  is expected to be returned for all input vector data points.

C  note that the index inputs ii,jj,kk and parameter inputs
C     xparam,yparam,zparam,hx,hxi,hy,hyi,hz,hzi are vectorized, and the
C     output array fval has a vector ** 1st dimension ** whose
C     size must be given as a separate argument

C  to use this routine in scalar mode, pass in ivec=ivecd=1

C---------------

      integer v

      REAL*8 sum
      REAL*8 sixth

      data sixth/0.166666666666666667d0/

C---------------

      z36th=sixth*sixth
      z216th=sixth*sixth*sixth

C  prepare useful parameters...

      do v=1,ivec
         i=ii(v)
         j=jj(v)
         k=kk(v)

C   ...in x direction

         xp=xparam(v)
         xpi=1.0d0-xp
         xp2=xp*xp
         xpi2=xpi*xpi

c        initialize values
         cx = 0.0
         cxi = 0.0
         hx2 = 0.0

         cxd = 0.0
         cxdi = 0.0

         if((ict(1).eq.1).or.(ict(3).eq.1).or.(ict(4).eq.1).or.
     >      (ict(6).eq.1).or.(ict(7).eq.1).or.(ict(10).eq.1)) then
            cx=xp*(xp2-1.0d0)
            cxi=xpi*(xpi2-1.0d0)
            hx2=hx(v)*hx(v)
         endif
         if((ict(2).eq.1).or.(ict(8).eq.1).or.(ict(9).eq.1)) then
            cxd=3.0d0*xp2-1.0d0
            cxdi=-3.0d0*xpi2+1.0d0
         endif

C   ...and in y direction

         yp=yparam(v)
         ypi=1.0d0-yp
         yp2=yp*yp
         ypi2=ypi*ypi

c        initialize values
         cy = 0.0
         cyi = 0.0
         hy2 = 0.0

         cyd = 0.0
         cydi = 0.0

         if((ict(1).eq.1).or.(ict(2).eq.1).or.(ict(4).eq.1).or.
     >      (ict(5).eq.1).or.(ict(7).eq.1).or.(ict(9).eq.1)) then
            cy=yp*(yp2-1.0d0)
            cyi=ypi*(ypi2-1.0d0)
            hy2=hy(v)*hy(v)
         endif
         if((ict(3).eq.1).or.(ict(8).eq.1).or.(ict(10).eq.1)) then
            cyd=3.0d0*yp2-1.0d0
            cydi=-3.0d0*ypi2+1.0d0
         endif

C   ...and in z direction

         zp=zparam(v)
         zpi=1.0d0-zp
         zp2=zp*zp
         zpi2=zpi*zpi

c        initialize values
         cz = 0.0
         czi = 0.0
         hz2 = 0.0

         czd = 0.0
         czdi = 0.0

         if((ict(1).eq.1).or.(ict(2).eq.1).or.(ict(3).eq.1).or.
     >      (ict(5).eq.1).or.(ict(6).eq.1).or.(ict(8).eq.1)) then
            cz=zp*(zp2-1.0d0)
            czi=zpi*(zpi2-1.0d0)
            hz2=hz(v)*hz(v)
         endif
         if((ict(4).eq.1).or.(ict(9).eq.1).or.(ict(10).eq.1)) then
            czd=3.0d0*zp2-1.0d0
            czdi=-3.0d0*zpi2+1.0d0
         endif

         iadr=0

C  get desired values:

         if(ict(1).eq.1) then

C  function value:

            iadr=iadr+1
            sum=(
     >         zpi*(
     >           xpi*(ypi*fin(0,i,j,k)  +yp*fin(0,i,j+1,k))+
     >            xp*(ypi*fin(0,i+1,j,k)+yp*fin(0,i+1,j+1,k)))
     >        +zp*(
     >           xpi*(ypi*fin(0,i,j,k+1)  +yp*fin(0,i,j+1,k+1))+
     >            xp*(ypi*fin(0,i+1,j,k+1)+yp*fin(0,i+1,j+1,k+1))))

            sum=sum+sixth*hx2*(
     >         zpi*(
     >           cxi*(ypi*fin(1,i,j,k)  +yp*fin(1,i,j+1,k))+
     >            cx*(ypi*fin(1,i+1,j,k)+yp*fin(1,i+1,j+1,k)))
     >        +zp*(
     >           cxi*(ypi*fin(1,i,j,k+1)  +yp*fin(1,i,j+1,k+1))+
     >            cx*(ypi*fin(1,i+1,j,k+1)+yp*fin(1,i+1,j+1,k+1))))

            sum=sum+sixth*hy2*(
     >         zpi*(
     >           xpi*(cyi*fin(2,i,j,k)  +cy*fin(2,i,j+1,k))+
     >            xp*(cyi*fin(2,i+1,j,k)+cy*fin(2,i+1,j+1,k)))
     >        +zp*(
     >           xpi*(cyi*fin(2,i,j,k+1)  +cy*fin(2,i,j+1,k+1))+
     >            xp*(cyi*fin(2,i+1,j,k+1)+cy*fin(2,i+1,j+1,k+1))))

            sum=sum+sixth*hz2*(
     >         czi*(
     >           xpi*(ypi*fin(3,i,j,k)  +yp*fin(3,i,j+1,k))+
     >            xp*(ypi*fin(3,i+1,j,k)+yp*fin(3,i+1,j+1,k)))
     >        +cz*(
     >           xpi*(ypi*fin(3,i,j,k+1)  +yp*fin(3,i,j+1,k+1))+
     >            xp*(ypi*fin(3,i+1,j,k+1)+yp*fin(3,i+1,j+1,k+1))))

            sum=sum+z36th*hx2*hy2*(
     >         zpi*(
     >           cxi*(cyi*fin(4,i,j,k)  +cy*fin(4,i,j+1,k))+
     >            cx*(cyi*fin(4,i+1,j,k)+cy*fin(4,i+1,j+1,k)))
     >        +zp*(
     >           cxi*(cyi*fin(4,i,j,k+1)  +cy*fin(4,i,j+1,k+1))+
     >            cx*(cyi*fin(4,i+1,j,k+1)+cy*fin(4,i+1,j+1,k+1))))

            sum=sum+z36th*hx2*hz2*(
     >         czi*(
     >           cxi*(ypi*fin(5,i,j,k)  +yp*fin(5,i,j+1,k))+
     >            cx*(ypi*fin(5,i+1,j,k)+yp*fin(5,i+1,j+1,k)))
     >        +cz*(
     >           cxi*(ypi*fin(5,i,j,k+1)  +yp*fin(5,i,j+1,k+1))+
     >            cx*(ypi*fin(5,i+1,j,k+1)+yp*fin(5,i+1,j+1,k+1))))

            sum=sum+z36th*hy2*hz2*(
     >         czi*(
     >           xpi*(cyi*fin(6,i,j,k)  +cy*fin(6,i,j+1,k))+
     >            xp*(cyi*fin(6,i+1,j,k)+cy*fin(6,i+1,j+1,k)))
     >        +cz*(
     >           xpi*(cyi*fin(6,i,j,k+1)  +cy*fin(6,i,j+1,k+1))+
     >            xp*(cyi*fin(6,i+1,j,k+1)+cy*fin(6,i+1,j+1,k+1))))

            sum=sum+z216th*hx2*hy2*hz2*(
     >         czi*(
     >           cxi*(cyi*fin(7,i,j,k)  +cy*fin(7,i,j+1,k))+
     >            cx*(cyi*fin(7,i+1,j,k)+cy*fin(7,i+1,j+1,k)))
     >        +cz*(
     >           cxi*(cyi*fin(7,i,j,k+1)  +cy*fin(7,i,j+1,k+1))+
     >            cx*(cyi*fin(7,i+1,j,k+1)+cy*fin(7,i+1,j+1,k+1))))

            fval(v,iadr)=sum
         endif

         if(ict(2).eq.1) then

C  df/dx:

            iadr=iadr+1

            sum=hxi(v)*(
     >         zpi*(
     >              -(ypi*fin(0,i,j,k)  +yp*fin(0,i,j+1,k))
     >              +(ypi*fin(0,i+1,j,k)+yp*fin(0,i+1,j+1,k)))
     >        +zp*(
     >              -(ypi*fin(0,i,j,k+1)  +yp*fin(0,i,j+1,k+1))
     >              +(ypi*fin(0,i+1,j,k+1)+yp*fin(0,i+1,j+1,k+1))))

            sum=sum+sixth*hx(v)*(
     >         zpi*(
     >           cxdi*(ypi*fin(1,i,j,k)  +yp*fin(1,i,j+1,k))+
     >            cxd*(ypi*fin(1,i+1,j,k)+yp*fin(1,i+1,j+1,k)))
     >        +zp*(
     >           cxdi*(ypi*fin(1,i,j,k+1)  +yp*fin(1,i,j+1,k+1))+
     >            cxd*(ypi*fin(1,i+1,j,k+1)+yp*fin(1,i+1,j+1,k+1))))

            sum=sum+sixth*hxi(v)*hy2*(
     >         zpi*(
     >              -(cyi*fin(2,i,j,k)  +cy*fin(2,i,j+1,k))
     >              +(cyi*fin(2,i+1,j,k)+cy*fin(2,i+1,j+1,k)))
     >        +zp*(
     >              -(cyi*fin(2,i,j,k+1)  +cy*fin(2,i,j+1,k+1))
     >              +(cyi*fin(2,i+1,j,k+1)+cy*fin(2,i+1,j+1,k+1))))

            sum=sum+sixth*hxi(v)*hz2*(
     >         czi*(
     >              -(ypi*fin(3,i,j,k)  +yp*fin(3,i,j+1,k))
     >              +(ypi*fin(3,i+1,j,k)+yp*fin(3,i+1,j+1,k)))
     >        +cz*(
     >              -(ypi*fin(3,i,j,k+1)  +yp*fin(3,i,j+1,k+1))
     >              +(ypi*fin(3,i+1,j,k+1)+yp*fin(3,i+1,j+1,k+1))))

            sum=sum+z36th*hx(v)*hy2*(
     >         zpi*(
     >           cxdi*(cyi*fin(4,i,j,k)  +cy*fin(4,i,j+1,k))+
     >            cxd*(cyi*fin(4,i+1,j,k)+cy*fin(4,i+1,j+1,k)))
     >        +zp*(
     >           cxdi*(cyi*fin(4,i,j,k+1)  +cy*fin(4,i,j+1,k+1))+
     >            cxd*(cyi*fin(4,i+1,j,k+1)+cy*fin(4,i+1,j+1,k+1))))

            sum=sum+z36th*hx(v)*hz2*(
     >         czi*(
     >           cxdi*(ypi*fin(5,i,j,k)  +yp*fin(5,i,j+1,k))+
     >            cxd*(ypi*fin(5,i+1,j,k)+yp*fin(5,i+1,j+1,k)))
     >        +cz*(
     >           cxdi*(ypi*fin(5,i,j,k+1)  +yp*fin(5,i,j+1,k+1))+
     >            cxd*(ypi*fin(5,i+1,j,k+1)+yp*fin(5,i+1,j+1,k+1))))

            sum=sum+z36th*hxi(v)*hy2*hz2*(
     >         czi*(
     >              -(cyi*fin(6,i,j,k)  +cy*fin(6,i,j+1,k))
     >              +(cyi*fin(6,i+1,j,k)+cy*fin(6,i+1,j+1,k)))
     >        +cz*(
     >              -(cyi*fin(6,i,j,k+1)  +cy*fin(6,i,j+1,k+1))
     >              +(cyi*fin(6,i+1,j,k+1)+cy*fin(6,i+1,j+1,k+1))))

            sum=sum+z216th*hx(v)*hy2*hz2*(
     >         czi*(
     >           cxdi*(cyi*fin(7,i,j,k)  +cy*fin(7,i,j+1,k))+
     >            cxd*(cyi*fin(7,i+1,j,k)+cy*fin(7,i+1,j+1,k)))
     >        +cz*(
     >           cxdi*(cyi*fin(7,i,j,k+1)  +cy*fin(7,i,j+1,k+1))+
     >            cxd*(cyi*fin(7,i+1,j,k+1)+cy*fin(7,i+1,j+1,k+1))))

            fval(v,iadr)=sum
         endif

         if(ict(3).eq.1) then

C  df/dy:

            iadr=iadr+1

            sum=hyi(v)*(
     >         zpi*(
     >           xpi*(-fin(0,i,j,k)  +fin(0,i,j+1,k))+
     >            xp*(-fin(0,i+1,j,k)+fin(0,i+1,j+1,k)))
     >        +zp*(
     >           xpi*(-fin(0,i,j,k+1)  +fin(0,i,j+1,k+1))+
     >            xp*(-fin(0,i+1,j,k+1)+fin(0,i+1,j+1,k+1))))

            sum=sum+sixth*hyi(v)*hx2*(
     >         zpi*(
     >           cxi*(-fin(1,i,j,k)  +fin(1,i,j+1,k))+
     >            cx*(-fin(1,i+1,j,k)+fin(1,i+1,j+1,k)))
     >        +zp*(
     >           cxi*(-fin(1,i,j,k+1)  +fin(1,i,j+1,k+1))+
     >            cx*(-fin(1,i+1,j,k+1)+fin(1,i+1,j+1,k+1))))

            sum=sum+sixth*hy(v)*(
     >         zpi*(
     >           xpi*(cydi*fin(2,i,j,k)  +cyd*fin(2,i,j+1,k))+
     >            xp*(cydi*fin(2,i+1,j,k)+cyd*fin(2,i+1,j+1,k)))
     >        +zp*(
     >           xpi*(cydi*fin(2,i,j,k+1)  +cyd*fin(2,i,j+1,k+1))+
     >            xp*(cydi*fin(2,i+1,j,k+1)+cyd*fin(2,i+1,j+1,k+1))))

            sum=sum+sixth*hyi(v)*hz2*(
     >         czi*(
     >           xpi*(-fin(3,i,j,k)  +fin(3,i,j+1,k))+
     >            xp*(-fin(3,i+1,j,k)+fin(3,i+1,j+1,k)))
     >        +cz*(
     >           xpi*(-fin(3,i,j,k+1)  +fin(3,i,j+1,k+1))+
     >            xp*(-fin(3,i+1,j,k+1)+fin(3,i+1,j+1,k+1))))

            sum=sum+z36th*hx2*hy(v)*(
     >         zpi*(
     >           cxi*(cydi*fin(4,i,j,k)  +cyd*fin(4,i,j+1,k))+
     >            cx*(cydi*fin(4,i+1,j,k)+cyd*fin(4,i+1,j+1,k)))
     >        +zp*(
     >           cxi*(cydi*fin(4,i,j,k+1)  +cyd*fin(4,i,j+1,k+1))+
     >            cx*(cydi*fin(4,i+1,j,k+1)+cyd*fin(4,i+1,j+1,k+1))))

            sum=sum+z36th*hyi(v)*hx2*hz2*(
     >         czi*(
     >           cxi*(-fin(5,i,j,k)  +fin(5,i,j+1,k))+
     >            cx*(-fin(5,i+1,j,k)+fin(5,i+1,j+1,k)))
     >        +cz*(
     >           cxi*(-fin(5,i,j,k+1)  +fin(5,i,j+1,k+1))+
     >            cx*(-fin(5,i+1,j,k+1)+fin(5,i+1,j+1,k+1))))

            sum=sum+z36th*hy(v)*hz2*(
     >         czi*(
     >           xpi*(cydi*fin(6,i,j,k)  +cyd*fin(6,i,j+1,k))+
     >            xp*(cydi*fin(6,i+1,j,k)+cyd*fin(6,i+1,j+1,k)))
     >        +cz*(
     >           xpi*(cydi*fin(6,i,j,k+1)  +cyd*fin(6,i,j+1,k+1))+
     >            xp*(cydi*fin(6,i+1,j,k+1)+cyd*fin(6,i+1,j+1,k+1))))

            sum=sum+z216th*hx2*hy(v)*hz2*(
     >         czi*(
     >           cxi*(cydi*fin(7,i,j,k)  +cyd*fin(7,i,j+1,k))+
     >            cx*(cydi*fin(7,i+1,j,k)+cyd*fin(7,i+1,j+1,k)))
     >        +cz*(
     >           cxi*(cydi*fin(7,i,j,k+1)  +cyd*fin(7,i,j+1,k+1))+
     >            cx*(cydi*fin(7,i+1,j,k+1)+cyd*fin(7,i+1,j+1,k+1))))

            fval(v,iadr)=sum
         endif

         if(ict(4).eq.1) then

C  df/dz:

            iadr=iadr+1

            sum=hzi(v)*(
     >           -(
     >           xpi*(ypi*fin(0,i,j,k)  +yp*fin(0,i,j+1,k))+
     >            xp*(ypi*fin(0,i+1,j,k)+yp*fin(0,i+1,j+1,k)))
     >           +(
     >           xpi*(ypi*fin(0,i,j,k+1)  +yp*fin(0,i,j+1,k+1))+
     >            xp*(ypi*fin(0,i+1,j,k+1)+yp*fin(0,i+1,j+1,k+1))))

            sum=sum+sixth*hx2*hzi(v)*(
     >           -(
     >           cxi*(ypi*fin(1,i,j,k)  +yp*fin(1,i,j+1,k))+
     >            cx*(ypi*fin(1,i+1,j,k)+yp*fin(1,i+1,j+1,k)))
     >           +(
     >           cxi*(ypi*fin(1,i,j,k+1)  +yp*fin(1,i,j+1,k+1))+
     >            cx*(ypi*fin(1,i+1,j,k+1)+yp*fin(1,i+1,j+1,k+1))))

            sum=sum+sixth*hy2*hzi(v)*(
     >           -(
     >           xpi*(cyi*fin(2,i,j,k)  +cy*fin(2,i,j+1,k))+
     >            xp*(cyi*fin(2,i+1,j,k)+cy*fin(2,i+1,j+1,k)))
     >           +(
     >           xpi*(cyi*fin(2,i,j,k+1)  +cy*fin(2,i,j+1,k+1))+
     >            xp*(cyi*fin(2,i+1,j,k+1)+cy*fin(2,i+1,j+1,k+1))))

            sum=sum+sixth*hz(v)*(
     >         czdi*(
     >           xpi*(ypi*fin(3,i,j,k)  +yp*fin(3,i,j+1,k))+
     >            xp*(ypi*fin(3,i+1,j,k)+yp*fin(3,i+1,j+1,k)))
     >        +czd*(
     >           xpi*(ypi*fin(3,i,j,k+1)  +yp*fin(3,i,j+1,k+1))+
     >            xp*(ypi*fin(3,i+1,j,k+1)+yp*fin(3,i+1,j+1,k+1))))

            sum=sum+z36th*hx2*hy2*hzi(v)*(
     >           -(
     >           cxi*(cyi*fin(4,i,j,k)  +cy*fin(4,i,j+1,k))+
     >            cx*(cyi*fin(4,i+1,j,k)+cy*fin(4,i+1,j+1,k)))
     >           +(
     >           cxi*(cyi*fin(4,i,j,k+1)  +cy*fin(4,i,j+1,k+1))+
     >            cx*(cyi*fin(4,i+1,j,k+1)+cy*fin(4,i+1,j+1,k+1))))

            sum=sum+z36th*hx2*hz(v)*(
     >         czdi*(
     >           cxi*(ypi*fin(5,i,j,k)  +yp*fin(5,i,j+1,k))+
     >            cx*(ypi*fin(5,i+1,j,k)+yp*fin(5,i+1,j+1,k)))
     >        +czd*(
     >           cxi*(ypi*fin(5,i,j,k+1)  +yp*fin(5,i,j+1,k+1))+
     >            cx*(ypi*fin(5,i+1,j,k+1)+yp*fin(5,i+1,j+1,k+1))))

            sum=sum+z36th*hy2*hz(v)*(
     >         czdi*(
     >           xpi*(cyi*fin(6,i,j,k)  +cy*fin(6,i,j+1,k))+
     >            xp*(cyi*fin(6,i+1,j,k)+cy*fin(6,i+1,j+1,k)))
     >        +czd*(
     >           xpi*(cyi*fin(6,i,j,k+1)  +cy*fin(6,i,j+1,k+1))+
     >            xp*(cyi*fin(6,i+1,j,k+1)+cy*fin(6,i+1,j+1,k+1))))

            sum=sum+z216th*hx2*hy2*hz(v)*(
     >         czdi*(
     >           cxi*(cyi*fin(7,i,j,k)  +cy*fin(7,i,j+1,k))+
     >            cx*(cyi*fin(7,i+1,j,k)+cy*fin(7,i+1,j+1,k)))
     >        +czd*(
     >           cxi*(cyi*fin(7,i,j,k+1)  +cy*fin(7,i,j+1,k+1))+
     >            cx*(cyi*fin(7,i+1,j,k+1)+cy*fin(7,i+1,j+1,k+1))))

            fval(v,iadr)=sum
         endif

         if(ict(5).eq.1) then

C  d2f/dx2:

            iadr=iadr+1

            sum=(
     >         zpi*(
     >           xpi*(ypi*fin(1,i,j,k)  +yp*fin(1,i,j+1,k))+
     >            xp*(ypi*fin(1,i+1,j,k)+yp*fin(1,i+1,j+1,k)))
     >        +zp*(
     >           xpi*(ypi*fin(1,i,j,k+1)  +yp*fin(1,i,j+1,k+1))+
     >            xp*(ypi*fin(1,i+1,j,k+1)+yp*fin(1,i+1,j+1,k+1))))

            sum=sum+sixth*hy2*(
     >         zpi*(
     >           xpi*(cyi*fin(4,i,j,k)  +cy*fin(4,i,j+1,k))+
     >            xp*(cyi*fin(4,i+1,j,k)+cy*fin(4,i+1,j+1,k)))
     >        +zp*(
     >           xpi*(cyi*fin(4,i,j,k+1)  +cy*fin(4,i,j+1,k+1))+
     >            xp*(cyi*fin(4,i+1,j,k+1)+cy*fin(4,i+1,j+1,k+1))))

            sum=sum+sixth*hz2*(
     >         czi*(
     >           xpi*(ypi*fin(5,i,j,k)  +yp*fin(5,i,j+1,k))+
     >            xp*(ypi*fin(5,i+1,j,k)+yp*fin(5,i+1,j+1,k)))
     >        +cz*(
     >           xpi*(ypi*fin(5,i,j,k+1)  +yp*fin(5,i,j+1,k+1))+
     >            xp*(ypi*fin(5,i+1,j,k+1)+yp*fin(5,i+1,j+1,k+1))))

            sum=sum+z36th*hy2*hz2*(
     >         czi*(
     >           xpi*(cyi*fin(7,i,j,k)  +cy*fin(7,i,j+1,k))+
     >            xp*(cyi*fin(7,i+1,j,k)+cy*fin(7,i+1,j+1,k)))
     >        +cz*(
     >           xpi*(cyi*fin(7,i,j,k+1)  +cy*fin(7,i,j+1,k+1))+
     >            xp*(cyi*fin(7,i+1,j,k+1)+cy*fin(7,i+1,j+1,k+1))))

            fval(v,iadr)=sum
         endif

         if(ict(6).eq.1) then

C  d2f/dy2:

            iadr=iadr+1

            sum=(
     >         zpi*(
     >           xpi*(ypi*fin(2,i,j,k)  +yp*fin(2,i,j+1,k))+
     >            xp*(ypi*fin(2,i+1,j,k)+yp*fin(2,i+1,j+1,k)))
     >        +zp*(
     >           xpi*(ypi*fin(2,i,j,k+1)  +yp*fin(2,i,j+1,k+1))+
     >            xp*(ypi*fin(2,i+1,j,k+1)+yp*fin(2,i+1,j+1,k+1))))

            sum=sum+sixth*hx2*(
     >         zpi*(
     >           cxi*(ypi*fin(4,i,j,k)  +yp*fin(4,i,j+1,k))+
     >            cx*(ypi*fin(4,i+1,j,k)+yp*fin(4,i+1,j+1,k)))
     >        +zp*(
     >           cxi*(ypi*fin(4,i,j,k+1)  +yp*fin(4,i,j+1,k+1))+
     >            cx*(ypi*fin(4,i+1,j,k+1)+yp*fin(4,i+1,j+1,k+1))))

            sum=sum+sixth*hz2*(
     >         czi*(
     >           xpi*(ypi*fin(6,i,j,k)  +yp*fin(6,i,j+1,k))+
     >            xp*(ypi*fin(6,i+1,j,k)+yp*fin(6,i+1,j+1,k)))
     >        +cz*(
     >           xpi*(ypi*fin(6,i,j,k+1)  +yp*fin(6,i,j+1,k+1))+
     >            xp*(ypi*fin(6,i+1,j,k+1)+yp*fin(6,i+1,j+1,k+1))))

            sum=sum+z36th*hx2*hz2*(
     >         czi*(
     >           cxi*(ypi*fin(7,i,j,k)  +yp*fin(7,i,j+1,k))+
     >            cx*(ypi*fin(7,i+1,j,k)+yp*fin(7,i+1,j+1,k)))
     >        +cz*(
     >           cxi*(ypi*fin(7,i,j,k+1)  +yp*fin(7,i,j+1,k+1))+
     >            cx*(ypi*fin(7,i+1,j,k+1)+yp*fin(7,i+1,j+1,k+1))))

            fval(v,iadr)=sum
         endif

         if(ict(7).eq.1) then

C  d2f/dz2:

            iadr=iadr+1

            sum=(
     >         zpi*(
     >           xpi*(ypi*fin(3,i,j,k)  +yp*fin(3,i,j+1,k))+
     >            xp*(ypi*fin(3,i+1,j,k)+yp*fin(3,i+1,j+1,k)))
     >        +zp*(
     >           xpi*(ypi*fin(3,i,j,k+1)  +yp*fin(3,i,j+1,k+1))+
     >            xp*(ypi*fin(3,i+1,j,k+1)+yp*fin(3,i+1,j+1,k+1))))

            sum=sum+sixth*hx2*(
     >         zpi*(
     >           cxi*(ypi*fin(5,i,j,k)  +yp*fin(5,i,j+1,k))+
     >            cx*(ypi*fin(5,i+1,j,k)+yp*fin(5,i+1,j+1,k)))
     >        +zp*(
     >           cxi*(ypi*fin(5,i,j,k+1)  +yp*fin(5,i,j+1,k+1))+
     >            cx*(ypi*fin(5,i+1,j,k+1)+yp*fin(5,i+1,j+1,k+1))))

            sum=sum+sixth*hy2*(
     >         zpi*(
     >           xpi*(cyi*fin(6,i,j,k)  +cy*fin(6,i,j+1,k))+
     >            xp*(cyi*fin(6,i+1,j,k)+cy*fin(6,i+1,j+1,k)))
     >        +zp*(
     >           xpi*(cyi*fin(6,i,j,k+1)  +cy*fin(6,i,j+1,k+1))+
     >            xp*(cyi*fin(6,i+1,j,k+1)+cy*fin(6,i+1,j+1,k+1))))

            sum=sum+z36th*hx2*hy2*(
     >         zpi*(
     >           cxi*(cyi*fin(7,i,j,k)  +cy*fin(7,i,j+1,k))+
     >            cx*(cyi*fin(7,i+1,j,k)+cy*fin(7,i+1,j+1,k)))
     >        +zp*(
     >           cxi*(cyi*fin(7,i,j,k+1)  +cy*fin(7,i,j+1,k+1))+
     >            cx*(cyi*fin(7,i+1,j,k+1)+cy*fin(7,i+1,j+1,k+1))))

            fval(v,iadr)=sum
         endif

         if(ict(8).eq.1) then

C  d2f/dxdy:

            iadr=iadr+1

            sum=hxi(v)*hyi(v)*(
     >         zpi*(
     >               (fin(0,i,j,k)  -fin(0,i,j+1,k))-
     >               (fin(0,i+1,j,k)-fin(0,i+1,j+1,k)))
     >        +zp*(
     >               (fin(0,i,j,k+1)  -fin(0,i,j+1,k+1))-
     >               (fin(0,i+1,j,k+1)-fin(0,i+1,j+1,k+1))))

            sum=sum+sixth*hyi(v)*hx(v)*(
     >         zpi*(
     >           cxdi*(-fin(1,i,j,k)  +fin(1,i,j+1,k))+
     >            cxd*(-fin(1,i+1,j,k)+fin(1,i+1,j+1,k)))
     >        +zp*(
     >           cxdi*(-fin(1,i,j,k+1)  +fin(1,i,j+1,k+1))+
     >            cxd*(-fin(1,i+1,j,k+1)+fin(1,i+1,j+1,k+1))))

            sum=sum+sixth*hxi(v)*hy(v)*(
     >         zpi*(
     >              -(cydi*fin(2,i,j,k)  +cyd*fin(2,i,j+1,k))
     >              +(cydi*fin(2,i+1,j,k)+cyd*fin(2,i+1,j+1,k)))
     >        +zp*(
     >              -(cydi*fin(2,i,j,k+1)  +cyd*fin(2,i,j+1,k+1))
     >              +(cydi*fin(2,i+1,j,k+1)+cyd*fin(2,i+1,j+1,k+1))))

            sum=sum+sixth*hxi(v)*hyi(v)*hz2*(
     >         czi*(
     >               (fin(3,i,j,k)  -fin(3,i,j+1,k))-
     >               (fin(3,i+1,j,k)-fin(3,i+1,j+1,k)))
     >        +cz*(
     >               (fin(3,i,j,k+1)  -fin(3,i,j+1,k+1))-
     >               (fin(3,i+1,j,k+1)-fin(3,i+1,j+1,k+1))))

            sum=sum+z36th*hx(v)*hy(v)*(
     >         zpi*(
     >           cxdi*(cydi*fin(4,i,j,k)  +cyd*fin(4,i,j+1,k))+
     >            cxd*(cydi*fin(4,i+1,j,k)+cyd*fin(4,i+1,j+1,k)))
     >        +zp*(
     >           cxdi*(cydi*fin(4,i,j,k+1)  +cyd*fin(4,i,j+1,k+1))+
     >            cxd*(cydi*fin(4,i+1,j,k+1)+cyd*fin(4,i+1,j+1,k+1))))

            sum=sum+z36th*hyi(v)*hx(v)*hz2*(
     >         czi*(
     >           cxdi*(-fin(5,i,j,k)  +fin(5,i,j+1,k))+
     >            cxd*(-fin(5,i+1,j,k)+fin(5,i+1,j+1,k)))
     >        +cz*(
     >           cxdi*(-fin(5,i,j,k+1)  +fin(5,i,j+1,k+1))+
     >            cxd*(-fin(5,i+1,j,k+1)+fin(5,i+1,j+1,k+1))))

            sum=sum+z36th*hxi(v)*hy(v)*hz2*(
     >         czi*(
     >               -(cydi*fin(6,i,j,k)  +cyd*fin(6,i,j+1,k))
     >               +(cydi*fin(6,i+1,j,k)+cyd*fin(6,i+1,j+1,k)))
     >        +cz*(
     >               -(cydi*fin(6,i,j,k+1)  +cyd*fin(6,i,j+1,k+1))
     >               +(cydi*fin(6,i+1,j,k+1)+cyd*fin(6,i+1,j+1,k+1))))

            sum=sum+z216th*hx(v)*hy(v)*hz2*(
     >         czi*(
     >           cxdi*(cydi*fin(7,i,j,k)  +cyd*fin(7,i,j+1,k))+
     >            cxd*(cydi*fin(7,i+1,j,k)+cyd*fin(7,i+1,j+1,k)))
     >        +cz*(
     >           cxdi*(cydi*fin(7,i,j,k+1)  +cyd*fin(7,i,j+1,k+1))+
     >            cxd*(cydi*fin(7,i+1,j,k+1)+cyd*fin(7,i+1,j+1,k+1))))

            fval(v,iadr)=sum
         endif

         if(ict(9).eq.1) then

C  d2f/dxdz:

            iadr=iadr+1

            sum=hxi(v)*hzi(v)*(
     >            (
     >              (ypi*fin(0,i,j,k)  +yp*fin(0,i,j+1,k)) -
     >              (ypi*fin(0,i+1,j,k)+yp*fin(0,i+1,j+1,k)))
     >           -(
     >              (ypi*fin(0,i,j,k+1)  +yp*fin(0,i,j+1,k+1)) -
     >              (ypi*fin(0,i+1,j,k+1)+yp*fin(0,i+1,j+1,k+1))))

            sum=sum+sixth*hx(v)*hzi(v)*(
     >           -(
     >           cxdi*(ypi*fin(1,i,j,k)  +yp*fin(1,i,j+1,k))+
     >            cxd*(ypi*fin(1,i+1,j,k)+yp*fin(1,i+1,j+1,k)))
     >           +(
     >           cxdi*(ypi*fin(1,i,j,k+1)  +yp*fin(1,i,j+1,k+1))+
     >            cxd*(ypi*fin(1,i+1,j,k+1)+yp*fin(1,i+1,j+1,k+1))))

            sum=sum+sixth*hxi(v)*hy2*hzi(v)*(
     >            (
     >              (cyi*fin(2,i,j,k)  +cy*fin(2,i,j+1,k)) -
     >              (cyi*fin(2,i+1,j,k)+cy*fin(2,i+1,j+1,k)))
     >           -(
     >              (cyi*fin(2,i,j,k+1)  +cy*fin(2,i,j+1,k+1)) -
     >              (cyi*fin(2,i+1,j,k+1)+cy*fin(2,i+1,j+1,k+1))))

            sum=sum+sixth*hxi(v)*hz(v)*(
     >         czdi*(
     >              -(ypi*fin(3,i,j,k)  +yp*fin(3,i,j+1,k))
     >              +(ypi*fin(3,i+1,j,k)+yp*fin(3,i+1,j+1,k)))
     >        +czd*(
     >              -(ypi*fin(3,i,j,k+1)  +yp*fin(3,i,j+1,k+1))
     >              +(ypi*fin(3,i+1,j,k+1)+yp*fin(3,i+1,j+1,k+1))))

            sum=sum+z36th*hx(v)*hy2*hzi(v)*(
     >           -(
     >           cxdi*(cyi*fin(4,i,j,k)  +cy*fin(4,i,j+1,k))+
     >            cxd*(cyi*fin(4,i+1,j,k)+cy*fin(4,i+1,j+1,k)))
     >           +(
     >           cxdi*(cyi*fin(4,i,j,k+1)  +cy*fin(4,i,j+1,k+1))+
     >            cxd*(cyi*fin(4,i+1,j,k+1)+cy*fin(4,i+1,j+1,k+1))))

            sum=sum+z36th*hx(v)*hz(v)*(
     >         czdi*(
     >           cxdi*(ypi*fin(5,i,j,k)  +yp*fin(5,i,j+1,k))+
     >            cxd*(ypi*fin(5,i+1,j,k)+yp*fin(5,i+1,j+1,k)))
     >        +czd*(
     >           cxdi*(ypi*fin(5,i,j,k+1)  +yp*fin(5,i,j+1,k+1))+
     >            cxd*(ypi*fin(5,i+1,j,k+1)+yp*fin(5,i+1,j+1,k+1))))

            sum=sum+z36th*hxi(v)*hy2*hz(v)*(
     >         czdi*(
     >              -(cyi*fin(6,i,j,k)  +cy*fin(6,i,j+1,k))
     >              +(cyi*fin(6,i+1,j,k)+cy*fin(6,i+1,j+1,k)))
     >        +czd*(
     >              -(cyi*fin(6,i,j,k+1)  +cy*fin(6,i,j+1,k+1))
     >              +(cyi*fin(6,i+1,j,k+1)+cy*fin(6,i+1,j+1,k+1))))

            sum=sum+z216th*hx(v)*hy2*hz(v)*(
     >         czdi*(
     >           cxdi*(cyi*fin(7,i,j,k)  +cy*fin(7,i,j+1,k))+
     >            cxd*(cyi*fin(7,i+1,j,k)+cy*fin(7,i+1,j+1,k)))
     >        +czd*(
     >           cxdi*(cyi*fin(7,i,j,k+1)  +cy*fin(7,i,j+1,k+1))+
     >            cxd*(cyi*fin(7,i+1,j,k+1)+cy*fin(7,i+1,j+1,k+1))))

            fval(v,iadr)=sum
         endif

         if(ict(10).eq.1) then

C  d2f/dydz:

            iadr=iadr+1

            sum=hyi(v)*hzi(v)*(
     >            (
     >           xpi*(fin(0,i,j,k)  -fin(0,i,j+1,k))+
     >            xp*(fin(0,i+1,j,k)-fin(0,i+1,j+1,k)))
     >           -(
     >           xpi*(fin(0,i,j,k+1)  -fin(0,i,j+1,k+1))+
     >            xp*(fin(0,i+1,j,k+1)-fin(0,i+1,j+1,k+1))))

            sum=sum+sixth*hyi(v)*hx2*hzi(v)*(
     >            (
     >           cxi*(fin(1,i,j,k)  -fin(1,i,j+1,k))+
     >            cx*(fin(1,i+1,j,k)-fin(1,i+1,j+1,k)))
     >           -(
     >           cxi*(fin(1,i,j,k+1)  -fin(1,i,j+1,k+1))+
     >            cx*(fin(1,i+1,j,k+1)-fin(1,i+1,j+1,k+1))))

            sum=sum+sixth*hy(v)*hzi(v)*(
     >           -(
     >           xpi*(cydi*fin(2,i,j,k)  +cyd*fin(2,i,j+1,k))+
     >            xp*(cydi*fin(2,i+1,j,k)+cyd*fin(2,i+1,j+1,k)))
     >           +(
     >           xpi*(cydi*fin(2,i,j,k+1)  +cyd*fin(2,i,j+1,k+1))+
     >            xp*(cydi*fin(2,i+1,j,k+1)+cyd*fin(2,i+1,j+1,k+1))))

            sum=sum+sixth*hyi(v)*hz(v)*(
     >         czdi*(
     >           xpi*(-fin(3,i,j,k)  +fin(3,i,j+1,k))+
     >            xp*(-fin(3,i+1,j,k)+fin(3,i+1,j+1,k)))
     >        +czd*(
     >           xpi*(-fin(3,i,j,k+1)  +fin(3,i,j+1,k+1))+
     >            xp*(-fin(3,i+1,j,k+1)+fin(3,i+1,j+1,k+1))))

            sum=sum+z36th*hx2*hy(v)*hzi(v)*(
     >           -(
     >           cxi*(cydi*fin(4,i,j,k)  +cyd*fin(4,i,j+1,k))+
     >            cx*(cydi*fin(4,i+1,j,k)+cyd*fin(4,i+1,j+1,k)))
     >           +(
     >           cxi*(cydi*fin(4,i,j,k+1)  +cyd*fin(4,i,j+1,k+1))+
     >            cx*(cydi*fin(4,i+1,j,k+1)+cyd*fin(4,i+1,j+1,k+1))))

            sum=sum+z36th*hyi(v)*hx2*hz(v)*(
     >         czdi*(
     >           cxi*(-fin(5,i,j,k)  +fin(5,i,j+1,k))+
     >            cx*(-fin(5,i+1,j,k)+fin(5,i+1,j+1,k)))
     >        +czd*(
     >           cxi*(-fin(5,i,j,k+1)  +fin(5,i,j+1,k+1))+
     >            cx*(-fin(5,i+1,j,k+1)+fin(5,i+1,j+1,k+1))))

            sum=sum+z36th*hy(v)*hz(v)*(
     >         czdi*(
     >           xpi*(cydi*fin(6,i,j,k)  +cyd*fin(6,i,j+1,k))+
     >            xp*(cydi*fin(6,i+1,j,k)+cyd*fin(6,i+1,j+1,k)))
     >        +czd*(
     >           xpi*(cydi*fin(6,i,j,k+1)  +cyd*fin(6,i,j+1,k+1))+
     >            xp*(cydi*fin(6,i+1,j,k+1)+cyd*fin(6,i+1,j+1,k+1))))

            sum=sum+z216th*hx2*hy(v)*hz(v)*(
     >         czdi*(
     >           cxi*(cydi*fin(7,i,j,k)  +cyd*fin(7,i,j+1,k))+
     >            cx*(cydi*fin(7,i+1,j,k)+cyd*fin(7,i+1,j+1,k)))
     >        +czd*(
     >           cxi*(cydi*fin(7,i,j,k+1)  +cyd*fin(7,i,j+1,k+1))+
     >            cx*(cydi*fin(7,i+1,j,k+1)+cyd*fin(7,i+1,j+1,k+1))))

            fval(v,iadr)=sum
         endif

      enddo                             ! vector loop

      return
      end
