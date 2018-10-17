      subroutine r8splinck(x,inx,ilinx,ztol,ier)

C  check if a grid is strictly ascending and if it is evenly spaced
C  to w/in ztol

!============
! idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
c     INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
      INTEGER inx,ix
!============
! idecl:  explicitize implicit REAL declarations:
      REAL*8 dxavg,zeps,zdiffx,zdiff
!============
      REAL*8 x(inx)                       ! input -- grid to check

      integer ilinx                     ! output -- =1 if evenly spaced =2 O.W.

      REAL*8 ztol                         ! input -- spacing check tolerance

      integer ier                       ! output -- =0 if OK

C  ier=1 is returned if x(1...inx) is NOT STRICTLY ASCENDING...

C-------------------------------

      ier=0
      ilinx=1
      if(inx.le.1) return

      dxavg=(x(inx)-x(1))/(inx-1)
      zeps=abs(ztol*dxavg)

      do ix=2,inx
         zdiffx=(x(ix)-x(ix-1))
         if(zdiffx.le.0.0d0) ier=2
         zdiff=zdiffx-dxavg
         if(abs(zdiff).gt.zeps) then
            ilinx=2
         endif
      enddo
 10   continue

      return
      end
