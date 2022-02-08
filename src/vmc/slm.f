      subroutine slm(l,rvec,r2,y,dy,ddy,ddy_lap,dlapy,ider)

c   l   |   1  2  3  4    5      6     7    8    9   10  11  12  13  14  15  16  17  18  19
c ------+-------------------------------------------------------------------------------------
c   y   |   1  x  y  z  3z2-r2  x2-y2  xy   xz   yz  xxx yyy zzz xxy xxz yyx yyz zzx zzy xyz
c

      use precision_kinds, only: dp
      implicit none

      integer :: i, j, l, ider
      real(dp) :: cd1, cd2, cf, cf2, cf3
      real(dp) :: cp, cs, ddy_lap
      real(dp) :: r2, y
      real(dp), dimension(3) :: rvec
      real(dp), dimension(3) :: dy
      real(dp), dimension(3,3) :: ddy
      real(dp), dimension(3) :: dlapy
      real(dp), parameter :: half = 0.5d0

c cs=1/sqrt(4*pi), cp=sqrt(3/(4*pi)), cd1=sqrt(5/(4*pi)), cd2=sqrt(15/(4*pi))
      data cs,cp,cd1,cd2/0.28209479d0,0.48860251d0,0.63078313d0,1.0925484d0/
c cf=sqrt(7/(4*pi)),cf2=cf*sqrt(5),cf3=cf*sqrt(15)
      data cf,cf2,cf3/0.746352665180231d0,1.66889529453114d0,2.89061144264055d0/

      ddy_lap=0.d0
      do i=1,3
        dy(i)=0.d0
        dlapy(i)=0.d0
        do j=1,3
          ddy(i,j)=0.d0
        enddo
      enddo

      if(ider.eq.3) goto 200

      go to (101,102,103,104,105,106,107,108,109,110
     &      ,111,112,113,114,115,116,117,118,119) l

      stop 'ylm: l.gt.19'

c s
 101  y=cs
      return

c p
 102  y=cp*rvec(1)
      dy(1)=cp
      return

 103  y=cp*rvec(2)
      dy(2)=cp
      return

 104  y=cp*rvec(3)
      dy(3)=cp
      return

c d
 105  y=half*cd1*(3.d0*rvec(3)*rvec(3)-r2)
      dy(1)=-cd1*rvec(1)
      dy(2)=-cd1*rvec(2)
      dy(3)= cd1*2.d0*rvec(3)
      return

 106  y=half*cd2*(rvec(1)*rvec(1)-rvec(2)*rvec(2))
      dy(1)= cd2*rvec(1)
      dy(2)=-cd2*rvec(2)
      return

 107  y=cd2*rvec(1)*rvec(2)
      dy(1)=cd2*rvec(2)
      dy(2)=cd2*rvec(1)
      return

 108  y=cd2*rvec(1)*rvec(3)
      dy(1)=cd2*rvec(3)
      dy(3)=cd2*rvec(1)
      return

 109  y=cd2*rvec(2)*rvec(3)
      dy(2)=cd2*rvec(3)
      dy(3)=cd2*rvec(2)
      return

c There are 10 f functions as in GAMESS
 110  y=cf*rvec(1)*rvec(1)*rvec(1)
      dy(1)=cf*3.d0*rvec(1)*rvec(1)
      ddy_lap=cf*6.d0*rvec(1)
      return

 111  y=cf*rvec(2)*rvec(2)*rvec(2)
      dy(2)=cf*3.d0*rvec(2)*rvec(2)
      ddy_lap=cf*6.d0*rvec(2)
      return

 112  y=cf*rvec(3)*rvec(3)*rvec(3)
      dy(3)=cf*3.d0*rvec(3)*rvec(3)
      ddy_lap=cf*6.d0*rvec(3)
      return

 113  y=cf2*rvec(1)*rvec(1)*rvec(2)
      dy(1)=cf2*2.d0*rvec(1)*rvec(2)
      dy(2)=cf2*rvec(1)*rvec(1)
      ddy_lap=cf2*2.d0*rvec(2)
      return

 114  y=cf2*rvec(1)*rvec(1)*rvec(3)
      dy(1)=cf2*2.d0*rvec(1)*rvec(3)
      dy(3)=cf2*rvec(1)*rvec(1)
      ddy_lap=cf2*2.d0*rvec(3)
      return

 115  y=cf2*rvec(2)*rvec(2)*rvec(1)
      dy(2)=cf2*2.d0*rvec(2)*rvec(1)
      dy(1)=cf2*rvec(2)*rvec(2)
      ddy_lap=cf2*2.d0*rvec(1)
      return

 116  y=cf2*rvec(2)*rvec(2)*rvec(3)
      dy(2)=cf2*2.d0*rvec(2)*rvec(3)
      dy(3)=cf2*rvec(2)*rvec(2)
      ddy_lap=cf2*2.d0*rvec(3)
      return

 117  y=cf2*rvec(3)*rvec(3)*rvec(1)
      dy(1)=cf2*rvec(3)*rvec(3)
      dy(3)=cf2*2.d0*rvec(3)*rvec(1)
      ddy_lap=cf2*2.d0*rvec(1)
      return

 118  y=cf2*rvec(3)*rvec(3)*rvec(2)
      dy(2)=cf2*rvec(3)*rvec(3)
      dy(3)=cf2*2.d0*rvec(3)*rvec(2)
      ddy_lap=cf2*2.d0*rvec(2)
      return

 119  y=cf3*rvec(1)*rvec(2)*rvec(3)
      dy(1)=cf3*rvec(2)*rvec(3)
      dy(2)=cf3*rvec(1)*rvec(3)
      dy(3)=cf3*rvec(1)*rvec(2)
      return

 200  go to (201,202,203,204,205,206,207,208,209,210
     &      ,211,212,213,214,215,216,217,218,219) l

      stop 'ylm: l.gt.19'

 201  y=cs
      return

 202  y=cp*rvec(1)
      dy(1)=cp
      return

 203  y=cp*rvec(2)
      dy(2)=cp
      return

 204  y=cp*rvec(3)
      dy(3)=cp
      return

 205  y=half*cd1*(3.d0*rvec(3)*rvec(3)-r2)
      dy(1)=-cd1*rvec(1)
      dy(2)=-cd1*rvec(2)
      dy(3)= cd1*2.d0*rvec(3)
      ddy(1,1)=-cd1
      ddy(2,2)=-cd1
      ddy(3,3)= 2.d0*cd1
      return

 206  y=half*cd2*(rvec(1)*rvec(1)-rvec(2)*rvec(2))
      dy(1)= cd2*rvec(1)
      dy(2)=-cd2*rvec(2)
      ddy(1,1)= cd2
      ddy(2,2)=-cd2
      return

 207  y=cd2*rvec(1)*rvec(2)
      dy(1)=cd2*rvec(2)
      dy(2)=cd2*rvec(1)
      ddy(1,2)=cd2
      ddy(2,1)=cd2
      return

 208  y=cd2*rvec(1)*rvec(3)
      dy(1)=cd2*rvec(3)
      dy(3)=cd2*rvec(1)
      ddy(1,3)=cd2
      ddy(3,1)=cd2
      return

 209  y=cd2*rvec(2)*rvec(3)
      dy(2)=cd2*rvec(3)
      dy(3)=cd2*rvec(2)
      ddy(2,3)=cd2
      ddy(3,2)=cd2
      return

c There are 10 f functions as in GAMESS
 210  y=cf*rvec(1)*rvec(1)*rvec(1)
      dy(1)=cf*3.d0*rvec(1)*rvec(1)
      ddy_lap=cf*6.d0*rvec(1)
      ddy(1,1)=ddy_lap
      dlapy(1)=cf*6.d0
      return

 211  y=cf*rvec(2)*rvec(2)*rvec(2)
      dy(2)=cf*3.d0*rvec(2)*rvec(2)
      ddy_lap=cf*6.d0*rvec(2)
      ddy(2,2)=ddy_lap
      dlapy(2)=cf*6.d0
      return

 212  y=cf*rvec(3)*rvec(3)*rvec(3)
      dy(3)=cf*3.d0*rvec(3)*rvec(3)
      ddy_lap=cf*6.d0*rvec(3)
      ddy(3,3)=ddy_lap
      dlapy(3)=cf*6.d0
      return

 213  y=cf2*rvec(1)*rvec(1)*rvec(2)
      dy(1)=cf2*2.d0*rvec(1)*rvec(2)
      dy(2)=cf2*rvec(1)*rvec(1)
      ddy_lap=cf2*2.d0*rvec(2)
      ddy(1,1)=ddy_lap
      ddy(1,2)=cf2*2.d0*rvec(1)
      ddy(2,1)=ddy(1,2)
      dlapy(2)=cf2*2.d0
      return

 214  y=cf2*rvec(1)*rvec(1)*rvec(3)
      dy(1)=cf2*2.d0*rvec(1)*rvec(3)
      dy(3)=cf2*rvec(1)*rvec(1)
      ddy_lap=cf2*2.d0*rvec(3)
      ddy(1,1)=ddy_lap
      ddy(1,3)=cf2*2.d0*rvec(1)
      ddy(3,1)=ddy(1,3)
      dlapy(3)=cf2*2.d0
      return

 215  y=cf2*rvec(2)*rvec(2)*rvec(1)
      dy(2)=cf2*2.d0*rvec(2)*rvec(1)
      dy(1)=cf2*rvec(2)*rvec(2)
      ddy_lap=cf2*2.d0*rvec(1)
      ddy(2,2)=ddy_lap
      ddy(1,2)=cf2*2.d0*rvec(2)
      ddy(2,1)=ddy(1,2)
      dlapy(1)=cf2*2.d0
      return

 216  y=cf2*rvec(2)*rvec(2)*rvec(3)
      dy(2)=cf2*2.d0*rvec(2)*rvec(3)
      dy(3)=cf2*rvec(2)*rvec(2)
      ddy_lap=cf2*2.d0*rvec(3)
      ddy(2,2)=ddy_lap
      ddy(2,3)=cf2*2.d0*rvec(2)
      ddy(3,2)=ddy(2,3)
      dlapy(3)=cf2*2.d0
      return

 217  y=cf2*rvec(3)*rvec(3)*rvec(1)
      dy(1)=cf2*rvec(3)*rvec(3)
      dy(3)=cf2*2.d0*rvec(3)*rvec(1)
      ddy_lap=cf2*2.d0*rvec(1)
      ddy(3,3)=ddy_lap
      ddy(1,3)=cf2*2.d0*rvec(3)
      ddy(3,1)=ddy(1,3)
      dlapy(1)=cf2*2.d0
      return

 218  y=cf2*rvec(3)*rvec(3)*rvec(2)
      dy(2)=cf2*rvec(3)*rvec(3)
      dy(3)=cf2*2.d0*rvec(3)*rvec(2)
      ddy_lap=cf2*2.d0*rvec(2)
      ddy(3,3)=ddy_lap
      ddy(2,3)=cf2*2.d0*rvec(3)
      ddy(3,2)=ddy(2,3)
      dlapy(2)=cf2*2.d0
      return

 219  y=cf3*rvec(1)*rvec(2)*rvec(3)
      dy(1)=cf3*rvec(2)*rvec(3)
      dy(2)=cf3*rvec(1)*rvec(3)
      dy(3)=cf3*rvec(1)*rvec(2)
      ddy(1,2)=cf3*rvec(3)
      ddy(1,3)=cf3*rvec(2)
      ddy(2,3)=cf3*rvec(1)
      ddy(2,1)=ddy(1,2)
      ddy(3,1)=ddy(1,3)
      ddy(3,2)=ddy(2,3)
      return

      end
