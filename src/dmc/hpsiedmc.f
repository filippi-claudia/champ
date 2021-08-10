      subroutine psiedmc(iel,iw,coord,psid,psij,iflag)
c Written by Claudia Filippi


      use const, only: nelec
      use config, only: xold_dmc

      use precision_kinds, only: dp
      implicit none

      integer :: i, ic, idum, iel, iflag
      integer :: iw
      real(dp) :: psid, psij
      real(dp), dimension(3) :: coord
      real(dp), dimension(3, nelec) :: x


      do 10 ic=1,3
      do 10 i=1,iel-1
  10    x(ic,i)=xold_dmc(ic,i,iw,1)

      do 20 ic=1,3
  20    x(ic,iel)=coord(ic)

      do 30 ic=1,3
      do 30 i=iel+1,nelec
  30    x(ic,i)=xold_dmc(ic,i,iw,1)

      idum=1
      call psie(iel,x,psid,psij,idum,iflag)

      return
      end
