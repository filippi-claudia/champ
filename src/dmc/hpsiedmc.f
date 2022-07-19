      module hpsiedmc
      contains
      subroutine psiedmc(iel,iw,coord,psid,psij,iflag)
c Written by Claudia Filippi


      use config, only: xold_dmc
      use hpsie, only: psie
      use mstates_mod, only: MSTATES

      use precision_kinds, only: dp
      use system, only: nelec
      implicit none

      integer :: i, ic, idum, iel, iflag
      integer :: iw
      real(dp) :: psid(*), psij
      real(dp), dimension(3) :: coord
      real(dp), dimension(3, nelec) :: x


      do ic=1,3
      do i=1,iel-1
        x(ic,i)=xold_dmc(ic,i,iw,1)
      enddo
      enddo

      do ic=1,3
        x(ic,iel)=coord(ic)
      enddo

      do ic=1,3
      do i=iel+1,nelec
        x(ic,i)=xold_dmc(ic,i,iw,1)
      enddo
      enddo

      idum=1
      call psie(iel,x,psid,psij,idum,iflag)

      return
      end
      end module
