      module jassav_mod
      contains
      subroutine jassav(iel,iflag)
c Written by Claudia Filippi

      use jastrow_update, only: d2ijo, d2o, fijo, fjo, fso, fsumo
      use velocity_jastrow, only: vj, vjn
      use jastrow_update, only: d2ijn, d2n, fijn, fjn, fsn, fsumn
      use system, only: nelec

      implicit none

      integer :: i, iel, iflag, j

      fsumo=fsumn
      do i=1,nelec
        fjo(1,i)=fjn(1,i)
        fjo(2,i)=fjn(2,i)
        fjo(3,i)=fjn(3,i)
      enddo

      do j=1,iel
        fso(iel,j)=fsn(iel,j)
      enddo

      do j=iel+1,nelec
        fso(j,iel)=fsn(j,iel)
      enddo

      do j=1,nelec
        fijo(1,iel,j)=fijn(1,iel,j)
        fijo(2,iel,j)=fijn(2,iel,j)
        fijo(3,iel,j)=fijn(3,iel,j)
        fijo(1,j,iel)=fijn(1,j,iel)
        fijo(2,j,iel)=fijn(2,j,iel)
        fijo(3,j,iel)=fijn(3,j,iel)
      enddo

      if(iflag.gt.0) then
        d2o=d2n
        do j=1,iel
          d2ijo(iel,j)=d2ijn(iel,j)
        enddo

        do j=iel+1,nelec
          d2ijo(j,iel)=d2ijn(j,iel)
        enddo
      endif

      do i=1,nelec
        vj(1,i)=vjn(1,i)
        vj(2,i)=vjn(2,i)
        vj(3,i)=vjn(3,i)
      enddo
      return
      end
      end module
