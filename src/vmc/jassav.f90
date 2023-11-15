module jassav_mod
contains
      subroutine jassav(iel,iflag)
! Written by Claudia Filippi

      use system, only: nelec
      use jastrow_update, only: d2ijo, d2o, fijo, fjo, fso, fsumo
      use velocity_jastrow, only: vj, vjn
      use jastrow_update, only: d2ijn, d2n, fijn, fjn, fsn, fsumn
      use vmc_mod, only: nwftypejas

      implicit none

      integer :: i, iel, iflag, j, k

      do k=1,nwftypejas
        fsumo(k)=fsumn(k)
        do i=1,nelec
          fjo(1,i,k)=fjn(1,i,k)
          fjo(2,i,k)=fjn(2,i,k)
          fjo(3,i,k)=fjn(3,i,k)
        enddo

        do j=1,iel
          fso(iel,j,k)=fsn(iel,j,k)
        enddo

        do j=iel+1,nelec
          fso(j,iel,k)=fsn(j,iel,k)
        enddo

        do j=1,nelec
          fijo(1,iel,j,k)=fijn(1,iel,j,k)
          fijo(2,iel,j,k)=fijn(2,iel,j,k)
          fijo(3,iel,j,k)=fijn(3,iel,j,k)
          fijo(1,j,iel,k)=fijn(1,j,iel,k)
          fijo(2,j,iel,k)=fijn(2,j,iel,k)
          fijo(3,j,iel,k)=fijn(3,j,iel,k)
        enddo

        if(iflag.gt.0) then
          d2o(k)=d2n(k)
          do j=1,iel
            d2ijo(iel,j,k)=d2ijn(iel,j,k)
          enddo

          do j=iel+1,nelec
            d2ijo(j,iel,k)=d2ijn(j,iel,k)
          enddo
        endif

        do i=1,nelec
          vj(1,i,k)=vjn(1,i,k)
          vj(2,i,k)=vjn(2,i,k)
          vj(3,i,k)=vjn(3,i,k)
        enddo
      enddo
      return
      end
end module
