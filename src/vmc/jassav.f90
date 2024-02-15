module jassav_mod
contains
      subroutine jassav(iel,iflag)
! Written by Claudia Filippi

      use contrl_per, only: iperiodic
      use ewald, only: cos_sum, cos_g, dcos_g
      use ewald, only: sin_sum, sin_g, dsin_g
      use ewald, only: cos_sum_new, cos_g_new, dcos_g_new
      use ewald, only: sin_sum_new, sin_g_new, dsin_g_new
      use jastrow_update, only: d2ijo, d2o, fijo, fjo, fso, fsumo
      use jastrow_update, only: d2ijn, d2n, fijn, fjn, fsn, fsumn
      use velocity_jastrow, only: vj, vjn
      use periodic, only: ngvec
      use system, only: nelec
      use vmc_mod, only: nwftypejas

      implicit none

      integer :: i, iel, iflag, j, k

      if(iperiodic.eq.1) then
        do i=1,ngvec
          cos_sum(i)=cos_sum_new(i)
          sin_sum(i)=sin_sum_new(i)
          cos_g(iel,i)=cos_g_new(iel,i)
          sin_g(iel,i)=sin_g_new(iel,i)
          do k=1,3
            dcos_g(k,iel,i)=dcos_g_new(k,iel,i)
            dsin_g(k,iel,i)=dsin_g_new(k,iel,i)
          enddo
        enddo
      endif

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
