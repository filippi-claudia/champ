module jastrowe_mod
contains
      subroutine jastrowe(iel,x,v,d2j,psij,iflag)
! Written by Claudia Filippi by modifying jastrow

      use system, only: nelec
      use precision_kinds, only: dp
      use jastrow4e_mod, only: jastrow4e
      use jastrow1e_mod, only: jastrow1e
      use jastrow_update, only: d2ijo, d2o, fijo, fjo, fso, fsumo
      use jastrow_update, only: d2ijn, d2n, fijn, fjn, fsn, fsumn
      use multiple_geo, only: iwf
      use vmc_mod, only: nwftypejas
      use jastrow, only: ijas, ijas_lr
      use contrl_per, only: iperiodic
      use ewald_breakup, only: jastrow_longrange

      implicit none

      integer :: i, iel, iflag
      real(dp), dimension(*) :: d2j
      real(dp), dimension(*) :: psij
      real(dp), dimension(3, *) :: x
      real(dp), dimension(3, nelec, *) :: v

      real(dp) :: psij_per, d2_per
      real(dp), dimension(3, nelec) :: v_per

      if(ijas.eq.1) then

         psij_per=0.d0
         d2_per=0.d0
         v_per=0.d0
      
         if(iperiodic.eq.1.and.ijas_lr.eq.1) call jastrow_longrange(iel,x,psij_per,d2_per,v_per,0)

         do iwf=1,nwftypejas
            call jastrow1e(iel,x,fjn(1,1,iwf),d2n(iwf),fsumn(iwf),fsn(1,1,iwf),fijn(1,1,1,iwf),d2ijn(1,1,iwf), &
                 fjo(1,1,iwf),d2o(iwf),fsumo(iwf),fso(1,1,iwf),fijo(1,1,1,iwf),d2ijo(1,1,iwf),iflag)
            do i=1,nelec
               v(1,i,iwf)=fjn(1,i,iwf)+v_per(1,i)
               v(2,i,iwf)=fjn(2,i,iwf)+v_per(2,i)
               v(3,i,iwf)=fjn(3,i,iwf)+v_per(3,i)
            enddo
            psij(iwf)=fsumn(iwf)+psij_per
            d2j(iwf)=d2n(iwf)+d2_per
        enddo
      else
         do iwf=1,nwftypejas
            call jastrow4e(iel,x,fjn(1,1,iwf),d2n(iwf),fsumn(iwf),fsn(1,1,iwf),fijn(1,1,1,iwf),d2ijn(1,1,iwf), &
                 fjo(1,1,iwf),d2o(iwf),fsumo(iwf),fso(1,1,iwf),fijo(1,1,1,iwf),d2ijo(1,1,iwf),iflag)
            do i=1,nelec
               v(1,i,iwf)=fjn(1,i,iwf)
               v(2,i,iwf)=fjn(2,i,iwf)
               v(3,i,iwf)=fjn(3,i,iwf)
            enddo
            psij(iwf)=fsumn(iwf)
            d2j(iwf)=d2n(iwf)
         enddo
      endif

      iwf=1

      return
      end
end module
