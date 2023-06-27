module deriv_jastrow_mod
contains
      subroutine deriv_jastrow(x,v,d2j,psij)
! Written by Claudia Filippi


      use precision_kinds, only: dp
      use jastrow_update, only: d2ijo, d2o, fijo, fjo, fso, fsumo      
      use system,  only: nelec
      use multiple_geo, only: iwf
      use vmc_mod, only: nwftypejas
      use derivjas, only: d2g, g, go, gvalue
      use optwf_control, only: ioptjas
      use deriv_jastrow4_mod, only: deriv_jastrow4      

      implicit none

      integer :: i
      real(dp), dimension(3, *) :: x
      real(dp), dimension(3, nelec, *) :: v
      real(dp), dimension(*) :: psij
      real(dp), dimension(*) :: d2j

      if (ioptjas.eq.0) return

      if(nwftypejas.gt.1) then
        do iwf=1,nwftypejas
          call deriv_jastrow4(x,fjo(1,1,iwf),d2o(iwf),fsumo(iwf), &
                              fso(1,1,iwf),fijo(1,1,1,iwf), &
                              d2ijo(1,1,iwf),g(1,1,1,iwf),go(1,1,1,iwf), &
                              d2g(1,iwf),gvalue(1,iwf))
          psij(iwf)=fsumo(iwf)
          d2j(iwf)=d2o(iwf)
          do i=1,nelec
            v(1,i,iwf)=fjo(1,i,iwf)
            v(2,i,iwf)=fjo(2,i,iwf)
            v(3,i,iwf)=fjo(3,i,iwf)
          enddo
        enddo
      else
        call deriv_jastrow4(x,fjo(1,1,1),d2o(1),fsumo(1),fso(1,1,1), &
                            fijo(1,1,1,1),d2ijo(1,1,1),g(1,1,1,1), &
                            go(1,1,1,1),d2g(1,1),gvalue(1,1))
        psij(1)=fsumo(1)
        d2j(1)=d2o(1)
         do i=1,nelec
          v(1,i,1)=fjo(1,i,1)
          v(2,i,1)=fjo(2,i,1)
          v(3,i,1)=fjo(3,i,1)
        enddo
      endif


      return
      end
end module
