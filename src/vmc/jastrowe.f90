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
      use contrl_file, only: ounit
      use jastrow4_mod, only: jastrow_factor4
      use distances_mod, only: distances
      use optwf_control, only: ioptjas

#if defined(TREXIO_FOUND) && defined(QMCKL_FOUND) 
      use jastrow_qmckl_mod, only: jastrow_qmckl, jastrowe_qmckl
      use qmckl_data
#endif

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
            !UNDO
#if defined(TREXIO_FOUND) && defined(QMCKL_FOUND) 
         !if (ioptjas.eq.0) then
           call jastrowe_qmckl(iel, x(:,iel),fjn(1,1,iwf),d2n(iwf),fsumn(iwf),2)

           fsumn(iwf)=fsumn(iwf)+fsumo(iwf)
           d2n(iwf)=d2n(iwf)+d2o(iwf)
           do i=1,nelec
            fjn(1,i,iwf)=fjn(1,i,iwf)+fjo(1,i,iwf)
            fjn(2,i,iwf)=fjn(2,i,iwf)+fjo(2,i,iwf)
            fjn(3,i,iwf)=fjn(3,i,iwf)+fjo(3,i,iwf)
           enddo
         !else 
         !   call jastrow4e(iel,x,fjn(1,1,iwf),d2n(iwf),fsumn(iwf),fsn(1,1,iwf),fijn(1,1,1,iwf),d2ijn(1,1,iwf), &
         !   fjo(1,1,iwf),d2o(iwf),fsumo(iwf),fso(1,1,iwf),fijo(1,1,1,iwf),d2ijo(1,1,iwf),iflag)
         !endif
           ! print*, "begin jastrowe_qmckl"
           ! print*, "iel", iel
           ! print*, "x", x(:,iel)
           ! print*, "fjn", fjn(:,:,iwf)
           ! print*, "d2n", d2n(iwf)
           ! print*, "fsumn", fsumn(iwf)
           ! print*, "end jastrowe_qmckl"
#else 
           call jastrow4e(iel,x,fjn(1,1,iwf),d2n(iwf),fsumn(iwf),fsn(1,1,iwf),fijn(1,1,1,iwf),d2ijn(1,1,iwf), &
                 fjo(1,1,iwf),d2o(iwf),fsumo(iwf),fso(1,1,iwf),fijo(1,1,1,iwf),d2ijo(1,1,iwf),iflag)
           
           ! print*, "begin jastrowe old"
           ! print*, "iel", iel
           ! print*, "x", x(:,iel)
           ! print*, "fjn", fjn(:,:,iwf)
           ! print*, "d2n", d2n(iwf)
           ! print*, "fsumn", fsumn(iwf)
           ! print*, "end jastrowe old"
           ! print*, "d2o", d2o(iwf)
#endif
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
