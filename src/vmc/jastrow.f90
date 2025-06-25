module jastrow_mod
contains
      subroutine jastrow_factor(x,v,d2j,psij,ifr)
! Written by Cyrus Umrigar

      use contrl_file, only: ounit
      use contrl_per, only: iperiodic
      use cuspmat4, only: nterms
      use derivjas, only: d2g, g, go, gvalue
      use deriv_jastrow4_mod, only: deriv_jastrow4
      use deriv_jastrow1_mod, only: deriv_jastrow1
      
      use ewald_breakup, only: jastrow_longrange
      use jastrow, only: ijas, ijas_lr
      use jastrow_update, only: d2ijo, d2o, fijo, fjo, fso, fsumo
      use jastrow1_mod, only: jastrow_factor1
      use jastrow4_mod, only: jastrow_factor4
      use multiple_geo, only: iwf, nforce
      use optwf_control, only: ioptjas
      use optwf_parms, only: nparmj
      use precision_kinds, only: dp
      use system, only: nelec, nctype
      use vmc_mod, only: nwftypejas
      use jastrow, only: norda, nordb, nordc
      use optwf_nparmj, only: nparmc

#if defined(TREXIO_FOUND) && defined(QMCKL_FOUND) 
      use jastrow_qmckl_mod, only: jastrow_qmckl
      use deriv_jastrow_qmckl_mod, only: deriv_jastrow4_qmckl
      
      use qmckl_data
#endif

      implicit none

      integer :: i, j, jwf, ifr, rc
      real(dp) :: psij_per, d2_per
      real(dp), dimension(3, *) :: x
      real(dp), dimension(3, nelec, nwftypejas) :: v
      real(dp), dimension(3, nelec) :: v_per
      real(dp), dimension(nwftypejas) :: psij
      real(dp), dimension(nwftypejas) :: d2j
      real(dp), dimension(norda+1,nctype) :: deriv_a
      real(dp), dimension(nordb+1) :: deriv_b
      real(dp), dimension(nterms,nctype) :: deriv_c

      real(dp), dimension(3,nelec,nparmj,1) :: g_champ
      real(dp), dimension(nparmj,1) :: d2g_champ
      integer :: idim, iel, iparm

      real(dp) :: a,b
!     real(dp) :: psij_per0
!     real(dp), dimension(3, nelec) :: v_per0

! add long-range Jastrow (fixed parameters and equal for all states/geometries)

      psij_per=0.d0
      d2_per=0.d0
      v_per=0.d0
      
!     if(iperiodic.eq.1) then
!       call jastrow_longrange(0,x,psij_per0,d2_per,v_per0,0)
!       b=0.d0
!       do i=1,nelec
!        do j=1,3
!         a=x(j,i)
!         x(j,i)=x(j,i)+1.d-6
!         call jastrow_longrange(0,x,psij_per,d2_per,v_per,0)
!         write(ounit,*)'grad ',v_per0(j,i),(psij_per-psij_per0)/1.d-6
!         b=b+(v_per(j,i)-v_per0(j,i))/1.d-6
!         x(j,i)=a
!        enddo
!       enddo
!       write(ounit,*)'lapl ',d2_per,b
!     endif

      if(iperiodic.eq.1.and.ijas_lr.eq.1) call jastrow_longrange(0,x,psij_per,d2_per,v_per,0)

!     psij_per=0.d0
!     d2_per=0.d0
!     v_per=0.d0

! keep an option for ifr 1 and ioptjas 0 so we don't play with iwf
      if (nforce.gt.1) then

         if(ifr.gt.1.or.ioptjas.eq.0) then
            if(ijas.eq.1) then
               call jastrow_factor1(x,fjo(1,1,1),d2o(1),fsumo(1),fso(1,1,1),fijo(1,1,1,1),d2ijo(1,1,1))
            else
               call jastrow_factor4(x,fjo(1,1,1),d2o(1),fsumo(1),fso(1,1,1),fijo(1,1,1,1),d2ijo(1,1,1))
            endif
         else
            if(ijas.eq.1) then
               call deriv_jastrow1(x,fjo(1,1,1),d2o(1),fsumo(1),fso(1,1,1), &
                    fijo(1,1,1,1),d2ijo(1,1,1),g(1,1,1,1), &
                    go(1,1,1,1),d2g(1,1),gvalue(1,1))
            else
               call deriv_jastrow4(x,fjo(1,1,1),d2o(1),fsumo(1),fso(1,1,1), &
                    fijo(1,1,1,1),d2ijo(1,1,1),g(1,1,1,1), &
                    go(1,1,1,1),d2g(1,1),gvalue(1,1))
            endif
         endif
         psij(1)=fsumo(1)+psij_per
         d2j(1)=d2o(1)+d2_per
         do i=1,nelec
           v(1,i,1)=fjo(1,i,1)+v_per(1,i)
           v(2,i,1)=fjo(2,i,1)+v_per(2,i)
           v(3,i,1)=fjo(3,i,1)+v_per(3,i)
         enddo

! nforce.eq.1 -> iwf not used in relation to nforce, maybe used for multiple jastrow
       else

        if(ioptjas.eq.0) then
          if(ijas.eq.1) then
            do jwf=1,nwftypejas
              iwf=jwf
              call jastrow_factor1(x,fjo(1,1,iwf),d2o(iwf),fsumo(iwf),fso(1,1,iwf), &
                   fijo(1,1,1,iwf),d2ijo(1,1,iwf))
            enddo
           else
            do jwf=1,nwftypejas
              iwf=jwf
!UNDO
#if defined(TREXIO_FOUND) && defined(QMCKL_FOUND) 

              call jastrow_qmckl(x,fjo(1,1,iwf),d2o(iwf),fsumo(iwf))
#else
              call jastrow_factor4(x,fjo(1,1,iwf),d2o(iwf),fsumo(iwf),fso(1,1,iwf), &
                   fijo(1,1,1,iwf),d2ijo(1,1,iwf))
#endif
            enddo
          endif
         elseif(ioptjas.gt.0) then
          if(ijas.eq.1) then
            do jwf=1,nwftypejas
              iwf=jwf
              call deriv_jastrow1(x,fjo(1,1,iwf),d2o(iwf),fsumo(iwf), &
                   fso(1,1,iwf),fijo(1,1,1,iwf), &
                   d2ijo(1,1,iwf),g(1,1,1,iwf),go(1,1,1,iwf), &
                   d2g(1,iwf),gvalue(1,iwf))
            enddo
           else
            do jwf=1,nwftypejas
              iwf=jwf

#if defined(TREXIO_FOUND) && defined(QMCKL_FOUND)
              call deriv_jastrow4_qmckl(x,fjo(:,:,iwf),d2o(iwf),fsumo(iwf),g(:,:,:,iwf),d2g(:,iwf),gvalue(:,iwf))
              ! print*, "begin deriv_jastrow4_qmckl"
              ! print*, "fjo", fjo(:,:,iwf)
              ! print*, "d2o", d2o(iwf)
              ! print*, "fsumo", fsumo(iwf)
              ! print*, "g", g(:,:,:,iwf)
              ! print*, "gvalue", gvalue(:,iwf)
              ! print*, "end deriv_jastrow4_qmckl"
#else
              call deriv_jastrow4(x,fjo(1,1,iwf),d2o(iwf),fsumo(iwf), &
                   fso(1,1,iwf),fijo(1,1,1,iwf), &
                   d2ijo(1,1,iwf),g(1,1,1,iwf),go(1,1,1,iwf), &
                   d2g(1,iwf),gvalue(1,iwf))
              
              ! print*, "begin deriv_jastrow4 old"
              ! print*, "fjo", fjo(:,:,iwf)
              ! print*, "d2o", d2o(iwf)
              ! print*, "fsumo", fsumo(iwf)
              ! print*, "g", g(:,:,:,iwf)
              ! print*, "gvalue", gvalue(:,iwf)
              ! print*, "end deriv_jastrow4 old"
#endif

            enddo
          endif
        endif
	do jwf=1,nwftypejas
          psij(jwf)=fsumo(jwf)+psij_per
          d2j(jwf)=d2o(jwf)+d2_per
	   do i=1,nelec
	     v(1,i,jwf)=fjo(1,i,jwf)+v_per(1,i)
	     v(2,i,jwf)=fjo(2,i,jwf)+v_per(2,i)
	     v(3,i,jwf)=fjo(3,i,jwf)+v_per(3,i)
	  enddo
        enddo

! endif nforce
      endif

      return
      end
end module
