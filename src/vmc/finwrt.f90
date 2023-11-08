module finwrt_mod
contains
      subroutine finwrt
! Written by Cyrus Umrigar, modified by Claudia Filippi
! routine to print out final results

      use ci000,   only: iciprt
      use contrl_file, only: errunit,ounit
      use contrl_per, only: iperiodic
      use control, only: mode
      use control_vmc, only: vmc_nblk,vmc_nstep
      use csfs,    only: nstates
      use est2cm,  only: ecm2,ecm21,pecm2,tpbcm2
      use estcum,  only: ecum,ecum1,iblk,pecum,tpbcum
      use estsig,  only: ecm21s,ecum1s
      use estsum,  only: acc
      use finwrt_more_mod, only: finwrt_more
      use forcewt, only: wcum
      use grdntspar, only: igrdtype,ngradnts
      use header,  only: title
      use inputflags, only: iqmmm
      use misc_grdnts, only: finwrt_diaghess_zmat,finwrt_grdnts_cart
      use misc_grdnts, only: finwrt_grdnts_zmat
      use multiple_geo, only: MFORCE,fcm2,fcum,nforce
      use multiple_states, only: efficiency_prt
      use optci_mod, only: optci_prt
      use optwf_corsam, only: energy,energy_err,force,force_err,sigma
      use precision_kinds, only: dp
      use properties_mod, only: prop_fin
      use qmmm_pot, only: qmmm_extpot_final
      use sa_check, only: energy_all,energy_err_all
      use sa_weights, only: weights
      use step,    only: suc,try
      use system,  only: ncent,nelec
      use tmpnode, only: distance_node_sum
      use vmc_mod, only: delri,nrad
!use contrl, only: nblk, nstep


      implicit none

      integer :: i, iciprt_sav, ifr, index, istate
      integer :: j
      real(dp) :: accept, delr, eerr, eerr1
      real(dp) :: eerr1s, eerr_p, efin, efin_p
      real(dp) :: ferr, ffin
      real(dp) :: passes, peerr, pefin, r2err
      real(dp) :: rtpass, sucsum
      real(dp) :: tcsq, term
      real(dp) :: tpberr, tpbfin, trysum = 0, x
      real(dp) :: x2
      real(dp), dimension(MFORCE) :: ffin_grdnts
      real(dp), dimension(MFORCE) :: ferr_grdnts
      real(dp), parameter :: one = 1.d0
      real(dp), parameter :: half = .5d0

      passes=dble(iblk*vmc_nstep)
      rtpass=dsqrt(passes)

! quantities not computed in acuest_write

      if(iperiodic.eq.0 .and. ncent.eq.1) then
        delr=one/delri
        term=one/(passes*delr)
        trysum=0
        sucsum=0
        do i=1,nrad
          trysum=trysum+try(i)
          sucsum=sucsum+suc(i)
        enddo
      endif

! quantities also computed in acuest_write

      if(index(mode,'mov1').eq.0) then
        accept=acc/passes
       else
        accept=acc/(passes*nelec)
      endif

      do ifr=1,nforce
        energy(ifr)=0
        energy_err(ifr)=0
        ffin_grdnts(ifr)=0
        ferr_grdnts(ifr)=0
        force(ifr)=0
        force_err(ifr)=0
      enddo

      do istate=1,nstates

      eerr1=err1(ecum1(istate),ecm21(istate),istate)
      eerr1s=err1(ecum1s(istate),ecm21s(istate),istate)

      efin=ecum(istate,1)/wcum(istate,1)
      pefin=pecum(istate)/wcum(istate,1)
      tpbfin=tpbcum(istate)/wcum(istate,1)

      eerr=err(ecum(istate,1),ecm2(istate,1),istate,1)
      peerr=err(pecum(istate),pecm2(istate),istate,1)
      tpberr=err(tpbcum(istate),tpbcm2(istate),istate,1)

      energy(1)=energy(1)+weights(istate)*efin

!     save the enegies (of all the states) of the last run for the check in lin_d and error
      energy_all(istate)=efin
      energy_err_all(istate)=eerr

! TMP
!     energy_err(1)=energy_err(1)+(weights(istate)*eerr)**2
      energy_err(1)=energy_err(1)+weights(istate)*eerr

! eerr1*rtpass differs from sigma in that eerr1 contains p*new+q*old,
! so eerr1 is a bit smaller than eerr1s. sigma is a property of the wave
! function only, whereas eerr1*rtpass depends on how quickly one evolves
! the system.
! In the calculation of T_corr, if one uses T_corr=(eerr/eerr1)^2, then
! T_corr=1 when nstep=1, whereas if one uses T_corr=(eerr/eerr1s)^2, then
! T_corr will be a bit < 1 when nstep=1. However, it makes sense to use
! the latter definition because p*new+q*old does reduce T_corr and that
! is precisely what is being reflected when we get T_corr < 1.
      if (eerr1s == 0.) then
        tcsq = 0.
      else
        tcsq=eerr/eerr1s
      endif
      sigma=eerr1s*rtpass

      if(istate.eq.1) then
        write(ounit,'(a12,2x,a20)') mode,title
        write(ounit,'(''results after '',f12.0,'' passes,  nstep, nblk ='',i6,2x,i6)') &
        passes, vmc_nstep,iblk
      endif
      if(nstates.gt.1) write(ounit,'(/,''State '',i4)') istate
      write(ounit,'(''physical variable'',t20,''average'',t34,''rms error'',t47,''rms er*rt(pass)'',t65,''sigma'',t74,''Tcor'')')

      write(ounit,'(''total E ='',t17,f12.7,'' +-'',f11.7,3f9.5,f8.2)') &
       efin,eerr,eerr*rtpass,eerr1*rtpass,sigma,tcsq*tcsq

      efin_p=efin
      eerr_p=eerr

      do ifr=2,nforce
        efin=ecum(istate,ifr)/wcum(istate,ifr)
        eerr=err(ecum(istate,ifr),ecm2(istate,ifr),istate,ifr)
        ffin=ecum(istate,1)/wcum(istate,1)-efin
        ferr=err(fcum(istate,ifr),fcm2(istate,ifr),istate,1)

! save energy, force, and, energy and force error for optimization
        energy(ifr)=energy(ifr)+weights(istate)*efin
! TMP
!       energy_err(ifr)=energy_err(ifr)+(weights(istate)*eerr)**2
        energy_err(ifr)=energy_err(ifr)+weights(istate)*eerr

        force(ifr)=force(ifr)+weights(istate)*ffin
! TMP
!       force_err(ifr)=force_err(ifr)+(weights(istate)*ferr)**2
        force_err(ifr)=force_err(ifr)+weights(istate)*ferr

! save forces and forces errors for calculations of energy gradients.
! Done by Omar Valsson 2008-12-01
        if(ngradnts.gt.0) then
          ffin_grdnts(ifr-1)=ffin
          ferr_grdnts(ifr-1)=ferr
        endif

        write(ounit,'(''total E ='',t17,f12.7,'' +-'',f11.7,f9.5)') efin,eerr,eerr*rtpass
        write(ounit,'(''force   ='',t17,e19.10,'' +-'',e16.8,f9.5)') ffin,ferr,ferr*rtpass
      enddo
      write(ounit,'(''potential E ='',t17,f12.7,'' +-'',f11.7,f9.5)') pefin,peerr,peerr*rtpass
      write(ounit,'(''pb kinetic E ='',t17,f12.7,'' +-'',f11.7,f9.5)') tpbfin,tpberr,tpberr*rtpass

      enddo

! TMP
!     do 250 ifr=1,nforce
!       energy_err(ifr)=sqrt(energy_err(ifr))
! 250   force_err(ifr)=sqrt(force_err(ifr))


      if(index(mode,'mov1').ne.0.and.iperiodic.eq.0.and.ncent.eq.1) then
        write(ounit,'(''acceptance ='',t17,2f12.7)') accept,sucsum/trysum
       else
        write(ounit,'(''acceptance ='',t17,2f12.7)') accept
      endif

      if(iqmmm.gt.0) call qmmm_extpot_final(nelec)

      iciprt_sav=iciprt
      iciprt=-1
      call optci_prt(wcum(1,1),iblk,ounit)
      iciprt=iciprt_sav

      call efficiency_prt(passes)

      call prop_fin(wcum(1,1),iblk,efin_p,eerr_p)

      call finwrt_more

      write(ounit,'(''distance from the nodes '',f10.5)') distance_node_sum/passes

      if(ngradnts.gt.0 .and. igrdtype.eq.1) call finwrt_grdnts_cart(ffin_grdnts,ferr_grdnts)
      if(ngradnts.gt.0 .and. igrdtype.eq.2) call finwrt_grdnts_zmat(ffin_grdnts,ferr_grdnts)
      if(ngradnts.gt.0 .and. igrdtype.eq.2) call finwrt_diaghess_zmat(ffin_grdnts,ferr_grdnts)

      return
contains
      elemental pure function err(x,x2,j,i)
        implicit none
        real(dp), intent(in) :: x, x2
        integer,  intent(in) :: j, i
        real(dp)             :: err
        err=dsqrt(abs(x2/wcum(j,i)-(x/wcum(j,i))**2)/iblk)
      end function
      elemental pure function err1(x,x2,j)
        implicit none
        real(dp), intent(in) :: x, x2
        integer,  intent(in) :: j
        real(dp)             :: err1
        err1=dsqrt(dabs(x2/wcum(j,1)-(x/wcum(j,1))**2)/passes)
      end function
      end
end module
