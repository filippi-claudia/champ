      module acuest_write_mod
      contains
      subroutine acuest_write(enow, nproc)
c Written by Claudia Filippi
c routine to write out estimators for energy etc.

      use contrl_file,    only: ounit, errunit
      use control, only: mode
      use control_vmc, only:  vmc_nstep
      use csfs, only: nstates
      use est2cm, only: ecm2, pecm2, tjfcm2, tpbcm2
      use estcum, only: ecum, iblk, pecum, tjfcum, tpbcum
      use estsum, only: acc
      use forcewt, only: wcum
      use mstates_mod, only: MSTATES
      use multiple_geo, only: MFORCE, nforce, fcm2, fcum
      use optci_mod, only: optci_prt
      use pcm_vmc, only: pcm_prt
      use precision_kinds, only: dp, i2b
      use properties_mod, only: prop_prt
      use system, only: nelec
      !use contrl, only: nstep
      implicit none

      integer :: i, ieerr, iferr, ifr, index
      integer :: ipeerr, istate, itpber
      integer :: j, nproc
      real(dp) :: acc_denom, accept, eave, eerr, err
      real(dp) :: fave, ferr, peave, peerr
      real(dp) :: tpbave, tpberr
      real(dp) :: x, x2
      real(dp), dimension(MSTATES, MFORCE) :: enow


c statement function for error calculation
      err(x,x2,j,i)=dsqrt(abs(x2/wcum(j,i)-(x/wcum(j,i))**2)/iblk)

c xave = current average value of x
c xerr = current error of x

c write out header first time

      if (iblk.eq.nproc) then
        if(nforce.gt.1) then
          write(ounit,'(t5,''enow'',t15,''eave'',t21,''(eerr )''
     &    ,t32,''peave'',t38,''(peerr)'',t48,''tpbave'',t54,''(tpberr)''
     &    ,t66,''fave'',t80,''(ferr)'',t93,''accept'',t101,''iter'')')
         else
          write(ounit,'(t5,''enow'',t15,''eave'',t21,''(eerr )''
     &    ,t32,''peave'',t38,''(peerr)'',t48,''tpbave'',t54,''(tpberr)''
     &    ,t67,''accept'',t79,''iter'')')
        endif
      endif

c write out current values of averages
      acc_denom=dfloat(vmc_nstep*iblk)
      if(index(mode,'mov1').eq.0) then
        accept=acc/acc_denom
       else
        accept=acc/(acc_denom*nelec)
      endif

      do ifr=1,nforce
        do istate=1,nstates
          eave=ecum(istate,ifr)/wcum(istate,ifr)
          if(iblk.eq.1) then
            eerr=0
          else
            eerr=err(ecum(istate,ifr),ecm2(istate,ifr),istate,ifr)
          endif

          ieerr=nint(100000*eerr)
          if(ifr.eq.1) then
            if(iblk.eq.1) then
              peerr=0.
              tpberr=0.
             else
              peerr=err(pecum(istate),pecm2(istate),istate,ifr)
              tpberr=err(tpbcum(istate),tpbcm2(istate),istate,ifr)
            endif
            peave=pecum(istate)/wcum(istate,ifr)
            tpbave=tpbcum(istate)/wcum(istate,ifr)

            ipeerr=nint(100000*peerr)
            itpber=nint(100000*tpberr)

            if(istate.eq.1) then
c           with single-state, fine, with multi gives write out errors, same for the else
c              write(ounit,'(f10.5,4(f10.5,''('',i5,'')''),25x,f10.5,i10)')
c     &        enow(1,1),eave,ieerr,peave,ipeerr,tpbave,itpber,tjfave,itjfer,accept,iblk*vmc_nstep

              if(nforce.gt.1) then
                write(ounit,'(f10.5,3(f10.5,''('',i5,'')''),25x,f10.5,i10)')
     &          enow(1,1),eave,ieerr,peave,ipeerr,tpbave,itpber,accept,iblk*vmc_nstep
              else
                write(ounit,'(f10.5,3(f10.5,''('',i5,'')''),1x,f10.5,i10)')
     &          enow(1,1),eave,ieerr,peave,ipeerr,tpbave,itpber,accept,iblk*vmc_nstep
c                write(ounit,'(f10.5,3(f10.5,a,i0,a),1x,f10.5,i10)')
c     &          enow(1,1),eave,"(", ieerr, ")", peave, "(", ipeerr , 
c     &          ")", tpbave, "(", itpber, ")", accept,iblk*vmc_nstep
              endif

              call prop_prt(wcum(1,ifr),iblk,ounit)
              call optci_prt(wcum(1,ifr),iblk,ounit)
c             call optorb_prt(wcum(1,ifr),eave,6)
c different meaning of last argument: 0 acuest, 1 finwrt
              call pcm_prt(wcum(1,ifr),iblk)

            else
              write(ounit,'(f10.5,3(f10.5,''('',i5,'')''))')
     &        enow(istate,1),eave,ieerr,peave,ipeerr,tpbave,itpber
            endif

          else
           fave=(ecum(istate,1)/wcum(istate,1)-ecum(istate,ifr)/wcum(istate,ifr))
     &    !/abs(deltot(ifr))
            ferr=err(fcum(istate,ifr),fcm2(istate,ifr),istate,1)!/abs(deltot(ifr))
            iferr=nint(1.0d9*ferr)
            write(ounit,'(f10.5,f10.5,''('',i5,'')'',51x,f14.9,''('',i9,'')'')
     &      ') enow(istate,ifr),eave,ieerr,fave,iferr
          endif
c endif for ifr.eq.1
        enddo
c enddo over istate, nstates
      enddo
c enddo over ifr, nforce

      return
      end
      end module
