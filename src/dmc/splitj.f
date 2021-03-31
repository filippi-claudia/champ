      subroutine splitj
c Written by Cyrus Umrigar

      use precision_kinds, only: dp
      use dmc_mod, only: MWALK
      use const, only: nelec
      use forcepar, only: nforce
      use age, only: iage
      use config, only: d2o, peo_dmc, psido_dmc, psijo_dmc, vold_dmc, xold_dmc
      use stats, only: nbrnch
      use force_dmc, only: nwprod
      use branch, only: eold, nwalk, pwt, wt
      use branch, only: wthist
      use jacobsave, only: ajacold
      use velratio, only: fratio, xdrifted

      implicit none

      integer :: i, ifr, ip, ipair, iunder
      integer :: iw, iw2, iwundr, j
      integer :: k, nwalk2
      real(dp) :: half, rannyu, two, wtsm, wtsm2
      real(dp) :: wttot, zero

      parameter (zero=0.d0,two=2.d0,half=.5d0)

      dimension iwundr(MWALK)

      iunder=0
      ipair=0
      wtsm=zero
      do 10 iw=1,nwalk
        wtsm=wtsm+wt(iw)
        if(wt(iw).lt.half) then
          if(wt(iw).eq.zero) then
            nbrnch=nbrnch+1
            iunder=iunder+1
            iwundr(iunder)=iw
           else
            if(ipair.eq.0) then
              ipair=1
              iw2=iw
             else
              nbrnch=nbrnch+1
              ipair=0
              iunder=iunder+1
              wttot=wt(iw)+wt(iw2)
              if(rannyu(0).gt.(wt(iw)/wttot)) then
                wt(iw2)=wttot
                iwundr(iunder)=iw
               else
                wt(iw)=wttot
                iwundr(iunder)=iw2
              endif
            endif
          endif
        endif
   10 continue

      nwalk2=nwalk
      do 20 iw=1,nwalk
        if(wt(iw).ge.two) then
          nbrnch=nbrnch+1
          if(iunder.gt.0) then
            iw2=iwundr(iunder)
            iunder=iunder-1
           else
            nwalk2=nwalk2+1
            iw2=nwalk2
            if(nwalk2.gt.MWALK) call fatal_error('SPLITJ: MWALK exceeded')
          endif
          wt(iw)=wt(iw)*half
          wt(iw2)=wt(iw)
          iage(iw2)=iage(iw)
          call splitjdet(iw,iw2)
          call splitjjas(iw,iw2)
c         call t_vpsp_splitj(iw,iw2)
          call prop_splitj(iw,iw2)
          call pcm_splitj(iw,iw2)
          call mmpol_splitj(iw,iw2)
          do 15 ifr=1,nforce
            ajacold(iw2,ifr)=ajacold(iw,ifr)
            eold(iw2,ifr)=eold(iw,ifr)
            psido_dmc(iw2,ifr)=psido_dmc(iw,ifr)
            psijo_dmc(iw2,ifr)=psijo_dmc(iw,ifr)
            peo_dmc(iw2,ifr)=peo_dmc(iw,ifr)
            d2o(iw2,ifr)=d2o(iw,ifr)
            pwt(iw2,ifr)=pwt(iw,ifr)
            fratio(iw2,ifr)=fratio(iw,ifr)
            do 12 ip=0,nwprod-1
   12         wthist(iw2,ip,ifr)=wthist(iw,ip,ifr)
            do 15 i=1,nelec
              do 15 k=1,3
                xdrifted(k,i,iw2,ifr)=xdrifted(k,i,iw,ifr)
                vold_dmc(k,i,iw2,ifr)=vold_dmc(k,i,iw,ifr)
   15           xold_dmc(k,i,iw2,ifr)=xold_dmc(k,i,iw,ifr)
        endif
   20 continue

      do 30 j=iunder,1,-1
        iw2=iwundr(j)
        iw=nwalk2
        nwalk2=nwalk2-1
        wt(iw2)=wt(iw)
        iage(iw2)=iage(iw)
        call splitjdet(iw,iw2)
        call splitjjas(iw,iw2)
c       call t_vpsp_splitj(iw,iw2)
        call prop_splitj(iw,iw2)
        call pcm_splitj(iw,iw2)
        call mmpol_splitj(iw,iw2)
        do 30 ifr=1,nforce
          ajacold(iw2,ifr)=ajacold(iw,ifr)
          eold(iw2,ifr)=eold(iw,ifr)
          psido_dmc(iw2,ifr)=psido_dmc(iw,ifr)
          psijo_dmc(iw2,ifr)=psijo_dmc(iw,ifr)
          peo_dmc(iw2,ifr)=peo_dmc(iw,ifr)
          d2o(iw2,ifr)=d2o(iw,ifr)
          pwt(iw2,ifr)=pwt(iw,ifr)
          fratio(iw2,ifr)=fratio(iw,ifr)
          do 25 ip=0,nwprod-1
   25       wthist(iw2,ip,ifr)=wthist(iw,ip,ifr)
          do 30 i=1,nelec
            do 30 k=1,3
              xdrifted(k,i,iw2,ifr)=xdrifted(k,i,iw,ifr)
              vold_dmc(k,i,iw2,ifr)=vold_dmc(k,i,iw,ifr)
   30         xold_dmc(k,i,iw2,ifr)=xold_dmc(k,i,iw,ifr)
      nwalk=nwalk2

      wtsm2=zero
      do 40 iw=1,nwalk
        wtsm2=wtsm2+wt(iw)
c       if(wt(iw).lt.half) write(11,'(i4,9d12.5)') iw,wt(iw),eold(iw)
c       if(wt(iw).gt.two) write(11,'(i4,9d12.5)') iw,wt(iw),eold(iw)
   40 continue
c     if(dabs(wtsm-wtsm2).gt.1.d-10) write(11,'(2f12.6)') wtsm,wtsm2

      return
      end
