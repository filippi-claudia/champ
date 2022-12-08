      module splitj_mod
      contains
      subroutine splitj
c Written by Cyrus Umrigar

      use dmc_mod, only: MWALK
      use system, only: nelec
      use multiple_geo, only: nforce, nwprod
      use age, only: iage
      use config, only: d2o, peo_dmc, psido_dmc, psijo_dmc, vold_dmc, xold_dmc
      use stats, only: nbrnch
      use branch, only: eold, nwalk, pwt, wt
      use branch, only: wthist
      use jacobsave, only: ajacold
      use velratio, only: fratio, xdrifted
      use precision_kinds, only: dp

      use mmpol_dmc,      only: mmpol_splitj
      use pcm_dmc,        only: pcm_splitj
      use prop_dmc,       only: prop_splitj
      use walksav_jas_mod,only: splitjjas
      use walksav_det_mod,only: splitjdet
      use error,          only: fatal_error
      use random_mod,     only: random_dp

      implicit none

      integer :: i, ifr, ip, ipair, iunder
      integer :: iw, iw2, j, k
      integer :: nwalk2
      integer, dimension(MWALK) :: iwundr
      real(dp) :: wtsm, wtsm2, wttot
      real(dp), parameter :: zero = 0.d0
      real(dp), parameter :: two = 2.d0
      real(dp), parameter :: half = .5d0



      iunder=0
      ipair=0
      wtsm=zero
      do iw=1,nwalk
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
              if(random_dp().gt.(wt(iw)/wttot)) then
                wt(iw2)=wttot
                iwundr(iunder)=iw
               else
                wt(iw)=wttot
                iwundr(iunder)=iw2
              endif
            endif
          endif
        endif
      enddo

      nwalk2=nwalk
      do iw=1,nwalk
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
          do ifr=1,nforce
            ajacold(iw2,ifr)=ajacold(iw,ifr)
            eold(iw2,ifr)=eold(iw,ifr)
            psido_dmc(iw2,ifr)=psido_dmc(iw,ifr)
            psijo_dmc(iw2,ifr)=psijo_dmc(iw,ifr)
            peo_dmc(iw2,ifr)=peo_dmc(iw,ifr)
            d2o(iw2,ifr)=d2o(iw,ifr)
            pwt(iw2,ifr)=pwt(iw,ifr)
            fratio(iw2,ifr)=fratio(iw,ifr)
            do ip=0,nwprod-1
              wthist(iw2,ip,ifr)=wthist(iw,ip,ifr)
            enddo
            do i=1,nelec
              do k=1,3
                xdrifted(k,i,iw2,ifr)=xdrifted(k,i,iw,ifr)
                vold_dmc(k,i,iw2,ifr)=vold_dmc(k,i,iw,ifr)
                xold_dmc(k,i,iw2,ifr)=xold_dmc(k,i,iw,ifr)
              enddo
            enddo
          enddo
        endif
      enddo

      do j=iunder,1,-1
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
        do ifr=1,nforce
          ajacold(iw2,ifr)=ajacold(iw,ifr)
          eold(iw2,ifr)=eold(iw,ifr)
          psido_dmc(iw2,ifr)=psido_dmc(iw,ifr)
          psijo_dmc(iw2,ifr)=psijo_dmc(iw,ifr)
          peo_dmc(iw2,ifr)=peo_dmc(iw,ifr)
          d2o(iw2,ifr)=d2o(iw,ifr)
          pwt(iw2,ifr)=pwt(iw,ifr)
          fratio(iw2,ifr)=fratio(iw,ifr)
          do ip=0,nwprod-1
            wthist(iw2,ip,ifr)=wthist(iw,ip,ifr)
          enddo
          do i=1,nelec
            do k=1,3
              xdrifted(k,i,iw2,ifr)=xdrifted(k,i,iw,ifr)
              vold_dmc(k,i,iw2,ifr)=vold_dmc(k,i,iw,ifr)
              xold_dmc(k,i,iw2,ifr)=xold_dmc(k,i,iw,ifr)
            enddo
          enddo
        enddo
      enddo
      nwalk=nwalk2

      wtsm2=zero
      do iw=1,nwalk
        wtsm2=wtsm2+wt(iw)
c       if(wt(iw).lt.half) write(11,'(i4,9d12.5)') iw,wt(iw),eold(iw)
c       if(wt(iw).gt.two) write(11,'(i4,9d12.5)') iw,wt(iw),eold(iw)
      enddo
c     if(dabs(wtsm-wtsm2).gt.1.d-10) write(11,'(2f12.6)') wtsm,wtsm2

      return
      end
      end module
