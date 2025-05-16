module splitj_mod
contains
      subroutine splitj
! Written by Cyrus Umrigar

      use age,     only: iage
      use branch,  only: eold,nwalk,pwt,wt,wthist
      use config,  only: d2o,peo_dmc,psido_dmc,psijo_dmc,vold_dmc
      use config,  only: xold_dmc
      use contrldmc, only: icut_e
      use dmc_mod, only: MWALK
      use error,   only: fatal_error
      use force_pth, only: PTH
      use fragments,  only: eloco_i, fratio_i, elocofrag, fratiofrag, nfrag
      use jacobsave, only: ajacold
      use mmpol_dmc, only: mmpol_splitj
      use multiple_geo, only: nforce,nwprod
      use m_force_analytic, only: iforce_analy
      use pcm_dmc, only: pcm_splitj
      use precision_kinds, only: dp
      use prop_dmc, only: prop_splitj
      use random_mod, only: random_dp
      use stats,   only: nbrnch
      use system,  only: nelec, ncent
      use vd_mod, only: esnake, ehist, deriv_eold, dmc_ivd
      use velratio, only: fratio,xdrifted
      use walksav_det_mod, only: splitjdet
      use walksav_jas_mod, only: splitjjas
      use pathak_mod, only: pold


      implicit none

      integer :: i, ifr, ip, ipair, iunder, ic, iph
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
!         call t_vpsp_splitj(iw,iw2)
          call prop_splitj(iw,iw2)
          call pcm_splitj(iw,iw2)
          call mmpol_splitj(iw,iw2)
          if(iforce_analy.eq.1) then
            if(dmc_ivd.gt.0) then
              do iph=1,PTH
                pold(iw2,iph)=pold(iw,iph)
                do ic=1,ncent
                  do k=1,3
                    esnake(k,ic,iw2,iph)=esnake(k,ic,iw,iph)
                    do ip=0,nwprod-1
                      ehist(k,ic,iw2,ip,iph)=ehist(k,ic,iw,ip,iph)
                    enddo
                  enddo
                enddo
              enddo
              do ic=1,ncent
                do k=1,3
                  deriv_eold(k,ic,iw2)=deriv_eold(k,ic,iw)
                enddo
              enddo
            else                
              do iph=1,PTH
                pold(iw2,iph)=pold(iw,iph)
              enddo
            endif
          endif
          do ifr=1,nforce
            ajacold(iw2,ifr)=ajacold(iw,ifr)
            eold(iw2,ifr)=eold(iw,ifr)
            psido_dmc(iw2,ifr)=psido_dmc(iw,ifr)
            psijo_dmc(iw2,ifr)=psijo_dmc(iw,ifr)
            peo_dmc(iw2,ifr)=peo_dmc(iw,ifr)
            d2o(iw2,ifr)=d2o(iw,ifr)
            pwt(iw2,ifr)=pwt(iw,ifr)
            fratio(iw2,ifr)=fratio(iw,ifr)
            
            if (nfrag.gt.1) then
              eloco_i(:,iw2,ifr)=eloco_i(:,iw,ifr)
              fratio_i(:,iw2,ifr)=fratio_i(:,iw,ifr)
            endif
            if (icut_e.lt.0) then
              elocofrag(:,iw2,ifr)=elocofrag(:,iw,ifr)
              fratiofrag(:,iw2,ifr)=fratiofrag(:,iw,ifr)
            endif

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
!       call t_vpsp_splitj(iw,iw2)
        call prop_splitj(iw,iw2)
        call pcm_splitj(iw,iw2)
        call mmpol_splitj(iw,iw2)

        if(iforce_analy.eq.1) then
          if(dmc_ivd.gt.0) then
            do iph=1,PTH
              do ic=1,ncent
                do k=1,3
                  esnake(k,ic,iw2,iph)=esnake(k,ic,iw,iph)
                  do ip=0,nwprod-1
                    ehist(k,ic,iw2,ip,iph)=ehist(k,ic,iw,ip,iph)
                  enddo
                enddo
             enddo
            enddo
            do ic=1,ncent
              do k=1,3
                deriv_eold(k,ic,iw2)=deriv_eold(k,ic,iw)
              enddo
            enddo
          else
            do iph=1,PTH
              pold(iw2,iph)=pold(iw,iph)
            enddo
          endif
        endif
        do ifr=1,nforce
          ajacold(iw2,ifr)=ajacold(iw,ifr)
          eold(iw2,ifr)=eold(iw,ifr)
          psido_dmc(iw2,ifr)=psido_dmc(iw,ifr)
          psijo_dmc(iw2,ifr)=psijo_dmc(iw,ifr)
          peo_dmc(iw2,ifr)=peo_dmc(iw,ifr)
          d2o(iw2,ifr)=d2o(iw,ifr)
          pwt(iw2,ifr)=pwt(iw,ifr)
          fratio(iw2,ifr)=fratio(iw,ifr)
          
          if (icut_e.lt.0) then
            eloco_i(:,iw2,ifr)=eloco_i(:,iw,ifr)
            fratio_i(:,iw2,ifr)=fratio_i(:,iw,ifr)
          endif
          if (nfrag.gt.1) then
            elocofrag(:,iw2,ifr)=elocofrag(:,iw,ifr)
            fratiofrag(:,iw2,ifr)=fratiofrag(:,iw,ifr)
          endif

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
!       if(wt(iw).lt.half) write(11,'(i4,9d12.5)') iw,wt(iw),eold(iw)
!       if(wt(iw).gt.two) write(11,'(i4,9d12.5)') iw,wt(iw),eold(iw)
      enddo
!     if(dabs(wtsm-wtsm2).gt.1.d-10) write(11,'(2f12.6)') wtsm,wtsm2

      return
      end
end module
