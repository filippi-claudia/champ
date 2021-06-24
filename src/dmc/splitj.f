      subroutine splitj
c Written by Cyrus Umrigar
      implicit real*8(a-h,o-z)
      include 'vmc.h'
      include 'dmc.h'
      include 'force.h'
      common /forcepar/ deltot(MFORCE),nforce,istrech
      common /forcest/ fgcum(MFORCE),fgcm2(MFORCE)
      common /force_dmc/ itausec,nwprod
      parameter (zero=0.d0,two=2.d0,half=.5d0)

      common /const/ pi,hb,etrial,delta,deltai,fbias,nelec,imetro,ipr
      common /config/ xold(3,MELEC,MWALK,MFORCE),vold(3,MELEC,MWALK,MFORCE),
     &psido(MWALK,MFORCE),psijo(MWALK,MFORCE),peo(MWALK,MFORCE),d2o(MWALK,MFORCE)
      common /velratio/ fratio(MWALK,MFORCE),xdrifted(3,MELEC,MWALK,MFORCE)
      common /age/ iage(MWALK),ioldest,ioldestmx
      common /stats/ dfus2ac,dfus2un,dr2ac,dr2un,acc,trymove,nacc,
     &nbrnch,nodecr
      common /branch/ wtgen(0:MFPRD1),ff(0:MFPRD1),eold(MWALK,MFORCE),
     &pwt(MWALK,MFORCE),wthist(MWALK,0:MFORCE_WT_PRD,MFORCE),
     &wt(MWALK),eigv,eest,wdsumo,wgdsumo,fprod,nwalk

      common /atom/ znuc(MCTYPE),cent(3,MCENT),pecent
     &,iwctype(MCENT),nctype,ncent
      common /derivanaly/ deriv_energy_sum(10,3,MCENT,PTH),deriv_energy_cum(10,3,MCENT,PTH),
     &energy_snake(3,MCENT,MWALK,PTH),energy_hist(3,MCENT,MWALK,0:MFORCE_WT_PRD,PTH),
     &deriv_energy_old(3,MCENT,MWALK),pathak_old(MWALK,PTH),eps_pathak(PTH),ipathak
      common /force_analy/ iforce_analy

      common /jacobsave/ ajacob,ajacold(MWALK,MFORCE)

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
          if(iforce_analy.eq.1) then
            if(ipathak.gt.0) then
              do 101 iph=1,ipathak
                pathak_old(iw2,iph)=pathak_old(iw,iph)
                do 101 ic=1,ncent
                  do 101 k=1,3
                    energy_snake(k,ic,iw2,iph)=energy_snake(k,ic,iw,iph)
                    do 101 ip=0,nwprod-1
  101                 energy_hist(k,ic,iw2,ip,iph)=energy_hist(k,ic,iw,ip,iph)
            else
              do 102 ic=1,ncent
                do 102 k=1,3
                  energy_snake(k,ic,iw2,1)=energy_snake(k,ic,iw,1)
                  do 102 ip=0,nwprod-1
  102               energy_hist(k,ic,iw2,ip,1)=energy_hist(k,ic,iw,ip,1)
            endif
            do 103 ic=1,ncent
              do 103 k=1,3
  103           deriv_energy_old(k,ic,iw2)=deriv_energy_old(k,ic,iw)
          endif
          do 15 ifr=1,nforce
            ajacold(iw2,ifr)=ajacold(iw,ifr)
            eold(iw2,ifr)=eold(iw,ifr)
            psido(iw2,ifr)=psido(iw,ifr)
            psijo(iw2,ifr)=psijo(iw,ifr)
            peo(iw2,ifr)=peo(iw,ifr)
            d2o(iw2,ifr)=d2o(iw,ifr)
            pwt(iw2,ifr)=pwt(iw,ifr)
            fratio(iw2,ifr)=fratio(iw,ifr)
            do 12 ip=0,nwprod-1
   12         wthist(iw2,ip,ifr)=wthist(iw,ip,ifr)
            do 15 i=1,nelec
              do 15 k=1,3
                xdrifted(k,i,iw2,ifr)=xdrifted(k,i,iw,ifr)
                vold(k,i,iw2,ifr)=vold(k,i,iw,ifr)
   15           xold(k,i,iw2,ifr)=xold(k,i,iw,ifr)
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
        if(iforce_analy.eq.1) then
          if(ipathak.gt.0) then
            do 201 iph=1,ipathak
              pathak_old(iw2,iph)=pathak_old(iw,iph)
              do 201 ic=1,ncent
                do 201 k=1,3
                  energy_snake(k,ic,iw2,iph)=energy_snake(k,ic,iw,iph)
                  do 201 ip=0,nwprod-1
201                 energy_hist(k,ic,iw2,ip,iph)=energy_hist(k,ic,iw,ip,iph)
          else
            do 202 ic=1,ncent
              do 202 k=1,3
                energy_snake(k,ic,iw2,1)=energy_snake(k,ic,iw,1)
                do 202 ip=0,nwprod-1
202               energy_hist(k,ic,iw2,ip,1)=energy_hist(k,ic,iw,ip,1)
          endif
          do 203 ic=1,ncent
            do 203 k=1,3
203           deriv_energy_old(k,ic,iw2)=deriv_energy_old(k,ic,iw)
        endif
        do 30 ifr=1,nforce
          ajacold(iw2,ifr)=ajacold(iw,ifr)
          eold(iw2,ifr)=eold(iw,ifr)
          psido(iw2,ifr)=psido(iw,ifr)
          psijo(iw2,ifr)=psijo(iw,ifr)
          peo(iw2,ifr)=peo(iw,ifr)
          d2o(iw2,ifr)=d2o(iw,ifr)
          pwt(iw2,ifr)=pwt(iw,ifr)
          fratio(iw2,ifr)=fratio(iw,ifr)
          do 25 ip=0,nwprod-1
   25       wthist(iw2,ip,ifr)=wthist(iw,ip,ifr)
          do 30 i=1,nelec
            do 30 k=1,3
              xdrifted(k,i,iw2,ifr)=xdrifted(k,i,iw,ifr)
              vold(k,i,iw2,ifr)=vold(k,i,iw,ifr)
   30         xold(k,i,iw2,ifr)=xold(k,i,iw,ifr)
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
