      subroutine init
c MPI version created by Claudia Filippi starting from serial version
c routine to accumulate estimators for energy etc.

      implicit real*8(a-h,o-z)
      include 'vmc.h'
      include 'dmc.h'
      include 'force.h'
      include 'basis.h'
      include 'pseudo.h'
      include 'mpif.h'
      parameter (zero=0.d0,one=1.d0)

      common /const/ pi,hb,etrial,delta,deltai,fbias,nelec,imetro,ipr
      common /contrl/ nstep,nblk,nblkeq,nconf,nconf_new,isite,idump,irstar
      common /iterat/ ipass,iblk
      common /atom/ znuc(MCTYPE),cent(3,MCENT),pecent
     &,iwctype(MCENT),nctype,ncent
      common /pseudo/ vps(MELEC,MCENT,MPS_L),vpso(MELEC,MCENT,MPS_L,MFORCE)
     &,lpot(MCTYPE),nloc
      common /qua/ xq0(MPS_QUAD),yq0(MPS_QUAD),zq0(MPS_QUAD)
     &,xq(MPS_QUAD),yq(MPS_QUAD),zq(MPS_QUAD),wq(MPS_QUAD),nquad
      common /casula/ t_vpsp(MCENT,MPS_QUAD,MELEC),icasula,i_vpsp
      common /config/ xold(3,MELEC,MWALK,MFORCE),vold(3,MELEC,MWALK,MFORCE),
     &psido(MWALK,MFORCE),psijo(MWALK,MFORCE),peo(MWALK,MFORCE),d2o(MWALK,MFORCE)
      common /branch/ wtgen(0:MFPRD1),ff(0:MFPRD1),eold(MWALK,MFORCE),
     &pwt(MWALK,MFORCE),wthist(MWALK,0:MFORCE_WT_PRD,MFORCE),
     &wt(MWALK),eigv,eest,wdsumo,wgdsumo,fprod,nwalk
      common /force_dmc/ itausec,nwprod
      common /forcepar/ deltot(MFORCE),nforce,istrech
      common /jacobsave/ ajacob,ajacold(MWALK,MFORCE)

      common /da_energy_now/ da_energy(3,MCENT),da_psi(3,MCENT)
      common /derivanaly/ deriv_esum(10,3,MCENT,PTH),deriv_ecum(10,3,MCENT,PTH),
     &esnake(3,MCENT,MWALK,PTH),ehist(3,MCENT,MWALK,0:MFORCE_WT_PRD,PTH),
     &deriv_eold(3,MCENT,MWALK),pold(MWALK,PTH),deriv_cm(3,MCENT,PTH),
     &deriv_cm2(3,MCENT,PTH),eps_pathak(PTH),ipathak
      common /force_analy/ iforce_analy

      character*12 mode
      common /contr3/ mode

      logical wid
      common /mpiconf/ idtask,nproc,wid

c Initialize various quantities at beginning of run
c the initial values of energy psi etc. are calculated here

      ipass=0

c set quadrature points
      if(nloc.gt.0) call gesqua (nquad,xq,yq,zq,wq)

c get nuclear potential energy
      call pot_nn(cent,znuc,iwctype,ncent,pecent)

      eigv=one
      eest=etrial
      nwalk=nconf
      fprod=one

      do 80 iw=1,nconf
        wt(iw)=one
        if(istrech.eq.0) then
          do 71 ifr=2,nforce
            do 71 ie=1,nelec
              do 71 k=1,3
   71           xold(k,ie,iw,ifr)=xold(k,ie,iw,1)
        endif
        do 72 ifr=1,nforce
          if(nforce.gt.1) then
            if(ifr.eq.1.or.istrech.eq.0) then
              call strech(xold(1,1,iw,1),xold(1,1,iw,ifr),ajacob,ifr,0)
               else
              call strech(xold(1,1,iw,1),xold(1,1,iw,ifr),ajacob,ifr,1)
            endif
           else
            ajacob=one
          endif
          ajacold(iw,ifr)=ajacob
          if(icasula.lt.0) i_vpsp=icasula
          call hpsi(xold(1,1,iw,ifr),psido(iw,ifr),psijo(iw,ifr),eold(iw,ifr),0,ifr)
          i_vpsp=0
          do 73 i=1,nelec
   73       call compute_determinante_grad(i,psido(iw,ifr),psido(iw,ifr),vold(1,i,iw,ifr),1)

          if(ifr.eq.1) then
            call walksav_det(iw)
            call walksav_jas(iw)
c           call t_vpsp_sav(iw)
            call t_vpsp_sav
            call prop_save_dmc(iw)
            call pcm_save(iw)
            call mmpol_save(iw)
            if(iforce_analy.eq.1) then
              if(ipathak.gt.0) then
                call nodes_distance(vold(1,1,iw,ifr),distance_node,1)
                do 74 iph=1,ipathak
                  call pathak(distance_node,pold(iw,iph),eps_pathak(iph))
                  do 74 ic=1,ncent
                    do 74 k=1,3
                      esnake(k,ic,iw,iph)=zero
                      do 74 ip=0,nwprod-1
   74                   ehist(k,ic,iw,ip,iph)=zero
              else
                do 76 ic=1,ncent
                  do 76 k=1,3
                    esnake(k,ic,iw,1)=zero
                    do 76 ip=0,nwprod-1
   76                 ehist(k,ic,iw,ip,1)=zero
              endif
              do 78 ic=1,ncent
                do 78 k=1,3
   78             deriv_eold(k,ic,iw)=da_energy(k,ic)
            endif
          endif
          pwt(iw,ifr)=0
          do 72 ip=0,nwprod-1
   72       wthist(iw,ip,ifr)=0
   80 continue

      if(mode.eq.'dmc_one_mpi2') nconf=nconf*nproc
      wdsumo=nconf
      wgdsumo=nconf
      do 70 i=0,MFPRD1
        wtgen(i)=nconf
   70   ff(i)=one

      call zerest

      return
      end
