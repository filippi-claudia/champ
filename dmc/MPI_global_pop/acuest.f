      subroutine acuest
c MPI version created by Claudia Filippi starting from serial version
c routine to accumulate estimators for energy etc.

      implicit real*8(a-h,o-z)
      include 'vmc.h'
      include 'dmc.h'
      include 'force.h'
      include 'pseudo.h'
      include 'basis.h'
      include 'mpif.h'
      parameter (zero=0.d0,one=1.d0)

      common /forcepar/ deltot(MFORCE),nforce,istrech
      common /forcest/ fgcum(MFORCE),fgcm2(MFORCE)
      common /force_dmc/ itausec,nwprod

      common /const/ pi,hb,etrial,delta,deltai,fbias,nelec,imetro,ipr
      common /contrl/ nstep,nblk,nblkeq,nconf,nconf_new,isite,idump,irstar
      common /contrldmc/ tau,rttau,taueff(MFORCE),tautot,nfprod,idmc,ipq
     &,itau_eff,iacc_rej,icross,icuspg,idiv_v,icut_br,icut_e
      common /iterat/ ipass,iblk
      common /atom/ znuc(MCTYPE),cent(3,MCENT),pecent
     &,iwctype(MCENT),nctype,ncent
      common /pseudo/ vps(MELEC,MCENT,MPS_L),vpso(MELEC,MCENT,MPS_L,MFORCE)
     &,lpot(MCTYPE),nloc
      common /qua/ xq0(MPS_QUAD),yq0(MPS_QUAD),zq0(MPS_QUAD)
     &,xq(MPS_QUAD),yq(MPS_QUAD),zq(MPS_QUAD),wq(MPS_QUAD),nquad
      common /config/ xold(3,MELEC,MWALK,MFORCE),vold(3,MELEC,MWALK,MFORCE),
     &psido(MWALK,MFORCE),psijo(MWALK,MFORCE),peo(MWALK,MFORCE),d2o(MWALK,MFORCE)
      common /stats/ dfus2ac,dfus2un,dr2ac,dr2un,acc,trymove,nacc,
     &nbrnch,nodecr
      common /estsum/ wsum,w_acc_sum,wfsum,wgsum(MFORCE),wg_acc_sum,wdsum,
     &wgdsum, wsum1(MFORCE),w_acc_sum1,wfsum1,wgsum1(MFORCE),wg_acc_sum1,
     &wdsum1, esum,efsum,egsum(MFORCE),esum1(MFORCE),efsum1,egsum1(MFORCE),
     &ei1sum,ei2sum,ei3sum, pesum(MFORCE),tpbsum(MFORCE),tjfsum(MFORCE),r2sum,
     &risum,tausum(MFORCE)
      common /estcum/ wcum,w_acc_cum,wfcum,wgcum(MFORCE),wg_acc_cum,wdcum,
     &wgdcum, wcum1,w_acc_cum1,wfcum1,wgcum1(MFORCE),wg_acc_cum1,
     &wdcum1, ecum,efcum,egcum(MFORCE),ecum1,efcum1,egcum1(MFORCE),
     &ei1cum,ei2cum,ei3cum, pecum(MFORCE),tpbcum(MFORCE),tjfcum(MFORCE),r2cum,
     &ricum,taucum(MFORCE)
      common /estcm2/ wcm2,wfcm2,wgcm2(MFORCE),wdcm2,wgdcm2, wcm21,
     &wfcm21,wgcm21(MFORCE),wdcm21, ecm2,efcm2,egcm2(MFORCE), ecm21,
     &efcm21,egcm21(MFORCE),ei1cm2,ei2cm2,ei3cm2, pecm2(MFORCE),tpbcm2(MFORCE),
     &tjfcm2(MFORCE),r2cm2,ricm2
      common /derivest/ derivsum(10,MFORCE),derivcum(10,MFORCE)
      common /step/try(nrad),suc(nrad),trunfb(nrad),rprob(nrad),
     &ekin(nrad),ekin2(nrad)
      common /denupdn/ rprobup(nrad),rprobdn(nrad)
      common /branch/ wtgen(0:MFPRD1),ff(0:MFPRD1),eold(MWALK,MFORCE),
     &pwt(MWALK,MFORCE),wthist(MWALK,0:MFORCE_WT_PRD,MFORCE),
     &wt(MWALK),eigv,eest,wdsumo,wgdsumo,fprod,nwalk
      common /age/ iage(MWALK),ioldest,ioldestmx
      common /jacobsave/ ajacob,ajacold(MWALK,MFORCE)

      common /casula/ t_vpsp(MCENT,MPS_QUAD,MELEC),icasula,i_vpsp

      logical wid
      common /mpiconf/ idtask,nproc,wid

      dimension egcollect(MFORCE),wgcollect(MFORCE),pecollect(MFORCE),
     +tpbcollect(MFORCE),tjfcollect(MFORCE),taucollect(MFORCE),
     +derivcollect(10,MFORCE)

      dimension ddx_ref(3)

c statement function for error calculation
      rn_eff(w,w2)=w**2/w2
      error(x,x2,w,w2)=dsqrt(max((x2/w-(x/w)**2)/(rn_eff(w,w2)-1),0.d0))
      errg(x,x2,i)=error(x,x2,wgcum(i),wgcm2(i))

c wt   = weight of configurations
c xsum = sum of values of x from dmc
c xnow = average of values of x from dmc
c xcum = accumulated sums of xnow
c xcm2 = accumulated sums of xnow**2
c xave = current average value of x
c xerr = current error of x

      iblk=iblk+1

      npass=iblk*nstep

      call mpi_reduce(pesum,pecollect,MFORCE
     &,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)
      call mpi_reduce(tpbsum,tpbcollect,MFORCE
     &,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)
      call mpi_reduce(tjfsum,tjfcollect,MFORCE
     &,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)
      call mpi_reduce(tausum,taucollect,MFORCE
     &,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)
      call mpi_reduce(derivsum,derivcollect,10*MFORCE
     &,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

      call mpi_allreduce(ioldest,ioldest_collect,1
     &,mpi_integer,mpi_max,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(ioldestmx,ioldestmx_collect,1
     &,mpi_integer,mpi_max,MPI_COMM_WORLD,ierr)

      ioldest=ioldest_collect
      ioldestmx=ioldestmx_collect

      call prop_reduce
      call pcm_reduce 
      call mmpol_reduce 

      call mpi_barrier(MPI_COMM_WORLD,ierr)

      if(.not.wid) goto 17

      wnow=wsum/nstep
      wfnow=wfsum/nstep
      enow=esum/wsum
      efnow=efsum/wfsum
      ei1now=wfsum/wdsum
      ei2now=wgsum(1)/wgdsum
      rinow=risum/wgsum(1)
      r2now=r2sum/wgsum(1)

      ei1cm2=ei1cm2+ei1now**2
      ei2cm2=ei2cm2+ei2now**2
      r2cm2=r2cm2+r2sum*r2now
      ricm2=ricm2+risum*rinow

      wdcum=wdcum+wdsum
      wgdcum=wgdcum+wgdsum
      ei1cum=ei1cum+ei1now
      ei2cum=ei2cum+ei2now
      r2cum=r2cum+r2sum
      ricum=ricum+risum

      wcm2=wcm2+wsum**2
      wfcm2=wfcm2+wfsum**2
      ecm2=ecm2+esum*enow
      efcm2=efcm2+efsum*efnow

      wcum=wcum+wsum
      wfcum=wfcum+wfsum
      ecum=ecum+esum
      efcum=efcum+efsum

      do 15 ifr=1,nforce

        pesum(ifr)=pecollect(ifr)
        tpbsum(ifr)=tpbcollect(ifr)
        tjfsum(ifr)=tjfcollect(ifr)
        tausum(ifr)=taucollect(ifr)
        do 13 k=1,3
  13      derivsum(k,ifr)=derivcollect(k,ifr)

        wgnow=wgsum(ifr)/nstep
        egnow=egsum(ifr)/wgsum(ifr)
        penow=pesum(ifr)/wgsum(ifr)
        tpbnow=tpbsum(ifr)/wgsum(ifr)
        tjfnow=tjfsum(ifr)/wgsum(ifr)

        wgcm2(ifr)=wgcm2(ifr)+wgsum(ifr)**2
        egcm2(ifr)=egcm2(ifr)+egsum(ifr)*egnow
        pecm2(ifr)=pecm2(ifr)+pesum(ifr)*penow
        tpbcm2(ifr)=tpbcm2(ifr)+tpbsum(ifr)*tpbnow
        tjfcm2(ifr)=tjfcm2(ifr)+tjfsum(ifr)*tjfnow

        wgcum(ifr)=wgcum(ifr)+wgsum(ifr)
        egcum(ifr)=egcum(ifr)+egsum(ifr)
        pecum(ifr)=pecum(ifr)+pesum(ifr)
        tpbcum(ifr)=tpbcum(ifr)+tpbsum(ifr)
        tjfcum(ifr)=tjfcum(ifr)+tjfsum(ifr)
        taucum(ifr)=taucum(ifr)+tausum(ifr)
        do 14 k=1,3
  14      derivcum(k,ifr)=derivcum(k,ifr)+derivsum(k,ifr)

        if(iblk.eq.1) then
          egerr=0
          peerr=0
          tpberr=0
          tjferr=0
         else
          egerr=errg(egcum(ifr),egcm2(ifr),ifr)
          peerr=errg(pecum(ifr),pecm2(ifr),ifr)
          tpberr=errg(tpbcum(ifr),tpbcm2(ifr),ifr)
          tjferr=errg(tjfcum(ifr),tjfcm2(ifr),ifr)
        endif

        egave=egcum(ifr)/wgcum(ifr)
        peave=pecum(ifr)/wgcum(ifr)
        tpbave=tpbcum(ifr)/wgcum(ifr)
        tjfave=tjfcum(ifr)/wgcum(ifr)

        if(ifr.gt.1) then
          fgcum(ifr)=fgcum(ifr)+wgsum(1)*(egnow-egsum(1)/wgsum(1))
          fgcm2(ifr)=fgcm2(ifr)+wgsum(1)*(egnow-egsum(1)/wgsum(1))**2
          fgave=egcum(1)/wgcum(1)-egcum(ifr)/wgcum(ifr)
          fgave=fgave/deltot(ifr)
          if(iblk.eq.1) then
            fgerr=0
            ifgerr=0
           else
            fgerr=errg(fgcum(ifr),fgcm2(ifr),1)
            fgerr=fgerr/abs(deltot(ifr))
            ifgerr=nint(100000* fgerr)
          endif
          egave1=egcum(1)/wgcum(1)
          derivtotave=-(derivcum(1,ifr)-derivcum(1,1)+derivcum(2,ifr)-derivcum(2,1)-egave1*(derivcum(3,ifr)-derivcum(3,1)))/wgcum(1)
         else
          call prop_cum(wgsum(ifr))
          call pcm_cum(wgsum(ifr))
          call mmpol_cum(wgsum(ifr))
        endif

c write out header first time

        if (iblk.eq.1.and.ifr.eq.1)
     &  write(6,'(t5,''egnow'',t15,''egave'',t21,''(egerr)'' ,t32
     &  ,''peave'',t38,''(peerr)'',t49,''tpbave'',t55,''(tpberr)'',t66
     &  ,''tjfave'',t72,''(tjferr)'',t83,''fgave'',t89,''(fgerr)'',
     &  t101,''npass'',t111,''wgsum'',t121,''ioldest'')')

c write out current values of averages etc.

        iegerr=nint(100000* egerr)
        ipeerr=nint(100000* peerr)
        itpber=nint(100000*tpberr)
        itjfer=nint(100000*tjferr)

        if(ifr.eq.1) then
          write(6,'(f10.5,4(f10.5,''('',i5,'')''),17x,3i10)')
     &    egsum(ifr)/wgsum(ifr),
     &    egave,iegerr,peave,ipeerr,tpbave,itpber,tjfave,itjfer,npass,
     &    nint(wgsum(ifr)),ioldest

          call prop_prt_dmc(iblk,0,wgcum,wgcm2)
          call pcm_prt(iblk,wgcum,wgcm2)
          call mmpol_prt(iblk,wgcum,wgcm2)
         else
          write(6,'(f10.5,5(f10.5,''('',i5,'')''),f10.5,10x,i10)')
     &    egsum(ifr)/wgsum(ifr),
     &    egave,iegerr,peave,ipeerr,tpbave,itpber,tjfave,itjfer,
     &    fgave,ifgerr,derivtotave,nint(wgsum(ifr))
        endif
   15 continue

c     call flush(6)

c zero out xsum variables for metrop

   17 wsum=zero
      wfsum=zero
      wdsum=zero
      wgdsum=zero
      esum=zero
      efsum=zero
      ei1sum=zero
      ei2sum=zero
      r2sum=zero
      risum=zero

      do 18 ifr=1,nforce
        egsum(ifr)=zero
        wgsum(ifr)=zero
        pesum(ifr)=zero
        tpbsum(ifr)=zero
        tjfsum(ifr)=zero
        tausum(ifr)=zero
        do 18 k=1,10
   18     derivsum(k,ifr)=zero

      call prop_init(1)
      call pcm_init(1)
      call mmpol_init(1)
      
      return

      entry acues1

c     call flush(11)
c statistical fluctuations without blocking
      wdsum1=wdsumo
      wgdsum1=wgdsumo

      call mpi_reduce(esum1(1),ecollect,1
     &,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)
      call mpi_reduce(wsum1(1),wcollect,1
     &,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)
      call mpi_reduce(efsum1,efcollect,1
     &,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)
      call mpi_reduce(wfsum1,wfcollect,1
     &,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

      call mpi_reduce(wgsum1,wgcollect,MFORCE
     &,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)
      call mpi_reduce(egsum1,egcollect,MFORCE
     &,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)
      call mpi_reduce(tausum(1),taublock,1
     &,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

      call mpi_barrier(MPI_COMM_WORLD,ierr)

      if(.not.wid) goto 23

      wsum1(1)=wcollect
      esum1(1)=ecollect
      efsum1=efcollect
      wfsum1=wfcollect
      do 20 ifr=1,nforce
        wgsum1(ifr)=wgcollect(ifr)
   20   egsum1(ifr)=egcollect(ifr)

      wcum1=wcum1+wsum1(1)
      wfcum1=wfcum1+wfsum1
      ecum1=ecum1+esum1(1)
      efcum1=efcum1+efsum1
      ei3cum=ei3cum+wfsum1/wdsum1

      wcm21=wcm21+wsum1(1)**2
      wfcm21=wfcm21+wfsum1**2
      ecm21=ecm21+esum1(1)**2/wsum1(1)
      efcm21=efcm21+efsum1**2/wfsum1
      ei3cm2=ei3cm2+(wfsum1/wdsum1)**2

      do 21 ifr=1,nforce
        wgcum1(ifr)=wgcum1(ifr)+wgsum1(ifr)
        egcum1(ifr)=egcum1(ifr)+egsum1(ifr)
        wgcm21(ifr)=wgcm21(ifr)+wgsum1(ifr)**2
   21   egcm21(ifr)=egcm21(ifr)+egsum1(ifr)**2/wgsum1(ifr)

c sum1 block averages
      wsum=wsum+wsum1(1)
      wfsum=wfsum+wfsum1
      wdsum=wdsum+wdsum1
      wgdsum=wgdsum+wgdsum1
      esum=esum+esum1(1)
      efsum=efsum+efsum1
      eisum=eisum+wfsum1/wdsum1
      do 22 ifr=1,nforce
        wgsum(ifr)=wgsum(ifr)+wgsum1(ifr)
   22   egsum(ifr)=egsum(ifr)+egsum1(ifr)

c Estimate eigenvalue of G from the energy
      ipmod=mod(ipass,nfprod)
      if(iabs(idmc).eq.1) then
        nfpro=min(nfprod,ipass)
        eigv=(wgsum1(1)/wtgen(ipmod))**(one/nfpro)
       else
        eest=(egcum(1)+egsum(1))/(wgcum(1)+wgsum(1))
        eigv=dexp((etrial-eest)*(taucum(1)+taublock)/
     &                          (wgcum(1)+wgsum(1)))
        if(ipr.ge.1) write(6,'(''eigv'',9f14.6)') eigv,eest,accavn,
     &  egcum(1),egsum(1),wgcum(1),wgsum(1),fprod
      endif

      wdsumo=wsum1(1)
      wgdsumo=wsum1(1)*fprod/ff(mod(ipass+1,nfprod))
      wtgen(ipmod)=wsum1(1)

   23 call mpi_bcast(eigv,1,mpi_double_precision,0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(eest,1,mpi_double_precision,0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(wdsumo,1,mpi_double_precision,0,MPI_COMM_WORLD
     &,ierr)

      call mpi_barrier(MPI_COMM_WORLD,ierr)

      call redistribute

      call mpi_barrier(MPI_COMM_WORLD,ierr)

c zero out step averages
      wfsum1=zero
      wdsum1=zero
      efsum1=zero
      do 24 ifr=1,nforce
        wsum1(ifr)=zero
        wgsum1(ifr)=zero
        esum1(ifr)=zero
   24   egsum1(ifr)=zero

      return

      entry zeres0
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
      wdsumo=nconf*nproc
      wgdsumo=nconf*nproc
      fprod=one
      do 70 i=0,MFPRD1
        wtgen(i)=nconf*nproc
   70   ff(i)=one

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
          endif
          pwt(iw,ifr)=one
          do 72 ip=0,nwprod-1
   72       wthist(iw,ip,ifr)=one
   80 continue

      nconf=nconf*nproc

      entry zerest
c entry point to zero out all averages etc. after equilibration runs

      iblk=0

c zero out estimators

      wcum1=zero
      wfcum1=zero
      wcum=zero
      wfcum=zero
      wdcum=zero
      wgdcum=zero
      ecum1=zero
      efcum1=zero
      ecum=zero
      efcum=zero
      ei1cum=zero
      ei2cum=zero
      ei3cum=zero
      r2cum=zero
      ricum=zero

      wcm21=zero
      wfcm21=zero
      wcm2=zero
      wfcm2=zero
      wdcm2=zero
      wgdcm2=zero
      ecm21=zero
      efcm21=zero
      ecm2=zero
      efcm2=zero
      ei1cm2=zero
      ei2cm2=zero
      ei3cm2=zero
      r2cm2=zero
      ricm2=zero

      wfsum1=zero
      wsum=zero
      wfsum=zero
      wdsum=zero
      wgdsum=zero
      efsum1=zero
      esum=zero
      efsum=zero
      ei1sum=zero
      ei2sum=zero
      ei3sum=zero
      r2sum=zero
      risum=zero


      do 85 ifr=1,nforce
        tausum(ifr)=zero
        taucum(ifr)=zero
        wgcum1(ifr)=zero
        wgcum(ifr)=zero
        egcum1(ifr)=zero
        egcum(ifr)=zero
        wgcm21(ifr)=zero
        wgcm2(ifr)=zero
        egcm21(ifr)=zero
        egcm2(ifr)=zero
        wsum1(ifr)=zero
        wgsum1(ifr)=zero
        wgsum(ifr)=zero
        esum1(ifr)=zero
        egsum1(ifr)=zero
        egsum(ifr)=zero
        pecum(ifr)=zero
        tpbcum(ifr)=zero
        tjfcum(ifr)=zero
        pecm2(ifr)=zero
        tpbcm2(ifr)=zero
        tjfcm2(ifr)=zero
        pesum(ifr)=zero
        tpbsum(ifr)=zero
        tjfsum(ifr)=zero
        fgcum(ifr)=zero
        fgcm2(ifr)=zero
        do 85 k=1,10
          derivsum(k,ifr)=zero
   85     derivcum(k,ifr)=zero

      nbrnch=0

      trymove=0
      acc=0
      nacc=0
      nodecr=0

c Zero out estimators for charge density of atom.
      do 90 i=1,nrad
        rprobup(i)=zero
        rprobdn(i)=zero
   90   rprob(i)=zero

      call prop_init(0)
      call pcm_init(0)
      call mmpol_init(0)

      return
      end
