      subroutine dumper
c Written by Cyrus Umrigar, modified by Claudia Filippi
c routine to pick up and dump everything needed to restart
c job where it left off
      implicit real*8(a-h,o-z)
      include 'vmc.h'
      include 'dmc.h'
      include 'force.h'
      include 'mstates.h'
      include 'pseudo.h'
      include 'basis.h'
      parameter (zero=0.d0,one=1.d0,small=1.d-6)

      common /forcepar/ deltot(MFORCE),nforce,istrech
      common /forcest/ fgcum(MFORCE),fgcm2(MFORCE)
      common /force_dmc/ itausec,nwprod
      common /const/ pi,hb,etrial,delta,deltai,fbias,nelec,imetro,ipr
      common /elec/ nup,ndn
      common /contrl/ nstep,nblk,nblkeq,nconf,nconf_new,isite,idump,irstar
      common /contrldmc/ tau,rttau,taueff(MFORCE),tautot,nfprod,idmc,ipq
     &,itau_eff,iacc_rej,icross,icuspg,idiv_v,icut_br,icut_e
      common /iterat/ ipass,iblk
      common /config/ xold(3,MELEC,MWALK,MFORCE),vold(3,MELEC,MWALK,MFORCE),
     &psido(MWALK,MFORCE),psijo(MWALK,MFORCE),peo(MWALK,MFORCE),d2o(MWALK,MFORCE)
      common /velratio/ fratio(MWALK,MFORCE),xdrifted(3,MELEC,MWALK,MFORCE)
      common /coefs/ coef(MBASIS,MORB,MWF),nbasis,norb
      common /atom/ znuc(MCTYPE),cent(3,MCENT),pecent
     &,iwctype(MCENT),nctype,ncent
      common /ghostatom/ newghostype,nghostcent
      common /pseudo/ vps(MELEC,MCENT,MPS_L),vpso(MELEC,MCENT,MPS_L,MFORCE)
     &,lpot(MCTYPE),nloc
      common /qua/ xq0(MPS_QUAD),yq0(MPS_QUAD),zq0(MPS_QUAD)
     &,xq(MPS_QUAD),yq(MPS_QUAD),zq(MPS_QUAD),wq(MPS_QUAD),nquad
      common /dets/ cdet(MDET,MSTATES,MWF),ndet
      common /jaspar1/ cjas1(MWF),cjas2(MWF)
      common /jaspar/ nspin1,nspin2,sspin,sspinn,is
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
      common /derivest/ derivsum(10,MFORCE),derivcum(10,MFORCE),derivcm2(MFORCE),
     &derivtotave_num_old(MFORCE)
      common /stats/ dfus2ac,dfus2un,dr2ac,dr2un,acc,trymove,nacc,
     &nbrnch,nodecr
      common /step/try(nrad),suc(nrad),trunfb(nrad),rprob(nrad),
     &ekin(nrad),ekin2(nrad)
      common /denupdn/ rprobup(nrad),rprobdn(nrad)
      common /branch/ wtgen(0:MFPRD1),ff(0:MFPRD1),eold(MWALK,MFORCE),
     &pwt(MWALK,MFORCE),wthist(MWALK,0:MFORCE_WT_PRD,MFORCE),
     &wt(MWALK),eigv,eest,wdsumo,wgdsumo,fprod,nwalk
      common /age/ iage(MWALK),ioldest,ioldestmx
      common /jacobsave/ ajacob,ajacold(MWALK,MFORCE)

      common /casula/ t_vpsp(MCENT,MPS_QUAD,MELEC),icasula,i_vpsp

      dimension irn(4)
      dimension coefx(MBASIS,MORB),zexx(MBASIS),centx(3,MCENT)
     &,znucx(MCENT),n1sx(MCENT),n2sx(MCENT),n2px(3,MCENT)
     &,n3sx(MCENT),n3px(3,MCENT),n3dzrx(MCENT),n3dx2x(MCENT)
     &,n3dxyx(MCENT),n3dxzx(MCENT),n3dyzx(MCENT),n4sx(MCENT)
     &,n4px(3,MCENT),nsax(MCENT),npax(3,MCENT),ndzrax(MCENT)
     &,ndx2ax(MCENT),ndxyax(MCENT),ndxzax(MCENT),ndyzax(MCENT)
     &,cdetx(MDET)

      dimension ddx_ref(3)

      if(nforce.gt.1) call strech(xold(1,1,1,1),xold(1,1,1,1),ajacob,1,0)

      call savern(irn)

      open(unit=10,status='unknown',form='unformatted',file='restart_dmc')
      rewind 10
      write(10) irn
      write(10) hb
      write(10) tau,rttau,taueff(1),tautot,idmc,nfprod
      write(10) nelec,nconf,nforce
      write(10) nwalk
      write(10) (wtgen(i),ff(i),i=0,nfprod),(wt(i),i=1,nwalk)
     &,eigv,eest,wdsumo,wgdsumo,fprod
      if(nforce.gt.1) write(10) (taueff(i),i=2,nforce)
      if(nforce.gt.1) write(10) nwprod
     &,((pwt(i,j),i=1,nwalk),j=1,nforce)
     &,(((wthist(i,l,j),i=1,nwalk),l=0,nwprod-1),j=1,nforce)
      write(10) (iage(i),i=1,nwalk),ioldest,ioldestmx
      write(10) (((xold(k,i,iw,1),k=1,3),i=1,nelec),iw=1,nwalk)
      write(10) ((fratio(iw,ifr),iw=1,nwalk),ifr=1,nforce)
      write(10) ((((xdrifted(k,i,iw,ifr),k=1,3),i=1,nelec),iw=1,nwalk),ifr=1,nforce)
      write(10) wcum,wfcum,(wgcum(i),i=1,nforce),wdcum,wgdcum,wcum1
     &,wfcum1,(wgcum1(i),i=1,nforce),wdcum1, ecum,efcum
     &,(egcum(i),i=1,nforce), ecum1,efcum1,(egcum1(i),i=1,nforce)
     &,ei1cum,ei2cum,ei3cum, (pecum(i),i=1,nforce)
     &,(tpbcum(i),i=1,nforce),(tjfcum(i),i=1,nforce),r2cum,ricum
     &,(taucum(i),i=1,nforce)
     &,((derivcum(k,i),k=1,3),i=1,nforce),(derivcm2(i),i=1,nforce)
     &,(derivtotave_num_old(i),i=1,nforce)
      write(10) ipass,iblk
      write(10) wcm2,wfcm2,(wgcm2(i),i=1,nforce),wdcm2,wgdcm2,wcm21
     &,wfcm21,(wgcm21(i),i=1,nforce),wdcm21, ecm2,efcm2
     &,(egcm2(i),i=1,nforce), ecm21,efcm21,(egcm21(i),i=1,nforce)
     &,ei1cm2,ei2cm2,ei3cm2, (pecm2(i),i=1,nforce)
     &,(tpbcm2(i),i=1,nforce),(tjfcm2(i),i=1,nforce),r2cm2,ricm2
      write(10) (fgcum(i),i=1,nforce),(fgcm2(i),i=1,nforce)
      write(10) (rprob(i),rprobup(i),rprobdn(i),i=1,nrad)
      write(10) dfus2ac,dfus2un,dr2ac,dr2un,acc
     &,trymove,nacc,nbrnch,nodecr
      call prop_dump(10)
      call pcm_dump(10)
      call mmpol_dump(10)
      if(nloc.gt.0) write(10) nquad,(xq(i),yq(i),zq(i),wq(i),i=1,nquad)

      write(10) ((coef(ib,i,1),ib=1,nbasis),i=1,norb)
      write(10) nbasis
      write(10) (zex(ib,1),ib=1,nbasis)
      write(10) nctype,ncent,newghostype,nghostcent,(iwctype(i),i=1,ncent+nghostcent)
      write(10) ((cent(k,ic),k=1,3),ic=1,ncent+nghostcent)
      write(10) pecent
      write(10) (znuc(i),i=1,nctype)
      write(10) (n1s(i),i=1,nctype)
      write(10) (n2s(i),i=1,nctype)
      write(10) ((n2p(k,i),k=1,3),i=1,nctype)
      write(10) (n3s(i),i=1,nctype)
      write(10) ((n3p(k,i),k=1,3),i=1,nctype)
      write(10) (n3dzr(i),i=1,nctype)
      write(10) (n3dx2(i),i=1,nctype)
      write(10) (n3dxy(i),i=1,nctype)
      write(10) (n3dxz(i),i=1,nctype)
      write(10) (n3dyz(i),i=1,nctype)
      write(10) (n4s(i),i=1,nctype)
      write(10) ((n4p(k,i),k=1,3),i=1,nctype)
      write(10) (nsa(i),i=1,nctype)
      write(10) ((npa(k,i),k=1,3),i=1,nctype)
      write(10) (ndzra(i),i=1,nctype)
      write(10) (ndx2a(i),i=1,nctype)
      write(10) (ndxya(i),i=1,nctype)
      write(10) (ndxza(i),i=1,nctype)
      write(10) (ndyza(i),i=1,nctype)
      write(10) (cdet(i,1,1),i=1,ndet)
      write(10) ndet,nup,ndn
      write(10) cjas1(1),cjas2(1)
      rewind 10

      write(6,'(1x,''successful dump to unit 10'')')
      close(10)
      return

c-------------------------------------------------------------------------
      entry startr

      write(6,'(1x,''attempting restart from unit 10'')')
      rewind 10
      read(10) irn
      call setrn(irn)
      read(10) hbx
      read(10) taux,rttau,taueff(1),tautot,idmc,nfprod
      read(10) nelecx,nconf,nforce
      if (dabs(hbx-hb).gt.small) call fatal_error('STARTR: hb')
      if (dabs(taux-tau).gt.small) call fatal_error('STARTR: tau')
      if (nelecx.ne.nelec) call fatal_error('STARTR: nelec')
      read(10) nwalk
      read(10) (wtgen(i),ff(i),i=0,nfprod),(wt(i),i=1,nwalk)
     &,eigv,eest,wdsumo,wgdsumo,fprod
      if(nforce.gt.1) read(10) (taueff(i),i=2,nforce)
      if(nforce.gt.1) read(10) nwprod
     &,((pwt(i,j),i=1,nwalk),j=1,nforce)
     &,(((wthist(i,l,j),i=1,nwalk),l=0,nwprod-1),j=1,nforce)
      read(10) (iage(i),i=1,nwalk),ioldest,ioldestmx
      read(10) (((xold(k,i,iw,1),k=1,3),i=1,nelec),iw=1,nwalk)
      read(10) ((fratio(iw,ifr),iw=1,nwalk),ifr=1,nforce)
      read(10) ((((xdrifted(k,i,iw,ifr),k=1,3),i=1,nelec),iw=1,nwalk),ifr=1,nforce)
      read(10) wcum,wfcum,(wgcum(i),i=1,nforce),wdcum,wgdcum,wcum1
     &,wfcum1,(wgcum1(i),i=1,nforce),wdcum1, ecum,efcum
     &,(egcum(i),i=1,nforce), ecum1,efcum1,(egcum1(i),i=1,nforce)
     &,ei1cum,ei2cum,ei3cum, (pecum(i),i=1,nforce)
     &,(tpbcum(i),i=1,nforce),(tjfcum(i),i=1,nforce),r2cum,ricum
     &,(taucum(i),i=1,nforce)
     &,((derivcum(k,i),k=1,3),i=1,nforce),(derivcm2(i),i=1,nforce)
     &,(derivtotave_num_old(i),i=1,nforce)
      read(10) ipass,iblk
      read(10) wcm2,wfcm2,(wgcm2(i),i=1,nforce),wdcm2,wgdcm2,wcm21
     &,wfcm21,(wgcm21(i),i=1,nforce),wdcm21, ecm2,efcm2
     &,(egcm2(i),i=1,nforce), ecm21,efcm21,(egcm21(i),i=1,nforce)
     &,ei1cm2,ei2cm2,ei3cm2, (pecm2(i),i=1,nforce)
     &,(tpbcm2(i),i=1,nforce),(tjfcm2(i),i=1,nforce),r2cm2,ricm2
      read(10) (fgcum(i),i=1,nforce),(fgcm2(i),i=1,nforce)
      read(10) (rprob(i),rprobup(i),rprobdn(i),i=1,nrad)
      read(10) dfus2ac,dfus2un,dr2ac,dr2un,acc
     &,trymove,nacc,nbrnch,nodecr
      call prop_rstrt(10)
      call pcm_rstrt(10)
      call mmpol_rstrt(10)
      if(nloc.gt.0) then
        read(10) nqx,(xq(i),yq(i),zq(i),wq(i),i=1,nqx)
        if(nqx.ne.nquad) call fatal_error('STARTR: nquad')
      endif

      read(10) ((coefx(ib,i),ib=1,nbasis),i=1,norb)
      read(10) nbasx
      do 10 j=1,norb
      do 10 i=1,nbasis
      if (dabs(coefx(i,j)-coef(i,j,1)).gt.small) call fatal_error('STARTR: coef')
   10 continue
      if (nbasx.ne.nbasis) call fatal_error('STARTR: nbasis')
      read(10) (zexx(ib),ib=1,nbasis)
      read(10) nctypex,ncentx,newghostypex,nghostcentx,(iwctype(i),i=1,ncentx+nghostcentx)
      read(10) ((centx(k,ic),k=1,3),ic=1,ncentx+nghostcentx)
      read(10) pecent
      read(10) (znucx(i),i=1,nctypex)
      read(10) (n1sx(i),i=1,nctypex)
      read(10) (n2sx(i),i=1,nctypex)
      read(10) ((n2px(k,i),k=1,3),i=1,nctypex)
      read(10) (n3sx(i),i=1,nctypex)
      read(10) ((n3px(k,i),k=1,3),i=1,nctypex)
      read(10) (n3dzrx(i),i=1,nctypex)
      read(10) (n3dx2x(i),i=1,nctypex)
      read(10) (n3dxyx(i),i=1,nctypex)
      read(10) (n3dxzx(i),i=1,nctypex)
      read(10) (n3dyzx(i),i=1,nctypex)
      read(10) (n4sx(i),i=1,nctypex)
      read(10) ((n4px(k,i),k=1,3),i=1,nctypex)
      read(10) (nsax(i),i=1,nctypex)
      read(10) ((npax(k,i),k=1,3),i=1,nctypex)
      read(10) (ndzrax(i),i=1,nctypex)
      read(10) (ndx2ax(i),i=1,nctypex)
      read(10) (ndxyax(i),i=1,nctypex)
      read(10) (ndxzax(i),i=1,nctypex)
      read(10) (ndyzax(i),i=1,nctypex)

      if (ncentx.ne.ncent) call fatal_error('STARTR: ncent')
      if (nctypex.ne.nctype) call fatal_error('STARTR: nctype')
      do 20 i=1,nbasis
      if (dabs(zexx(i)-zex(i,1)).gt.small) call fatal_error('STARTR: zex')
   20 continue
      do 30 i=1,ncent+nghostcent
      do 30 k=1,3
      if (dabs(cent(k,i)-centx(k,i)).gt.small) call fatal_error('STARTR: cent')
   30 continue
      do 40 i=1,nctype
      if (dabs(znucx(i)-znuc(i)).gt.small) call fatal_error('STARTR: znuc')
      if (n1s(i).ne.n1sx(i)) call fatal_error('STARTR: n1s')
      if (n2s(i).ne.n2sx(i)) call fatal_error('STARTR: n2s')
      if (n3s(i).ne.n3sx(i)) call fatal_error('STARTR: n3s')
      if (n3dzr(i).ne.n3dzrx(i)) call fatal_error('STARTR: n3dzrx')
      if (n3dx2(i).ne.n3dx2x(i)) call fatal_error('STARTR: n3dx2x')
      if (n3dxy(i).ne.n3dxyx(i)) call fatal_error('STARTR: n3dxy')
      if (n3dxz(i).ne.n3dxzx(i)) call fatal_error('STARTR: n3dxzx')
      if (n3dyz(i).ne.n3dyzx(i)) call fatal_error('STARTR: n3dyz')
      if (n4s(i).ne.n4sx(i)) call fatal_error('STARTR: n4s')
      if (nsa(i).ne.nsax(i)) call fatal_error('STARTR: nsa')
      if (ndzra(i).ne.ndzrax(i)) call fatal_error('STARTR: ndzra')
      if (ndx2a(i).ne.ndx2ax(i)) call fatal_error('STARTR: ndx2a')
      if (ndxya(i).ne.ndxyax(i)) call fatal_error('STARTR: ndxya')
      if (ndxza(i).ne.ndxzax(i)) call fatal_error('STARTR: ndxza')
      if (ndyza(i).ne.ndyzax(i)) call fatal_error('STARTR: ndyza')
      do 40 k=1,3
      if (n2p(k,i).ne.n2px(k,i)) call fatal_error('STARTR: n2p')
      if (n3p(k,i).ne.n3px(k,i)) call fatal_error('STARTR: n3p')
      if (n4p(k,i).ne.n4px(k,i)) call fatal_error('STARTR: n4p')
      if (npa(k,i).ne.npax(k,i)) call fatal_error('STARTR: npa')
   40 continue
      read(10) (cdetx(i),i=1,ndet)
      read(10) ndetx,nupx,ndnx
      do 50 i=1,ndet
      if (dabs(cdetx(i)-cdet(i,1,1)).gt.small) call fatal_error('STARTR: cdet')
   50 continue
      if (ndetx.ne.ndet) call fatal_error('STARTR: ndet')
      if (nupx.ne.nup) call fatal_error('STARTR: nup')
      if (ndnx.ne.ndn) call fatal_error('STARTR: ndn')
      read(10) cjas1x,cjas2x
      if (dabs(cjas1x-cjas1(1)).gt.small) call fatal_error('STARTR: cjas1')
      if (dabs(cjas2x-cjas2(1)).gt.small) call fatal_error('STARTR: cjas2')
      write(6,'(1x,''succesful read from unit 10'')')
      write(6,'(t5,''egnow'',t15,''egave'',t21
     &,''(egerr)'' ,t32,''peave'',t38,''(peerr)'',t49,''tpbave'',t55
     &,''(tpberr)'' ,t66,''tjfave'',t72,''(tjferr)'',t83,''npass'',t93
     &,''wgsum'',t103 ,''ioldest'')')

      do 70 iw=1,nwalk
        if(istrech.eq.0) then
          do 60 ifr=2,nforce
            do 60 ie=1,nelec
              do 60 k=1,3
   60           xold(k,ie,iw,ifr)=xold(k,ie,iw,1)
        endif
        do 70 ifr=1,nforce
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
          call hpsi(xold(1,1,iw,ifr),psido(iw,ifr),psijo(iw,ifr),eold(iw,ifr),ifr)
          i_vpsp=0
          do 65 i=1,nelec
   65       call compute_determinante_grad(i,psido(iw,ifr),psido(iw,ifr),vold(1,i,iw,ifr),1)
          if(ifr.eq.1) then
            call walksav_det(iw)
            call walksav_jas(iw)
c           call t_vpsp_sav(iw)
            call t_vpsp_sav
            call prop_save_dmc(iw)
            call pcm_save(iw)
            call mmpol_save(iw)
          endif
   70 continue

c zero out xsum variables

      wsum=zero
      wfsum=zero
      wdsum=zero
      wgdsum=zero
      esum=zero
      efsum=zero
      ei1sum=zero
      ei2sum=zero
      r2sum=zero
      risum=zero

      do 75 ifr=1,nforce
        egsum(ifr)=zero
        wgsum(ifr)=zero
        pesum(ifr)=zero
        tpbsum(ifr)=zero
        tjfsum(ifr)=zero
        tausum(ifr)=zero
        do 75 k=1,3
   75     derivsum(k,ifr)=zero

      call prop_init(1)
      call pcm_init(1)
      call mmpol_init(1)

      if(ipr.ge.-2) then
        open(unit=11,file='walkalize',status='old')
        do 80 i=1,2000000000
   80     read(11,fmt=*,end=90)
      endif
   90 backspace 11
      backspace 11

      return
      end
