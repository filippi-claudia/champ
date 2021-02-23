      subroutine startr_gpop

      use vmc_mod, only: MELEC, MORB, MBASIS, MDET, MCENT, MCTYPE, MCTYP3X,
     &NSPLIN, nrad, MORDJ, MORDJ1, MMAT_DIM, MMAT_DIM2, MMAT_DIM20,
     &radmax, delri, NEQSX, MTERMS, MCENT3, NCOEF, MEXCIT
      use dmc_mod, only: MWALK, MFPROD, MFPRD1, MPATH
      use basis, only: zex, betaq, n1s, n2s, n2p, n3s, n3p, n3dzr, n3dx2, n3dxy, n3dxz, n3dyz,
     & n4s, n4p, n4fxxx, n4fyyy, n4fzzz, n4fxxy, n4fxxz, n4fyyx, n4fyyz,
     & n4fzzx, n4fzzy, n4fxyz, nsa, npa, ndzra, ndz2a, ndxya, ndxza, ndyza, ndx2a

      use const, only: delta, deltai, etrial, fbias, hb, imetro, ipr, nelec, pi
      use forcest, only: fgcm2, fgcum
      use forcepar, only: deltot, istrech, nforce
      use age, only: iage, ioldest, ioldestmx
      use contrldmc, only: iacc_rej, icross, icuspg, icut_br, icut_e, idiv_v, idmc, ipq,
     &itau_eff, nfprod, rttau, tau, taueff, tautot
      use atom, only: cent, iwctype, ncent, nctype, pecent, znuc

      use iterat, only: iblk, ipass
      use config, only: d2o, peo_dmc, psido_dmc, psijo_dmc, vold_dmc, xold_dmc

      use stats, only: acc, dfus2ac, dfus2un, dr2ac, dr2un, nacc, nbrnch, nodecr, trymove

      use estsum, only: efsum, efsum1, egsum, egsum1, ei1sum, ei2sum, ei3sum, esum1_dmc, esum_dmc,
     &pesum_dmc, r2sum, risum, tausum, tjfsum_dmc, tpbsum_dmc, w_acc_sum, w_acc_sum1, wdsum,
     &wdsum1, wfsum, wfsum1, wg_acc_sum, wg_acc_sum1, wgdsum, wgsum, wgsum1, wsum1, wsum_dmc
      use estcum, only: ecum1_dmc, ecum_dmc, efcum, efcum1, egcum, egcum1, ei1cum, ei2cum,
     &ei3cum, pecum_dmc, r2cum_dmc, ricum, taucum, tjfcum_dmc, tpbcum_dmc, w_acc_cum, w_acc_cum1,
     &wcum1, wcum_dmc, wdcum, wdcum1, wfcum, wfcum1, wg_acc_cum, wg_acc_cum1, wgcum, wgcum1,
     &wgdcum
      use force_dmc, only: itausec, nwprod
      use est2cm, only: ecm21_dmc, ecm2_dmc, efcm2, efcm21, egcm2, egcm21, ei1cm2, ei2cm2,
     &ei3cm2, pecm2_dmc, r2cm2_dmc, ricm2, tjfcm_dmc, tpbcm2_dmc, wcm2, wcm21, wdcm2, wdcm21,
     &wfcm2, wfcm21, wgcm2, wgcm21, wgdcm2
      use derivest, only: derivcm2, derivcum, derivsum, derivtotave_num_old
      use step, only: ekin, ekin2, rprob, suc, trunfb, try
      use mpiconf, only: idtask, nproc, wid, NPROCX
      use denupdn, only: rprobdn, rprobup
      use force_mod, only: MFORCE, MFORCE_WT_PRD, MWF
      use mstates_mod, only: MSTATES, MDETCSFX
      use pseudo_mod, only: MPS_L, MPS_QUAD, MPS_GRID, MGAUSS

      use qua, only: nquad, wq, xq, xq0, yq, yq0, zq, zq0

      use branch, only: eest, eigv, eold, ff, fprod, nwalk, pwt, wdsumo, wgdsumo, wt, wtgen,
     &wthist
      use casula, only: i_vpsp, icasula, t_vpsp
      implicit real*8(a-h,o-z)




      include 'mpif.h'

      parameter (zero=0.d0,one=1.d0)
      parameter (small=1.e-6)

      common /contrl/ nstep,nblk,nblkeq,nconf,nconf_new,isite,idump,irstar
      common /velratio/ fratio(MWALK,MFORCE)
      common /coefs/ coef(MBASIS,MORB,MWF),nbasis,norb
      common /ghostatom/ newghostype,nghostcent
      common /pseudo/ vps(MELEC,MCENT,MPS_L),vpso(MELEC,MCENT,MPS_L,MFORCE)
     &,lpot(MCTYPE),nloc
      common /dets/ cdet(MDET,MSTATES,MWF),ndet
      common /elec/ nup,ndn
      common /jaspar1/ cjas1(MWF),cjas2(MWF)
      common /jaspar/ nspin1,nspin2,sspin,sspinn,is
      common /jacobsave/ ajacob,ajacold(MWALK,MFORCE)


      character*13 filename

      dimension irn(4,0:NPROCX),istatus(MPI_STATUS_SIZE)
      dimension irn_tmp(4,0:NPROCX)

      dimension coefx(MBASIS,MORB),zexx(MBASIS),centx(3,MCENT)
     &,znucx(MCENT),n1sx(MCENT),n2sx(MCENT),n2px(3,MCENT)
     &,n3sx(MCENT),n3px(3,MCENT),n3dzrx(MCENT),n3dx2x(MCENT)
     &,n3dxyx(MCENT),n3dxzx(MCENT),n3dyzx(MCENT),n4sx(MCENT)
     &,n4px(3,MCENT),nsax(MCENT),npax(3,MCENT),ndzrax(MCENT)
     &,ndx2ax(MCENT),ndxyax(MCENT),ndxzax(MCENT),ndyzax(MCENT)
     &,cdetx(MDET)

      write(6,'(1x,''attempting restart from unit 10'')')
      rewind 10
      read(10) nprock
      read(10) nfprod,(ff(i),i=0,nfprod),fprod,eigv,eest,wdsumo
     &,ioldest,ioldestmx
      if(nprock.ne.nproc) call fatal_error('STARTR: different num procs')
      do 4 id=0,idtask
        read(10) nwalk
        read(10) (wt(i),i=1,nwalk),(iage(i),i=1,nwalk)
        read(10) (((xold_dmc(ic,i,iw,1),ic=1,3),i=1,nelec),iw=1,nwalk)
        read(10) nforce,((fratio(iw,ifr),iw=1,nwalk),ifr=1,nforce)
    4   if(nloc.gt.0)
     &  read(10) nquad,(xq(i),yq(i),zq(i),wq(i),i=1,nquad)
      do 5 id=idtask+1,nproc-1
        read(10) nwalk_id
        read(10) (wt_id,i=1,nwalk_id),(iage_id,i=1,nwalk_id)
        read(10) (xold_dmc_id,i=1,3*nelec*nwalk_id)
        read(10) nf_id,((fratio_id,iw=1,nwalk_id),ifr=1,nf_id)
    5   if(nloc.gt.0)
     &  read(10) nq_id,(xq_id,yq_id,zq_id,wq_id,i=1,nquad)
c     if(nforce.gt.1) read(10) nwprod
c    &,((pwt(i,j),i=1,nwalk),j=1,nforce)
c    &,(((wthist(i,l,j),i=1,nwalk),l=0,nwprod-1),j=1,nforce)
      read(10) (wgcum(i),egcum(i),pecum_dmc(i),tpbcum_dmc(i),tjfcum_dmc(i),
     &wgcm2(i),egcm2(i),pecm2_dmc(i),tpbcm2_dmc(i),tjfcm_dmc(i),taucum(i),
     &i=1,nforce)
      read(10) ((irn(i,j),i=1,4),j=0,nproc-1)
      call setrn(irn(1,idtask))
      read(10) hbx
      read(10) taux,rttau,idmc
      read(10) nelecx,nconf
      if (dabs(hbx-hb).gt.small) call fatal_error('STARTR: hb')
      if (dabs(taux-tau).gt.small) call fatal_error('STARTR: tau')
      if (nelecx.ne.nelec) call fatal_error('STARTR: nelec')
      read(10) (wtgen(i),i=0,nfprod),wgdsumo
      read(10) wcum_dmc,wfcum,wdcum,wgdcum,wcum1
     &,wfcum1,(wgcum1(i),i=1,nforce),wdcum1, ecum_dmc,efcum
     &,ecum1_dmc,efcum1,(egcum1(i),i=1,nforce)
     &,ei1cum,ei2cum,ei3cum, r2cum_dmc,ricum
      read(10) ipass,iblk
      read(10) wcm2,wfcm2,wdcm2,wgdcm2,wcm21
     &,wfcm21,(wgcm21(i),i=1,nforce),wdcm21, ecm2_dmc,efcm2
     &,ecm21_dmc,efcm21,(egcm21(i),i=1,nforce)
     &,ei1cm2,ei2cm2,ei3cm2,r2cm2_dmc,ricm2
      read(10) (fgcum(i),i=1,nforce),(fgcm2(i),i=1,nforce)
     &,((derivcum(k,i),k=1,3),i=1,nforce)
      read(10) (rprob(i),rprobup(i),rprobdn(i),i=1,nrad)
      read(10) dfus2ac,dfus2un,dr2ac,dr2un,acc
     &,trymove,nacc,nbrnch,nodecr
      if(.not.wid) then
        acc=0
        nacc=0
        trymove=0
        nodecr=0
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
      read(10) ((n2px(ic,i),ic=1,3),i=1,nctypex)
      read(10) (n3sx(i),i=1,nctypex)
      read(10) ((n3px(ic,i),ic=1,3),i=1,nctypex)
      read(10) (n3dzrx(i),i=1,nctypex)
      read(10) (n3dx2x(i),i=1,nctypex)
      read(10) (n3dxyx(i),i=1,nctypex)
      read(10) (n3dxzx(i),i=1,nctypex)
      read(10) (n3dyzx(i),i=1,nctypex)
      read(10) (n4sx(i),i=1,nctypex)
      read(10) ((n4px(ic,i),ic=1,3),i=1,nctypex)
      read(10) (nsax(i),i=1,nctypex)
      read(10) ((npax(ic,i),ic=1,3),i=1,nctypex)
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
      do 40 ic=1,3
      if (n2p(ic,i).ne.n2px(ic,i)) call fatal_error('STARTR: n2p')
      if (n3p(ic,i).ne.n3px(ic,i)) call fatal_error('STARTR: n3p')
      if (n4p(ic,i).ne.n4px(ic,i)) call fatal_error('STARTR: n4p')
      if (npa(ic,i).ne.npax(ic,i)) call fatal_error('STARTR: npa')
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
   60           xold_dmc(k,ie,iw,ifr)=xold_dmc(k,ie,iw,1)
        endif
        do 70 ifr=1,nforce
          if(nforce.gt.1) then
            if(ifr.eq.1.or.istrech.eq.0) then
              call strech(xold_dmc(1,1,iw,1),xold_dmc(1,1,iw,ifr),ajacob,ifr,0)
               else
              call strech(xold_dmc(1,1,iw,1),xold_dmc(1,1,iw,ifr),ajacob,ifr,1)
            endif
           else
            ajacob=one
          endif
          ajacold(iw,ifr)=ajacob
          if(icasula.lt.0) i_vpsp=icasula
          call hpsi(xold_dmc(1,1,iw,ifr),psido_dmc(iw,ifr),psijo_dmc(iw,ifr),eold(iw,ifr),0,ifr)
          i_vpsp=0
          do 65 i=1,nelec
   65       call compute_determinante_grad(i,psido_dmc(iw,ifr),psido_dmc(iw,ifr),vold_dmc(1,i,iw,ifr),1)
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

c zero out xsum variables for metrop

      wsum_dmc=zero
      wfsum=zero
      wdsum=zero
      wgdsum=zero
      esum_dmc=zero
      efsum=zero
      ei1sum=zero
      ei2sum=zero
      r2sum=zero
      risum=zero

      do 80 ifr=1,nforce
        egsum(ifr)=zero
        wgsum(ifr)=zero
        pesum_dmc(ifr)=zero
        tpbsum_dmc(ifr)=zero
        tjfsum_dmc(ifr)=zero
        tausum(ifr)=zero
        do 80 k=1,3
   80     derivsum(k,ifr)=zero

      call prop_init(1)
      call pcm_init(1)
      call mmpol_init(1)

      if(ipr.ge.-2) then
        if(idtask.le.9) then
          write(filename,'(''walkalize.'',i1)') idtask
         elseif(idtask.le.99) then
          write(filename,'(''walkalize.'',i2)') idtask
         elseif(idtask.le.999) then
          write(filename,'(''walkalize.'',i3)') idtask
         else
          call fatal_error('STARTR: idtask > 999')
        endif
        open(unit=11,file=filename,status='old')
        do 90 i=1,2000000000
   90     read(11,fmt=*,end=100)
      endif
  100 backspace 11
      backspace 11

      return
      end
