      subroutine dumper_more
c Written by Cyrus Umrigar, modified by Claudia Filippi
c routine to pick up and dump everything needed to restart
c job where it left off
      use atom, only: znuc, cent, pecent, iwctype, nctype, ncent
      use ghostatom, only: newghostype, nghostcent
      use const, only: pi, hb, etrial, delta, deltai, fbias, nelec, imetro, ipr
      use config, only: delttn, enew, eold, nearestn, nearesto, pen, peo, psi2n, psi2o,
     &psido, psijo, rminn, rminno, rmino, rminon, rvminn, rvminno, rvmino, rvminon, tjfn, tjfo,
     &vnew, vold, xnew, xold
      use jaspar1, only: cjas1, cjas2
      implicit real*8(a-h,o-z)




      include 'vmc.h'
      include 'force.h'
      include 'mstates.h'
      include 'pseudo.h'
      include 'basis.h'
      include 'optorb.h'
      include 'optci.h'

      parameter(half=0.5d0,small=1.d-6)

      common /const2/ deltar,deltat
      common /contrl/ nstep,nblk,nblkeq,nconf,nconf_new,isite,idump,irstar
      common /coefs/ coef(MBASIS,MORB,MWF),nbasis,norb
      common /dets/ cdet(MDET,MSTATES,MWF),ndet
      common /elec/ nup,ndn

      common /estsum/ esum1(MSTATES),esum(MSTATES,MFORCE),pesum(MSTATES),tpbsum(MSTATES),tjfsum(MSTATES),r2sum,acc
      common /estcum/ ecum1(MSTATES),ecum(MSTATES,MFORCE),pecum(MSTATES),tpbcum(MSTATES),tjfcum(MSTATES),r2cum,iblk
      common /est2cm/ ecm21(MSTATES),ecm2(MSTATES,MFORCE),pecm2(MSTATES),tpbcm2(MSTATES),tjfcm2(MSTATES),r2cm2
      common /estsig/ ecum1s(MSTATES),ecm21s(MSTATES)

      common /forcepar/ deltot(MFORCE),nforce,istrech
      common /forcest/ fcum(MSTATES,MFORCE),fcm2(MSTATES,MFORCE)
      common /forcewt/ wsum(MSTATES,MFORCE),wcum(MSTATES,MFORCE)

      common /stats/ rejmax
      common /step/try(nrad),suc(nrad),trunfb(nrad),rprob(nrad),
     &ekin(nrad),ekin2(nrad)
      common /denupdn/ rprobup(nrad),rprobdn(nrad)
      common /wfsec/ iwftype(MFORCE),iwf,nwftype

      common /csfs/ ccsf(MDET,MSTATES,MWF),cxdet(MDET*MDETCSFX)
     &,icxdet(MDET*MDETCSFX),iadet(MDET),ibdet(MDET),ncsf,nstates

      common /optwf_contrl/ ioptjas,ioptorb,ioptci,nparm

      dimension coefx(MBASIS,MORB),zexx(MBASIS),centx(3,MCENT)
     &,znucx(MCTYPE),n1sx(MCTYPE),n2sx(MCTYPE),n2px(3,MCTYPE)
     &,n3sx(MCTYPE),n3px(3,MCTYPE),n3dzrx(MCTYPE),n3dx2x(MCTYPE)
     &,n3dxyx(MCTYPE),n3dxzx(MCTYPE),n3dyzx(MCTYPE),n4sx(MCTYPE)
     &,n4px(3,MCTYPE),nsax(MCTYPE),npax(3,MCTYPE),ndzrax(MCTYPE)
     &,ndx2ax(MCTYPE),ndxyax(MCTYPE),ndxzax(MCTYPE),ndyzax(MCTYPE)
     &,cdetx(MDET)
      dimension xstrech(3,MELEC)
      dimension d2(MSTATES)

      write(10) delta,deltar,deltat

      write(10) nstep,iblk
      do 1 istate=1,nstates
        write(10) ecum1(istate),(ecum(istate,i),i=1,nforce),pecum(istate),tpbcum(istate),tjfcum(istate),r2cum,acc
        write(10) ecm21(istate),(ecm2(istate,i),i=1,nforce),pecm2(istate),tpbcm2(istate),tjfcm2(istate),r2cm2
        if(nforce.gt.1) then
          write(10) (wcum(istate,i),fcum(istate,i),fcm2(istate,i),i=1,nforce)
         else
          write(10) wcum(istate,1)
        endif
  1     write(10) ecum1s(istate),ecm21s(istate)
      write(10) (try(i),suc(i),trunfb(i),rprob(i),
     &rprobup(i),rprobdn(i),ekin(i),ekin2(i),i=1,nrad)
      call optorb_dump(10)
      call optci_dump(10)
      call prop_dump(10)
      call efficiency_dump(10)
      call pcm_dump(10)
      call force_analy_dump(10)
      write(10) rejmax

      write(10) nbasis,norb
      write(10) ((coef(ib,i,1),ib=1,nbasis),i=1,norb)
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

      call optjas_dump(10)
      call optx_jas_orb_dump(10)
      call optx_jas_ci_dump(10)
      call optx_orb_ci_dump(10)

      rewind 10

      write(6,'(1x,''successful dump to unit 10'')')

      return

c-----------------------------------------------------------------------
      entry startr_more

      read(10) deltax,deltarx,deltatx
      if (dabs(deltax-delta).gt.small) call fatal_error('STARTR: delta')
      if (dabs(deltarx-deltar).gt.small) call fatal_error('STARTR: deltar')
      if (dabs(deltatx-deltat).gt.small) call fatal_error('STARTR: deltat')

      read(10) nstepx,iblk
      if (nstepx.ne.nstep) call fatal_error('STARTR: nstep')
      do 2 istate=1,nstates
        read(10) ecum1(istate),(ecum(istate,i),i=1,nforce),pecum(istate),tpbcum(istate),tjfcum(istate),r2cum,acc
        read(10) ecm21(istate),(ecm2(istate,i),i=1,nforce),pecm2(istate),tpbcm2(istate),tjfcm2(istate),r2cm2
        if(nforce.gt.1) then
          read(10) (wcum(istate,i),fcum(istate,i),fcm2(istate,i),i=1,nforce)
         else
          read(10) wcum(istate,1)
        endif
  2     read(10) ecum1s(istate),ecm21s(istate)
      read(10) (try(i),suc(i),trunfb(i),rprob(i),
     &rprobup(i),rprobdn(i),ekin(i),ekin2(i),i=1,nrad)

      call prop_rstrt(10)
      call optorb_rstrt(10)
      call optci_rstrt(10)
      call efficiency_rstrt(10)
      call pcm_rstrt(10)
      call force_analy_rstrt(10)
      read(10) rejmax

      read(10) nbasx,norbx
      if (nbasx.ne.nbasis) call fatal_error('STARTR: nbasis')
      if (norbx.ne.norb) call fatal_error('STARTR: norb')
      read(10) ((coefx(ib,i),ib=1,nbasis),i=1,norb)
      read(10) (zexx(ib),ib=1,nbasis)
      read(10) nctypex,ncentx,newghostypex,nghostcentx,(iwctype(i),i=1,ncentx+nghostcentx)
      if (ncentx.ne.ncent) call fatal_error('STARTR: ncent')
      if (nctypex.ne.nctype) call fatal_error('STARTR: nctype')
      read(10) ((centx(k,ic),k=1,3),ic=1,ncentx+nghostcentx)
      read(10) pecx
      read(10) (znucx(i),i=1,nctype)
      read(10) (n1sx(i),i=1,nctype)
      read(10) (n2sx(i),i=1,nctype)
      read(10) ((n2px(k,i),k=1,3),i=1,nctype)
      read(10) (n3sx(i),i=1,nctype)
      read(10) ((n3px(k,i),k=1,3),i=1,nctype)
      read(10) (n3dzrx(i),i=1,nctype)
      read(10) (n3dx2x(i),i=1,nctype)
      read(10) (n3dxyx(i),i=1,nctype)
      read(10) (n3dxzx(i),i=1,nctype)
      read(10) (n3dyzx(i),i=1,nctype)
      read(10) (n4sx(i),i=1,nctype)
      read(10) ((n4px(k,i),k=1,3),i=1,nctype)
      read(10) (nsax(i),i=1,nctype)
      read(10) ((npax(k,i),k=1,3),i=1,nctype)
      read(10) (ndzrax(i),i=1,nctype)
      read(10) (ndx2ax(i),i=1,nctype)
      read(10) (ndxyax(i),i=1,nctype)
      read(10) (ndxzax(i),i=1,nctype)
      read(10) (ndyzax(i),i=1,nctype)
      do 10 j=1,norb
        do 10 i=1,nbasis
   10     if (dabs(coefx(i,j)-coef(i,j,1)).gt.small) call fatal_error('STARTR: coef')
      do 20 i=1,nbasis
        if (dabs(zexx(i)-zex(i,1)).gt.small) call fatal_error('STARTR: zex')
   20 continue
      do 30 i=1,ncent+nghostcent
        do 30 k=1,3
          if (dabs(cent(k,i)-centx(k,i)).gt.small) call fatal_error('STARTR: cent')
   30 continue
      if (pecx.ne.pecent) call fatal_error('STARTR: pec')
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

      call optjas_rstrt(10)
      call optx_jas_orb_rstrt(10)
      call optx_jas_ci_rstrt(10)
      call optx_orb_ci_rstrt(10)

      write(6,'(1x,''succesful read from unit 10'')')
      write(6,'(t5,''enow'',t15,''eave'',t25,''eerr'',t35,''peave'',
     &t45,''peerr'',t55,''tpbave'',t65,''tpberr'',t75,''tjfave'',
     &t85,''tjferr'',t95,''accept'',t105,''iter'')')

      if(nforce.gt.1) then
        call setup_force
       else
        nwftype=1
        iwftype(1)=1
      endif

c loop over secondary config
      do 60 ifr=2,nforce
c set n- and e-coord and n-n potential
        call strech(xold,xstrech,ajacob,ifr,1)
        call hpsi(xstrech,psido,psijo,eold(1,ifr),0,ifr)
        do 60 istate=1,nforce
   60     psi2o(istate,ifr)=2*(dlog(dabs(psido(istate)))+psijo)+dlog(ajacob)

c primary config
c set n-coord and n-n potential
      if(nforce.gt.1) call strech(xold,xstrech,ajacob,1,0)
      call hpsi(xold,psido,psijo,eold(1,1),0,1)
      do 65 istate=1,nforce
        psi2o(istate,1)=2*(dlog(dabs(psido(istate)))+psijo)
        tjfo(istate)=d2(istate)
   65   tjfo(istate)=-tjfo(istate)*half*hb

      if(iguiding.gt.0) then
        call determinant_psig(psido,psidg)
c rewrite psi2o if you are sampling guiding
        psi2o(1,1)=2*(dlog(dabs(psidg))+psijo)
      endif

      if(node_cutoff.gt.0) then
        do 83 jel=1,nelec
   83     call compute_determinante_grad(jel,psido,psido,vold(1,jel),1)
        call nodes_distance(vold,distance_node,1)
        rnorm_nodes=rnorm_nodes_num(distance_node,eps_node_cutoff)/distance_node

        psi2o(1,1)=psi2o(1,1)+2*dlog(rnorm_nodes)
      endif

      if(ioptorb.gt.0) ns_current=0

      call prop_save
      call optjas_save
      call optci_save
      call optorb_save

      do 70 i=1,nelec
        do 70 k=1,3
   70     xnew(k,i)=xold(k,i)

      do 100 i=1,nelec
        rmino(i)=99.d9
        do 90 j=1,ncent
          dist=0
          do 80 k=1,3
   80       dist=dist+(xold(k,i)-cent(k,j))**2
          if(dist.lt.rmino(i)) then
            rmino(i)=dist
            nearesto(i)=j
          endif
   90   continue
        rmino(i)=dsqrt(rmino(i))
        do 100 k=1,3
  100     rvmino(k,i)=xold(k,i)-cent(k,nearesto(i))

      do 130 istate=1,nstates
        do 110 ifr=1,nforce
          esum(istate,ifr)=0
  110     wsum(istate,ifr)=0
        pesum(istate)=0
        tpbsum(istate)=0
  130   tjfsum(istate)=0
      r2sum=0

      return
      end
