      subroutine dumper_more
c Written by Cyrus Umrigar, modified by Claudia Filippi
c routine to pick up and dump everything needed to restart
c job where it left off
      use vmc_mod, only: MELEC, MORB, MBASIS, MDET, MCENT, MCTYPE
      use vmc_mod, only: nrad
      use atom, only: znuc, cent, pecent, iwctype, nctype, ncent
      use mstates_mod, only: MSTATES
      use ghostatom, only: newghostype, nghostcent
      use const, only: hb, delta, nelec
      use config, only: eold, nearesto, psi2o
      use config, only: psido, psijo, rmino, rvmino, tjfo
      use config, only: vold, xnew, xold
      use jaspar1, only: cjas1, cjas2
      use csfs, only: nstates
      use denupdn, only: rprobdn, rprobup
      use dets, only: cdet, ndet
      use elec, only: ndn, nup
      use est2cm, only: ecm2, ecm21, pecm2, r2cm2, tjfcm2, tpbcm2
      use estcum, only: ecum, ecum1, iblk, pecum, r2cum, tjfcum, tpbcum
      use estsig, only: ecm21s, ecum1s
      use estsum, only: acc, esum, pesum, r2sum, tjfsum, tpbsum
      use forcepar, only: nforce
      use forcest, only: fcm2, fcum
      use forcewt, only: wcum, wsum
      use optwf_contrl, only: ioptorb
      use stats, only: rejmax
      use step, only: ekin, ekin2, rprob, suc, trunfb, try
      use wfsec, only: iwftype, nwftype
      use coefs, only: coef, nbasis, norb
      use const2, only: deltar, deltat
      use contrl, only: nstep
      use basis, only: zex, n1s, n2s, n2p, n3s, n3p, n3dzr, n3dx2, n3dxy, n3dxz, n3dyz
      use basis, only: n4s, n4p
      use basis, only: nsa, npa, ndzra, ndxya, ndxza, ndyza, ndx2a
      use mstates_ctrl, only: iguiding
      use inputflags, only: node_cutoff, eps_node_cutoff

      ! I'm 50% sure it's needed
      ! it was in master as part of the include optorb.h 
      use optorb_cblock, only: ns_current
      implicit real*8(a-h,o-z)

      parameter(half=0.5d0,small=1.d-6)

      dimension coefx(MBASIS,MORB),zexx(MBASIS),centx(3,MCENT)
     &,znucx(MCTYPE),n1sx(MCTYPE),n2sx(MCTYPE),n2px(3,MCTYPE)
     &,n3sx(MCTYPE),n3px(3,MCTYPE),n3dzrx(MCTYPE),n3dx2x(MCTYPE)
     &,n3dxyx(MCTYPE),n3dxzx(MCTYPE),n3dyzx(MCTYPE),n4sx(MCTYPE)
     &,n4px(3,MCTYPE),nsax(MCTYPE),npax(3,MCTYPE),ndzrax(MCTYPE)
     &,ndx2ax(MCTYPE),ndxyax(MCTYPE),ndxzax(MCTYPE),ndyzax(MCTYPE)
     &,cdetx(MDET)
      dimension xstrech(3,MELEC)
      dimension d2(MSTATES)

c     RLPB
      kstate=1

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
      write(10) ((coef(ib,i,kstate,1),ib=1,nbasis),i=1,norb)
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

c     RLPB
      kstate=1

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
   10     if (dabs(coefx(i,j)-coef(i,j,kstate,1)).gt.small) call fatal_error('STARTR: coef')
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
