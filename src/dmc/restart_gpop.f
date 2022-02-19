      module restart_gpop
      contains
      subroutine startr_gpop

      use vmc_mod, only: norb_tot
      use vmc_mod, only: nrad
      use basis, only: zex, n1s, n2s, n2p, n3s, n3p, n3dzr, n3dx2, n3dxy, n3dxz, n3dyz
      use basis, only: n4s, n4p
      use basis, only: nsa, npa, ndzra, ndxya, ndxza, ndyza, ndx2a
      use const, only: hb, ipr, nelec
      use forcest, only: fgcm2, fgcum
      use forcepar, only: istrech, nforce
      use age, only: iage, ioldest, ioldestmx
      use contrldmc, only: idmc
      use contrldmc, only: nfprod, rttau, tau
      use atom, only: cent, iwctype, ncent, nctype, pecent, znuc
      use atom, only: ncent_tot
      use estcum, only: iblk, ipass
      use config, only: psido_dmc, psijo_dmc, vold_dmc, xold_dmc
      use stats, only: acc, dfus2ac, dfus2un, dr2ac, dr2un, nacc, nbrnch, nodecr, trymove
      use estsum, only: efsum, egsum, ei1sum, ei2sum, esum_dmc
      use estsum, only: pesum_dmc, r2sum, risum, tausum, tjfsum_dmc, tpbsum_dmc, wdsum
      use estsum, only: wfsum, wgdsum, wgsum, wsum_dmc
      use estcum, only: ecum1_dmc, ecum_dmc, efcum, efcum1, egcum, egcum1, ei1cum, ei2cum
      use estcum, only: ei3cum, pecum_dmc, r2cum_dmc, ricum, taucum, tjfcum_dmc, tpbcum_dmc
      use estcum, only: wcum1, wcum_dmc, wdcum, wdcum1, wfcum, wfcum1, wgcum, wgcum1
      use estcum, only: wgdcum
      use est2cm, only: ecm21_dmc, ecm2_dmc, efcm2, efcm21, egcm2, egcm21, ei1cm2, ei2cm2
      use est2cm, only: ei3cm2, pecm2_dmc, r2cm2_dmc, ricm2, tjfcm_dmc, tpbcm2_dmc, wcm2, wcm21, wdcm2, wdcm21
      use est2cm, only: wfcm2, wfcm21, wgcm2, wgcm21, wgdcm2
      use derivest, only: derivcum, derivsum
      use step, only: rprob
      use mpiconf, only: idtask, nproc, wid, NPROCX
      use denupdn, only: rprobdn, rprobup
      use qua, only: nquad, wq, xq, yq, zq
      use branch, only: eest, eigv, eold, ff, fprod, nwalk, wdsumo, wgdsumo, wt, wtgen
      use casula, only: i_vpsp, icasula
      use jacobsave, only: ajacob, ajacold
      use pseudo, only: nloc
      use dets, only: cdet, ndet
      use elec, only: ndn, nup
      use coefs, only: coef, nbasis, norb
      use ghostatom, only: nghostcent
      use velratio, only: fratio
!      use contrl, only: nconf
      use control_dmc, only: dmc_nconf
      use mpi
      use contrl_file,    only: ounit
      use precision_kinds, only: dp

      use error,           only: fatal_error
      use mmpol,           only: mmpol_init
      use mmpol_dmc,       only: mmpol_save
      use pcm_dmc,         only: pcm_save
      use prop_dmc,        only: prop_save_dmc
      use pcm_mod,         only: pcm_init
      use properties_mod,  only: prop_init
      use nonloc_grid_mod, only: t_vpsp_sav
      use rannyu_mod,      only: setrn
      use strech_mod,      only: strech
      use hpsi_mod,        only: hpsi
      use determinante_mod,only: compute_determinante_grad
      use walksav_det_mod, only: walksav_det
      use walksav_jas_mod, only: walksav_jas
      implicit none

      integer :: i, iage_id, ib, ic, id
      integer :: ie, ifr, iw, j
      integer :: k, nbasx, ncentx, nctypex
      integer :: ndetx, ndnx, nelecx, newghostypex
      integer :: nf_id, nghostcentx, nprock, nq_id
      integer :: num, nupx, nwalk_id
      integer, dimension(4, 0:NPROCX) :: irn
      integer, dimension(ncent_tot) :: n1sx
      integer, dimension(ncent_tot) :: n2sx
      integer, dimension(3, ncent_tot) :: n2px
      integer, dimension(ncent_tot) :: n3sx
      integer, dimension(3, ncent_tot) :: n3px
      integer, dimension(ncent_tot) :: n3dzrx
      integer, dimension(ncent_tot) :: n3dx2x
      integer, dimension(ncent_tot) :: n3dxyx
      integer, dimension(ncent_tot) :: n3dxzx
      integer, dimension(ncent_tot) :: n3dyzx
      integer, dimension(ncent_tot) :: n4sx
      integer, dimension(3, ncent_tot) :: n4px
      integer, dimension(ncent_tot) :: nsax
      integer, dimension(3, ncent_tot) :: npax
      integer, dimension(ncent_tot) :: ndzrax
      integer, dimension(ncent_tot) :: ndx2ax
      integer, dimension(ncent_tot) :: ndxyax
      integer, dimension(ncent_tot) :: ndxzax
      integer, dimension(ncent_tot) :: ndyzax
      real(dp) :: different, fmt
      real(dp) :: fratio_id, hbx, taux, wq_id
      real(dp) :: wt_id, xold_dmc_id, xq_id, yq_id
      real(dp) :: zq_id
      real(dp), dimension(nbasis, norb_tot) :: coefx
      real(dp), dimension(nbasis) :: zexx
      real(dp), dimension(3, ncent_tot) :: centx
      real(dp), dimension(ncent_tot) :: znucx
      real(dp), dimension(ndet) :: cdetx
      real(dp), parameter :: zero = 0.d0
      real(dp), parameter :: one = 1.d0
      real(dp), parameter :: small = 1.e-6




      character*13 filename



      write(ounit,'(1x,''attempting restart from unit 10'')')
      rewind 10
      read(10) nprock
      read(10) nfprod,(ff(i),i=0,nfprod),fprod,eigv,eest,wdsumo
     &,ioldest,ioldestmx
      if(nprock.ne.nproc) call fatal_error('STARTR: different num procs')
      do id=0,idtask
        read(10) nwalk
        read(10) (wt(i),i=1,nwalk),(iage(i),i=1,nwalk)
        read(10) (((xold_dmc(ic,i,iw,1),ic=1,3),i=1,nelec),iw=1,nwalk)
        read(10) nforce,((fratio(iw,ifr),iw=1,nwalk),ifr=1,nforce)
        if(nloc.gt.0)
     &  read(10) nquad,(xq(i),yq(i),zq(i),wq(i),i=1,nquad)
      enddo
      do id=idtask+1,nproc-1
        read(10) nwalk_id
        read(10) (wt_id,i=1,nwalk_id),(iage_id,i=1,nwalk_id)
        read(10) (xold_dmc_id,i=1,3*nelec*nwalk_id)
        read(10) nf_id,((fratio_id,iw=1,nwalk_id),ifr=1,nf_id)
        if(nloc.gt.0)
     &  read(10) nq_id,(xq_id,yq_id,zq_id,wq_id,i=1,nquad)
      enddo
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
      read(10) nelecx,dmc_nconf
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
      do j=1,norb
      do i=1,nbasis
      if (dabs(coefx(i,j)-coef(i,j,1)).gt.small) call fatal_error('STARTR: coef')
      enddo
      enddo
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
      do i=1,nbasis
      if (dabs(zexx(i)-zex(i,1)).gt.small) call fatal_error('STARTR: zex')
      enddo
      do i=1,ncent+nghostcent
      do k=1,3
      if (dabs(cent(k,i)-centx(k,i)).gt.small) call fatal_error('STARTR: cent')
      enddo
      enddo
      do i=1,nctype
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
      do ic=1,3
      if (n2p(ic,i).ne.n2px(ic,i)) call fatal_error('STARTR: n2p')
      if (n3p(ic,i).ne.n3px(ic,i)) call fatal_error('STARTR: n3p')
      if (n4p(ic,i).ne.n4px(ic,i)) call fatal_error('STARTR: n4p')
      if (npa(ic,i).ne.npax(ic,i)) call fatal_error('STARTR: npa')
      enddo
      enddo
      read(10) (cdetx(i),i=1,ndet)
      read(10) ndetx,nupx,ndnx
      do i=1,ndet
      if (dabs(cdetx(i)-cdet(i,1,1)).gt.small) call fatal_error('STARTR: cdet')
      enddo
      if (ndetx.ne.ndet) call fatal_error('STARTR: ndet')
      if (nupx.ne.nup) call fatal_error('STARTR: nup')
      if (ndnx.ne.ndn) call fatal_error('STARTR: ndn')
      write(ounit,'(1x,''succesful read from unit 10'')')
      write(ounit,'(t5,''egnow'',t15,''egave'',t21
     &,''(egerr)'' ,t32,''peave'',t38,''(peerr)'',t49,''tpbave'',t55
     &,''(tpberr)'' ,t66,''tjfave'',t72,''(tjferr)'',t83,''npass'',t93
     &,''wgsum'',t103 ,''ioldest'')')

      do iw=1,nwalk
        if(istrech.eq.0) then
          do ifr=2,nforce
            do ie=1,nelec
              do k=1,3
                xold_dmc(k,ie,iw,ifr)=xold_dmc(k,ie,iw,1)
              enddo
            enddo
          enddo
        endif
        do ifr=1,nforce
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
          do i=1,nelec
            call compute_determinante_grad(i,psido_dmc(iw,ifr),psido_dmc(iw,ifr),vold_dmc(1,i,iw,ifr),1)
          enddo
          if(ifr.eq.1) then
            call walksav_det(iw)
            call walksav_jas(iw)
c           call t_vpsp_sav(iw)
            call t_vpsp_sav
            call prop_save_dmc(iw)
            call pcm_save(iw)
            call mmpol_save(iw)
          endif
        enddo
      enddo

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

      do ifr=1,nforce
        egsum(ifr)=zero
        wgsum(ifr)=zero
        pesum_dmc(ifr)=zero
        tpbsum_dmc(ifr)=zero
        tjfsum_dmc(ifr)=zero
        tausum(ifr)=zero
        do k=1,3
          derivsum(k,ifr)=zero
        enddo
      enddo

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
        do i=1,2000000000
          read(11,fmt=*,end=100)
        enddo
      endif
  100 backspace 11
      backspace 11

      return
      end
      end module
