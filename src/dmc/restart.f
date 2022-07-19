      module restart
      contains
      subroutine startr

      use vmc_mod, only: norb_tot
      use vmc_mod, only: nrad
      use basis, only: zex
      use basis, only: ns, npx, npy, npz, ndxx, ndxy, ndxz, ndyy, ndyz, ndzz
      use basis, only: nfxxx, nfxxy, nfxxz, nfxyy, nfxyz, nfxzz, nfyyy, nfyyz, nfyzz, nfzzz
      use basis, only: ngxxxx, ngxxxy, ngxxxz, ngxxyy, ngxxyz, ngxxzz, ngxyyy, ngxyyz
      use basis, only: ngxyzz, ngxzzz, ngyyyy, ngyyyz, ngyyzz, ngyzzz, ngzzzz
      use const, only: hb
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
      use derivest, only: derivcm2, derivcum, derivsum, derivtotave_num_old
      use step, only: rprob
      use mpiconf, only: idtask, nproc, wid
      use denupdn, only: rprobdn, rprobup
      use mpiblk, only: iblk_proc
      use qua, only: nquad, wq, xq, yq, zq
      use branch, only: eest, eigv, eold, ff, fprod, nwalk, wdsumo, wgdsumo, wt, wtgen
      use casula, only: i_vpsp, icasula
      use jacobsave, only: ajacob, ajacold
      use pseudo, only: nloc
      use dets, only: cdet, ndet
      use coefs, only: coef, nbasis, norb
      use ghostatom, only: nghostcent
      use velratio, only: fratio
!      use contrl, only: nconf
      use control_dmc, only: dmc_nconf
      use mpi
      use contrl_file, only: ounit
      use precision_kinds, only: dp
      
      use restart_gpop, only: startr_gpop
      use error, only: fatal_error
      use rannyu_mod, only: setrn
      use mmpol, only: mmpol_init, mmpol_rstrt
      use pcm_mod, only: pcm_init, pcm_rstrt
      use properties_mod, only: prop_init, prop_rstrt
      use mmpol_dmc, only: mmpol_save
      use pcm_dmc, only: pcm_save
      use prop_dmc, only: prop_save_dmc
      use nonloc_grid_mod, only: t_vpsp_sav
      use walksav_det_mod, only: walksav_det
      use walksav_jas_mod, only: walksav_jas
      use determinante_mod, only: compute_determinante_grad
      use hpsi_mod, only: hpsi
      use strech_mod, only: strech
      use control, only: ipr
      use control, only: mode
      use system, only: nelec
      use system, only: nup
      use system, only: ndn
      implicit none

      integer :: i, iage_id, ib, ic, id
      integer :: ie, ifr, ioldest_id, ioldestmx_id
      integer :: iw, j, k, n1_id
      integer :: n2_id, nbasx, ncentx, nctypex
      integer :: ndetx, ndnx, nelecx, newghostypex
      integer :: nghostcentx, nprock, nq_id, num
      integer :: nupx, nwalk_id
      integer, dimension(4, 0:nproc) :: irn
      integer, dimension(ncent_tot) :: nsx
      integer, dimension(ncent_tot) :: npxx
      integer, dimension(ncent_tot) :: npyx
      integer, dimension(ncent_tot) :: npzx
      integer, dimension(ncent_tot) :: ndxxx
      integer, dimension(ncent_tot) :: ndxyx
      integer, dimension(ncent_tot) :: ndxzx
      integer, dimension(ncent_tot) :: ndyyx
      integer, dimension(ncent_tot) :: ndyzx
      integer, dimension(ncent_tot) :: ndzzx
      real(dp) :: different, eest_id
      real(dp) :: eigv_id, ff_id, fmt, fprod_id
      real(dp) :: fratio_id, hbx, taux, wdsumo_id
      real(dp) :: wq_id, wt_id, xold_dmc_id, xq_id
      real(dp) :: yq_id, zq_id
      real(dp), dimension(nbasis, norb_tot) :: coefx
      real(dp), dimension(nbasis) :: zexx
      real(dp), dimension(3, ncent_tot) :: centx
      real(dp), dimension(ncent_tot) :: znucx
      real(dp), dimension(ndet) :: cdetx
      real(dp), parameter :: zero = 0.d0
      real(dp), parameter :: one = 1.d0
      real(dp), parameter :: small = 1.e-6





      character*13 filename



      if(mode.eq.'dmc_one_mpi2') then
        call startr_gpop
        return
      endif

      write(ounit,'(1x,''attempting restart from unit 10'')')
      rewind 10
      read(10) nprock
      if(nprock.ne.nproc) call fatal_error('STARTR: different num procs')
      do id=0,idtask
        read(10) nwalk
        read(10) (((xold_dmc(ic,i,iw,1),ic=1,3),i=1,nelec),iw=1,nwalk)
        read(10) nfprod,(ff(i),i=0,nfprod),(wt(i),i=1,nwalk),fprod
     &  ,eigv,eest,wdsumo
        read(10) (iage(i),i=1,nwalk),ioldest,ioldestmx
        read(10) nforce,((fratio(iw,ifr),iw=1,nwalk),ifr=1,nforce)
c       read(10) (wgcum(i),egcum(i),pecum_dmc(i),tpbcum_dmc(i),tjfcum_dmc(i),
c    &  wgcm2(i),egcm2(i),pecm2_dmc(i),tpbcm2_dmc(i),tjfcm_dmc(i),taucum(i),
c    &  i=1,nforce)
        if(nloc.gt.0)
     &  read(10) nquad,(xq(i),yq(i),zq(i),wq(i),i=1,nquad)
      enddo
      do id=idtask+1,nproc-1
        read(10) nwalk_id
        read(10) (xold_dmc_id,i=1,3*nelec*nwalk_id)
        read(10) n1_id,(ff_id,i=0,n1_id),(wt_id,i=1,nwalk_id),fprod_id
     &  ,eigv_id,eest_id,wdsumo_id
        read(10) (iage_id,i=1,nwalk_id),ioldest_id,ioldestmx_id
        read(10) n2_id,((fratio_id,iw=1,nwalk_id),ifr=1,n2_id)
c       read(10) (wgcum_id,egcum_id,pecum_dmc_id,tpbcum_dmc_id,tjfcum_dmc_id,
c    &  wgcm2_id,egcm2_id,pecm2_dmc_id,tpbcm2_dmc_id,tjfcm_dmc_id,taucum_id,
c    &  i=1,nforce)
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
      read(10) ipass,iblk,iblk_proc
      read(10) wcm2,wfcm2,wdcm2,wgdcm2,wcm21
     &,wfcm21,(wgcm21(i),i=1,nforce),wdcm21, ecm2_dmc,efcm2
     &,ecm21_dmc,efcm21,(egcm21(i),i=1,nforce)
     &,ei1cm2,ei2cm2,ei3cm2,r2cm2_dmc,ricm2
      read(10) (fgcum(i),i=1,nforce),(fgcm2(i),i=1,nforce)
     &,((derivcum(k,i),k=1,3),i=1,nforce),(derivcm2(i),i=1,nforce)
     &,(derivtotave_num_old(i),i=1,nforce)
      read(10) (rprob(i),rprobup(i),rprobdn(i),i=1,nrad)
      read(10) dfus2ac,dfus2un,dr2ac,dr2un,acc
     &,trymove,nacc,nbrnch,nodecr
      if(.not.wid) then
        acc=0
        nacc=0
        trymove=0
        nodecr=0
      endif
      call prop_rstrt(10)
      call pcm_rstrt(10)
      call mmpol_rstrt(10)
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
      read(10) (nsx(i),i=1,nctype)
      read(10) (npxx(i),i=1,nctype)
      read(10) (npyx(i),i=1,nctype)
      read(10) (npzx(i),i=1,nctype)
      read(10) (ndxxx(i),i=1,nctype)
      read(10) (ndxyx(i),i=1,nctype)
      read(10) (ndxzx(i),i=1,nctype)
      read(10) (ndyyx(i),i=1,nctype)
      read(10) (ndyzx(i),i=1,nctype)
      read(10) (ndzzx(i),i=1,nctype)

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
      if (ns(i).ne.nsx(i)) call fatal_error('STARTR: ns')
      if (npx(i).ne.npxx(i)) call fatal_error('STARTR: npx')
      if (npy(i).ne.npyx(i)) call fatal_error('STARTR: npy')
      if (npz(i).ne.npzx(i)) call fatal_error('STARTR: npz')
      if (ndxx(i).ne.ndxxx(i)) call fatal_error('STARTR: ndxx')
      if (ndxy(i).ne.ndxyx(i)) call fatal_error('STARTR: ndxy')
      if (ndxz(i).ne.ndxzx(i)) call fatal_error('STARTR: ndxz')
      if (ndyy(i).ne.ndyyx(i)) call fatal_error('STARTR: ndyy')
      if (ndyz(i).ne.ndyzx(i)) call fatal_error('STARTR: ndyz')
      if (ndzz(i).ne.ndzzx(i)) call fatal_error('STARTR: ndzz')
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
