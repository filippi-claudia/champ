      module dumper_more_mod
      contains
      subroutine dumper_more
c Written by Cyrus Umrigar, modified by Claudia Filippi
c routine to pick up and dump everything needed to restart
c job where it left off
      use vmc_mod, only: norb_tot
      use vmc_mod, only: nrad
      use system, only: znuc, cent, pecent, iwctype, nctype, ncent, ncent_tot, nctype_tot
      use mstates_mod, only: MSTATES
      use system, only: newghostype, nghostcent
      use const, only: hb, delta
      use config, only: eold, nearesto, psi2o
      use config, only: psido, psijo, rmino, rvmino, tjfo
      use config, only: vold, xnew, xold
      use csfs, only: nstates
      use denupdn, only: rprobdn, rprobup
      use dets, only: cdet, ndet
      use est2cm, only: ecm2, ecm21, pecm2, r2cm2, tjfcm2, tpbcm2
      use estcum, only: ecum, ecum1, iblk, pecum, r2cum, tjfcum, tpbcum
      use estsig, only: ecm21s, ecum1s
      use estsum, only: acc, esum, pesum, r2sum, tjfsum, tpbsum
      use forcest, only: fcm2, fcum
      use forcewt, only: wcum, wsum
      use optwf_contrl, only: ioptorb
      use stats, only: rejmax
      use step, only: ekin, ekin2, rprob, suc, trunfb, try
      use wfsec, only: iwftype, nwftype
      use coefs, only: coef, nbasis, norb
      use const2, only: deltar, deltat
!      use contrl, only: nstep
      use control_vmc, only: vmc_nstep
      use basis, only: zex
      use basis, only: ns, npx, npy, npz, ndxx, ndxy, ndxz, ndyy, ndyz, ndzz
      use basis, only: nfxxx, nfxxy, nfxxz, nfxyy, nfxyz, nfxzz, nfyyy, nfyyz, nfyzz, nfzzz
      use basis, only: ngxxxx, ngxxxy, ngxxxz, ngxxyy, ngxxyz, ngxxzz, ngxyyy, ngxyyz
      use basis, only: ngxyzz, ngxzzz, ngyyyy, ngyyyz, ngyyzz, ngyzzz, ngzzzz
      use mstates_ctrl, only: iguiding
      use inputflags, only: node_cutoff, eps_node_cutoff
      use contrl_file, only: ounit, errunit
      ! I'm 50% sure it's needed
      ! it was in master as part of the include optorb.h
      use optorb_cblock, only: ns_current
      use precision_kinds, only: dp

      use optorb_f_mod, only: optorb_save, optorb_rstrt
      use optci_mod, only: optci_save, optci_rstrt
      use optjas_mod, only: optjas_dump, optjas_rstrt, optjas_save
      use prop_vmc, only: prop_save
      use pcm_mod, only: pcm_rstrt, pcm_dump
      use nodes_distance_mod, only: nodes_distance
      use determinante_mod, only: compute_determinante_grad
      use error, only: fatal_error
      use optorb_f_mod, only: optorb_dump
      use optci_mod, only: optci_dump
      use nodes_distance_mod, only: rnorm_nodes_num
      use determinant_psig_mod, only: determinant_psig
      use hpsi_mod, only: hpsi
      use strech_mod, only: strech, setup_force
      use optx_orb_ci, only: optx_orb_ci_rstrt, optx_orb_ci_dump
      use optx_jas_ci, only: optx_jas_ci_rstrt, optx_jas_ci_dump
      use optx_jas_orb, only: optx_jas_orb_rstrt, optx_jas_orb_dump
      use force_analytic, only: force_analy_rstrt, force_analy_dump
      use multiple_states, only: efficiency_rstrt, efficiency_dump
      use properties_mod, only: prop_rstrt, prop_dump
      use system, only: nelec
      use system, only: nup
      use system, only: ndn
      use multiple_geo, only: nforce
      
      implicit none

      integer :: i, ib, ic, ifr, istate
      integer :: j, jel, k, nbasx
      integer :: ncentx, nctypex, ndetx, ndnx
      integer :: newghostypex, nghostcentx, norbx, nstepx
      integer :: nupx
      integer, dimension(nctype_tot) :: nsx
      integer, dimension(nctype_tot) :: npxx
      integer, dimension(nctype_tot) :: npyx
      integer, dimension(nctype_tot) :: npzx
      integer, dimension(nctype_tot) :: ndxxx
      integer, dimension(nctype_tot) :: ndxyx
      integer, dimension(nctype_tot) :: ndxzx
      integer, dimension(nctype_tot) :: ndyyx
      integer, dimension(nctype_tot) :: ndyzx
      integer, dimension(nctype_tot) :: ndzzx

      real(dp) :: ajacob, deltarx
      real(dp) :: deltatx, deltax, dist, distance_node
      real(dp) :: pecx, psidg, rnorm_nodes
      real(dp), dimension(nbasis,norb_tot) :: coefx
      real(dp), dimension(nbasis) :: zexx
      real(dp), dimension(3,ncent_tot) :: centx
      real(dp), dimension(nctype_tot) :: znucx
      real(dp), dimension(ndet) :: cdetx
      real(dp), dimension(3,nelec) :: xstrech
      real(dp), dimension(MSTATES) :: d2
      real(dp), parameter :: half = 0.5d0
      real(dp), parameter :: small = 1.d-6





      write(10) delta,deltar,deltat

      write(10) vmc_nstep,iblk
      do istate=1,nstates
        write(10) ecum1(istate),(ecum(istate,i),i=1,nforce),pecum(istate),tpbcum(istate),tjfcum(istate),r2cum,acc
        write(10) ecm21(istate),(ecm2(istate,i),i=1,nforce),pecm2(istate),tpbcm2(istate),tjfcm2(istate),r2cm2
        if(nforce.gt.1) then
          write(10) (wcum(istate,i),fcum(istate,i),fcm2(istate,i),i=1,nforce)
         else
          write(10) wcum(istate,1)
        endif
        write(10) ecum1s(istate),ecm21s(istate)
      enddo
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
      write(10) (ns(i),i=1,nctype)
      write(10) (npx(i),i=1,nctype)
      write(10) (npy(i),i=1,nctype)
      write(10) (npz(i),i=1,nctype)
      write(10) (ndxx(i),i=1,nctype)
      write(10) (ndxy(i),i=1,nctype)
      write(10) (ndxz(i),i=1,nctype)
      write(10) (ndyy(i),i=1,nctype)
      write(10) (ndyz(i),i=1,nctype)
      write(10) (ndzz(i),i=1,nctype)
      write(10) (cdet(i,1,1),i=1,ndet)
      write(10) ndet,nup,ndn

      call optjas_dump(10)
      call optx_jas_orb_dump(10)
      call optx_jas_ci_dump(10)
      call optx_orb_ci_dump(10)

      rewind 10

      write(ounit,'(1x,''successful dump to unit 10'')')

      return

c-----------------------------------------------------------------------
      entry startr_more

      read(10) deltax,deltarx,deltatx
      if (dabs(deltax-delta).gt.small) call fatal_error('STARTR: delta')
      if (dabs(deltarx-deltar).gt.small) call fatal_error('STARTR: deltar')
      if (dabs(deltatx-deltat).gt.small) call fatal_error('STARTR: deltat')

      read(10) nstepx,iblk
      if (nstepx.ne.vmc_nstep) call fatal_error('STARTR: nstep')
      do istate=1,nstates
        read(10) ecum1(istate),(ecum(istate,i),i=1,nforce),pecum(istate),tpbcum(istate),tjfcum(istate),r2cum,acc
        read(10) ecm21(istate),(ecm2(istate,i),i=1,nforce),pecm2(istate),tpbcm2(istate),tjfcm2(istate),r2cm2
        if(nforce.gt.1) then
          read(10) (wcum(istate,i),fcum(istate,i),fcm2(istate,i),i=1,nforce)
         else
          read(10) wcum(istate,1)
        endif
        read(10) ecum1s(istate),ecm21s(istate)
      enddo
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
      do j=1,norb
        do i=1,nbasis
          if (dabs(coefx(i,j)-coef(i,j,1)).gt.small) call fatal_error('STARTR: coef')
        enddo
      enddo
      do i=1,nbasis
        if (dabs(zexx(i)-zex(i,1)).gt.small) call fatal_error('STARTR: zex')
      enddo
      do i=1,ncent+nghostcent
        do k=1,3
          if (dabs(cent(k,i)-centx(k,i)).gt.small) call fatal_error('STARTR: cent')
        enddo
      enddo
      if (pecx.ne.pecent) call fatal_error('STARTR: pec')
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

      call optjas_rstrt(10)
      call optx_jas_orb_rstrt(10)
      call optx_jas_ci_rstrt(10)
      call optx_orb_ci_rstrt(10)

      write(ounit,'(1x,''succesful read from unit 10'')')
      write(ounit,'(t5,''enow'',t15,''eave'',t25,''eerr'',t35,''peave'',
     &t45,''peerr'',t55,''tpbave'',t65,''tpberr'',t75,''tjfave'',
     &t85,''tjferr'',t95,''accept'',t105,''iter'')')

      if(nforce.gt.1) then
        call setup_force
       else
        nwftype=1
        iwftype(1)=1
      endif

c loop over secondary config
      do ifr=2,nforce
c set n- and e-coord and n-n potential
        call strech(xold,xstrech,ajacob,ifr,1)
        call hpsi(xstrech,psido,psijo,eold(1,ifr),0,ifr)
        do istate=1,nforce
          psi2o(istate,ifr)=2*(dlog(dabs(psido(istate)))+psijo)+dlog(ajacob)
        enddo
      enddo

c primary config
c set n-coord and n-n potential
      if(nforce.gt.1) call strech(xold,xstrech,ajacob,1,0)
      call hpsi(xold,psido,psijo,eold(1,1),0,1)
      do istate=1,nforce
        psi2o(istate,1)=2*(dlog(dabs(psido(istate)))+psijo)
        tjfo(istate)=d2(istate)
        tjfo(istate)=-tjfo(istate)*half*hb
      enddo

      if(iguiding.gt.0) then
        call determinant_psig(psido,psidg)
c rewrite psi2o if you are sampling guiding
        psi2o(1,1)=2*(dlog(dabs(psidg))+psijo)
      endif

      if(node_cutoff.gt.0) then
        do jel=1,nelec
          call compute_determinante_grad(jel,psido(1),psido,vold(1,jel),1)
        enddo
        call nodes_distance(vold,distance_node,1)
        rnorm_nodes=rnorm_nodes_num(distance_node,eps_node_cutoff)/distance_node

        psi2o(1,1)=psi2o(1,1)+2*dlog(rnorm_nodes)
      endif

      if(ioptorb.gt.0) ns_current=0

      call prop_save
      call optjas_save
      call optci_save
      call optorb_save

      do i=1,nelec
        do k=1,3
          xnew(k,i)=xold(k,i)
        enddo
      enddo

      do i=1,nelec
        rmino(i)=99.d9
        do j=1,ncent
          dist=0
          do k=1,3
            dist=dist+(xold(k,i)-cent(k,j))**2
          enddo
          if(dist.lt.rmino(i)) then
            rmino(i)=dist
            nearesto(i)=j
          endif
        enddo
        rmino(i)=dsqrt(rmino(i))
        do k=1,3
          rvmino(k,i)=xold(k,i)-cent(k,nearesto(i))
        enddo
      enddo

      do istate=1,nstates
        do ifr=1,nforce
          esum(istate,ifr)=0
          wsum(istate,ifr)=0
        enddo
        pesum(istate)=0
        tpbsum(istate)=0
        tjfsum(istate)=0
      enddo
      r2sum=0

      return
      end
      end module
