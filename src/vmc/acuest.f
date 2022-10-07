      module acuest_mod
      contains
      subroutine acuest
c Written by Cyrus Umrigar, modified by Claudia Filippi
c routine to accumulate estimators for energy etc.

      use precision_kinds, only: dp
      use force_mod, only: MFORCE
      use vmc_mod, only: nrad
      use system, only: znuc, cent, pecent, iwctype, ncent
      use mstates_mod, only: MSTATES
      use const, only: nelec, ipr
      use config, only: eold, nearesto, psi2o
      use config, only: psido, psijo, rmino, rvmino
      use config, only: vold, xold
      use csfs, only: nstates
      use denupdn, only: rprobdn, rprobup
      use est2cm, only: ecm2, ecm21, pecm2, r2cm2, tjfcm2, tpbcm2
      use estcum, only: ecum, ecum1, iblk, pecum, r2cum, tjfcum, tpbcum
      use estpsi, only: apsi, aref, detref
      use estsig, only: ecm21s, ecum1s
      use estsum, only: acc, esum, esum1, pesum, r2sum, tjfsum, tpbsum
      use forcepar, only: nforce
      use forcest, only: fcm2, fcum
      use forcewt, only: wcum, wsum
      use multidet, only: kref
      use optwf_contrl, only: ioptorb
      use step, only: ekin, ekin2, rprob, suc, trunfb, try
      use pseudo, only: nloc
      use qua, only: nquad, wq, xq, yq, zq
      use mstates_ctrl, only: iguiding

      use optorb_cblock, only: ns_current
      use distance_mod, only: rshift, r_en, rvec_en
      use multislater, only: detiab
      use distance_mod, only: rshift, r_en, rvec_en
      use inputflags, only: node_cutoff, eps_node_cutoff
      use contrl_file, only: ounit
      
      use distances_mod, only: distances
      use force_analytic, only: force_analy_save
      use optorb_f_mod, only: optorb_save
      use optci_mod, only: optci_save
      use optjas_mod,   only: optjas_save
      use mmpol_vmc,    only: mmpol_save, mmpol_cum
      use pcm_vmc,      only: pcm_save, pcm_cum
      use prop_vmc,     only: prop_save
      use nodes_distance_mod, only: nodes_distance
      use determinante_mod, only: compute_determinante_grad
      use determinant_psig_mod,  only: determinant_psig
      use hpsi_mod, only: hpsi
      use strech_mod, only: strech
      use pot, only: pot_nn
      use multiple_states, only: efficiency_init
      use force_analytic, only: force_analy_init, force_analy_cum
      use properties_mod, only: prop_init, prop_cum
      use pcm_mod, only: pcm_init
      use mmpol, only: mmpol_init
      use optci_mod, only: optci_init, optci_cum
      use optorb_f_mod, only: optorb_init, optorb_cum
      use optjas_mod, only: optjas_init, optjas_cum
      use optx_orb_ci, only: optx_orb_ci_init
      use optx_jas_ci, only: optx_jas_ci_init
      use optx_jas_orb, only: optx_jas_orb_init
      use rotqua_mod, only: gesqua
      use acuest_reduce_mod, only: acuest_reduce, acues1_reduce
      use nodes_distance_mod, only: rnorm_nodes_num


      implicit none

      integer :: i, ic, ifr, istate, jel
      integer :: k
      real(dp) :: ajacob, distance_node, penow
      real(dp) :: psidg, r2now, rnorm_nodes, tjfnow
      real(dp) :: tpbnow
      real(dp), dimension(3,nelec) :: xstrech
      real(dp), dimension(MSTATES) :: wtg
      real(dp), parameter :: half = .5d0



      real(dp), dimension(:,:), allocatable :: enow

c xsum = sum of values of x from metrop
c xnow = average of values of x from metrop
c xcum = accumulated sums of xnow
c xcm2 = accumulated sums of xnow**2

      allocate(enow(MSTATES, MFORCE))


c collect cumulative averages

      do ifr=1,nforce
        do istate=1,nstates

        enow(istate,ifr)=esum(istate,ifr)/wsum(istate,ifr)
        wcum(istate,ifr)=wcum(istate,ifr)+wsum(istate,ifr)
        ecum(istate,ifr)=ecum(istate,ifr)+esum(istate,ifr)
        ecm2(istate,ifr)=ecm2(istate,ifr)+esum(istate,ifr)*enow(istate,ifr)

        if(ifr.eq.1) then
          penow=pesum(istate)/wsum(istate,ifr)
          tpbnow=tpbsum(istate)/wsum(istate,ifr)
          tjfnow=tjfsum(istate)/wsum(istate,ifr)
          r2now=r2sum/(wsum(istate,ifr)*nelec)

          pecm2(istate)=pecm2(istate)+pesum(istate)*penow
          tpbcm2(istate)=tpbcm2(istate)+tpbsum(istate)*tpbnow
          tjfcm2(istate)=tjfcm2(istate)+tjfsum(istate)*tjfnow
          r2cm2=r2cm2+r2sum*r2now/nelec

          pecum(istate)=pecum(istate)+pesum(istate)
          tpbcum(istate)=tpbcum(istate)+tpbsum(istate)
          tjfcum(istate)=tjfcum(istate)+tjfsum(istate)
          r2cum=r2cum+r2sum/nelec

         else
          fcum(istate,ifr)=fcum(istate,ifr)+wsum(istate,1)*(enow(istate,ifr)-esum(istate,1)/wsum(istate,1))
          fcm2(istate,ifr)=fcm2(istate,ifr)+wsum(istate,1)*(enow(istate,ifr)-esum(istate,1)/wsum(istate,1))**2
        endif
        enddo
      enddo

c only called for ifr=1
      call optjas_cum(wsum(1,1),enow(1,1))
      call optorb_cum(wsum(1,1),esum(1,1))
      call optci_cum(wsum(1,1))
      call prop_cum(wsum(1,1))
      call pcm_cum(wsum(1,1))
      call mmpol_cum(wsum(1,1))
      call force_analy_cum(wsum(1,1),ecum(1,1)/wcum(1,1),wcum(1,1))

c zero out xsum variables for metrop

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

      call prop_init(1)
      call optorb_init(1)
      call optci_init(1)
      call pcm_init(1)
      call mmpol_init(1)
      call force_analy_init(1)


      call acuest_reduce(enow)
      if(allocated(enow)) deallocate(enow)

      return
c-----------------------------------------------------------------------
      entry acues1(wtg)
c statistical fluctuations without blocking
      do istate=1,nstates
        ecum1(istate)=ecum1(istate)+esum1(istate)*wtg(istate)
        ecm21(istate)=ecm21(istate)+esum1(istate)**2*wtg(istate)
        esum1(istate)=0

        apsi(istate)=apsi(istate)+dabs(psido(istate))
      enddo

      aref=aref+dabs(detiab(kref,1)*detiab(kref,2))

      detref(1)=detref(1)+dlog10(dabs(detiab(kref,1)))
      detref(2)=detref(2)+dlog10(dabs(detiab(kref,2)))

      call acues1_reduce

      return
c-----------------------------------------------------------------------
      entry acusig(wtg)
c sigma evaluation
      do istate=1,nstates
        ecum1s(istate)=ecum1s(istate)+esum1(istate)*wtg(istate)
        ecm21s(istate)=ecm21s(istate)+esum1(istate)**2*wtg(istate)
        esum1(istate)=0
      enddo
      return
c-----------------------------------------------------------------------
      entry zerest

c entry point to zero out all averages etc.
c the initial values of energy psi etc. is also calculated here
c although that really only needs to be done before the equil. blocks.

      iblk=0

c set quadrature points
      if(nloc.gt.0) call gesqua(nquad,xq,yq,zq,wq)

c zero out estimators
      acc=0
      do istate=1,nstates
        pecum(istate)=0
        tpbcum(istate)=0
        tjfcum(istate)=0
        ecum1(istate)=0
        ecum1s(istate)=0

        pecm2(istate)=0
        tpbcm2(istate)=0
        tjfcm2(istate)=0
        ecm21(istate)=0
        ecm21s(istate)=0

        pesum(istate)=0
        tpbsum(istate)=0
        tjfsum(istate)=0

        apsi(istate)=0
      enddo

      detref(1)=0
      detref(2)=0

      aref=0

      r2cm2=0
      r2cum=0
      r2sum=0

      call optjas_init
      call optci_init(0)
      call optorb_init(0)
      call optx_jas_orb_init
      call optx_jas_ci_init
      call optx_orb_ci_init

      call prop_init(0)
      call pcm_init(0)
      call mmpol_init(0)
      call force_analy_init(0)
      call efficiency_init

      do ifr=1,nforce
        do istate=1,nstates
          ecum(istate,ifr)=0
          ecm2(istate,ifr)=0
          wcum(istate,ifr)=0
          esum(istate,ifr)=0
          wsum(istate,ifr)=0
          fcum(istate,ifr)=0
          fcm2(istate,ifr)=0
        enddo
      enddo

c Zero out estimators for acceptance, force-bias trun., kin. en. and density
      do i=1,nrad
        try(i)=0
        suc(i)=0
        trunfb(i)=0
        ekin(i)=0
        ekin2(i)=0
        rprobup(i)=0
        rprobdn(i)=0
        rprob(i)=0
      enddo

c get nuclear potential energy
      call pot_nn(cent,znuc,iwctype,ncent,pecent)

c get wavefunction etc. at initial point

c secondary configs
c set n- and e-coords and n-n potentials before getting wavefn. etc.
      do ifr=2,nforce
        call strech(xold,xstrech,ajacob,ifr,1)
        call hpsi(xstrech,psido,psijo,eold(1,ifr),0,ifr)
        do istate=1,nstates
          psi2o(istate,ifr)=2*(dlog(dabs(psido(istate)))+psijo)+dlog(ajacob)
        enddo
      enddo

c primary config
c set n- and e-coords and n-n potentials before getting wavefn. etc.
      if(nforce.gt.1) call strech(xold,xstrech,ajacob,1,0)
      call hpsi(xold,psido,psijo,eold(1,1),0,1)

      do istate=1,nstates
        psi2o(istate,1)=2*(dlog(dabs(psido(istate)))+psijo)
      enddo

      if(iguiding.gt.0) then
        call determinant_psig(psido,psidg)
c rewrite psi2o if you are sampling guiding
        psi2o(1,1)=2*(dlog(dabs(psidg))+psijo)
      endif

      if(ipr.gt.1) then
        write(ounit,'(''psid, psidg='',2d12.4)') psido(1),psidg
        write(ounit,'(''psid2o='',f9.4)') psi2o(1,1)
      endif

      if(node_cutoff.gt.0) then
        do jel=1,nelec
          call compute_determinante_grad(jel,psido(1),psido,vold(1,jel),1)
        enddo
        call nodes_distance(vold,distance_node,1)
        rnorm_nodes=rnorm_nodes_num(distance_node,eps_node_cutoff)/distance_node

        psi2o(1,1)=psi2o(1,1)+2*dlog(rnorm_nodes)

        if(ipr.gt.1) then
          write(ounit,'(''distance_node='',d12.4)') distance_node
          write(ounit,'(''rnorm_nodes='',d12.4)') rnorm_nodes
          write(ounit,'(''psid2o_ncut='',f9.4)') psi2o(1,1)
          do i=1,nelec
              write(ounit,'(''vd'',3e20.10)') (vold(k,i),k=1,3)
          enddo
        endif
      endif

      if(ioptorb.gt.0) ns_current=0

      call prop_save
      call pcm_save
      call mmpol_save
      call optjas_save
      call optci_save
      call optorb_save
      call force_analy_save

c get interparticle distances
      call distances(0,xold)

      do i=1,nelec
        rmino(i)=99.d9
        do ic=1,ncent
          if(r_en(i,ic).lt.rmino(i)) then
            rmino(i)=r_en(i,ic)
            nearesto(i)=ic
          endif
        enddo
        do  k=1,3
          rvmino(k,i)=rvec_en(k,i,nearesto(i))
        enddo
      enddo

      return
      end
      end module
