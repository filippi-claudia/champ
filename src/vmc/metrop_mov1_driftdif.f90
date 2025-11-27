      module metrop_mov1_driftdif
      contains

      subroutine metrop1(ipass,irun)
      use acuest_mod, only: acues1,acusig
      use config, only: eold,psi2n,psi2o, psido,psijo
      use config, only: vnew,vold,xnew,xold
      use constants, only: pi
      use contrl_file, only: ounit
      use control, only: ipr
      use csfs, only: nstates
      use determinant_psig_mod, only: determinant_psig
      use determinante_mod, only: compute_determinante_grad
      use detsav_mod, only: detsav
      use distances_mod, only: distancese_restore
      use estsum,  only: acc,esum,esum1,pesum,tpbsum
      use force_analytic, only: force_analy_sum
      use force_pth, only: PTH
      use forcewt, only: wsum
      use gauss_mod, only: gauss
      use hpsi_mod, only: hpsi
      use hpsie, only: psie
      use inputflags, only: eps_node_cutoff,node_cutoff
      use jassav_mod, only: jassav
      use m_force_analytic, only: iforce_analy
      use metropolis, only: vmc_tau
      use mmpol, only: mmpol_efield
      use mmpol_cntrl, only: ich_mmpol
      use mmpol_vmc, only: mmpol_sum
      use mstates_ctrl, only: iguiding
      use mstates_mod, only: MSTATES
      use multideterminant_mod, only: update_ymat
      use multiple_geo, only: nforce
      use multiple_states, only: efficiency_sample
      use nodes_distance_mod, only: nodes_distance,rnorm_nodes_num
      use optci_mod, only: optci_sum
      use optjas_mod, only: optjas_sum
      use optorb_f_mod, only: check_orbitals,check_orbitals_reset
      use optorb_f_mod, only: optorb_sum
      use optwf_handle_wf, only: optwf_store
      use optx_jas_ci, only: optx_jas_ci_sum
      use optx_jas_orb, only: optx_jas_orb_sum
      use optx_orb_ci, only: optx_orb_ci_sum
      use pathak_mod, only: ipathak, eps_pathak, pnew, pathak
      use pcm_cntrl, only: ichpol
      use pcm_mod, only: qpcm_efield
      use pcm_vmc, only: pcm_sum
      use precision_kinds, only: dp
      use prop_vmc, only: prop_sum
      use random_mod, only: random_dp
      use strech_mod, only: strech
      use system,  only: nelec,nup
      use tmpnode, only: distance_node_sum
      use vmc_mod, only: nwftypejas, stoj

#if defined(TREXIO_FOUND) && defined(QMCKL_FOUND) 
      use qmckl_data
#endif

      implicit none

      integer :: i, iab, ic, iel, iflag_dn, iflag_up, ijas1, rc
      integer :: ifr, ii, ipass, irun, istate
      integer :: j, jel, k, iph
      real(dp) :: ajacob, distance_node, dmin1
      real(dp) :: p, psi2g, q
      real(dp) :: rnorm_nodes, vavvt, v2new, v2old, wstro
      real(dp) :: rttau, tau, drift, dfus2o, dfus2n, dfus, dx
      real(dp), dimension(3,nelec) :: xstrech
      real(dp), dimension(3) :: xbac
      real(dp), dimension(MSTATES) :: ekino
      real(dp), dimension(MSTATES) :: psidn
      real(dp), dimension(nwftypejas) :: psijn
      real(dp), dimension(MSTATES) :: wtg
      real(dp), dimension(MSTATES) :: wtg_sqrt 
      real(dp), dimension(MSTATES) :: zero_array = 0.0_dp
      real(dp), parameter :: adrift = 0.5d0
      real(dp), parameter :: zero = 0.d0
      real(dp), parameter :: one = 1.d0
      real(dp), parameter :: two = 2.d0

      tau=vmc_tau
      rttau=sqrt(tau)

      call check_orbitals

      ijas1=stoj(1)
      do i=1,nelec

        if(i.le.nup) then
          iab=1
          iflag_up=2
          iflag_dn=3
         else
          iab=2
          iflag_up=3
          iflag_dn=2
        endif

        if(iguiding.eq.0) then
          psi2g=1.d0
        else
          call determinant_psig(psido,psijo,psi2g)
        endif
        call compute_determinante_grad(i,psi2g,psido,psijo,vold(1,i),1)

        v2old=vold(1,i)**2+vold(2,i)**2+vold(3,i)**2
        vavvt=(dsqrt(one+two*adrift*v2old*tau)-one)/(adrift*v2old)

        dfus2o=zero
        do k=1,3
          drift=vavvt*vold(k,i)*tau
          dfus=gauss()*rttau
          dx=drift+dfus
          dfus2o=dfus2o+dfus**2
          xnew(k,i)=xold(k,i)+dx
        enddo

! calculate psi at new configuration
        iel=i

        call psie(iel,xnew,psidn,psijn,ipass,0)

        if(iguiding.eq.0) then
          psi2g=1.d0
        else
          call determinant_psig(psidn,psijn,psi2g)
        endif

        if(psi2g.eq.0.d0) then
          p=zero
          q=one
          goto 208
        endif

        call compute_determinante_grad(iel,psi2g,psidn,psijn,vnew(1,iel),0)

        if(ipr.gt.1) then
          write(ounit,'(''psidn,psi2g ='',3d12.4)') psidn(1),psi2g
          write(ounit,'(''psidn,psi2g ='',3d12.4)') psijn(1),psi2g
        endif

        psi2n(1)=dlog(dabs(psi2g))+2*(dlog(dabs(psidn(1)))+psijn(ijas1))

        if(node_cutoff.ne.0) then
          do jel=1,nup
            if(jel.ne.iel) call compute_determinante_grad(jel,psi2g,psidn,psijn,vnew(1,jel),iflag_up)
          enddo
  
          do jel=nup+1,nelec
            if(jel.ne.iel) call compute_determinante_grad(jel,psi2g,psidn,psijn,vnew(1,jel),iflag_dn)
          enddo
  
          call nodes_distance(vnew,distance_node,0)
          rnorm_nodes=rnorm_nodes_num(distance_node,eps_node_cutoff)/distance_node
  
          psi2n(1)=psi2n(1)+2*dlog(rnorm_nodes)
  
          if(ipr.gt.1) then
            write(ounit,'(''distance_node='',d12.4)') distance_node
            write(ounit,'(''rnorm_nodes='',d12.4)') rnorm_nodes
            write(ounit,'(''psid2n_ncut='',f9.4)') psi2n(1)
          endif
        endif

        v2new=vnew(1,iel)**2+vnew(2,iel)**2+vnew(3,iel)**2
        vavvt=(dsqrt(one+two*adrift*v2new*tau)-one)/(adrift*v2new)

! calculate probability for reverse transition
        dfus2n=zero
        do k=1,3
          drift=vavvt*vnew(k,iel)*tau
          xbac(k)=xnew(k,i)+drift
          dfus=xbac(k)-xold(k,iel)
          dfus2n=dfus2n+dfus**2
        enddo

! p is the probability of accepting new move
        p=exp(psi2n(1)-psi2o(1,1))*exp(-0.5*(dfus2n-dfus2o)/tau)

        if(ipr.ge.1) then
          write(ounit,'(''vold='',9d12.4)') (vold(ic,iel),ic=1,3)
          write(ounit,'(''vnew='',9d12.4)') (vnew(ic,iel),ic=1,3)
          write(ounit,'(''dfus2o,dfus2n,psi2n,psi2o,p'',9d12.4)') &
                        dfus2o,dfus2n,psi2n(1),psi2o(1,1),p
          if(dabs(vnew(1,i))+dabs(vnew(1,i))+dabs(vnew(1,i)).gt.10d+8) then
            do ii=1,i-1
              write(ounit,*) (xold(k,ii),k=1,3)
            enddo
            write(ounit,*) (xnew(k,i),k=1,3)
            do ii=i+1,nelec
              write(ounit,*) (xold(k,ii),k=1,3)
            enddo
          endif
        endif

      p=dmin1(one,p)
      q=one-p

  208 continue

! accept new move with probability p
! Note when one electron moves the velocity on all electrons change.
      if (random_dp().lt.p) then
#if defined(TREXIO_FOUND) && defined(QMCKL_FOUND) 
        if (use_qmckl_jastrow) then
          rc = qmckl_get_jastrow_champ_single_accept(qmckl_ctx(qmckl_no_ctx))
        endif
#endif
        psi2o(1,1)=psi2n(1)
        do ic=1,3
          xold(ic,i)=xnew(ic,i)
        enddo
        if(node_cutoff.gt.0) then
          do ic=1,3
            do ii=1,nelec
              vold(ic,ii)=vnew(ic,ii)
            enddo
          enddo
        endif
        do istate=1,nstates
          psido(istate)=psidn(istate)
          psijo(stoj(istate))=psijn(stoj(istate))
        enddo
        acc=acc+one
        call jassav(i,0)
        call detsav(i,0)
        if(ipr.ge.1) write(ounit,*)'METROP ACCEPT'
       else
        if(ipr.ge.1) write(ounit,*)'METROP REJECT'
        do ic=1,3
          xnew(ic,i)=xold(ic,i)
        enddo
        call distancese_restore(i)
      endif

      call update_ymat(i)

      enddo

! loop over secondary configurations
      do ifr=2,nforce
        call strech(xold,xstrech,ajacob,ifr,1)
        call hpsi(xstrech,psido(1),psijo,ekino,eold(1,ifr),ipass,ifr)
        do istate=1,nstates
          j=stoj(istate)
          psi2o(istate,ifr)=2*(dlog(dabs(psido(istate)))+psijo(j))+dlog(ajacob)
        enddo
      enddo

      call check_orbitals_reset

! primary configuration
      if(nforce.gt.1) call strech(xold,xstrech,ajacob,1,0)
      call hpsi(xold,psido(1),psijo,ekino,eold(1,1),ipass,1)
      do istate=1,nstates
         j=stoj(istate)
         psi2o(istate,1)=2*(dlog(dabs(psido(istate)))+psijo(j))
      enddo

! psig is relative to the first state
      if(iguiding.eq.0) then
        psi2g=1.d0
       else
        call determinant_psig(psido,psijo,psi2g)
      endif

      if(ipr.gt.1) then
        write(ounit,'(''psid,psi2g ='',2d12.4)') psido(1),psi2g
      endif

      rnorm_nodes=1.d0
      if(node_cutoff.gt.0) then
        do jel=1,nelec
          call compute_determinante_grad(jel,psi2g,psido(1),psijo,vold(1,jel),1)
        enddo
        call nodes_distance(vold,distance_node,1)
        rnorm_nodes=rnorm_nodes_num(distance_node,eps_node_cutoff)/distance_node
        psi2g=psi2g*rnorm_nodes*rnorm_nodes
        if(ipr.gt.1) then
          write(ounit,'(''distance_node='',d12.4)') distance_node
          write(ounit,'(''rnorm_nodes='',d12.4)') rnorm_nodes
          write(ounit,'(''psi2g_ncut='',d12.4)') psi2g
        endif
        distance_node_sum=distance_node_sum+distance_node
      endif

      if(iforce_analy.eq.1) then
         if (ipathak.gt.0) then
            do jel=1,nelec
              call compute_determinante_grad(jel,psi2g,psido(1),psijo,vold(1,jel),1)
            enddo
            call nodes_distance(vold,distance_node,1)
            do iph=1,PTH
               call pathak(distance_node,pnew(iph),eps_pathak(iph))
            enddo
            distance_node_sum=distance_node_sum+distance_node
         endif
      endif

      do istate=1,nstates
        j=stoj(istate)
        wtg_sqrt(istate)=psido(istate)/psido(1)*exp(psijo(j)-psijo(ijas1))/sqrt(psi2g)
        wtg(istate)=wtg_sqrt(istate)*wtg_sqrt(istate)

! form expected values of e, pe, etc.
        esum1(istate)=eold(istate,1)
        wsum(istate,1)=wsum(istate,1)+wtg(istate)
        esum(istate,1)=esum(istate,1)+eold(istate,1)*wtg(istate)
        pesum(istate)=pesum(istate)+(eold(istate,1)-ekino(istate))*wtg(istate)
        tpbsum(istate)=tpbsum(istate)+ekino(istate)*wtg(istate)
      enddo

      if(ipr.gt.1) write(ounit,'(''energy reweighted '',d12.4)') eold(1,1)*wtg(1)

! normal component efield on cavity surface to compute a new set of polarization charges
      if(ichpol.eq.1) call qpcm_efield(nelec,xold)
! efield dovuto agli elettroni sui siti dei dipoli
      if(ich_mmpol.eq.1) call mmpol_efield(nelec,xold)

! use 'new' not 'old' value
      call pcm_sum(wtg(1),0.d0)
      call mmpol_sum(wtg(1),0.d0)
      call prop_sum(wtg(1),0.d0)

      call force_analy_sum(wtg(1),0.d0,eold(1,1),0.0d0)

      call optjas_sum(wtg,zero_array,eold(1,1),eold(1,1),0)
      call optorb_sum(wtg,zero_array,eold(1,1),eold(1,1),0)
      call optci_sum(wtg(1),0.d0,eold(1,1),eold(1,1))

      call optx_jas_orb_sum(wtg(1),zero_array,0)
      call optx_jas_ci_sum(wtg(1),0.d0,eold(1,1),eold(1,1))
      call optx_orb_ci_sum(wtg(1),0.d0)

      if(irun.eq.1) call optwf_store(ipass,wtg,wtg_sqrt,psido,eold(1,1))

      call efficiency_sample(ipass,psido,psijo,psi2g)

      call acues1(wtg)

      do istate=1,nstates
        esum1(istate)=eold(istate,1)
      enddo
      call acusig(wtg)

      do ifr=2,nforce
        do istate=1,nstates
          wstro=exp(psi2o(istate,ifr)-psi2o(istate,1))*wtg(istate)
          esum(istate,ifr)=esum(istate,ifr)+eold(istate,ifr)*wstro
          wsum(istate,ifr)=wsum(istate,ifr)+wstro
        enddo
      enddo

! rewrite psi2o for next metropolis step if you are sampling guiding
      if(iguiding.gt.0) psi2o(1,1)=dlog(dabs(psi2g))+2*(dlog(dabs(psido(1)))+psijo(ijas1))

      if(node_cutoff.gt.0) then
        psi2o(1,1)=psi2o(1,1)+2*dlog(rnorm_nodes)
      endif
      return

      end
      end module
