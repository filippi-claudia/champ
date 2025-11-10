      module dmc_ae_mov1
      contains

      subroutine dmc_ae(lpass,irun)

      use age,     only: iage,ioldest,ioldestmx
      use assignment_mod, only: assign_elecs
      use averages, only: average
      use branch,  only: eest,esigma,eigv,eold,ff,fprod,nwalk,pwt,wdsumo
      use branch,  only: wgdsumo,wt,wthist
      use branching, only: calculate_fratio, calculate_reweight
      use config,  only: psido_dmc,psijo_dmc,vold_dmc,xold_dmc
      use const,   only: etrial,esigmatrial
      use constants, only: hb
      use contrl_file, only: ounit
      use contrldmc, only: iacc_rej,icross,icut_br,icut_e,idmc,ipq,limit_wt_dmc
      use contrldmc, only: nfprod,rttau,tau, ibranching_c
      use control, only: ipr
      use control_dmc, only: dmc_irstar,dmc_nconf
      use da_energy_now, only: da_energy,da_psi
      use derivest, only: derivsum
      use determinante_mod, only: compute_determinante_grad
      use detsav_mod, only: detsav
      use distances_mod, only: distances,distancese_restore
      use distance_mod, only: r_en,rvec_en
      use estcum,  only: ipass
      use estsum,  only: efsum1,egsum1,esum1_dmc,pesum_dmc
      use estsum,  only: tausum,tpbsum_dmc,wfsum1,wgsum1,wsum1
      use force_analytic, only: force_analy_sum, force_analy_save
      use force_analytic, only: force_analy_vd
      use force_pth, only: PTH
      use fragments, only: eloc_i, eloco_i, esum_i, eest_i, v2_i, vav2_i, fration_i, fratio_i
      use fragments,  only: nfrag, ifragelec, ifragcent, elocfrag, elocofrag, esumfrag, eestfrag
      use fragments,  only: v2frag, vav2frag, fratiofrag, frationfrag, potnnfrag, ibranching_cfrag, sqrt_nelecfrag
      use fragments,  only: dfus2acfrag, dfus2unfrag, tauefffrag, etrialfrag, egsum1frag
      use gauss_mod, only: gauss
      use general, only: write_walkalize
      use hpsi_mod, only: hpsi
      use hpsiedmc, only: psiedmc
      use inputflags, only: eps_node_cutoff,icircular
      use inputflags, only: node_cutoff
      use jacobsave, only: ajacob,ajacold
      use jassav_mod, only: jassav
      use mmpol_dmc, only: mmpol_save,mmpol_sum
      use multideterminant_mod, only: update_ymat
      use multideterminant_tmove_mod, only: multideterminant_tmove
      use multiple_geo, only: istrech,itausec,nforce,nwprod
      use m_force_analytic, only: iforce_analy
      use mpiconf, only: idtask,mpiconf_init,nproc,wid
      use nodes_distance_mod, only: nodes_distance,rnorm_nodes_num
      use nonloc_grid_mod, only: nonloc_grid,t_vpsp_get,t_vpsp_sav
      use nonloc_pot_mod, only: nonloc_pot
      use optci_mod, only: optci_sum
      use optjas_mod, only: optjas_sum
      use optorb_f_mod, only: optorb_sum
      use optwf_handle_wf, only: optwf_store
      use optwf_parms, only: nparmj
      use optx_jas_ci, only: optx_jas_ci_sum
      use optx_jas_orb, only: optx_jas_orb_sum
      use optx_orb_ci, only: optx_orb_ci_sum
      use pathak_mod, only: ipathak, eps_pathak, pold, pnew, pathak
      use pcm_dmc, only: pcm_save,pcm_sum
      use precision_kinds, only: dp
      use prop_dmc, only: prop_save_dmc,prop_sum_dmc
      use random_mod, only: random_dp
      use splitj_mod, only: splitj
      use stats, only: acc,dfus2ac,dfus2un,nacc,nodecr
      use stats, only: trymove
      use strech_mod, only: strech
      use system,  only: cent,iwctype, znuc, nelec,nup, ncent
      use velratio, only: fratio
      use vd_mod, only: dmc_ivd
      use vmc_mod, only: nbjx
      use walksav_det_mod, only: walksav_det,walkstrdet
      use walksav_jas_mod, only: walksav_jas,walkstrjas

#if defined(TREXIO_FOUND) && defined(QMCKL_FOUND) 
      use qmckl_data
#endif

      implicit none

      integer :: i, iaccept, iel, ic, iph, rc
      integer :: iflag_dn, iflag_up, ifr, ii
      integer :: imove, imove_dn, imove_up, ipmod, ipmod2, iw, irun
      integer :: iwmod, iwnuc, j, jel, k, lpass
      integer :: ncall
      integer, dimension(nelec) :: itryo
      integer, dimension(nelec) :: itryn
      integer, dimension(nelec) :: iacc_elec
      real(dp) :: adrift, costht, d2n, den, deo, dfus, dfusb
      real(dp) :: dfus2a, dfus2b, dfus2n, dfus2o, distance_node, distance_node_ratio2
      real(dp) :: dmin1, dr2, drifdif, drifdifgfunc
      real(dp) :: drifdifr, drifdifs, drift, dwt
      real(dp) :: ecuto, ecutn, fnormn, fnormo
      real(dp) :: driftr, dx, e_cutoff, dwt_cutoff, ekino(1), enew(1)
      real(dp) :: ewtn, ewto, expon, ffi
      real(dp) :: ffn, fration, fratio_aux, ginv, hafzr2
      real(dp) :: p, phi, pen, pgaus, pp, psi2savo
      real(dp) :: psidn(1), psijn(1), q, qgaus
      real(dp) :: ren2, ren2mn, rmino, rminn
      real(dp) :: rnorm_nodes, rnorm_nodes_new
      real(dp) :: rnorm_nodes_old, ro, rtry, taunow, tauprim, tratio
      real(dp) :: term, volda, voldp, voldr, vnewa, vnewp, vnewr, v2new, v2old, v2
      real(dp) :: v2sumn, v2sumo, vavfac, vav2sumn, vav2sumo
      real(dp) :: vavvn, vavvo, vavvt, wtg(1), wtg_sqrt(1)
      real(dp) :: wtg_derivsum1, wtnow
      real(dp) :: pi, sqrt_pi_o2, sintht, sqrt_nelec
      real(dp) :: t_norm, t_norm_new, zeta, zprime, xprime

      real(dp) :: pe
      real(dp), dimension(2,nbjx) :: vpsp_det
      real(dp), dimension(nparmj,nbjx) :: dvpsp_dj

      real(dp), dimension(3) :: rvminn, rvmino
      real(dp), dimension(3) :: xaxis, zaxis
      real(dp), dimension(3) :: xbac
      real(dp), dimension(3) :: xnew
      real(dp), dimension(3) :: x_tmove_old
      real(dp), dimension(3, nelec) :: xstrech
      real(dp), dimension(3, nelec) :: vnew
      real(dp), dimension(nelec) :: unacp
      real(dp), dimension(10, 3, ncent) :: deriv_esum
      real(dp), dimension(nelec) :: deo_i, den_i
      real(dp), dimension(nelec) :: esum_i1, fratio_aux_i
      real(dp), dimension(nfrag) :: deofrag, denfrag
      real(dp), dimension(nfrag) :: esumfrag1, fratio_auxfrag
      real(dp), dimension(nfrag) :: ewtofrag, ewtnfrag
      real(dp), parameter :: zero = 0.d0
      real(dp), parameter :: one = 1.d0
      real(dp), parameter :: two = 2.d0
      real(dp), parameter :: half = .5d0
      real(dp), parameter :: eps = 1.d-10
      real(dp), parameter :: huge = 1.d+100
      real(dp), parameter :: adrift0 = 0.1d0
      real(dp), parameter :: small = 1.d-10
      real(dp), parameter :: zero_1d(1) = (/0.d0/)

      data ncall /0/

      esum_i1 = zero
      esumfrag1 = zero

      sqrt_pi_o2 = 0.88622692545d0
      sqrt_nelec = dsqrt(dble(nelec))

      pi=4.d0*datan(1.d0)
      term=(sqrt(two*pi*tau))**3/pi

      eps_node_cutoff=eps_node_cutoff*sqrt(tau)
      e_cutoff=0.2d0*sqrt(nelec/tau)

      if(idmc.lt.0) then
        expon=1
        dwt=1
      endif
      ! Undo products
      ipmod=mod(ipass,nfprod)
      ipmod2=mod(ipass+1,nfprod)
      ginv=min(1.d0,tau)
      ffn=eigv*(wdsumo/dmc_nconf)**ginv
      ffi=one/ffn
      fprod=fprod*ffn/ff(ipmod)
      ff(ipmod)=ffn

      ! Undo weights
      iwmod=mod(ipass,nwprod)

      ! Store (well behaved velocity/velocity)
      if(idmc.gt.0.and.ncall.eq.0.and.dmc_irstar.eq.0) then
        do ifr=1,nforce
          do iw=1,nwalk
            tratio = one
            ! Todo Joris - Think of a way to not pass in the fragments if they are not needed. 
            call calculate_fratio(icut_e, adrift0, tratio, taunow, sqrt_nelec, &
                                  xold_dmc(:,:,iw,ifr), vold_dmc(:,:,iw,ifr), &
                                  eest, eold(iw,ifr), fratio(iw,ifr), &
                                  eest_i, eloco_i(:,iw,ifr), fratio_i(:,iw,ifr), &
                                  eestfrag, elocofrag(:,iw,ifr), fratiofrag(:,iw,ifr))
          enddo
        enddo
        ncall=ncall+1
      endif

      imove=0
      ioldest=0
      do iw=1,nwalk
        ! Loop over primary walker
#if defined(TREXIO_FOUND) && defined(QMCKL_FOUND)        
        rc = qmckl_set_point(qmckl_ctx(qmckl_no_ctx), 'N', nelec*1_8, xold_dmc(:,:,iw,1), nelec*3_8)
#endif

        call distances(0,xold_dmc(1,1,iw,1))
        ! Set nuclear coordinates and n-n potential (0 flag = no strech e-coord)
        if(nforce.gt.1) &
        call strech(xold_dmc(1,1,iw,1),xold_dmc(1,1,iw,1),ajacob,1,0)

        call walkstrdet(iw)
        call walkstrjas(iw)

        ! Sample Green function for forward move
        dfus2ac=zero
        dfus2un=zero
        drifdif=zero
        iaccept=0

        if (nfrag.gt.1) then
          dfus2acfrag = zero
          dfus2unfrag = zero
        endif

        dwt=1

        !     to initilialize pp
        pp=1.d0
        do i=1,nelec

          if(i.le.nup) then
            iflag_up=2
            iflag_dn=3
           else
            iflag_up=3
            iflag_dn=2
          endif
          call compute_determinante_grad(i,psido_dmc(iw,1),psido_dmc(iw,1),psijo_dmc(iw,1),vold_dmc(1,i,iw,1),1)

! Find nearest nucleus, vector from that nucleus to electron, component of velocity in that direction
          ren2mn=huge
          do ic=1,ncent
            ren2=(xold_dmc(1,i,iw,1)-cent(1,ic))**2 &
                +(xold_dmc(2,i,iw,1)-cent(2,ic))**2 &
                +(xold_dmc(3,i,iw,1)-cent(3,ic))**2
            if(ren2.lt.ren2mn) then
              ren2mn=ren2
              iwnuc=ic
            endif
          enddo
          rmino=zero
          voldr=zero
          v2old=zero
          do k=1,3
            rvmino(k)=xold_dmc(k,i,iw,1)-cent(k,iwnuc)
            rmino=rmino+rvmino(k)**2
            voldr=voldr+vold_dmc(k,i,iw,1)*rvmino(k)
            v2old=v2old+vold_dmc(k,i,iw,1)**2
          enddo
          rmino=sqrt(rmino)
          voldr=voldr/rmino
          volda=sqrt(v2old)
          zeta=dsqrt(one/tau+znuc(iwctype(iwnuc))**2)

! Place zaxis along direction from nearest nucleus to electron and
! x-axis along direction of angular component of velocity.
! Calculate the velocity in the phi direction
          voldp=zero
          do k=1,3
            zaxis(k)=rvmino(k)/rmino
            xaxis(k)=vold_dmc(k,i,iw,1)-voldr*zaxis(k)
            voldp=voldp+xaxis(k)**2
          enddo
          voldp=sqrt(voldp)
          if(voldp.lt.eps) then
            xaxis(1)=eps*(one-zaxis(1)**2)
            xaxis(2)=eps*(-zaxis(1)*zaxis(2))
            xaxis(3)=eps*(-zaxis(1)*zaxis(3))
            voldp=eps*dsqrt(one+eps-zaxis(1)**2)
          endif
          do k=1,3
            xaxis(k)=xaxis(k)/voldp
          enddo

! Use more accurate formula for the drift
          hafzr2=(half*znuc(iwctype(iwnuc))*rmino)**2
          adrift=(half*(1+eps+voldr/volda))+adrift0*hafzr2/(1+hafzr2)

! Tau primary -> tratio=one
          vavvt=(dsqrt(one+two*adrift*v2old*tau)-one)/(adrift*v2old)

          driftr=vavvt*voldr
          rtry=rmino+driftr

! Prob. of sampling exponential rather than gaussian is
! half*derfc(rtry/dsqrt(two*tau)) = half*(one+derf(-rtry/dsqrt(two*tau)))
! We use both expressions because under AIX the former is faster if rtry>0
! and the latter is faster if rtry<0.
! Note that if adrift is always set to 1 then it may be better to use
! vavvt rather than tau since the max drift distance is dsqrt(2*tau/adrift),
! so both the position of the gaussian and its width are prop to dsqrt(tau)
! if tau is used in derfc, and so qgaus does not tend to 1 for large tau.
! However, we are using a variable adrift that can be very small and then
! using tau is the better choice.

          if(rtry.gt.zero) then
            qgaus=half*derfc(rtry/dsqrt(two*tau))

! Calculate drifted x and y coordinates in local coordinate system centered
! on nearest nucleus
            xprime=vavvt*voldp*rtry/(half*(rmino+rtry))
            zprime=rtry

! Convert back to original coordinate system
            do k=1,3
              xnew(k)=cent(k,iwnuc)+xaxis(k)*xprime+zaxis(k)*zprime
            enddo
           else
            qgaus=half*(one+derf(-rtry/dsqrt(two*tau)))
            rtry=zero
            do k=1,3
              xnew(k)=cent(k,iwnuc)
            enddo
          endif
          pgaus=one-qgaus

          if(ipr.ge.1) write(6,'(''xnewdr'',2i4,9f8.5)') iw,i,(xnew(k),k=1,3)

! Do the diffusion
! Sample gaussian with prob pgaus, exponential with prob. qgaus
          dr2=zero
          dfus2a=zero
          if(random_dp().lt.pgaus) then
            dfus2b=zero
            do k=1,3
              drift=xnew(k)-xold_dmc(k,i,iw,1)
              dfus=gauss()*rttau
              dx=drift+dfus
              dr2=dr2+dx**2
              dfus2a=dfus2a+dfus**2
              xnew(k)=xnew(k)+dfus
              dfus2b=dfus2b+(xnew(k)-cent(k,iwnuc))**2
            enddo
            dfusb=sqrt(dfus2b)
           else
            dfusb=(-half/zeta)*dlog(random_dp()*random_dp()*random_dp())
            costht=two*(random_dp()-half)
            sintht=sqrt(one-costht*costht)
            phi=two*pi*random_dp()
            do k=1,3
              drift=xnew(k)-xold_dmc(k,i,iw,1)
              if(k.eq.1) then
                dfus=dfusb*sintht*dcos(phi)
               elseif(k.eq.2) then
                dfus=dfusb*sintht*dsin(phi)
               else
                dfus=dfusb*costht
              endif
              dx=drift+dfus
              dr2=dr2+dx**2
              dfus2a=dfus2a+(cent(k,iwnuc)+dfus-xnew(k))**2
              xnew(k)=cent(k,iwnuc)+dfus
            enddo
          endif
          dfus2o=dfus2a
          fnormo=(pgaus+qgaus*term*(zeta**3)*dexp(-two*zeta*dfusb+half*dfus2a/tau))

          if(ipr.ge.1) then
            write(ounit,'(''xold_dmc'',2i4,9f8.5)') iw,i,(xold_dmc(k,i,iw,1),k=1,3)
            write(ounit,'(''vold_dmc'',2i4,9f8.5)') iw,i,(vold_dmc(k,i,iw,1),k=1,3)
            write(ounit,'(''psido_dmc'',2i4,9f8.5)') iw,i,psido_dmc(iw,1)
            write(ounit,'(''xnewdr'',2i4,9f8.5)') iw,i,(xnew(k),k=1,3)
          endif

          ! calculate psi and velocity at new configuration
          call psiedmc(i,iw,xnew,psidn,psijn,0)

          call compute_determinante_grad(i,psidn(1),psidn,psijn,vnew(1,i),0)

          distance_node_ratio2=1.d0
          if(node_cutoff.gt.0) then
            do jel=1,nup
              if(jel.ne.i) call compute_determinante_grad(jel,psidn(1),psidn,psijn,vnew(1,jel),iflag_up)
            enddo

            do jel=nup+1,nelec
              if(jel.ne.i) call compute_determinante_grad(jel,psidn(1),psidn,psijn,vnew(1,jel),iflag_dn)
            enddo

            call nodes_distance(vold_dmc(1,1,iw,1),distance_node,1)
            rnorm_nodes_old=rnorm_nodes_num(distance_node,eps_node_cutoff)/distance_node

            call nodes_distance(vnew,distance_node,0)
            rnorm_nodes_new=rnorm_nodes_num(distance_node,eps_node_cutoff)/distance_node
            distance_node_ratio2=(rnorm_nodes_new/rnorm_nodes_old)**2
          endif

          ! Check for node crossings
          if(psidn(1)*psido_dmc(iw,1).le.zero) then
            nodecr=nodecr+1
            if(icross.le.0) then
              p=zero
              goto 160
            endif
          endif

          ! Calculate Green function for the reverse move

! Find the nearest nucleus, vector from that nucleus to electron, component of velocity in that direction
          ren2mn=huge
          do ic=1,ncent
            ren2=(xnew(1)-cent(1,ic))**2 &
                +(xnew(2)-cent(2,ic))**2 &
                +(xnew(3)-cent(3,ic))**2
            if(ren2.lt.ren2mn) then
              ren2mn=ren2
              iwnuc=ic
            endif
          enddo
          rminn=zero
          vnewr=zero
          v2new=zero
          do k=1,3
            rvminn(k)=xnew(k)-cent(k,iwnuc)
            rminn=rminn+rvminn(k)**2
            vnewr=vnewr+vnew(k,i)*rvminn(k)
            v2new=v2new+vnew(k,i)**2
          enddo
          rminn=sqrt(rminn)
          vnewr=vnewr/rminn
          vnewa=sqrt(v2new)
          zeta=dsqrt(one/tau+znuc(iwctype(iwnuc))**2)

! Place zaxis along direction from nearest nucleus to electron and
! x-axis along direction of angular component of velocity.
! Calculate the velocity in the phi direction
          vnewp=zero
          do k=1,3
            zaxis(k)=rvminn(k)/rminn
            xaxis(k)=vnew(k,i)-vnewr*zaxis(k)
            vnewp=vnewp+xaxis(k)**2
          enddo
          vnewp=sqrt(vnewp)
          if(vnewp.lt.eps) then
            xaxis(1)=eps*(one-zaxis(1)**2)
            xaxis(2)=eps*(-zaxis(1)*zaxis(2))
            xaxis(3)=eps*(-zaxis(1)*zaxis(3))
            vnewp=eps*dsqrt(one+eps-zaxis(1)**2)
          endif
          do k=1,3
            xaxis(k)=xaxis(k)/vnewp
          enddo

! Use more accurate formula for the drift
          hafzr2=(half*znuc(iwctype(iwnuc))*rminn)**2
          adrift=(half*(1+eps+vnewr/vnewa))+adrift0*hafzr2/(1+hafzr2)

          vavvt=(dsqrt(one+two*adrift*v2new*tau)-one)/(adrift*v2new)

          driftr=vavvt*vnewr
          rtry=rminn+driftr
          dfus2a=zero
          dfus2b=zero
          if(rtry.gt.zero) then
            qgaus=half*derfc(rtry/dsqrt(two*tau))

! Calculate drifted x and y coordinates in local coordinate system centered
! on nearest nucleus
            xprime=vavvt*vnewp*rtry/(half*(rminn+rtry))
            zprime=rtry

! Convert back to original coordinate system
            do k=1,3
              xbac(k)=cent(k,iwnuc)+xaxis(k)*xprime+zaxis(k)*zprime
              dfus2b=dfus2b+(cent(k,iwnuc)-xold_dmc(k,i,iw,1))**2
              dfus=xbac(k)-xold_dmc(k,i,iw,1)
              dfus2a=dfus2a+dfus**2
            enddo
           else
            qgaus=half*(one+derf(-rtry/dsqrt(two*tau)))
            rtry=zero
            do k=1,3
              xbac(k)=cent(k,iwnuc)
              dfus=xbac(k)-xold_dmc(k,i,iw,1)
              dfus2b=dfus2b+dfus**2
            enddo
            dfus2a=dfus2b
          endif
          dfus2n=dfus2a
          pgaus=one-qgaus
          dfusb=sqrt(dfus2b)

          fnormn=pgaus+qgaus*term*(zeta**3)*dexp(-two*zeta*dfusb+half*dfus2a/tau)

          if(ipr.ge.1) then
            write(ounit,'(''xold_dmc'',9f10.6)')(xold_dmc(k,i,iw,1),k=1,3), &
            (xnew(k),k=1,3), (xbac(k),k=1,3)
            write(ounit,'(''dfus2o'',9f10.6)')dfus2o,dfus2n, &
            psido_dmc(iw,1),psidn,psijo_dmc(iw,1),psijn(1),fnormo,fnormn
          endif

          p=(psidn(1)/psido_dmc(iw,1))**2*exp(2*(psijn(1)-psijo_dmc(iw,1)))* &
          exp((dfus2o-dfus2n)/(two*tau))*fnormn/fnormo*distance_node_ratio2

          if(ipr.ge.1) write(6,'(''p'',11f10.6)') &
          p,(psidn/psido_dmc(iw,1))**2*exp(2*(psijn(1)-psijo_dmc(iw,1))), &
          exp((dfus2o-dfus2n)/(two*tau)),psidn,psido_dmc(iw,1), &
          psijn(1),psijo_dmc(iw,1),dfus2o,dfus2n,fnormo,fnormn

          !          if(ipr.ge.1) write(ounit,'(''parts p'',11f10.6)')
          !     &         psidn(1), psijn, psido_dmc(iw,1), dfus2o, dfus2n, distance_node_ratio2

          ! Way to cure persistent configurations; not needed if itau_eff <=0; in practice never needed
          if(iage(iw).gt.50) p=p*1.1d0**(iage(iw)-50)

          pp=pp*p
          p=dmin1(one,p)
          160     q=one-p

          acc=acc+p
          trymove=trymove+1
          dfus2ac=dfus2ac+p*dfus2o
          dfus2un=dfus2un+dfus2o
        
          if (nfrag.gt.1) then
            dfus2acfrag(ifragelec(i)) = dfus2acfrag(ifragelec(i)) + p * dfus2o
            dfus2unfrag(ifragelec(i)) = dfus2unfrag(ifragelec(i)) + dfus2o
          endif

          ! If we are using weights rather than accept/reject
          if(iacc_rej.le.0) then
            p=one
            q=zero
          endif

          iacc_elec(i)=0
          if(random_dp().lt.p) then
#if defined(TREXIO_FOUND) && defined(QMCKL_FOUND) 
            rc = qmckl_get_jastrow_champ_single_accept(qmckl_ctx(qmckl_no_ctx))
#endif
            iaccept=1
            nacc=nacc+1
            iacc_elec(i)=1
            if(ipq.le.0) p=one

            iage(iw)=0
            do k=1,3
              drifdif=drifdif+(xold_dmc(k,i,iw,1)-xnew(k))**2
              xold_dmc(k,i,iw,1)=xnew(k)
            enddo
            psido_dmc(iw,1)=psidn(1)
            psijo_dmc(iw,1)=psijn(1)
            call jassav(i,0)
            call detsav(i,0)

           else
            if(ipq.le.0) p=zero
            call distancese_restore(i)
          endif
          q=one-p

          ! Calculate moments of r and save rejection probability for primary walk
          unacp(i)=q

          call update_ymat(i)

        enddo

        ! Effective tau for branching
        tauprim=tau*dfus2ac/dfus2un
        tauefffrag = tau * dfus2acfrag / dfus2unfrag
      
        do ifr=1,nforce

          if(ifr.eq.1) then
            ! Primary configuration
            drifdifr=one
            if(nforce.gt.1) &
            call strech(xold_dmc(1,1,iw,1),xold_dmc(1,1,iw,1),ajacob,1,0)
            if(nfrag.gt.1) call assign_elecs(xold_dmc(:,:,iw,1),2)
            call hpsi(xold_dmc(1,1,iw,1),psidn(1),psijn,ekino,enew,ipass,1)
            !if (abs(enew(1)- sum(eloc_i(:))) > 1e-9) print *, 'hpsi', enew(1)- sum(eloc_i(:))

            if(irun.eq.1) then
               wtg_sqrt(1)=dsqrt(wtg(1))
               call optwf_store(lpass,wtg(1),wtg_sqrt(1),psidn(1),enew(1))
            endif


            call walksav_det(iw)
            call walksav_jas(iw)
            rnorm_nodes=1.d0
            if(node_cutoff.gt.0) then
              call nodes_distance(vold_dmc(1,1,iw,1),distance_node,1)
              rnorm_nodes=rnorm_nodes_num(distance_node,eps_node_cutoff)/distance_node
            endif
           else
            ! Secondary configuration
            if(istrech.eq.0) then
              call strech(xold_dmc(1,1,iw,ifr),xold_dmc(1,1,iw,ifr),ajacob,ifr,0)
              drifdifr=one
              ! No streched positions for electrons
              do i=1,nelec
                do k=1,3
                  xold_dmc(k,i,iw,ifr)=xold_dmc(k,i,iw,1)
                enddo
              enddo
              ajacold(iw,ifr)=one
             else
              ! Compute streched electronic positions for all nucleus displacement
              call strech(xold_dmc(1,1,iw,1),xstrech,ajacob,ifr,1)
              drifdifs=zero
              do i=1,nelec
                do k=1,3
                  drifdifs=drifdifs+(xstrech(k,i)-xold_dmc(k,i,iw,ifr))**2
                  xold_dmc(k,i,iw,ifr)=xstrech(k,i)
                enddo
              enddo
              ajacold(iw,ifr)=ajacob
              if(drifdif.eq.0.d0) then
                drifdifr=one
               else
                drifdifr=drifdifs/drifdif
              endif
            endif
            call hpsi(xold_dmc(1,1,iw,ifr),psidn,psijn,ekino,enew,ipass,ifr)
          endif

          do i=1,nelec
              call compute_determinante_grad(i,psidn(1),psidn,psijn,vold_dmc(1,i,iw,ifr),1)
          enddo

          tratio=one
          if(ifr.gt.1.and.itausec.eq.1) tratio=drifdifr

          taunow=tauprim*drifdifr

          if(idmc.gt.0) then

          call calculate_fratio(icut_e, adrift0, tratio, taunow, sqrt_nelec, &
               xold_dmc(:,:,iw,ifr), vold_dmc(:,:,iw,ifr), eest, enew(1), fration, &
               eest_i, eloc_i(:), fration_i(:), eestfrag, elocfrag(:), frationfrag(:))
          
          if(ipr.ge.1) write(ounit,'(''wt'',9f10.5)') wt(iw),etrial,eest
          
          call calculate_reweight(idmc, icut_e, icut_br, taunow, e_cutoff, &
               etrial, eest, eold(iw,ifr), enew(1), fratio(iw,ifr), fration, &
               eest_i, eloco_i(:,iw,ifr), eloc_i, fratio_i(:,iw,ifr), fration_i, &
               etrialfrag, eestfrag, elocofrag(:,iw,ifr), elocfrag, fratiofrag(:,iw,ifr), frationfrag, tauefffrag, dwt)

          if(ipr.ge.1) write(ounit,*) 'dwt', taunow, fratio(iw,ifr), fration, etrial, eest, eold(iw,ifr), enew(1), dwt
          !if (abs(dwt - 1) > 0.1) then 
          !  print*, 'dwt:', dwt
          !  print*, 'ewto:', ewto, 'ewtn:', ewtn
          !  print*, 'deo:', deo_i(:)
          !  print*, 'den:', den_i(:)
          !  print*, ' '
          !endif

          ! Limit the weights for LA
          if(limit_wt_dmc.gt.0) then
            dwt_cutoff=exp((etrial-eest+limit_wt_dmc*esigma/sqrt(rttau))*taunow)
            if(dwt.gt.dwt_cutoff) dwt=dwt_cutoff
          endif

          endif

          ! If we are using weights rather than accept/reject
          if(iacc_rej.eq.0) dwt=dwt*pp

          ! Exercise population control if dmc or vmc with weights
          if(idmc.gt.0.or.iacc_rej.eq.0) dwt=dwt*ffi

          ! Set weights and product of weights over last nwprod steps
          if(ifr.eq.1) then

            wt(iw)=wt(iw)*dwt
            wtnow=wt(iw)
            pwt(iw,ifr)=pwt(iw,ifr)+log(dwt)-wthist(iw,iwmod,ifr)
            wthist(iw,iwmod,ifr)=dlog(dwt)

           elseif(ifr.gt.1) then

            pwt(iw,ifr)=pwt(iw,ifr)+dlog(dwt)-wthist(iw,iwmod,ifr)
            wthist(iw,iwmod,ifr)=dlog(dwt)
            wtnow=wt(iw)*dexp(pwt(iw,ifr)-pwt(iw,1))

          endif

          wtnow=wtnow/rnorm_nodes**2

          if(ipr.ge.1)write(ounit,'(''eold,enew,wt'',9f10.5)') &
          eold(iw,ifr),enew,wtnow

          if(idmc.gt.0) then
            wtg=wtnow*fprod
           else
            wtg=wtnow
          endif
          tausum(ifr)=tausum(ifr)+wtg(1)*taunow

          if(ipr.gt.5.and.dabs((enew(1)-etrial)/etrial).gt.0.2d+0) then
           write(18,'(i6,f8.2,2d10.2,(8f8.4))') ipass,  &
            enew(1)-etrial,psidn,psijn(1),(xnew(ii),ii=1,3)
          endif

          if(ipr.gt.5.and.wt(iw).gt.3) write(18,'(i6,i4,3f8.2,30f8.4)') ipass,iw, &
            wt(iw),enew(1)-etrial,eold(iw,ifr)-etrial,(xnew(ii),ii=1,3)

          if(iforce_analy.eq.1) then
            if (ipathak.gt.0) then
              call nodes_distance(vold_dmc(1,1,iw,ifr),distance_node,1)
              do iph=1,PTH
                call pathak(distance_node,pnew(iph),eps_pathak(iph))
              enddo
            endif
            if (dmc_ivd.gt.0) then
              if (icut_e.eq.3) then 
                fratio_aux = max(fration*taunow, 1e-9)
                call force_analy_vd(enew(1), sqrt_pi_o2 * derf(fratio_aux)/(fratio_aux), taunow, iw, iwmod)
              else
                call force_analy_vd(ecutn, ecuto, e_cutoff, iw, iwmod)
              endif
            endif

            do iph=1,PTH
              do ic=1,ncent
                do k=1,3
                  if (ipathak.gt.0) then
                    derivsum(1,k,ic,iph)=derivsum(1,k,ic,iph)+wtg(1)*da_energy(k,ic)*pnew(iph)
                    derivsum(2,k,ic,iph)=derivsum(2,k,ic,iph)+wtg(1)*enew(1)*da_psi(k,ic)*pnew(iph)
                    derivsum(3,k,ic,iph)=derivsum(3,k,ic,iph)+wtg(1)*da_psi(k,ic)*pnew(iph)
                  else
                    derivsum(1,k,ic,iph)=derivsum(1,k,ic,iph)+wtg(1)*da_energy(k,ic)
                    derivsum(2,k,ic,iph)=derivsum(2,k,ic,iph)+wtg(1)*enew(1)*da_psi(k,ic)
                    derivsum(3,k,ic,iph)=derivsum(3,k,ic,iph)+wtg(1)*da_psi(k,ic)
                  endif
                enddo
              enddo
            enddo
          endif
          
          ! call update_reweight
          eold(iw,ifr)=enew(1)
          !if (abs(sum(eloc_i(:)))<1e-6) print *, 'eloc zero'
          if (icut_e.lt.0) then
            eloco_i(:,iw,ifr) = eloc_i(:)
            fratio_i(:,iw,ifr)=fration_i(:)
          end if 
          
          if (nfrag.gt.1) then 
            elocofrag(:,iw,ifr) = elocfrag(:)
            fratiofrag(:,iw,ifr)=frationfrag(:)
          endif

          psido_dmc(iw,ifr)=psidn(1)
          psijo_dmc(iw,ifr)=psijn(1)
          fratio(iw,ifr)=fration
          
          call prop_save_dmc(iw)
          call pcm_save(iw)
          call mmpol_save(iw)
          call force_analy_save

          if(ifr.eq.1) then
            if(iaccept.eq.0) then
              iage(iw)=iage(iw)+1
              ioldest=max(ioldest,iage(iw))
              ioldestmx=max(ioldestmx,iage(iw))
            endif

            psi2savo=2*(dlog(dabs(psido_dmc(iw,1)))+psijo_dmc(iw,1))

            wsum1(ifr)=wsum1(ifr)+wtnow
            esum1_dmc(ifr)=esum1_dmc(ifr)+wtnow*eold(iw,ifr)
            if (icut_e.lt.0) then
              esum_i1(:) = esum_i1(:) + wtnow*eloc_i(:)
            endif
            if (nfrag.gt.1) then
              esumfrag1(:) = esumfrag1(:) + wtnow*elocfrag(:)
            endif

            !print *, 'esum', esum1_dmc(1),  sum(esum_i1(:))!, esum1_dmc(1) - sum(esum_i(:))
            pesum_dmc(ifr)=pesum_dmc(ifr)+wtg(1)*(eold(iw,ifr)-ekino(1))
            tpbsum_dmc(ifr)=tpbsum_dmc(ifr)+wtg(1)*ekino(1)

            call prop_sum_dmc(0.d0,wtg(1),iw)
            call pcm_sum(0.d0,wtg(1),iw)
            call mmpol_sum(0.d0,wtg(1),iw)
            call force_analy_sum(wtg(1),0.d0,eold(iw,1),0.0d0)

            call optjas_sum(wtg,zero_1d,eold(iw,1),eold(iw,1),0)
            call optorb_sum(wtg,zero_1d,eold(iw,1),eold(iw,1),0)
            call optci_sum(wtg(1),0.d0,eold(iw,1),eold(iw,1))

            call optx_jas_orb_sum(wtg,zero_1d,0)
            call optx_jas_ci_sum(wtg(1),0.d0,eold(iw,1),eold(iw,1))
            call optx_orb_ci_sum(wtg(1),0.d0)

          else

            ro=ajacold(iw,ifr)*psido_dmc(iw,ifr)**2*exp(2*psijo_dmc(iw,ifr)-psi2savo)

            wsum1(ifr)=wsum1(ifr)+wtnow*ro
            esum1_dmc(ifr)=esum1_dmc(ifr)+wtnow*eold(iw,ifr)*ro
            pesum_dmc(ifr)=pesum_dmc(ifr)+wtg(1)*(eold(iw,ifr)-ekino(1))*ro
            tpbsum_dmc(ifr)=tpbsum_dmc(ifr)+wtg(1)*ekino(1)*ro

            wtg=wt(iw)*fprod/rnorm_nodes**2
            wtg_derivsum1=wtg(1)
          endif
        enddo

      enddo

      if(ipr.gt.5.and.wsum1(1).gt.1.1d0*dmc_nconf) write(18,'(i6,9d12.4)') ipass,ffn,fprod,fprod/ff(ipmod2),wsum1(1),wgdsumo

      if(idmc.gt.0.or.iacc_rej.eq.0) then
        wfsum1=wsum1(1)*ffn
        efsum1=esum1_dmc(1)*ffn
      endif
      do ifr=1,nforce
        if(idmc.gt.0.or.iacc_rej.eq.0) then
          wgsum1(ifr)=wsum1(ifr)*fprod
          egsum1(ifr)=esum1_dmc(ifr)*fprod
         else
          wgsum1(ifr)=wsum1(ifr)
          egsum1(ifr)=esum1_dmc(ifr)
        endif
      enddo

      if (icut_e.lt.0) then
        if (idmc.gt.0.or.iacc_rej.eq.0) then
          esum_i(:) = esum_i(:) + esum_i1(:)*fprod
        endif
      endif
      if (nfrag.gt.1) then
        if (idmc.gt.0.or.iacc_rej.eq.0) then
          egsum1frag(:)=esumfrag1(:) * fprod
          esumfrag(:) = esumfrag(:) + esumfrag1(:)*fprod
        else 
          egsum1frag(:)=esumfrag1(:)
        endif
      endif

      if(idmc.gt.0) call splitj
      if(write_walkalize) write(11,'(i8,f9.6,f12.5,f11.6,i5,f11.5)') ipass,ffn, &
      wsum1(1),esum1_dmc(1)/wsum1(1),nwalk

      return
      endsubroutine

      end module
