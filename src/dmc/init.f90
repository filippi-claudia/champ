module init_mod
contains
      subroutine init
! MPI version created by Claudia Filippi starting from serial version
! routine to accumulate estimators for energy etc.

      use system, only: ncent
      use assignment_mod, only: assign_elecs
      use branch,  only: eest,esigma,eigv,eold,ff,fprod,nwalk,pwt,wdsumo
      use branch,  only: wgdsumo,wt,wtgen,wthist
      use casula,  only: i_vpsp,icasula
      use config,  only: psido_dmc,psijo_dmc,vold_dmc,xold_dmc
      use const,   only: etrial,esigmatrial
      use contrldmc, only: icut_e
      use control, only: mode
      use control_dmc, only: dmc_nconf
      use da_energy_now, only: da_energy
      use determinante_mod, only: compute_determinante_grad
      use dmc_mod, only: MFPRD1
      use estcum,  only: ipass
      use ewald, only: cos_n_sum, sin_n_sum
      use force_pth, only: PTH
      use force_analytic, only: force_analy_save
      use fragments,  only: eest_i, eloco_i, eloc_i ! electron fragments
      use fragments,  only: eestfrag, elocofrag, elocfrag, nfrag, ifragelec, sqrt_nelecfrag, etrialfrag ! fragments
      use hpsi_mod, only: hpsi
      use jacobsave, only: ajacob,ajacold
      use mmpol_dmc, only: mmpol_save
      use mpi
      use mpiconf, only: nproc
      use mpitimer, only: elapsed_time
      use multiple_geo, only: istrech,nforce,nwprod,pecent
      use m_force_analytic, only: iforce_analy
      use nodes_distance_mod,   only: nodes_distance
      use nonloc_grid_mod, only: t_vpsp_sav
      use pcm_dmc, only: pcm_save
      use pot, only: pot_nn
      use precision_kinds, only: dp
      use prop_dmc, only: prop_save_dmc
      use pseudo, only: nloc
      use qua,only: nquad,wq,xq,yq,zq
      use rotqua_mod, only: gesqua
      use strech_mod, only: strech
      use system, only: cent,iwctype,ncent,nelec,znuc
      use vd_mod, only: deriv_eold, esnake, ehist, dmc_ivd
      use walksav_det_mod, only: walksav_det
      use walksav_jas_mod, only: walksav_jas
      use zerest_mod, only: zerest
      use pathak_mod, only: init_eps_pathak, pathak
      use pathak_mod, only: ipathak, eps_pathak, pold
!      use contrl, only: nconf

      implicit none

      integer :: i, ie, ifr, ip, iw, ic
      integer :: k, iph

      real(dp) :: ekino(1)
      real(dp), parameter :: zero = 0.d0
      real(dp), parameter :: one = 1.d0
      real(dp) :: distance_node

! Initialize various quantities at beginning of run
! the initial values of energy psi etc. are calculated here

      ipass=0

! set quadrature points
      if(nloc.gt.0) call gesqua (nquad,xq,yq,zq,wq)

      
! get nuclear potential energy
      call pot_nn(cent,znuc,iwctype,ncent,pecent,cos_n_sum,sin_n_sum)

      eestfrag = 0.d0
      eigv=one
      eest=etrial

      if (nfrag.gt.1) then 
        call assign_elecs(xold_dmc(:,:,1,1), 2)
        sqrt_nelecfrag = 0
        do i = 1, nelec
          sqrt_nelecfrag(ifragelec(i)) = sqrt_nelecfrag(ifragelec(i)) + 1
          ! Inital value for eest maybe think of something better
          eestfrag(ifragelec(i)) = etrialfrag(ifragelec(i))
          !eestfrag(ifragelec(i)) = eestfrag(ifragelec(i)) + etrial/nelec
          !etrialfrag(ifragelec(i)) = etrialfrag(ifragelec(i)) + etrial/nelec
        end do
        sqrt_nelecfrag = sqrt(sqrt_nelecfrag)
      end if
      if (icut_e.lt.1) then
        do i = 1, nelec
          eest_i(i) = etrial/nelec
        enddo
      end if
      
      esigma=esigmatrial
      nwalk=dmc_nconf
      fprod=one

      do iw=1,dmc_nconf
        wt(iw)=one
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
          
          call hpsi(xold_dmc(1,1,iw,ifr),psido_dmc(iw,ifr),psijo_dmc(iw,ifr),ekino,eold(iw,ifr),0,ifr)
          if (icut_e.lt.0) then
            eloco_i(:,iw,ifr) = eloc_i(:)
          endif
          if (nfrag.gt.1) then
            elocofrag(:,iw,ifr) = elocfrag(:)
          endif

          i_vpsp=0
          do i=1,nelec !STU check psijo_dmc, should be one state so should be ok?
            call compute_determinante_grad(i,psido_dmc(iw,ifr),psido_dmc(iw,ifr),psijo_dmc(iw,ifr),vold_dmc(1,i,iw,ifr),1)
          enddo

          if(ifr.eq.1) then
            call walksav_det(iw)
            call walksav_jas(iw)
!           call t_vpsp_sav(iw)
            call t_vpsp_sav
            call prop_save_dmc(iw)
            call pcm_save(iw)
            call mmpol_save(iw)
            if(iforce_analy.eq.1) then
              call force_analy_save
              if(ipathak.gt.0) then
                if (iw.eq.1) call init_eps_pathak()
                call nodes_distance(vold_dmc(1,1,iw,ifr), distance_node, 1)
              endif
              if(dmc_ivd.gt.0) then
                do iph=1,PTH
                  do ic=1,ncent
                    do k=1,3
                      esnake(k,ic,iw,iph)=zero
                      do ip=0,nwprod-1
                        ehist(k,ic,iw,ip,iph)=zero
                      enddo
                    enddo
                  enddo       
                enddo
                do ic=1,ncent
                  do k=1,3
                    deriv_eold(k,ic,iw)=da_energy(k,ic)
                  enddo
                enddo           
              endif
              do iph=1,PTH
                if(ipathak.gt.0) call pathak(distance_node,pold(iw,iph),eps_pathak(iph))
              enddo
            endif
          endif
          pwt(iw,ifr)=0
          do ip=0,nwprod-1
            wthist(iw,ip,ifr)=0
          enddo
        enddo
      enddo


      if(mode.eq.'dmc_one_mpi2') dmc_nconf=dmc_nconf*nproc
      wdsumo=dmc_nconf
      wgdsumo=dmc_nconf
      do i=0,MFPRD1
        wtgen(i)=dmc_nconf
        ff(i)=one
      enddo

      call zerest

      return
      end
end module
