      module init_mod
      contains
      subroutine init
c MPI version created by Claudia Filippi starting from serial version
c routine to accumulate estimators for energy etc.

      use dmc_mod, only: MFPRD1
      use const, only: etrial
      use forcepar, only: istrech, nforce
      use system, only: cent, iwctype, ncent, pecent, znuc
      use estcum, only: ipass
      use config, only: psido_dmc, psijo_dmc, vold_dmc, xold_dmc
      use force_dmc, only: nwprod
      use mpiconf, only: nproc
      use qua, only: nquad, wq, xq, yq, zq
      use branch, only: eest, eigv, eold, ff, fprod, nwalk, pwt, wdsumo, wgdsumo, wt, wtgen
      use branch, only: wthist
      use casula, only: i_vpsp, icasula
      use jacobsave, only: ajacob, ajacold
      use pseudo, only: nloc
!      use contrl, only: nconf
      use control_dmc, only: dmc_nconf
      use mpi

      use precision_kinds, only: dp

      use rotqua_mod, only: gesqua
      use pot, only: pot_nn
      use zerest_mod, only: zerest
      use strech_mod, only: strech
      use mmpol_dmc, only: mmpol_save
      use pcm_dmc, only: pcm_save
      use prop_dmc, only: prop_save_dmc
      use nonloc_grid_mod, only: t_vpsp_sav
      use hpsi_mod, only: hpsi
      use walksav_det_mod, only: walksav_det
      use walksav_jas_mod, only: walksav_jas
      use determinante_mod, only: compute_determinante_grad
      use control, only: mode
      use system, only: nelec
      implicit none

      integer :: i, ie, ifr, ip, iw
      integer :: k

      real(dp), parameter :: zero = 0.d0
      real(dp), parameter :: one = 1.d0


c Initialize various quantities at beginning of run
c the initial values of energy psi etc. are calculated here

      ipass=0

c set quadrature points
      if(nloc.gt.0) call gesqua (nquad,xq,yq,zq,wq)

c get nuclear potential energy
      call pot_nn(cent,znuc,iwctype,ncent,pecent)

      eigv=one
      eest=etrial
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
