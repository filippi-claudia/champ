      module optwf_mix_sa
      contains
      subroutine optwf_mix

      use sr_mod, only: mparm
      use csfs, only: nstates
      use mstates_mod, only: MSTATES
      use optwf_corsam, only: energy, energy_err
      use sa_check, only: energy_all, energy_err_all
!      use contrl, only: nblk, nblk_max, nblk_ci
      use control_vmc, only: vmc_nblk, vmc_nblk_max, vmc_nblk_ci
      use orbval, only: nadorb
      use force_analy, only: iforce_analy, alfgeo
      use mstates_ctrl, only: iguiding
      use method_opt, only: method
      use optwf_sr_mod, only: sr
      use optwf_corsam, only: sigma
      use force_analy, only: iforce_analy
      use optwf_contrl, only: ioptci, ioptjas, ioptorb, nparm
      use optwf_contrl, only: iroot_geo
      use optwf_contrl, only: dparm_norm_min, nopt_iter, micro_iter_sr
      use optwf_contrl, only: sr_tau, sr_adiag, sr_eps
      use optwf_contrl, only: nvec, nvecx, alin_adiag, alin_eps
      use precision_kinds, only: dp
      use contrl_file,    only: ounit

      use error, only: fatal_error
      use optwf_handle_wf,only: save_nparms, write_wf, restore_wf
      use optwf_handle_wf,only: set_nparms, save_wf, compute_parameters
      use optwf_handle_wf,only: test_solution_parm, save_ci_best
      use optwf_handle_wf,only: restore_ci_best, set_nparms_tot
      use optgeo_lib, only: write_geometry, compute_positions
      use optwf_lin_dav_extra, only: select_ci_root
      use optwf_lin_dav_more,  only: lin_d
      use sr_more, only: dscal
      implicit none
      interface
      subroutine qmc
      end subroutine
      end interface

      integer :: i, iflag, iforce_analy_sav, iguiding_sav, inc_nblk
      integer :: ioptci_sav, ioptjas_sav, ioptorb_sav, iqmc_again
      integer :: iqmc_check, istate, istate0, iter
      integer :: miter, nblk_sav, nparmci, nstates_sav, nadorb_sav
      integer, dimension(MSTATES) :: i_deltap
      integer, dimension(5) :: index_min_energy
      integer, dimension(5,MSTATES) :: index_more
      real(dp) :: adiag, alin_adiag_sav, denergy, denergy_err, diffene
      real(dp) :: dparm_norm, energy_err_sav, energy_sav, errdiff
      real(dp) :: sigma_sav, sr_adiag_sav
      real(dp), dimension(mparm*MSTATES) :: deltap
      real(dp), dimension(mparm*MSTATES,5) :: deltap_more
      real(dp), dimension(MSTATES) :: energy_old
      real(dp), dimension(MSTATES) :: energy_err_old
      real(dp), dimension(6,MSTATES) :: energy_davidson
      real(dp), dimension(mparm) :: deltap_new
      character*20 method_sav

      save method_sav

      if(method.ne.'mix_n') return

      if(nstates.eq.1.or.ioptci.eq.0) call fatal_error('OPTWF_MIX: no need mix_n')

      call set_nparms_tot

      if(nparm.gt.mparm)call fatal_error('SR_OPTWF: nparmtot gt mparm')

      write(ounit,'(''Starting dparm_norm_min'',g12.4)') dparm_norm_min
      write(ounit,'(/,''SR adiag: '',f10.5)') sr_adiag
      write(ounit,'(''SR tau:   '',f10.5)') sr_tau
      write(ounit,'(''SR eps:   '',f10.5)') sr_eps


      write(ounit,'(/,''LIN_D adiag: '',f10.5)') alin_adiag
      write(ounit,'(''LIN_D ethr:  '',f10.5)') alin_eps
      write(ounit,'(''LIN_D nvec:  '',i4)') nvec
      write(ounit,'(''LIN_D nvecx: '',i4)') nvecx

      if(nstates.gt.1.and.nvec.lt.nstates) call fatal_error('SR_OPTWF: nvec < nstates')

      inc_nblk=0

      sr_adiag_sav=sr_adiag
      alin_adiag_sav=alin_adiag
      nblk_sav=vmc_nblk

      nstates_sav=nstates
      iforce_analy_sav=iforce_analy

      ioptjas_sav=ioptjas
      ioptorb_sav=ioptorb
      ioptci_sav=ioptci

      nadorb_sav=nadorb

      call save_nparms
      call write_geometry(0)
      call save_wf
      call save_ci_best

c do iteration
      do iter=1,nopt_iter
        write(ounit,'(/,''Optimization iteration'',i5,'' of'',i5)')iter,nopt_iter

        iforce_analy=0

c do micro_iteration
        do miter=1,micro_iter_sr

          if(micro_iter_sr.gt.1) write(ounit,'(/,''Micro iteration'',i5,'' of'',i5)')miter,micro_iter_sr

          method='sr_n'
          ioptci=0
          ioptorb=ioptorb_sav
          ioptjas=ioptjas_sav
          nblk_sav=vmc_nblk
          if(miter.eq.micro_iter_sr) then
            method='lin_d'
            ioptci=ioptci_sav

            vmc_nblk=vmc_nblk_ci
            write(ounit,'(''NBLOCK changed from '',i7, '' to '',i7)') nblk_sav,vmc_nblk

            ioptorb=0
            ioptjas=0
          endif
          call set_nparms
          if(iter.eq.1.and. miter.eq.micro_iter_sr) nparmci=nparm

c if the last step was a davidson then save the old energy before recomputing it (for the check)

          if(miter.eq.1.and.iter.gt.1) then
             do istate=1,nstates
               energy_old(istate)=energy_all(istate)
               energy_err_old(istate)=energy_err_all(istate)
             enddo
          endif

          iqmc_check=0

   5      call qmc
          vmc_nblk=nblk_sav
          write(ounit,'(/,''Completed sampling'')')

          if(miter.eq.1 .and. iter.gt.1) then
            if(iqmc_check.eq.0) then
              do istate=1,nstates
                i_deltap(istate)=0
              enddo
            endif

            if(iqmc_check.lt.2) then
              iqmc_again=0
              do istate=1,nstates
                diffene=abs(energy_all(istate)-energy_old(istate))
                errdiff=sqrt(energy_err_all(istate)**2+energy_err_old(istate)**2)

                if(diffene.ge.10*errdiff)then
                  i_deltap(istate)=i_deltap(istate)+1
                  istate0=(istate-1)*nparmci+1

                  if(i_deltap(istate).gt.5) call fatal_error('OPTWF_MIX: only 5 more deltap stored per state')
                  call change_ci(deltap_more(istate0,i_deltap(istate)),istate)
                  write(ounit,*) istate0,i_deltap(istate)

                  write(ounit,'(''STATE, N OVERLAP, ENERGY OLD, ENERGY NEW,10*ERRDIFF '',2i3,3f12.5)')
     &            istate,index_more(i_deltap(istate),istate),energy_old(istate),energy_all(istate),10*errdiff
                  iqmc_again=1
                endif

              enddo
              if(iqmc_again.gt.0) then
                iqmc_check=iqmc_check+1
                write(ounit,'(''Use new set of CI coefficients'',i4)')
                go to 5
               else
                call save_ci_best
                write(ounit,'(''Save current CI coefficients as best'')')
              endif
             else
              call restore_ci_best
              write(ounit,'(''Restore CI cofficients to previous iteration'')')

              call qmc
            endif
          endif

   6      continue

          if(method.eq.'sr_n') then
            call sr_hs(nparm,sr_adiag)
            call sr(nparm,deltap,sr_adiag,sr_eps,i)
            call dscal(nparm,-sr_tau,deltap,1)
            adiag=sr_adiag
           else
            call lin_d(nparm,nvec,nvecx,deltap,deltap_more,index_more,alin_adiag,alin_eps)
            if(nstates.eq.1) call dscal(nparm,-1.d0,deltap,1)
            adiag=alin_adiag
          endif
          call test_solution_parm(nparm,deltap,dparm_norm,dparm_norm_min,adiag,iflag)
          write(ounit,'(''Norm of parm variation '',g12.5)') dparm_norm
          if(iflag.ne.0) then
            write(ounit,'(''Warning: dparm_norm>1'')')
            adiag=10*adiag
            write(ounit,'(''adiag increased to '',f10.5)') adiag

            sr_adiag=adiag
            alin_adiag=adiag
            go to 6
           else
            sr_adiag=sr_adiag_sav
            alin_adiag=alin_adiag_sav
          endif

          call compute_parameters(deltap,iflag,1)
          call write_wf(1,iter)

          if(miter.eq.micro_iter_sr-1) call save_ci_best
          call save_wf
        enddo
c enddo micro_iteration

        if(iforce_analy_sav.eq.1) then

          call select_ci_root(iroot_geo)

          iforce_analy=iforce_analy_sav
          iguiding_sav=iguiding

          nstates=1
          iguiding=0

          ioptjas=0
          ioptorb=0
          ioptci=0

          call set_nparms
          call qmc

          call compute_positions
          call write_geometry(iter)

          nstates=nstates_sav
          iguiding=iguiding_sav

          ioptjas=ioptjas_sav
          ioptorb=ioptorb_sav
          ioptci=ioptci_sav

          call restore_wf(1)
        endif

        if(iter.ge.2) then
          denergy=energy(1)-energy_sav
          denergy_err=sqrt(energy_err(1)**2+energy_err_sav**2)
c         call check_length_run_sr(iter,inc_nblk,vmc_nblk,vmc_nblk_max,denergy,denergy_err,energy_err_sav,energy_tol)
          vmc_nblk=vmc_nblk*1.2
          vmc_nblk=min(vmc_nblk,vmc_nblk_max)
c         if(-denergy.gt.3*denergy_err) alfgeo=alfgeo/1.2
        endif
        write(ounit,'(''nblk = '',i6)') vmc_nblk
        write(ounit,'(''alfgeo = '',f10.4)') alfgeo

        energy_sav=energy(1)
        energy_err_sav=energy_err(1)
        sigma_sav=sigma
      enddo
c enddo iteration

      write(ounit,'(/,''Check last iteration'')')

      ioptjas=0
      ioptorb=0
      ioptci=0
      iforce_analy=0

      call set_nparms

      call qmc

      nadorb=nadorb_sav
      call write_wf(1,-1)
      call write_geometry(-1)

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine change_ci(dparm_new,istate)

      use csfs, only: ccsf, cxdet, iadet, ibdet, icxdet, ncsf

      use dets, only: cdet, ndet
      use precision_kinds, only: dp
      use contrl_file,    only: ounit

      implicit none

      integer :: icsf, idet, istate, j, jx
      real(dp), dimension(*) :: dparm_new

c      write(ounit,*) "COPUTING NEW CI, ccsf(1,state,1)", ccsf(1,istate,1), dparm_new(1)
c update the ci coef

       if(ncsf.eq.0) then
         do idet=1,ndet
          cdet(idet,istate,1)=dparm_new(idet)
         enddo
       else
         do idet=1,ndet
           cdet(idet,istate,1)=0.d0
         enddo
         do icsf=1,ncsf
           do j=iadet(icsf),ibdet(icsf)
              jx=icxdet(j)
              cdet(jx,istate,1)=cdet(jx,istate,1)+dparm_new(icsf)*cxdet(j)
           enddo
           ccsf(icsf,istate,1)=dparm_new(icsf)
         enddo
       endif
c      write(ounit,*) "COPUTING NEW CI, ccsf(1,istate,1)", ccsf(1,istate,1), dparm_new(1)

      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      end module
