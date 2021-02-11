      subroutine optwf_mix

      use sr_mod, only: MPARM, MVEC, nvecx
      use csfs, only: nstates
      use mstates_mod, only: MSTATES
      use optwf_corsam, only: energy, energy_err, force
      use sa_check, only: energy_all, energy_err_all
      use contrl, only: nblk, nblk_max, nblk_ci
      use force_analy, only: iforce_analy, alfgeo
      use mstates_ctrl, only: iguiding
      use method_opt, only: method
      use optwf_sr_mod, only: sr
      use force_analy, only: iforce_analy
      use optwf_contrl, only: ioptci, ioptjas, ioptorb, nparm
      use optwf_contrl, only: iroot_geo
      use optwf_contrl, only: energy_tol, dparm_norm_min, nopt_iter, micro_iter_sr
      use optwf_contrl, only: sr_tau , sr_adiag, sr_eps 
      use optwf_contrl, only: nvec, nvecx, alin_adiag, alin_eps

      implicit real*8(a-h,o-z)

      character*20 method_sav

      
      dimension deltap(MPARM*MSTATES),deltap_more(MPARM*MSTATES,5)
      dimension energy_old(MSTATES), energy_err_old(MSTATES), i_deltap(MSTATES), energy_davidson(6,MSTATES) 
      dimension index_min_energy(5), deltap_new(MPARM)
      dimension index_more(5,MSTATES)
      
      save method_sav

      if(method.ne.'mix_n') return

      if(nstates.eq.1.or.ioptci.eq.0) call fatal_error('OPTWF_MIX: no need mix_n')

      call set_nparms_tot

      if(nparm.gt.MPARM)call fatal_error('SR_OPTWF: nparmtot gt MPARM')

      write(6,'(''Starting dparm_norm_min'',g12.4)') dparm_norm_min
      write(6,'(/,''SR adiag: '',f10.5)') sr_adiag
      write(6,'(''SR tau:   '',f10.5)') sr_tau
      write(6,'(''SR eps:   '',f10.5)') sr_eps

      ! if(nvecx.gt.MVEC) call fatal_error('SR_OPTWF: nvecx > MVEC')

      write(6,'(/,''LIN_D adiag: '',f10.5)') alin_adiag
      write(6,'(''LIN_D ethr:  '',f10.5)') alin_eps
      write(6,'(''LIN_D nvec:  '',i4)') nvec
      write(6,'(''LIN_D nvecx: '',i4)') nvecx

      if(nstates.gt.1.and.nvec.lt.nstates) call fatal_error('SR_OPTWF: nvec < nstates')

      inc_nblk=0

      sr_adiag_sav=sr_adiag
      alin_adiag_sav=alin_adiag
      nblk_sav=nblk

      nstates_sav=nstates
      iforce_analy_sav=iforce_analy

      ioptjas_sav=ioptjas
      ioptorb_sav=ioptorb
      ioptci_sav=ioptci

      call save_nparms
      call write_geometry(0)
      call save_wf
      call save_ci_best

c do iteration
      do iter=1,nopt_iter
        write(6,'(/,''Optimization iteration'',i5,'' of'',i5)')iter,nopt_iter

        iforce_analy=0

c do micro_iteration
        do miter=1,micro_iter_sr

          if(micro_iter_sr.gt.1) write(6,'(/,''Micro iteration'',i5,'' of'',i5)')miter,micro_iter_sr

          method='sr_n'
          ioptci=0
          ioptorb=ioptorb_sav
          ioptjas=ioptjas_sav
          nblk_sav=nblk
          if(miter.eq.micro_iter_sr) then
            method='lin_d'
            ioptci=ioptci_sav

            nblk=nblk_ci
            write(6,'(''NBLOCK changed from '',i7, '' to '',i7)') nblk_sav,nblk

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
          nblk=nblk_sav
          write(6,'(/,''Completed sampling'')')

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
                  write(6,*) istate0,i_deltap(istate)

                  write(6,'(''STATE, N OVERLAP, ENERGY OLD, ENERGY NEW,10*ERRDIFF '',2i3,3f12.5)') 
     &            istate,index_more(i_deltap(istate),istate),energy_old(istate),energy_all(istate),10*errdiff
                  iqmc_again=1
                endif

              enddo   
              if(iqmc_again.gt.0) then
                iqmc_check=iqmc_check+1
                write(6,'(''Use new set of CI coefficients'',i4)')
                go to 5
               else
                call save_ci_best
                write(6,'(''Save current CI coefficients as best'')')
              endif
             else
              call restore_ci_best
              write(6,'(''Restore CI cofficients to previous iteration'')')

              call qmc
            endif   
          endif

   6      continue

          if(method.eq.'sr_n') then
            call sr(nparm,deltap,sr_adiag,sr_eps,i)
            call dscal(nparm,-sr_tau,deltap,1)
            adiag=sr_adiag
           else
            call lin_d(nparm,nvec,nvecx,deltap,deltap_more,index_more,alin_adiag,alin_eps)
            if(nstates.eq.1) call dscal(nparm,-1.d0,deltap,1)
            adiag=alin_adiag
          endif
          call test_solution_parm(nparm,deltap,dparm_norm,dparm_norm_min,adiag,iflag)
          write(6,'(''Norm of parm variation '',g12.5)') dparm_norm
          if(iflag.ne.0) then
            write(6,'(''Warning: dparm_norm>1'')')
            adiag=10*adiag
            write(6,'(''adiag increased to '',f10.5)') adiag

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
c         call check_length_run_sr(iter,inc_nblk,nblk,nblk_max,denergy,denergy_err,energy_err_sav,energy_tol)
          nblk=nblk*1.2
          nblk=min(nblk,nblk_max)
c         if(-denergy.gt.3*denergy_err) alfgeo=alfgeo/1.2
        endif
        write(6,'(''nblk = '',i6)') nblk
        write(6,'(''alfgeo = '',f10.4)') alfgeo

        energy_sav=energy(1)
        energy_err_sav=energy_err(1)
        sigma_sav=sigma
      enddo
c enddo iteration

      write(6,'(/,''Check last iteration'')')

      ioptjas=0
      ioptorb=0
      ioptci=0
      iforce_analy=0

      call set_nparms

      call qmc

      call write_wf(1,-1)
      call write_geometry(-1)

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine change_ci(dparm_new,istate)
      use csfs, only: ccsf, cxdet, iadet, ibdet, icxdet, ncsf

      use dets, only: cdet, ndet
      implicit real*8(a-h,o-z)




      dimension dparm_new(*)

c      write(6,*) "COPUTING NEW CI, ccsf(1,state,1)", ccsf(1,istate,1), dparm_new(1)
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
c      write(6,*) "COPUTING NEW CI, ccsf(1,istate,1)", ccsf(1,istate,1), dparm_new(1)

      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
