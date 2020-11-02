      subroutine optwf_mix_nogeo

      implicit real*8 (a-h,o-z)
      include 'vmc.h'
      include 'force.h'
      include 'mstates.h'
      include 'sr.h'

      character*20 method_sav

      common /const/ pi,hb,etrial,delta,deltai,fbias,nelec,imetro,ipr
      common /contrl/ nstep,nblk,nblkeq,nconf_old,nconf_new,isite,idump,irstar
      common /optwf_corsam/ add_diag(MFORCE),energy(MFORCE),energy_err(MFORCE),force(MFORCE),force_err(MFORCE),sigma
      common /optwf_contrl/ ioptjas,ioptorb,ioptci,nparm
      common /optwf_func/ omega,omega_hes,ifunc_omega
      common /sa_check/ energy_all(MSTATES),energy_err_all(MSTATES) 
      common /force_analy/ iforce_analy,iuse_zmat,alfgeo
      common /atom/ znuc(MCTYPE),cent(3,MCENT),pecent
     &,iwctype(MCENT),nctype,ncent
      common /csfs/ ccsf(MDET,MSTATES,MWF),cxdet(MDET*MDETCSFX)
     &,icxdet(MDET*MDETCSFX),iadet(MDET),ibdet(MDET),ncsf,nstates
      dimension grad(MPARM*MSTATES),grad_more(MPARM*MSTATES,5),index_more(5,MSTATES)
      dimension deltap(MPARM*MSTATES),deltap_more(MPARM*MSTATES,5)
      dimension energy_old(MSTATES), energy_err_old(MSTATES), i_deltap(MSTATES), energy_davidson(6,MSTATES) 
      dimension index_min_energy(5), deltap_new(MPARM)

      if(nstates.eq.1.or.ioptci.eq.0) call fatal_error('OPTWF_MIX: no need mix_n')

      call set_nparms_tot

      if(nparm.gt.MPARM)call fatal_error('SR_OPTWF: nparmtot gt MPARM')

      call p2gtid('optwf:nopt_iter',nopt_iter,6,1)
      if(nopt_iter==0) return

      call p2gtid('optwf:nblk_max',nblk_max,nblk,1)
      call p2gtfd('optwf:energy_tol',energy_tol,1.d-3,1)

      call p2gtfd('optwf:dparm_norm_min',dparm_norm_min,1.0d0,1)
      write(6,'(''Starting dparm_norm_min'',g12.4)') dparm_norm_min

      call p2gtfd('optwf:sr_tau',sr_tau,0.02,1)
      call p2gtfd('optwf:sr_adiag',sr_adiag,0.01,1)
      call p2gtfd('optwf:sr_eps',sr_eps,0.001,1)

      write(6,'(/,''SR adiag: '',f10.5)') sr_adiag
      write(6,'(''SR tau:   '',f10.5)') sr_tau
      write(6,'(''SR eps:   '',f10.5)') sr_eps

      call p2gtid('optwf:lin_nvec',nvec,5,1)
      call p2gtid('optwf:lin_nvecx',nvecx,MVEC,1)
      call p2gtfd('optwf:lin_adiag',alin_adiag,0.01,1)
      call p2gtfd('optwf:lin_eps',alin_eps,0.001,1)

      if(nvecx.gt.MVEC) call fatal_error('SR_OPTWF: nvecx > MVEC')
      write(6,'(/,''LIN_D adiag: '',f10.5)') alin_adiag
      write(6,'(''LIN_D ethr:  '',f10.5)') alin_eps
      write(6,'(''LIN_D nvec:  '',i4)') nvec
      write(6,'(''LIN_D nvecx: '',i4)') nvecx

      if(nstates.gt.1.and.nvec.lt.nstates) call fatal_error('SR_OPTWF: nvec < nstates')

      call p2gtid('optwf:micro_iter_sr',micro_iter_sr,1,1)

      if(iforce_analy.gt.0) call p2gtid('optgeo:iroot_geo',iroot_geo,0,0)

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
            call p2gtid('optwf:nblk_ci',nblk_ci,nblk,0)
            nblk=nblk_ci
            write(6,*) "NBLOCK CI", "NBLOCK SAV", nblk, nblk_sav
            ioptorb=0
            ioptjas=0
          endif
          call set_nparms
          if(iter.eq.1.and. miter.eq.micro_iter_sr) nparmci=nparm


   5      call qmc 
          nblk=nblk_sav
          write(6,'(/,''Completed sampling'')')

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

      enddo
c enddo iteration

      return
      end

