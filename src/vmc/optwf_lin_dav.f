      subroutine optwf_lin_d


      use sr_mod, only: MPARM, MOBS, MCONF, MVEC, nvecx
      use force_mod, only: MFORCE, MFORCE_WT_PRD, MWF
      use vmc_mod, only: MELEC, MORB, MBASIS, MDET, MCENT, MCTYPE, MCTYP3X
      use vmc_mod, only: NSPLIN, nrad, MORDJ, MORDJ1, MMAT_DIM, MMAT_DIM2, MMAT_DIM20
      use vmc_mod, only: radmax, delri
      use vmc_mod, only: NEQSX, MTERMS
      use vmc_mod, only: MCENT3, NCOEF, MEXCIT
      use csfs, only: nstates
      use mstates_mod, only: MSTATES
      use optwf_contrl, only: ioptci, ioptjas, ioptorb, nparm
      use optwf_corsam, only: energy, energy_err, force
      use optwf_func, only: ifunc_omega, omega
      use contrl, only: nblk
      use force_analy, only: iforce_analy, alfgeo

      use method_opt, only: method

      implicit real*8(a-h,o-z)



      character*20 method_sav

      dimension grad(MPARM*nstates),grad_more(MPARM*nstates,5),index_more(5,nstates)

      if(method .ne.'lin_d')return

      call set_nparms_tot

      if(nparm.gt.MPARM)call fatal_error('OPTWF_LIN_D: nparmtot gt MPARM')

      call p2gtid('optwf:nopt_iter',nopt_iter,6,1)
      call p2gtid('optwf:nblk_max',nblk_max,nblk,1)
      call p2gtfd('optwf:energy_tol',energy_tol,1.d-3,1)

      call p2gtfd('optwf:dparm_norm_min',dparm_norm_min,1.0d0,1)
      write(6,'(''Starting dparm_norm_min'',g12.4)') dparm_norm_min

      call p2gtid('optwf:lin_nvec',nvec,5,1)
      call p2gtid('optwf:lin_nvecx',nvecx,MVEC,1)
      call p2gtfd('optwf:lin_adiag',alin_adiag,0.01,1)
      call p2gtfd('optwf:lin_eps',alin_eps,0.001,1)


      call p2gtid('optwf:func_omega',ifunc_omega,0,1)
      if(ifunc_omega.gt.0) then
       call p2gtfd('optwf:omega',omega0,0.d0,1)
       call p2gtid('optwf:n_omegaf',n_omegaf,nopt_iter,1)
       call p2gtid('optwf:n_omegat',n_omegat,0,1)
       if(n_omegaf+n_omegat.gt.nopt_iter) call fatal_error('OPTWF_LIN_D: n_omegaf+n_omegat > nopt_iter')
       omega=omega0
       write(6,'(/,''LIN_D ifunc_omega: '',i3)') ifunc_omega
       write(6,'(''LIN_D omega: '',f10.5)') omega
       write(6,'(''LIN_D n_omegaf: '',i4)') n_omegaf
       write(6,'(''LIN_D n_omegat: '',i4)') n_omegat
      endif

      call p2gtid('optwf:micro_iter_sr',micro_iter_sr,1,1)

      ! if(nvecx.gt.MVEC) call fatal_error('SR_OPTWF: nvecx > MVEC')
      write(6,'(/,''LIN_D adiag: '',f10.5)') alin_adiag
      write(6,'(''LIN_D ethr:  '',f10.5)') alin_eps
      write(6,'(''LIN_D nvec:  '',i4)') nvec
      write(6,'(''LIN_D nvecx: '',i4)') nvecx

      if(nstates.gt.1.and.nvec.lt.nstates) call fatal_error('OPTWF_LIN_D: nvec < nstates')

      inc_nblk=0
      
      alin_adiag_sav=alin_adiag

      nstates_sav=nstates
      iforce_analy_sav=iforce_analy

      ioptjas_sav=ioptjas
      ioptorb_sav=ioptorb
      ioptci_sav=ioptci
      call save_nparms

      call write_geometry(0)

c do iteration
      do iter=1,nopt_iter
        write(6,'(/,''Optimization iteration'',i5,'' of'',i5)')iter,nopt_iter

        iforce_analy=0

        if(ifunc_omega.gt.0) then
          if(iter.gt.n_omegaf) then
            alpha_omega=dfloat(n_omegaf+n_omegat-iter)/n_omegat
            omega=alpha_omega*omega0+(1.d0-alpha_omega)*(energy_sav-sigma_sav)
            if(ifunc_omega.eq.2) omega=alpha_omega*omega0+(1.d0-alpha_omega)*energy_sav
          endif
          if(iter.gt.n_omegaf+n_omegat) then
            omega=energy_sav-sigma_sav
            if(ifunc_omega.eq.2) omega=energy_sav
          endif
          write(6,'(''LIN_D omega: '',f10.5)') omega
        endif

c do micro_iteration
        do miter=1,micro_iter_sr

          if(micro_iter_sr.gt.1) write(6,'(/,''Micro iteration'',i5,'' of'',i5)')miter,micro_iter_sr

          if(miter.eq.micro_iter_sr) iforce_analy=iforce_analy_sav

c        efin_old = efin define efin_old as the energy before

          call qmc

          write(6,'(/,''Completed sampling'')')

   6      continue

          call lin_d(nparm,nvec,nvecx,grad,grad_more,index_more,alin_adiag,alin_eps)
          if(nstates.eq.1) call dscal(nparm,-1.d0,grad,1)

          if(method.eq.'lin_d'.and.ioptorb+ioptjas.gt.0) then
            adiag=alin_adiag
            call test_solution_parm(nparm,grad,dparm_norm,dparm_norm_min,adiag,iflag)
            write(6,'(''Norm of parm variation '',g12.5)') dparm_norm
            if(iflag.ne.0) then
              write(6,'(''Warning: dparm_norm>1'')')
              adiag=10*adiag
              write(6,'(''adiag increased to '',f10.5)') adiag

              alin_adiag=adiag
              go to 6
             else
              alin_adiag=alin_adiag_sav
            endif
          endif


c Here I should save the old parameters 
 
          call compute_parameters(grad,iflag,1)
          call write_wf(1,iter)

          call save_wf

          if(iforce_analy.gt.0) then
            call compute_positions
            call write_geometry(iter)
          endif
        enddo
c enddo micro_iteration
 
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
      subroutine h_psi_lin_d(ndim,nvec,psi,hpsi )
      use sr_mod, only: MPARM, MOBS, MCONF, MVEC
      use optwf_func, only: ifunc_omega
      implicit real*8(a-h,o-z)




      dimension psi(MPARM,*),hpsi(MPARM,*)

      if(ifunc_omega.eq.0) then
        call h_psi_energymin(ndim,nvec,psi,hpsi )
       elseif(ifunc_omega.le.2) then
        call h_psi_varmin(ndim,nvec,psi,hpsi )
       else
        call h_psi_omegamin(ndim,nvec,psi,hpsi )
      endif
      
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine s_psi_lin_d(ndim,nvec,psi,spsi )
      use sr_mod, only: MPARM, MOBS, MCONF, MVEC
      use optwf_func, only: ifunc_omega
      implicit real*8(a-h,o-z)


 

      dimension psi(MPARM,*),spsi(MPARM,*)

      if(ifunc_omega.le.2) then
        call s_psi_energymin(ndim,nvec,psi,spsi )
       else
        call s_psi_omegamin(ndim,nvec,psi,spsi )
      endif
      
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine select_ci_root(iroot)
      use force_mod, only: MFORCE, MFORCE_WT_PRD, MWF
      use vmc_mod, only: MELEC, MORB, MBASIS, MDET, MCENT, MCTYPE, MCTYP3X
      use vmc_mod, only: NSPLIN, nrad, MORDJ, MORDJ1, MMAT_DIM, MMAT_DIM2, MMAT_DIM20
      use vmc_mod, only: radmax, delri
      use vmc_mod, only: NEQSX, MTERMS
      use vmc_mod, only: MCENT3, NCOEF, MEXCIT
      use csfs, only: ccsf, ncsf

      use dets, only: cdet, ndet
      implicit real*8(a-h,o-z)




      do 30 i=1,ndet
   30   cdet(i,1,1)=cdet(i,iroot,1)

      do 40 icsf=1,ncsf
   40   ccsf(icsf,1,iadiag)=ccsf(icsf,iroot,1)

      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
