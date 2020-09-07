      subroutine optwf_sr
      use vmc, only: MELEC, MORB, MBASIS, MDET, MCENT, MCTYPE, MCTYP3X
      use vmc, only: NSPLIN, nrad, MORDJ, MORDJ1, MMAT_DIM, MMAT_DIM2, MMAT_DIM20
      use vmc, only: radmax, delri
      use vmc, only: NEQSX, MTERMS
      use vmc, only: MCENT3, NCOEF, MEXCIT
      use optwf_contrl, only: ioptci, ioptjas, ioptorb, nparm
      use mstates_mod, only: MSTATES
      use optwf_corsam, only: energy, energy_err, force
      use optwf_func, only: ifunc_omega, omega, omega_hes
      use contrl, only: nblk
      use force_analy, only: iforce_analy, alfgeo

      use method_opt, only: method

      implicit real*8(a-h,o-z)


      include 'force.h'
      include 'sr.h'

      dimension grad(MPARM*MSTATES)

      if(method.ne.'sr_n') return

      call set_nparms_tot

      if(nparm.gt.MPARM)call fatal_error('SR_OPTWF: nparmtot gt MPARM')

      call p2gtid('optwf:nopt_iter',nopt_iter,6,1)
      call p2gtid('optwf:nblk_max',nblk_max,nblk,1)
      call p2gtfd('optwf:energy_tol',energy_tol,1.d-3,1)

      call p2gtfd('optwf:dparm_norm_min',dparm_norm_min,1.0d0,1)
      write(6,'(''Starting dparm_norm_min'',g12.4)') dparm_norm_min

      call p2gtfd('optwf:sr_tau',sr_tau,0.02,1)
      call p2gtfd('optwf:sr_adiag',sr_adiag,0.01,1)
      call p2gtfd('optwf:sr_eps',sr_eps,0.001,1)
      call p2gtid('optwf:sr_rescale',i_sr_rescale,0,1)

      call p2gtid('optwf:func_omega',ifunc_omega,0,1)
      if(ifunc_omega.gt.0) then
       call p2gtfd('optwf:omega',omega0,0.d0,1)
       call p2gtid('optwf:n_omegaf',n_omegaf,nopt_iter,1)
       call p2gtid('optwf:n_omegat',n_omegat,0,1)
       if(n_omegaf+n_omegat.gt.nopt_iter) call fatal_error('SR_OPTWF: n_omegaf+n_omegat > nopt_iter')
       omega=omega0
       write(6,'(/,''SR ifunc_omega: '',i3)') ifunc_omega
       write(6,'(''SR omega: '',f10.5)') omega
       write(6,'(''SR n_omegaf: '',i4)') n_omegaf
       write(6,'(''SR n_omegat: '',i4)') n_omegat
      endif

      call p2gtid('optwf:micro_iter_sr',micro_iter_sr,1,1)

      call p2gtid('optgeo:izvzb',izvzb,0,1)

      write(6,'(/,''SR adiag: '',f10.5)') sr_adiag
      write(6,'(''SR tau:   '',f10.5)') sr_tau
      write(6,'(''SR eps:   '',f10.5)') sr_eps

      inc_nblk=0
      
      sr_adiag_sav=sr_adiag

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
          omega_hes=energy_sav
          if(iter.gt.n_omegaf) then
            alpha_omega=dfloat(n_omegaf+n_omegat-iter)/n_omegat
            omega=alpha_omega*omega0+(1.d0-alpha_omega)*(energy_sav-sigma_sav)
            if(ifunc_omega.eq.1.or.ifunc_omega.eq.2) omega=alpha_omega*omega0+(1.d0-alpha_omega)*energy_sav
          endif
          if(iter.gt.n_omegaf+n_omegat) then
            omega=energy_sav-sigma_sav
            if(ifunc_omega.eq.1.or.ifunc_omega.eq.2) omega=energy_sav
          endif
          write(6,'(''SR omega: '',f10.5)') omega
        endif

c do micro_iteration
        do miter=1,micro_iter_sr

          if(micro_iter_sr.gt.1) write(6,'(/,''Micro iteration'',i5,'' of'',i5)')miter,micro_iter_sr

          if(miter.eq.micro_iter_sr) iforce_analy=iforce_analy_sav

          call qmc

          write(6,'(/,''Completed sampling'')')

   6      continue

          call sr(nparm,grad,sr_adiag,sr_eps,i)
          call dscal(nparm,-sr_tau,grad,1)

          adiag=sr_adiag
          call test_solution_parm(nparm,grad,dparm_norm,dparm_norm_min,adiag,iflag)
          write(6,'(''Norm of parm variation '',d12.5)') dparm_norm
          if(iflag.ne.0) then
            write(6,'(''Warning: dparm_norm>1'')')
            adiag=10*adiag
            write(6,'(''adiag increased to '',f10.5)') adiag

            sr_adiag=adiag
            go to 6
           else
            sr_adiag=sr_adiag_sav
          endif
 
          call compute_parameters(grad,iflag,1)
          call write_wf(1,iter)

          call save_wf

          if(iforce_analy.gt.0) then

            if(izvzb.gt.0) call forces_zvzb(nparm)

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

      subroutine sr(nparm,deltap,sr_adiag,sr_eps,i)
c solve S*deltap=h_sr (call in optwf)

      use sr_mat_n, only: h_sr
      implicit real*8(a-h,o-z)

      include 'sr.h'


      dimension deltap(nparm)

      call sr_hs(nparm,sr_adiag)

      imax=nparm          ! max n. iterations conjugate gradients
      imod=50             ! inv. freq. of calc. r=b-Ax vs. r=r-\alpha q (see pcg)
      do i=1,nparm   
       deltap(i)=0.d0     ! initial guess of solution
      enddo
      call pcg(nparm,h_sr,deltap,i,imax,imod,sr_eps)
      write(6,*) 'CG iter ',i

      call sr_rescale_deltap(nparm,deltap)

      return              ! deltap
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine check_length_run_sr(iter,increase_nblk,nblk,nblk_max,denergy,denergy_err,energy_err_sav,energy_tol)

      implicit real*8(a-h,o-z)

c Increase nblk if near convergence to value needed to get desired statistical error
      increase_nblk=increase_nblk+1

c Increase if subsequent energies are within errorbar
      if(dabs(denergy).lt.3*denergy_err.and.energy_tol.gt.0.d0) then
        nblk_new=nblk*max(1.d0,(energy_err_sav/energy_tol)**2)
        nblk_new=min(nblk_new,nblk_max)
        if(nblk_new.gt.nblk) then
          increase_nblk=0
          nblk=nblk_new
          write(6,'(''nblk reset to'',i8,9d12.4)') nblk,dabs(denergy),energy_tol
        endif
      endif
c Always increase nblk by a factor of 2 every other iteration
      if(increase_nblk.eq.2.and.nblk.lt.nblk_max) then
        increase_nblk=0
        nbkl=1.2*nblk
        nblk=min(nblk,nblk_max)
        write(6,'(''nblk reset to'',i8,9d12.4)') nblk
      endif

      return
      end
