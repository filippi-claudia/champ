      subroutine optwf_sr_nogeo

      implicit real*8 (a-h,o-z)
      include 'vmc.h'
      include 'force.h'
      include 'mstates.h'
      include 'sr.h'

      common /const/ pi,hb,etrial,delta,deltai,fbias,nelec,imetro,ipr
      common /contrl/ nstep,nblk,nblkeq,nconf_old,nconf_new,isite,idump,irstar
      common /optwf_corsam/ add_diag(MFORCE),energy(MFORCE),energy_err(MFORCE),force(MFORCE),force_err(MFORCE),sigma
      common /optwf_contrl/ ioptjas,ioptorb,ioptci,nparm
      common /optwf_func/ omega,omega_hes,ifunc_omega

      common /force_analy/ iforce_analy,iuse_zmat,alfgeo
      common /atom/ znuc(MCTYPE),cent(3,MCENT),pecent
     &,iwctype(MCENT),nctype,ncent
      common /csfs/ ccsf(MDET,MSTATES,MWF),cxdet(MDET*MDETCSFX)
     &,icxdet(MDET*MDETCSFX),iadet(MDET),ibdet(MDET),ncsf,nstates

      dimension grad(MPARM*MSTATES)

      if(method.ne.'sr_n') return

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

      call p2gtid('optgeo:izvzb',izvzb,0,1)

      write(6,'(/,''SR adiag: '',f10.5)') sr_adiag
      write(6,'(''SR tau:   '',f10.5)') sr_tau
      write(6,'(''SR eps:   '',f10.5)') sr_eps

      inc_nblk=0
      
      sr_adiag_sav=sr_adiag


      ioptjas_sav=ioptjas
      ioptorb_sav=ioptorb
      ioptci_sav=ioptci
      call save_nparms

      call write_geometry(0)

c do iteration
      do iter=1,nopt_iter
        write(6,'(/,''Optimization iteration'',i5,'' of'',i5)')iter,nopt_iter

        iforce_analy=0

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
 
          call write_wf(1,iter)

          call save_wf

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

      return
      end

