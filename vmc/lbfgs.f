      subroutine lbfgs_optwf

      implicit real*8 (a-h,o-z)
      character*20 dl_alg

      include 'vmc.h'
      include 'force.h'
      include 'mstates.h'
      include 'sr.h'

      common /const/ pi,hb,etrial,delta,deltai,fbias,nelec,imetro,ipr
      common /contrl/ nstep,nblk,nblkeq,nconf_old,nconf_new,isite,idump,irstar
      common /optwf_corsam/ add_diag(MFORCE),energy(MFORCE),energy_err(MFORCE),force(MFORCE),force_err(MFORCE)
      common /optwf_contrl/ ioptjas,ioptorb,ioptci,nparm

      common /force_analy/ iforce_analy,iuse_zmat,alfgeo
      common /atom/ znuc(MCTYPE),cent(3,MCENT),pecent
     &,iwctype(MCENT),nctype,ncent
      common /csfs/ ccsf(MDET,MSTATES,MWF),cxdet(MDET*MDETCSFX)
     &,icxdet(MDET*MDETCSFX),iadet(MDET),ibdet(MDET),ncsf,nstates

      common /sr_mat_n/ sr_o(MPARM,MCONF),sr_ho(MPARM,MCONF),obs(MOBS,MSTATES),s_diag(MPARM,MSTATES)
     &,s_ii_inv(MPARM),h_sr(MPARM),wtg(MCONF,MSTATES),elocal(MCONF,MSTATES),jfj,jefj,jhfj,nconf

c vector of wave function parameters
      dimension deltap(MPARM), parameters(MPARM), diag(MPARM), workspace(MPARM*11 + 10)

      call p2gtid('optwf:idl_flag',idl_flag,0,1)

      if(method.ne.'sr_n')return

      write(6,'(''Started lbfgs optimization'')')

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
      call p2gtid('optwf:idl_flag',idl_flag,0,1)

      write(6,'(/,''SR adiag: '',f10.5)') sr_adiag
      write(6,'(''SR tau:   '',f10.5)') sr_tau
      write(6,'(''SR eps:   '',f10.5)') sr_eps
      write(6,'(''DL flag:   '',I10)') idl_flag

      inc_nblk=0
c Initialize DL vectors to zero
      do i=1,nparm
        parameters(i) = 0.d0
        diag(i) = 0.d0
      enddo

      call save_nparms

      call fetch_parameters(parameters)

c do iteration
      do iter=1,nopt_iter
        write(6,'(/,''LBFGS Optimization iteration'',i5,'' of'',i5)')iter,nopt_iter

        call qmc

        write(6,'(/,''Completed sampling'')')

   6    continue

        call lbfgs_more(iter, nparm, deltap, parameters, energy(1), diag, workspace)

c historically, we input -deltap in compute_parameters, so we multiply actual deltap by -1
        call dscal(nparm,-1.d0,deltap,1)

        call test_solution_parm(nparm,deltap,dparm_norm,dparm_norm_min,sr_adiag,iflag)
        write(6,'(''Norm of parm variation '',g12.5)') dparm_norm
        if(iflag.ne.0) then
          write(6,'(''Warning: dparm_norm>1'')')
          stop
c         sr_adiag=10*sr_adiag
c         write(6,'(''sr_adiag increased to '',f10.5)') sr_adiag
c         go to 6
        endif

        call compute_parameters(deltap,iflag,1)
        call write_wf(1,iter)
 
        call save_wf

        if(iter.ge.2) then
          denergy=energy(1)-energy_sav
          denergy_err=sqrt(energy_err(1)**2+energy_err_sav**2)
c         call check_length_run_sr(iter,inc_nblk,nblk,nblk_max,denergy,denergy_err,energy_err_sav,energy_tol)
          nblk=nblk*1.2
          nblk=min(nblk,nblk_max)
        endif
        write(6,'(''nblk = '',i6)') nblk

        energy_sav=energy(1)
        energy_err_sav=energy_err(1)
      enddo
c enddo iteration

      write(6,'(/,''Check last iteration'')')

      ioptjas=0
      ioptorb=0
      ioptci=0

      call set_nparms_tot

      call qmc

      call write_wf(1,-1)

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
