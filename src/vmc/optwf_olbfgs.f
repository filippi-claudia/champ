      subroutine optwf_olbfgs
      use olbfgs, only: initialize_olbfgs
      use optwf_contrl, only: ioptci, ioptjas, ioptorb, nparm
      use optwf_corsam, only: energy, energy_err
      use control_vmc, only: vmc_nblk, vmc_nblk_max
      use optwf_contrl, only: idl_flag, ilbfgs_flag
      use optwf_contrl, only: sr_tau, sr_adiag, sr_eps
      use optwf_contrl, only: dparm_norm_min, nopt_iter
      use method_opt, only: method
      use contrl_file,    only: ounit
      use precision_kinds, only: dp
      implicit none

      integer :: i, iflag, ilbfgs_m, inc_nblk, iter
      real(dp) :: denergy, denergy_err, dparm_norm, energy_err_sav, energy_sav
      real(dp), dimension(nparm) :: deltap
      real(dp), dimension(nparm) :: parameters

      character*20 dl_alg

c vector of wave function parameters


      if(method.ne.'sr_n'.or.ilbfgs_flag.eq.0)return

      write(ounit,'(''Started oLBFGS optimization'')')

      call set_nparms_tot


      write(ounit,'(''Starting dparm_norm_min'',g12.4)') dparm_norm_min

      write(ounit,'(/,''SR adiag: '',f10.5)') sr_adiag
      write(ounit,'(''SR tau:   '',f10.5)') sr_tau
      write(ounit,'(''SR eps:   '',f10.5)') sr_eps
      write(ounit,'(''DL flag:   '',I10)') idl_flag
      write(ounit,'(''LBFGS flag:   '',I10)') ilbfgs_flag

      inc_nblk=0
c Initialize DL vectors to zero
      do i=1,nparm
      enddo

      call save_nparms

      call fetch_parameters(parameters)

      ! initialize olbfgs
      call initialize_olbfgs(nparm, ilbfgs_m)

c do iteration
      do iter=1,nopt_iter
        write(ounit,'(/,''OLBFGS Optimization iteration'',i5,'' of'',i5)')iter,nopt_iter

        call qmc

        write(ounit,'(/,''Completed sampling'')')

   6    continue

        call olbfgs_more(iter, nparm, deltap, parameters)

c historically, we input -deltap in compute_parameters, so we multiply actual deltap by -1
        call dscal(nparm,-1.d0,deltap,1)

        call test_solution_parm(nparm,deltap,dparm_norm,dparm_norm_min,sr_adiag,iflag)
        write(ounit,'(''Norm of parm variation '',g12.5)') dparm_norm
        if(iflag.ne.0) then
          write(ounit,'(''Warning: dparm_norm>1'')')
          stop
c         sr_adiag=10*sr_adiag
c         write(ounit,'(''sr_adiag increased to '',f10.5)') sr_adiag
c         go to 6
        endif

        call compute_parameters(deltap,iflag,1)
        call write_wf(1,iter)

        call save_wf

        if(iter.ge.2) then
          denergy=energy(1)-energy_sav
          denergy_err=sqrt(energy_err(1)**2+energy_err_sav**2)
c         call check_length_run_sr(iter,inc_nblk,nblk,nblk_max,denergy,denergy_err,energy_err_sav,energy_tol)
          vmc_nblk=vmc_nblk*1.2
          vmc_nblk=min(vmc_nblk,vmc_nblk_max)
        endif
        write(ounit,'(''nblk = '',i6)') vmc_nblk

        energy_sav=energy(1)
        energy_err_sav=energy_err(1)
      enddo
c enddo iteration

      write(ounit,'(/,''Check last iteration'')')

      ioptjas=0
      ioptorb=0
      ioptci=0

      call set_nparms_tot

      call qmc

      call write_wf(1,-1)

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
