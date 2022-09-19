!------------------------------------------------------------------------------
!        Optimization routine using stochastic reconfiguration
!------------------------------------------------------------------------------
!> @author
!> Saverio Moroni, Claudia Filippi
!
! DESCRIPTION:
!> Opitmize the wave function parameters using the SR method
!
! URL           : https://github.com/filippi-claudia/champ
!---------------------------------------------------------------------------

module optwf_sr_mod

    use precision_kinds, only: dp
    use optwf_contrl, only: ioptci, ioptjas, ioptorb
    use force_analy, only: iforce_analy
!    use contrl, only: nblk_max
    use control_vmc, only: vmc_nblk_max
    use optwf_contrl, only: energy_tol, nopt_iter, micro_iter_sr, dparm_norm_min
    use optwf_contrl, only: sr_tau , sr_adiag, sr_eps
    use orbval, only: nadorb
    use contrl_file,    only: ounit
    use mpitimer,    only: elapsed_time
    use error, only: fatal_error
    use optwf_handle_wf, only : set_nparms_tot, save_nparms, test_solution_parm
    use optwf_handle_wf, only : compute_parameters, write_wf, save_wf
    use optwf_handle_wf, only : set_nparms
    use optgeo_lib, only: write_geometry, compute_positions
    use sr_mod, only: izvzb, i_sr_rescale
    use fssd

    real(dp) :: omega0
    integer :: n_omegaf, n_omegat
    integer ::  i_func_omega
    real(dp) :: sr_adiag_sav

    integer :: ioptjas_sav, ioptorb_sav, ioptci_sav, iforce_analy_sav
    real(dp), dimension(:), allocatable :: deltap

    private
    public :: optwf_sr, sr, sr_hs
    save
    interface !interface to LAPACK
! -- LAPACK routine (version 3.1) --
!    Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!    November 2006
      SUBROUTINE DGETRF( M, N, A, LDA, IPIV, INFO )
        INTEGER            INFO, LDA, M, N
        INTEGER            IPIV( * )
        DOUBLE PRECISION   A( LDA, * )
      END SUBROUTINE
      SUBROUTINE DGETRI( N, A, LDA, IPIV, WORK, LWORK, INFO )
        INTEGER            INFO, LDA, LWORK, N
        INTEGER            IPIV( * )
        DOUBLE PRECISION   A( LDA, * ), WORK( * )
      END SUBROUTINE
        SUBROUTINE dscal(N,DA,DX,INCX)
! -- Reference BLAS level1 routine --
          DOUBLE PRECISION DA
          INTEGER INCX,N
          DOUBLE PRECISION DX(*)
        END SUBROUTINE
    end interface

    interface ! Let linking decide between dmc/vmc
      subroutine qmc
      end subroutine
    end interface
contains

    subroutine optwf_sr

        use atom, only: cent, ncent
        use sr_mod, only: mparm
        use optwf_contrl, only: ioptci, ioptjas, ioptorb, nparm
        use mstates_mod, only: MSTATES
        use optwf_corsam, only: energy, energy_err
        use optwf_func, only: ifunc_omega, omega0, n_omegaf, n_omegat, omega_hes
        !use contrl, only: nblk
        use control_vmc, only: vmc_nblk
        use force_analy, only: alfgeo
        use optwf_contrl, only: nparm
        use method_opt, only: method
        use orbval, only: nadorb
        use contrl_file,    only: ounit
        use fssd, only: errbar_fssd, div_fssd, norm_fssd

        implicit none

        real(dp) :: adiag, denergy, alpha_omega, denergy_err, dparm_norm
        real(dp) :: energy_sav, energy_err_sav, omega, sigma, sigma_sav, nadorb_sav
        real(dp) :: error_bar_fssd
        real(dp), dimension(:, :, :), allocatable :: cent_history
        real(dp), dimension(:, :), allocatable :: cent_ave
        integer :: i, iflag, iter, miter
        integer :: keep_track, stage_track, max_conv, config_conv, max_conf, Na, Nave, Nmin, j, k, ic, n_stage, max_stage
        logical :: check_config_conv

        sigma_sav = 0.0
        energy_sav= 0.0
        allocate (deltap(mparm*MSTATES))

        nadorb_sav = nadorb

        if (method .ne. 'sr_n') return

        call set_nparms_tot

        if (nparm .gt. mparm) call fatal_error('SR_OPTWF: nparmtot gt mparm')

        write (ounit, '(''Starting dparm_norm_min'',g12.4)') dparm_norm_min

        if (ifunc_omega .gt. 0) then
           if (n_omegaf + n_omegat .gt. nopt_iter) call fatal_error('SR_OPTWF: n_omegaf+n_omegat > nopt_iter')
           omega = omega0
           write (ounit, '(/,''SR ifunc_omega: '',i3)') ifunc_omega
           write (ounit, '(''SR omega: '',f10.5)') omega
           write (ounit, '(''SR n_omegaf: '',i4)') n_omegaf
           write (ounit, '(''SR n_omegat: '',i4)') n_omegat
        endif

        write (ounit, '(/,''SR adiag: '',f10.5)') sr_adiag
        write (ounit, '(''SR tau:   '',f10.5)') sr_tau
        write (ounit, '(''SR eps:   '',f10.5)') sr_eps


        call save_params()

        call save_nparms

        call write_geometry(0)

        ! do set stage
        call allocate_fssd()
        d_fssd = 0.0

        keep_track = 0
        error_bar_fssd = errbar_fssd
        !vmc_nblk = 1 / (error_bar_fssd * error_bar_fssd)
        !error_bar_fssd = 1 / sqrt(vmc_nblk)

        max_conv = 10
        max_conf = 10
        Na = 5
        Nave = 10
        Nmin = 2 * Na + Nave
        max_stage = 2

        allocate(cent_history(max_conf * max_conv * nopt_iter, 3, ncent))
        allocate(cent_ave(3, ncent))

        cent_history = 0.0

        check_config_conv = .false.
        config_conv = 1
        n_stage = 1

        if (ifssd.gt.0) then
           if (div_fssd.gt.0.and.norm_fssd.gt.0) then
              write (ounit, '(/,''You are running the SR with Shiwei algorithm in the standard fashion'')')
           elseif (div_fssd.le.0.and.norm_fssd.gt.0) then
              write (ounit, '(/,''You are running the SR with Shiwei algorithm in the modified version where we do not normalize the displacement direction'')')
           else
              write (ounit, '(/,''You are running the SR with Shiwei algorithm in the modified version where we do not normalize the displacement direction and else the displacement direction numerator is not divided by alpha + 1'')')
           endif
        else
           write (ounit, '(/,''You are running the SR with steepest descent algorithm'')')
        endif

        do while (n_stage.le.max_stage)
           write (ounit, '(/,''Starting stage '', i5)') n_stage
           stage_track = 0
           check_config_conv = .false.

           if (iset.gt.0) then
              d_fssd = 0.0
              cent_history = 0.0
              cent_ave = 0.0
           endif

           if (ifssd.gt.0) then
              do while (stage_track.le.Nmin)
                 write (ounit, '(/,''Prepatory steps ... '')') 
                 call macro_iteration(energy_sav, alpha_omega, omega, sigma_sav, adiag, dparm_norm, iflag, cent_history, keep_track, stage_track, error_bar_fssd, config_conv)
              enddo
           endif

           do while (.not. check_config_conv)
              cent_ave = 0.0
              if (ifssd .gt. 0) then 
                 write (ounit, '(/,''Config_conv is '', i5)') config_conv

                 write (ounit, '(/,''Starting macro_iteration with alfgeo, error_bar_fssd, keep_track '', d10.3, '' '', d10.3, i5)') alfgeo, error_bar_fssd, keep_track

              endif

              call macro_iteration(energy_sav, alpha_omega, omega, sigma_sav, adiag, dparm_norm, iflag, cent_history, keep_track, stage_track, error_bar_fssd, config_conv)

              if(ifssd.gt.0) then
                 do j = stage_track - Nave + 1, stage_track
                    write (ounit, '(/,''Configuration at '', i5)') j 
                    do ic = 1,ncent
                       write(ounit,'(3f10.5)') (cent_history(j, k, ic),k=1,3)
                       do k = 1,3
                          cent_ave(k, ic) = cent_ave(k, ic) + cent_history(j, k, ic) / Nave
                       enddo
                    enddo

                    write (ounit, '(/,''Average configuration at '', i5)') j
                    do ic = 1,ncent
                       write(ounit,'(3f10.5)') (cent_ave(k,ic),k=1,3)
                    enddo
                 enddo

                 write (ounit, *) 'Average configuration is: '
                 do ic=1,ncent
                    write(ounit,'(3f10.5)') (cent_ave(k,ic),k=1,3)
                 enddo
              endif

              if (ifssd.gt.0) then ! .and. stage_track.gt.Nmin) then

                 write (ounit, '(/,''Checking convergence ...'')') 
                 call config_check(stage_track, max_conv, cent_history, error_bar_fssd, cent_ave, check_config_conv, Nave)

                 config_conv = config_conv + 1
              else
                 check_config_conv = .true.
              endif
           enddo
           if (ifssd.gt.0) then 
              n_stage = n_stage + 1

              vmc_nblk = vmc_nblk * 4
              alfgeo = alfgeo / 2.d0
           else
              n_stage = 2 * max_stage
           endif
        enddo

        if (ifssd .gt. 0) write (ounit, '(/,''Calculation terminated with alfgeo, error_bar_fssd, keep_track, config_conv '', d10.3, '' '', d10.3, i5, i5)') alfgeo, error_bar_fssd, keep_track, config_conv

        deallocate(cent_ave)
        call deallocate_fssd()

        write (ounit, '(/,''Check last iteration'')')

        ioptjas = 0
        ioptorb = 0
        ioptci = 0
        iforce_analy = 0

        call set_nparms

        call qmc

        nadorb_sav=nadorb
        call write_wf(1, -1)
        call write_geometry(-1)

        deallocate (deltap)
        call elapsed_time( "Last iteration of QMC")
        return
      end subroutine optwf_sr

      subroutine macro_iteration(energy_sav, alpha_omega, omega, sigma_sav, adiag, dparm_norm, iflag, cent_history, keep_track, stage_track, error_bar_fssd, config_conv)
        use optwf_func, only: ifunc_omega, omega0, n_omegaf, n_omegat, omega_hes
        use control_vmc, only: vmc_nblk
        use optwf_corsam, only: energy, energy_err
        use force_analy, only: alfgeo
        use fssd, only: ifssd

        integer :: iter, iflag, keep_track, stage_track, config_conv
        real(dp) :: energy_sav, alpha_omega, omega, sigma_sav, adiag, dparm_norm, error_bar_fssd, denergy, energy_err_sav, denergy_err, sigma
        real(dp), dimension(:, :, :) :: cent_history                

        ! do iteration
        do iter = 1, nopt_iter
           write (ounit, '(/,''Optimization iteration'',i5,'' of'',i5)') iter, nopt_iter

           iforce_analy = 0

           if (ifunc_omega .gt. 0) then
              omega_hes = energy_sav
              if (iter .gt. n_omegaf) then
                 alpha_omega = dfloat(n_omegaf + n_omegat - iter)/n_omegat
                 omega = alpha_omega*omega0 + (1.d0 - alpha_omega)*(energy_sav - sigma_sav)
                 if (ifunc_omega .eq. 1 .or. ifunc_omega .eq. 2) omega = alpha_omega*omega0 + (1.d0 - alpha_omega)*energy_sav
              endif
              if (iter .gt. n_omegaf + n_omegat) then
                 omega = energy_sav - sigma_sav
                 if (ifunc_omega .eq. 1 .or. ifunc_omega .eq. 2) omega = energy_sav
              endif
              write (ounit, '(''SR omega: '',f10.5)') omega
           endif

           call micro_iteration(adiag, dparm_norm, iflag, iter, cent_history, keep_track, stage_track, config_conv)             

           if (iter .ge. 2) then
              denergy = energy(1) - energy_sav
              denergy_err = sqrt(energy_err(1)**2 + energy_err_sav**2)

              if (ifssd.le.0) then
                 vmc_nblk = vmc_nblk*1.2
                 vmc_nblk = min(vmc_nblk, vmc_nblk_max)
              endif

           endif
           write (ounit, '(''nblk = '',i6)') vmc_nblk
           write (ounit, '(''alfgeo = '',f10.4)') alfgeo

           energy_sav = energy(1)
           energy_err_sav = energy_err(1)
           sigma_sav = sigma
           call elapsed_time( "CG iteration ", iter )
        enddo
        ! enddo iteration
    end subroutine macro_iteration

    subroutine micro_iteration(adiag, dparm_norm, iflag, iter, cent_history, keep_track, stage_track, config_conv)

      use optwf_contrl, only: nparm
      use optwf_corsam, only: energy
      use atom, only: cent, ncent

      integer :: i, miter, n_sr, iflag, iter, k, ic, keep_track, stage_track, j, config_conv
      real(dp) :: adiag, dparm_norm
      real(dp), dimension(:, :, :) :: cent_history

      ! do micro_iteration
      do miter = 1, micro_iter_sr

         if (micro_iter_sr .gt. 1) write (ounit, '(/,''Micro iteration'',i5,'' of'',i5)') miter, micro_iter_sr

         if (miter .eq. micro_iter_sr) iforce_analy = iforce_analy_sav

         call qmc

         write (ounit, '(/,''Completed sampling'')')

6        continue

         call sr(nparm, deltap, sr_adiag, sr_eps, i)
         call dscal(nparm, -sr_tau, deltap, 1)

         adiag = sr_adiag
         call test_solution_parm(nparm, deltap, dparm_norm, dparm_norm_min, adiag, iflag)
         write (ounit, '(''Norm of parm variation '',d12.5)') dparm_norm
         if (iflag .ne. 0) then
            write (ounit, '(''Warning: dparm_norm>1'')')
            adiag = 10*adiag
            write (ounit, '(''adiag increased to '',f10.5)') adiag

            sr_adiag = adiag
            go to 6
         else
            sr_adiag = sr_adiag_sav
         endif

         call compute_parameters(deltap, iflag, 1)
         call write_wf(1, iter)

         call save_wf

         keep_track = keep_track + 1
         stage_track = stage_track + 1

         if (iforce_analy .gt. 0) then

            if (izvzb .gt. 0) call forces_zvzb(nparm)

            call compute_positions
            call write_geometry(iter)

            write  (ounit, '(''keep_track is '',i5, '' for config_conv '', i5,'', iter '',i5,'',  miter '',i5)') keep_track, config_conv, iter, miter
            do ic = 1,ncent
               do k = 1,3
                  cent_history(stage_track/micro_iter_sr, k, ic) = cent(k, ic)
               enddo
            enddo

         endif

         iforce_analy = 0

         call elapsed_time("CG micro iteration", miter)
      enddo
      ! enddo micro_iteration
      
    end subroutine micro_iteration

    subroutine config_check(keep_track, max_conv, cent_history, error_bar_fssd, cent_ave, check_config_conv, Nave)

      use atom, only: cent, ncent
      use force_analy, only: alfgeo
      use fssd, only: d_fssd

      real(dp), dimension(:), allocatable :: cent_dist, R_conv
      real(dp), dimension(:, :) :: cent_ave
      real(dp), dimension(:, :, :) :: cent_history
      real(dp) :: ave_dist_a, ave_dist_b, stderr_dist_a, stderr_dist_b, error_bar_fssd, R_th, R_max, infinity
      integer :: last_track, iter, miter, track, ic, k, max_conv, n_history, j, N_a, keep_track, Nave
      logical :: check_config_conv

      last_track = keep_track - Nave
      allocate(cent_dist(last_track/micro_iter_sr))

      track = 0

      do iter = 1,last_track/micro_iter_sr
         track = track + 1
         cent_dist(track) = 0.0

         do ic = 1,ncent
            do k = 1,3
               cent_dist(track) = cent_dist(track) + (cent_history(iter, k, ic) - cent_ave(k, ic))**2
               !if (track .eq. keep_track/3) then
               !   write (ounit, '(f10.5)') (cent_history(iter, k, ic) - cent_ave(k, ic))**2
               !endif
            enddo
         enddo
         cent_dist(track) = sqrt(cent_dist(track))
         write (ounit, '(''Distance from reference configuration is '', f10.5)') cent_dist(track)
      enddo

      allocate(R_conv(track))

      R_conv = 0.0
      N_a = 5 ! modify also in optwf_sr
      R_th = 5.0
      R_max = 0.0

      if (track > 2 * N_a + 2) then 
         do j = N_a,track-N_a

            ave_dist_a = 0.0
            ave_dist_b = 0.0
            stderr_dist_a = 0.0
            stderr_dist_b = 0.0

            do k = 1,j
               ave_dist_a = ave_dist_a + cent_dist(k) / j
            enddo

            do k = j+1,track
               ave_dist_b = ave_dist_b + cent_dist(k) / (track - j)
            enddo

            do k = 1,j
               stderr_dist_a = stderr_dist_a + (cent_dist(k) - ave_dist_a)**2 / j
            enddo

            do k = j+1,track
               stderr_dist_b = stderr_dist_b + (cent_dist(k) - ave_dist_b)**2 / (track - j)
            enddo

            stderr_dist_a = sqrt(stderr_dist_a)
            stderr_dist_b = sqrt(stderr_dist_b)
            write (ounit, '(''stderr(D) = '',f10.2,'' before '',i5, '' and stderr(D) = '',f10.2,'' after '',i5)') stderr_dist_a, j, stderr_dist_b, j

            R_conv(j) = stderr_dist_a / stderr_dist_b

            write (ounit, '(''R = '',f10.2,'' at step '',i5)') R_conv(j), j

            !if (j .ge. 2) then
            !  if (R_conv(j) .gt. R_max) R_max = R_conv(j)
            !else
            !  R_max = R_conv(j)
            !endif

            infinity = 1.d0 / 0.d0

            if (R_conv(j) .ge. R_th .and. R_conv(j) .lt. infinity) then
               write (ounit, *) "Convergence is reached for the configuration!"
               write (ounit, '(''R '', f10.5)')  R_conv(j)
               if (.not.check_config_conv) then
                  cent_ave = 0.0
                  do iter = j,keep_track
                     write (ounit, '(/,''Averaging '', i5)') iter
                     do ic = 1,ncent
                        write(ounit,'(3f10.5)') (cent_history(iter,k,ic),k=1,3)
                        do k = 1,3
                           cent_ave(k, ic) = cent_ave(k, ic) + cent_history(iter, k, ic) / (keep_track - j + 1)
                        enddo
                     enddo
                  enddo
                  write (ounit, *) 'Reference configuration is: '
                  !write (ounit,*) 'CENT'
                  do ic=1,ncent
                     write(ounit,'(3f10.5)') (cent_ave(k,ic),k=1,3)
                  enddo
                  check_config_conv = .true.
               endif

            endif

         enddo
      endif

      !if(R_max .ge. R_th .and. R_max .lt. infinity) then
      !   write (ounit, *) "Convergence is reached for the configuration!"
      !   write (ounit, '(''R_max '', f10.5)')  R_max
      !   write (ounit, *) 'Reference configuration is: '
      !   !write (ounit,*) 'CENT'
      !   do ic=1,ncent
      !      write(ounit,'(3f10.5)') (cent_ave(k,ic),k=1,3)
      !   enddo
      !    check_config_conv = .true.
      ! endif

      do ic = 1,ncent
         do k = 1,3
            cent(k, ic) = cent_ave(k, ic)
         enddo
      enddo

      !write (ounit, *) 'Average configuration after checking convergence is: '
      !do ic=1,ncent
      !   write(ounit,'(3f10.5)') (cent(k,ic),k=1,3)
      !enddo

      if (check_config_conv .eq. .false.) then
         write(ounit, *) 'Running another FSSD'

         !write(ounit, *) 'Setting parameters for next SET'
         !error_bar_fssd = error_bar_fssd / 2.d0
         !alfgeo = alfgeo / 2.d0
         !d_fssd = 0.0
      endif

      deallocate(cent_dist)
      deallocate(R_conv)
    end subroutine config_check

    subroutine save_params()

        sr_adiag_sav = sr_adiag
        iforce_analy_sav = iforce_analy
        ioptjas_sav = ioptjas
        ioptorb_sav = ioptorb
        ioptci_sav = ioptci

    end subroutine save_params

    subroutine sr(nparm, deltap, sr_adiag, sr_eps, i)

        ! solve S*deltap=h_sr (call in optwf)
        use sr_more, only: pcg
        use sr_mat_n, only: h_sr
        use contrl_file,    only: ounit
        implicit none

        integer, intent(in) :: nparm
        real(dp), dimension(:), intent(inout) :: deltap
        real(dp), intent(in) :: sr_adiag
        real(dp), intent(in) :: sr_eps
        integer, intent(inout) :: i

        integer :: imax, imod

        call sr_hs(nparm, sr_adiag)

        imax = nparm          ! max n. iterations conjugate gradients
        imod = 50             ! inv. freq. of calc. r=b-Ax vs. r=r-alpha q (see pcg)
        do i = 1, nparm
            deltap(i) = 0.d0     ! initial guess of solution
        enddo
        call pcg(nparm, h_sr, deltap, i, imax, imod, sr_eps)
        !write (ounit, *) 'CG iter ', i

        call sr_rescale_deltap(nparm, deltap)
        return

    end subroutine sr

    subroutine check_length_run_sr(iter, increase_nblk, nblk, nblk_max, denergy, denergy_err, energy_err_sav, energy_tol)
        use contrl_file,    only: ounit
        implicit none

        integer :: iter, increase_nblk, nblk, nblk_max, nblk_new, nbkl
        real(dp) :: energy_err_sav, denergy, denergy_err, energy_tol

        ! Increase nblk if near convergence to value needed to get desired statistical error
        increase_nblk = increase_nblk + 1

        ! Increase if subsequent energies are within errorbar
        if (dabs(denergy) .lt. 3*denergy_err .and. energy_tol .gt. 0.d0) then
            nblk_new = nblk*max(1.d0, (energy_err_sav/energy_tol)**2)
            nblk_new = min(nblk_new, nblk_max)
            if (nblk_new .gt. nblk) then
                increase_nblk = 0
                nblk = nblk_new
                write(ounit, '(''nblk reset to'',i8,9d12.4)') nblk, dabs(denergy), energy_tol
            endif
        endif

        ! Always increase nblk by a factor of 2 every other iteration
        if (increase_nblk .eq. 2 .and. nblk .lt. nblk_max) then
            increase_nblk = 0
            nbkl = 1.2*nblk
            nblk = min(nblk, nblk_max)
            write(ounit, '(''nblk reset to'',i8,9d12.4)') nblk
        endif

        return
    end

    subroutine sr_hs(nparm, sr_adiag)
        ! <elo>, <o_i>, <elo o_i>, <o_i o_i>; s_diag, s_ii_inv, h_sr

        use mpi
        use sr_mod, only: mobs
        use csfs, only: nstates
        use mstates_mod, only: MSTATES
        use mpiconf, only: idtask
        use optwf_func, only: ifunc_omega, omega
        use sa_weights, only: weights
        use sr_index, only: jelo, jelo2, jelohfj
        use sr_mat_n, only: elocal, h_sr, jefj, jfj, jhfj, nconf_n, s_diag, s_ii_inv, sr_ho
        use sr_mat_n, only: sr_o, wtg, obs_tot
        use optorb_cblock, only: norbterm
        use method_opt, only: method
        use contrl_file,    only: ounit
        implicit none

        real(dp), DIMENSION(:, :), allocatable :: obs
        real(dp), DIMENSION(:), allocatable :: obs_wtg
        real(dp), DIMENSION(:), allocatable :: obs_wtg_tot
        integer :: i, j, k, kk, ish, istate, nparm, iconf
        integer :: nstates_eff, jwtg
        integer :: jfifj, jfhfj, n_obs, nparm_jasci, ier, nreduce
        real(dp) :: sr_adiag, var, wts, aux, den, dum1, dum2, smax
        real(dp), parameter :: eps_eigval = 0.d0 ! in the original implementation, it is not used

        allocate (obs_wtg(MSTATES))
        allocate (obs_wtg_tot(MSTATES))
        allocate (obs(MOBS, MSTATES))

        nstates_eff = nstates
        if (method .eq. 'lin_d') nstates_eff = 1

        jwtg = 1
        jelo = 2
        n_obs = 2
        jfj = n_obs + 1
        n_obs = n_obs + nparm
        jefj = n_obs + 1
        n_obs = n_obs + nparm
        jfifj = n_obs + 1
        n_obs = n_obs + nparm

        jhfj = n_obs + 1
        n_obs = n_obs + nparm
        jfhfj = n_obs + 1
        n_obs = n_obs + nparm

        ! for omega functional
        jelo2 = n_obs + 1
        n_obs = n_obs + 1
        jelohfj = n_obs + 1
        n_obs = n_obs + nparm

        if (n_obs .gt. mobs) call fatal_error('SR_HS LIN: n_obs > mobs)')

        do k = 1, nparm
            h_sr(k) = 0.d0
            s_ii_inv(k) = 0.d0
        enddo

        nparm_jasci = max(nparm - norbterm, 0)

        do istate = 1, nstates
            obs(jwtg, istate) = 0.d0
            do iconf = 1, nconf_n
                obs(jwtg, istate) = obs(jwtg, istate) + wtg(iconf, istate)
            enddo
            obs_wtg(istate) = obs(jwtg, istate)
        enddo

        call MPI_REDUCE(obs_wtg, obs_wtg_tot, nstates, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ier)
        do istate = 1, nstates
            obs_tot(jwtg, istate) = obs_wtg_tot(istate)
        enddo

        do istate = 1, nstates_eff
            do i = 2, n_obs
                obs(i, istate) = 0.d0
            enddo

            ish = (istate - 1)*norbterm
            do iconf = 1, nconf_n
                obs(jelo, istate) = obs(jelo, istate) + elocal(iconf, istate)*wtg(iconf, istate)
                do i = 1, nparm_jasci
                    obs(jfj + i - 1, istate) = obs(jfj + i - 1, istate) + sr_o(i, iconf)*wtg(iconf, istate)
                    obs(jefj + i - 1, istate) = obs(jefj + i - 1, istate) + elocal(iconf, istate)*sr_o(i, iconf)*wtg(iconf, istate)
                    obs(jfifj + i - 1, istate) = obs(jfifj + i - 1, istate) + sr_o(i, iconf)*sr_o(i, iconf)*wtg(iconf, istate)
                enddo
                do i = nparm_jasci + 1, nparm
                    obs(jfj + i - 1, istate) = obs(jfj + i - 1, istate) + sr_o(ish + i, iconf)*wtg(iconf, istate)
                    obs(jefj + i - 1, istate) = obs(jefj + i - 1, istate) + &
                                    elocal(iconf, istate)*sr_o(ish + i, iconf)*wtg(iconf, istate)
                    obs(jfifj + i - 1, istate) = obs(jfifj + i - 1, istate) + &
                                    sr_o(ish + i, iconf)*sr_o(ish + i, iconf)*wtg(iconf, istate)
                enddo
            enddo

            call MPI_REDUCE(obs(1, istate), obs_tot(1, istate), n_obs, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ier)
        enddo

        if (idtask .eq. 0) then
            do istate = 1, nstates_eff
                wts = weights(istate)
                if (method .eq. 'lin_d') wts = 1.d0

                do i = 2, n_obs
                    obs_tot(i, istate) = obs_tot(i, istate)/obs_tot(1, istate)
                enddo

                do k = 1, nparm
                    aux = obs_tot(jfifj + k - 1, istate) - obs_tot(jfj + k - 1, istate)*obs_tot(jfj + k - 1, istate)
                    s_diag(k, istate) = aux*sr_adiag
                    s_ii_inv(k) = s_ii_inv(k) + wts*(aux + s_diag(k, istate))
                    h_sr(k) = h_sr(k) - 2*wts*(obs_tot(jefj + k - 1, istate) - obs_tot(jfj + k - 1, istate)*obs_tot(jelo, istate))
                enddo
            enddo

            smax = 0.d0
            do k = 1, nparm
                if (s_ii_inv(k) .gt. smax) smax = s_ii_inv(k)
            enddo
            write (ounit, '(''max S diagonal element '',t41,f16.8)') smax

            kk = 0
            do k = 1, nparm
                if (s_ii_inv(k)/smax .gt. eps_eigval) then
                    kk = kk + 1
                    s_ii_inv(k) = 1.d0/s_ii_inv(k)
                else
                    s_ii_inv(k) = 0.d0
                endif
            enddo
            write (ounit, '(''nparm, non-zero S diag'',t41,2i5)') nparm, kk

        endif

        if (method .eq. 'sr_n' .and. i_sr_rescale .eq. 0 .and. izvzb .eq. 0 .and. ifunc_omega .eq. 0) return

        if (method .ne. 'sr_n') then
            s_diag(1, 1) = sr_adiag !!!

            do k = 1, nparm
                h_sr(k) = -0.5d0*h_sr(k)
            enddo
        elseif (ifunc_omega .ne. 0) then
            s_diag(1, 1) = sr_adiag !!!
        endif

        if (n_obs .gt. mobs) call fatal_error('SR_HS LIN: n_obs > mobs)')

        do i = jhfj, n_obs
            obs(i, 1) = 0.d0
        enddo
        do iconf = 1, nconf_n
            obs(jelo2, 1) = obs(jelo2, 1) + elocal(iconf, 1)*elocal(iconf, 1)*wtg(iconf, 1)
            do i = 1, nparm
                obs(jhfj + i - 1, 1) = obs(jhfj + i - 1, 1) + sr_ho(i, iconf)*wtg(iconf, 1)
                obs(jfhfj + i - 1, 1) = obs(jfhfj + i - 1, 1) + sr_o(i, iconf)*sr_ho(i, iconf)*wtg(iconf, 1)
                obs(jelohfj + i - 1, 1) = obs(jelohfj + i - 1, 1) + elocal(iconf, 1)*sr_ho(i, iconf)*wtg(iconf, 1)
            enddo
        enddo

        nreduce = n_obs - jhfj + 1
        call MPI_REDUCE(obs(jhfj, 1), obs_tot(jhfj, 1), nreduce, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, j)

        if (idtask .eq. 0) then
            do i = jhfj, n_obs
                obs_tot(i, 1) = obs_tot(i, 1)/obs_tot(1, 1)
            enddo

            if (ifunc_omega .eq. 1) then
                ! variance
                var = obs_tot(jelo2, 1) - obs_tot(jelo, 1)**2
                do k = 1, nparm
                h_sr(k) = -2*(obs_tot(jelohfj + k - 1, 1) - (obs_tot(jhfj + k - 1, 1) - obs_tot(jefj + k - 1, 1))*obs_tot(jelo, 1)&
                                  - obs_tot(jfj + k - 1, 1)*obs_tot(jelo2, 1) &
                                  - 2*obs_tot(jelo, 1)*(obs_tot(jefj + k - 1, 1) - obs_tot(jfj + k - 1, 1)*obs_tot(jelo, 1)))
                enddo
            elseif (ifunc_omega .eq. 2) then
                ! variance with fixed average energy (omega)
                var = omega*omega + obs_tot(jelo2, 1) - 2*omega*obs_tot(jelo, 1)
                dum1 = -2
                do k = 1, nparm
                    h_sr(k) = dum1*(omega*omega*obs_tot(jfj + k - 1, 1) + obs_tot(jelohfj + k - 1, 1) &
                                    - omega*(obs_tot(jhfj + k - 1, 1) + obs_tot(jefj + k - 1, 1)) &
                                    - var*obs_tot(jfj + k - 1, 1))
                    ! adding a term which intergrates to zero
                    !    &     -(obs_tot(jelo,1)-omega)*(obs_tot(jhfj+k-1,1)-obs_tot(jefj+k-1,1)))
                enddo

            elseif (ifunc_omega .eq. 3 .and. method .eq. 'sr_n') then
                !  Neuscamman's functional
                den = omega*omega + obs_tot(jelo2, 1) - 2*omega*obs_tot(jelo, 1)
                dum1 = -2/den
                dum2 = (omega - obs_tot(jelo, 1))/den
                do k = 1, nparm
                    h_sr(k) = dum1*(omega*obs_tot(jfj + k - 1, 1) - obs_tot(jefj + k - 1, 1) &
                                    - dum2*(omega*omega*obs_tot(jfj + k - 1, 1) + obs_tot(jelohfj + k - 1, 1) &
                                            - omega*(obs_tot(jhfj + k - 1, 1) + obs_tot(jefj + k - 1, 1))))
                enddo
            endif

        endif

        deallocate (obs_wtg)
        deallocate (obs_wtg_tot)
        deallocate (obs)

        return
    end

    subroutine sr_rescale_deltap(nparm, deltap)

        use mpi
        use mpiconf, only: idtask
        use sr_mat_n, only: jefj, jfj, jhfj
        use sr_mat_n, only: obs_tot
        use sr_index, only: jelo, jelo2, jelohfj
        use optwf_contrl, only: sr_tau
        use contrl_file,    only: ounit

        implicit none

        integer, intent(in)                     :: nparm
        real(dp) :: de, top, bot
        real(dp), dimension(:), intent(inout)   :: deltap
        integer :: i, j, jfifj, jwtg, jfhfj, n_obs

        if (i_sr_rescale .eq. 0) return

        jwtg = 1
        jelo = 2
        n_obs = 2
        jfj = n_obs + 1
        n_obs = n_obs + nparm
        jefj = n_obs + 1
        n_obs = n_obs + nparm
        jfifj = n_obs + 1
        n_obs = n_obs + nparm

        jhfj = n_obs + 1
        n_obs = n_obs + nparm
        jfhfj = n_obs + 1
        n_obs = n_obs + nparm

        jelo2 = n_obs + 1
        n_obs = n_obs + 1
        jelohfj = n_obs + 1
        n_obs = n_obs + nparm

        if (idtask .eq. 0) then
        do i = 1, nparm
            de=obs_tot(jfhfj + i - 1, 1)/obs_tot(jfifj + i - 1, 1)-obs_tot(jelo, 1)
!           write(ounit, *) 'CIAO', obs_tot(jfhfj + i - 1, 1)/obs_tot(jfifj + i - 1, 1), obs_tot(jelo, 1), &
!               obs_tot(jfhfj + i - 1, 1)/obs_tot(jfifj + i - 1, 1) - obs_tot(jelo, 1)

            top = obs_tot(jfhfj + i - 1, 1) + obs_tot(jelo, 1)*obs_tot(jfj + i - 1, 1)*obs_tot(jfj + i - 1, 1) &
                - (obs_tot(jefj + i - 1, 1) + obs_tot(jhfj + i - 1, 1))*obs_tot(jfj + i - 1, 1)
            bot = obs_tot(jfifj + i - 1, 1) - obs_tot(jfj + i - 1, 1)*obs_tot(jfj + i - 1, 1)
            de = top/bot - obs_tot(jelo, 1)

            de = de + 1.0

            write(ounit,*) "CIAO", de

            if(de .le. 0.1) then
                deltap(i) = deltap(i)*0.01/sr_tau
            else
                deltap(i) = deltap(i)/de
            endif
        enddo
        endif

        call MPI_BCAST(deltap, nparm, MPI_REAL8, 0, MPI_COMM_WORLD, j)

        return
    end subroutine sr_rescale_deltap

    subroutine forces_zvzb(nparm)

        use mpi
        use sr_mod, only: mparm
        use atom, only: ncent
        use force_fin, only: da_energy_ave
        use force_mat_n, only: force_o
        use mpiconf, only: idtask
        use sr_mat_n, only: elocal, jefj, jfj, jhfj, nconf_n, obs_tot, sr_ho
        use sr_mat_n, only: sr_o, wtg
        use sr_index, only: jelo
        use contrl_file,    only: ounit
        use error, only: fatal_error

        implicit none

        integer, parameter :: MTEST = 1500
        real(dp), dimension(:, :), allocatable :: cloc
        real(dp), dimension(:, :), allocatable :: c
        real(dp), dimension(:), allocatable :: oloc
        real(dp), dimension(:), allocatable :: o
        real(dp), dimension(:), allocatable :: p
        real(dp), dimension(:), allocatable :: tmp
        real(dp), dimension(:), allocatable :: work
        integer, dimension(:), allocatable :: ipvt
        integer :: nparm, i, j, k, jfhfj, ia, icent, iparm, ish, n_obs, l, info
        integer :: jparm, jwtg, jfifj
        real(dp) :: dum, energy_tot, force_tmp, wtoti

        wtoti = 0.0
        allocate (cloc(MTEST, MTEST))
        allocate (c(MTEST, MTEST))
        allocate (oloc(mparm))
        allocate (o(mparm))
        allocate (p(mparm))
        allocate (tmp(mparm))
        allocate (work(MTEST))
        allocate (ipvt(MTEST))

        if (nparm .gt. MTEST) stop 'mparm>MTEST'

        jwtg = 1
        jelo = 2
        n_obs = 2
        jfj = n_obs + 1
        n_obs = n_obs + nparm
        jefj = n_obs + 1
        n_obs = n_obs + nparm
        jfifj = n_obs + 1
        n_obs = n_obs + nparm

        jhfj = n_obs + 1
        n_obs = n_obs + nparm
        jfhfj = n_obs + 1
        n_obs = n_obs + nparm

        do i = 1, nparm
            do j = i, nparm
                cloc(i, j) = 0.d0
            enddo
        enddo

        do l = 1, nconf_n
            do i = 1, nparm
                tmp(i) = (sr_ho(i, l) - elocal(l, 1)*sr_o(i, l))*sqrt(wtg(l, 1))
            enddo

            do k = 1, nparm
                do j = k, nparm
                    cloc(k, j) = cloc(k, j) + tmp(k)*tmp(j)
                enddo
            enddo
        enddo

        call MPI_REDUCE(cloc, c, MTEST*nparm, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, i)

        if (idtask .eq. 0) then

            wtoti = 1.d0/obs_tot(1, 1)
            do i = 1, nparm
                dum = (obs_tot(jhfj + i - 1, 1) - obs_tot(jefj + i - 1, 1))
                c(i, i) = c(i, i)*wtoti - dum*dum
                do j = i + 1, nparm
                    c(i, j) = c(i, j)*wtoti - dum*(obs_tot(jhfj + j - 1, 1) - obs_tot(jefj + j - 1, 1))
                    c(j, i) = c(i, j)
                enddo
            enddo

            call dgetrf(nparm, nparm, c, MTEST, ipvt, info)
            if (info .gt. 0) then
                write(ounit,'(''optwf_sr.f'')')
                write(ounit, '(''MATINV: u(k,k)=0 with k= '',i5)') info
                call fatal_error('MATINV: info ne 0 in dgetrf')
            endif
            call dgetri(nparm, c, MTEST, ipvt, work, MTEST, info)

            do iparm = 1, nparm
                tmp(iparm) = obs_tot(jhfj + iparm - 1, 1) - obs_tot(jefj + iparm - 1, 1)
            enddo

        endif

        energy_tot = obs_tot(2, 1)

        call MPI_BCAST(energy_tot, 1, MPI_REAL8, 0, MPI_COMM_WORLD, j)

        ia = 0
        ish = 3*ncent
        do icent = 1, ncent
            write(ounit, '(''FORCE before'',i4,3e15.7)') icent, (da_energy_ave(k, icent), k=1, 3)
            do k = 1, 3
                ia = ia + 1

                do i = 1, nparm
                    oloc(i) = 0.d0
                    do l = 1, nconf_n
                        oloc(i) = oloc(i) + (sr_ho(i, l) - elocal(l, 1)*sr_o(i, l)) &
                                  *(force_o(ia + ish, l) - 2*energy_tot*force_o(ia, l))*wtg(l, 1)
                    enddo
                enddo

                call MPI_REDUCE(oloc, o, nparm, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, i)

                if (idtask .eq. 0) then

                    do i = 1, nparm
                        o(i) = o(i)*wtoti - (obs_tot(jhfj + i - 1, 1) - obs_tot(jefj + i - 1, 1))*da_energy_ave(k, icent)
                    enddo

                    do iparm = 1, nparm
                        p(iparm) = 0.d0
                        do jparm = 1, nparm
                            p(iparm) = p(iparm) + c(iparm, jparm)*o(jparm)
                        enddo
                        p(iparm) = -0.5*p(iparm)
                    enddo

                    force_tmp = da_energy_ave(k, icent)
                    do iparm = 1, nparm
                        force_tmp = force_tmp + p(iparm)*tmp(iparm)
                    enddo
                    da_energy_ave(k, icent) = force_tmp

                endif
            enddo
            write(ounit, '(''FORCE after '',i4,3e15.7)') icent, (da_energy_ave(k, icent), k=1, 3)
        enddo

        return
    end subroutine forces_zvzb

end module optwf_sr_mod
