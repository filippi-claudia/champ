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

    real(dp) :: omega0
    integer :: n_omegaf, n_omegat
    integer ::  i_func_omega
    real(dp) :: sr_adiag_sav

    integer :: ioptjas_sav, ioptorb_sav, ioptci_sav, iforce_analy_sav
    real(dp), dimension(:,:), allocatable :: deltap

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

        use sr_mod, only: mparm
        use sr_mat_n, only: sr_lambda, sr_state
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

        use csfs, only: nstates

        implicit none

        real(dp) :: adiag, denergy, alpha_omega, denergy_err, dparm_norm
        real(dp) :: energy_sav, energy_err_sav, omega, sigma, sigma_sav, nadorb_sav
        integer :: i, iflag, iter, miter, istate, jstate
        integer :: iflagin

        sigma_sav = 0.0
        energy_sav= 0.0
        !allocate (deltap(mparm*MSTATES))
        allocate (deltap(mparm,MSTATES))

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
        if (nstates .gt. 1) then
            write (ounit, '(''SR lambda:   '',f10.5)') & 
                  ((istate,jstate,sr_lambda(istate,jstate),istate=1,nstates-1),jstate=istate+1,nstates)
        endif

        call save_params()

        call save_nparms

        call write_geometry(0)

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

            ! do micro_iteration
            do miter = 1, micro_iter_sr

                if (micro_iter_sr .gt. 1) write (ounit, '(/,''Micro iteration'',i5,'' of'',i5)') miter, micro_iter_sr

                if (miter .eq. micro_iter_sr) iforce_analy = iforce_analy_sav

                call qmc

                write (ounit, '(/,''Completed sampling'')')

6               continue

                !call sr(nparm, deltap, sr_adiag, sr_eps, i)
                !call dscal(nparm, -sr_tau, deltap, 1)

                !new sr part
                if (nstates .eq. 1) then
                    call sr_hs(nparm, sr_adiag)
                elseif (nstates .gt. 1) then
                    call compute_gradients_sr_ortho(nparm, sr_adiag)
                endif
                
                ! maybe distribute pcg over nodes at some point, reduce by MSTATE
                adiag = sr_adiag
                iflagin=0
                do istate = 1, nstates
                    sr_state = istate
                    call sr(nparm, deltap(:, istate), sr_adiag, sr_eps, i)
                    call dscal(nparm, -sr_tau, deltap(:, istate), 1)
                    call test_solution_parm(nparm, deltap(:, istate), dparm_norm, dparm_norm_min, adiag, iflag)
                    write (ounit, '(''Norm of parm variation '',d12.5)') dparm_norm
                    if (iflag .ne. 0) iflagin = 1
                enddo

                !call test_solution_parm(nparm, deltap, dparm_norm, dparm_norm_min, adiag, iflag)
                !write (ounit, '(''Norm of parm variation '',d12.5)') dparm_norm
                if (iflag .ne. 0) then
                    write (ounit, '(''Warning: dparm_norm>1'')')
                    adiag = 10*adiag
                    write (ounit, '(''adiag increased to '',f10.5)') adiag

                    sr_adiag = adiag
                    go to 6
                else
                    sr_adiag = sr_adiag_sav
                endif
                
                ! if statement? do we only use compute norm lin in multistate?
                if(nstates.gt.1) call compute_norm_lin(nparm,-deltap)
                call compute_parameters(deltap, iflag, 1)
                call write_wf(1, iter)

                call save_wf

                if (iforce_analy .gt. 0) then

                    if (izvzb .gt. 0) call forces_zvzb(nparm)

                    call compute_positions
                    call write_geometry(iter)
                endif
                call elapsed_time("CG micro iteration", miter)
            enddo
            ! enddo micro_iteration

            if (iter .ge. 2) then
                denergy = energy(1) - energy_sav
                denergy_err = sqrt(energy_err(1)**2 + energy_err_sav**2)

                vmc_nblk = vmc_nblk*1.2
                vmc_nblk = min(vmc_nblk, vmc_nblk_max)

            endif
            write (ounit, '(''nblk = '',i6)') vmc_nblk
            write (ounit, '(''alfgeo = '',f10.4)') alfgeo

            energy_sav = energy(1)
            energy_err_sav = energy_err(1)
            sigma_sav = sigma
            call elapsed_time( "CG iteration ", iter )
        enddo
        ! enddo iteration

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
    end

    subroutine compute_norm_lin(nparm,deltap)  ! check these modules and rest of subroutine
        use mpi
        use mpiconf, only: idtask
        use sr_mod, only: mobs, mparm
        use mstates_mod, only: MSTATES
        use csfs, only: nstates, anormo
        use sr_mat_n, only: elocal, h_sr, jefj, jfj, jhfj, nconf_n, s_diag, s_ii_inv, sr_ho
        use sr_mat_n, only: sr_o, wtg, obs_tot
        use estcum, only: iblk
        use control_vmc, only: vmc_nstep
        use optorb_cblock, only: norbterm
        use contrl_file, only: ounit

        implicit none

        integer, intent(in) :: nparm
        real(dp), dimension(mparm,MSTATES), intent(in) :: deltap
        real(dp) :: passes, tmp, tmp1
        real(dp), dimension(mobs,MSTATES) :: obs_norm
        real(dp), dimension(mobs,MSTATES) :: obs_norm_tot

        integer :: jwtg, jwfj, jsqfj, n_obs, istate, iconf, i, ier

        jwtg = 1
        jwfj = 1
        jsqfj = 2
        n_obs = 2

        obs_norm = 0.0d0
        passes = dfloat(iblk*vmc_nstep)

        do istate = 1, nstates
           do iconf = 1, nconf_n
              tmp = 0.0d0
              do i = 1, nparm
                 tmp = tmp + deltap(i, istate)*sr_o(iconf, i, istate) ! check for bugs, was a '.' not a ','
              enddo
              obs_norm(jwfj, istate) = obs_norm(jwfj, istate) + tmp*wtg(iconf, istate)
              obs_norm(jsqfj,istate) = obs_norm(jsqfj, istate) + tmp*tmp*wtg(iconf, istate)
           enddo
           obs_norm(jwfj, istate) = 2.0d0*obs_norm(jwfj, istate)
           call MPI_REDUCE(obs_norm(1, istate), obs_norm_tot(1, istate),&
                   n_obs, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ier)
        enddo

        if(idtask.eq.0) then
           do istate = 1, nstates
              tmp = obs_tot(jwtg, istate) + obs_norm_tot(jwfj, istate)
              if(istate.eq.1) tmp1 = tmp
              anormo(istate) = obs_tot(jwtg, istate) + obs_norm_tot(jwfj, istate)&
                      + obs_norm_tot(jsqfj, istate)
              write(ounit, '(''NORMS'',i3,10f20.8)') istate, obs_tot(jwtg,istate)/passes, anormo(istate)/passes,&
                                 obs_tot(jwtg, istate)/obs_tot(jwtg, 1), anormo(istate)/anormo(1),tmp/tmp1
           enddo
           ! Ramon was using the anormo
           !anormo=anormo/passes
           ! Claudia tries this
           do istate = 2, nstates
             anormo(istate) = anormo(istate)/anormo(1)
           enddo
           anormo(1) = 1.d0
        endif
        ! TMP
        ! anormo=1.d0
        call MPI_BCAST(anormo(1), nstates, MPI_REAL8, 0, MPI_COMM_WORLD, i)
        !call bcast(anormo(1)) ! check this bcast, we set this to 1 already, do we need to bcast now?

    end subroutine

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
        use sr_mat_n, only: h_sr, sr_state
        use contrl_file,    only: ounit
        implicit none

        integer, intent(in) :: nparm
        real(dp), dimension(:), intent(inout) :: deltap
        real(dp), intent(in) :: sr_adiag
        real(dp), intent(in) :: sr_eps
        integer, intent(inout) :: i

        integer :: imax, imod

        !call sr_hs(nparm, sr_adiag) ! moved outside
        
        imax = nparm          ! max n. iterations conjugate gradients
        imod = 50             ! inv. freq. of calc. r=b-Ax vs. r=r-alpha q (see pcg)
        do i = 1, nparm
            deltap(i) = 0.d0     ! initial guess of solution
        enddo
        
        call pcg(nparm, h_sr(1:nparm,sr_state), deltap, i, imax, imod, sr_eps)
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
        allocate (obs(mobs, MSTATES))

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
            h_sr(k,1) = 0.d0
            s_ii_inv(k,1) = 0.d0
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
                    obs(jfj + i - 1, istate) = obs(jfj + i - 1, istate) + sr_o(i, iconf, 1)*wtg(iconf, istate)
                    obs(jefj + i - 1, istate) = obs(jefj + i - 1, istate) + elocal(iconf, istate)*sr_o(i, iconf, 1)*wtg(iconf, istate)
                    obs(jfifj + i - 1, istate) = obs(jfifj + i - 1, istate) + sr_o(i, iconf, 1)*sr_o(i, iconf, 1)*wtg(iconf, istate)
                enddo
                do i = nparm_jasci + 1, nparm
                    obs(jfj + i - 1, istate) = obs(jfj + i - 1, istate) + sr_o(ish + i, iconf, 1)*wtg(iconf, istate)
                    obs(jefj + i - 1, istate) = obs(jefj + i - 1, istate) + &
                                    elocal(iconf, istate)*sr_o(ish + i, iconf, 1)*wtg(iconf, istate)
                    obs(jfifj + i - 1, istate) = obs(jfifj + i - 1, istate) + &
                                    sr_o(ish + i, iconf, 1)*sr_o(ish + i, iconf, 1)*wtg(iconf, istate)
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
                    s_ii_inv(k,1) = s_ii_inv(k,1) + wts*(aux + s_diag(k, istate))
                    h_sr(k,1) = h_sr(k,1) - 2*wts*(obs_tot(jefj + k - 1, istate) - obs_tot(jfj + k - 1, istate)*obs_tot(jelo, istate))
                enddo
            enddo

            smax = 0.d0
            do k = 1, nparm
                if (s_ii_inv(k,1) .gt. smax) smax = s_ii_inv(k,1)
            enddo
            write (ounit, '(''max S diagonal element '',t41,f16.8)') smax

            kk = 0
            do k = 1, nparm
                if (s_ii_inv(k,1)/smax .gt. eps_eigval) then
                    kk = kk + 1
                    s_ii_inv(k,1) = 1.d0/s_ii_inv(k,1)
                else
                    s_ii_inv(k,1) = 0.d0
                endif
            enddo
            write (ounit, '(''nparm, non-zero S diag'',t41,2i5)') nparm, kk

        endif
        
        if (method .eq. 'sr_n' .and. i_sr_rescale .eq. 0 .and. izvzb .eq. 0 .and. ifunc_omega .eq. 0) return
        
        if (method .ne. 'sr_n') then
            s_diag(1, 1) = sr_adiag !!!

            do k = 1, nparm
                h_sr(k,1) = -0.5d0*h_sr(k,1)
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
                obs(jfhfj + i - 1, 1) = obs(jfhfj + i - 1, 1) + sr_o(i, iconf, 1)*sr_ho(i, iconf)*wtg(iconf, 1)
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
                h_sr(k,1) = -2*(obs_tot(jelohfj + k - 1, 1) - (obs_tot(jhfj + k - 1, 1) - obs_tot(jefj + k - 1, 1))*obs_tot(jelo, 1)&
                                  - obs_tot(jfj + k - 1, 1)*obs_tot(jelo2, 1) &
                                  - 2*obs_tot(jelo, 1)*(obs_tot(jefj + k - 1, 1) - obs_tot(jfj + k - 1, 1)*obs_tot(jelo, 1)))
                enddo
            elseif (ifunc_omega .eq. 2) then
                ! variance with fixed average energy (omega)
                var = omega*omega + obs_tot(jelo2, 1) - 2*omega*obs_tot(jelo, 1)
                dum1 = -2
                do k = 1, nparm
                    h_sr(k,1) = dum1*(omega*omega*obs_tot(jfj + k - 1, 1) + obs_tot(jelohfj + k - 1, 1) &
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
                    h_sr(k,1) = dum1*(omega*obs_tot(jfj + k - 1, 1) - obs_tot(jefj + k - 1, 1) &
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

    subroutine compute_gradients_sr_ortho(nparm,sr_adiag) ! check subroutine
        use mpi
        use sr_mod, only: mobs
        use csfs, only: nstates
        use mstates_mod, only: MSTATES
        use mpiconf, only: idtask
        use optwf_func, only: ifunc_omega, omega
        use sr_index, only: jelo, jelo2, jelohfj
        use sr_mat_n, only: elocal, h_sr, jefj, jfj, jhfj, nconf_n, s_diag, s_ii_inv, sr_ho
        use sr_mat_n, only: sr_o, wtg, obs_tot, h_sr_penalty, sr_lambda
        use optorb_cblock, only: norbterm
        use method_opt, only: method

        implicit none

        integer, intent(in) :: nparm

        real(dp), DIMENSION(:, :), allocatable :: obs
        real(dp) :: sr_adiag, dum, aux2, smax, penalty
        real(dp), parameter :: eps_eigval=1.d-14

        integer :: jwtg, jfifj, jfjsi, n_obs, istate, jstate, iconf, i, ier, k, kk, n_obs_reduce


        jwtg=1
        n_obs=nstates+1
        jelo=n_obs+1
        n_obs=n_obs+1
        jfj=n_obs+1
        n_obs=n_obs+nparm
        jefj=n_obs+1
        n_obs=n_obs+nparm
        jfifj=n_obs+1
        n_obs=n_obs+nparm

        if(nstates.gt.1) then
           jfjsi=n_obs+1
           n_obs=n_obs+nparm
        endif

        if(n_obs.gt.mobs) call fatal_error('SR_HS LIN: n_obs_ortho BS)')

    ! initialize
        do istate=1,nstates
           do i=1,n_obs
              obs(i,istate) = 0.0d0
           enddo
        enddo

        do istate=1,nstates

           do iconf=1,nconf_n

! <psi(istate)^2/psig^2>
              obs(jwtg,istate)=obs(jwtg,istate)+wtg(iconf,istate)

! <(psi(istate)/psig)*(*psi(jstate)/psig)>
              do jstate=1,istate-1
                 obs(jwtg+jstate,istate)=obs(jwtg+jstate,istate)&
                         +sr_o(nparm+2,iconf,jstate)*sr_o(nparm+2,iconf,istate)
              enddo

! <eloc(istate)*psi(istate)^2/psig^2>
              obs(jelo,istate)=obs(jelo,istate)&
                      +elocal(iconf,istate)*wtg(iconf,istate)

              do i=1,nparm
! <psi_i(istate)/psi(istate)*psi(istate)^2/psig^2>
                 obs(jfj+i-1,istate)=obs(jfj+i-1,istate)&
                         +sr_o(i,iconf,istate)*wtg(iconf,istate)

! <eloc(istate)*psi_i(istate)/psi(istate)*psi(istate)^2/psig^2>
                 obs(jefj+i-1,istate)=obs(jefj+i-1,istate)&
                         +elocal(iconf,istate)*sr_o(i,iconf,istate)*wtg(iconf,istate)

! <(psi_i(istate)/psi(istate))^2*psi(istate)^2/psig^2>
                 obs(jfifj+i-1,istate)=obs(jfifj+i-1,istate)&
                         +sr_o(i,iconf,istate)*sr_o(i,iconf,istate)*wtg(iconf,istate)
              enddo
           enddo

! reduce what computed so far
           n_obs_reduce=n_obs
           if(nstates.gt.1) n_obs_reduce=n_obs-nparm

           call MPI_REDUCE(obs(1,istate),obs_tot(1,istate),&
                   n_obs_reduce,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ier)

! bcast overlaps only
           call MPI_BCAST(obs_tot(1,istate),nstates,MPI_REAL8,0,MPI_COMM_WORLD,i)

! symmetrize overlaps
           do jstate=1,istate-1
             obs_tot(jwtg+istate,jstate)=obs_tot(jwtg+jstate,istate)
           enddo

! close loop over istate: overlaps computed among all states
        enddo
! loop over the states

! compute part of gradient of cost function sum_jstate^{istate-1} chi_jstate,istate
! sum_jstate^{istate-1} <(psi_i(istate)/psi(istate)*(psi(istate)/psig)*(psi(jstate/psig)>
!                      *<(psi(istate)/psig)*(psi(jstate)/psig)>_tot/<psi(jstate)^2/psig^2>_tot
        if(nstates.gt.1) then

           do istate=1,nstates

           !call cpu_time(dstart_time)
              do jstate=1,istate-1
                 !dum=obs_tot(jwtg+jstate,istate)/obs_tot(jwtg,jstate)
                 dum=obs_tot(jwtg+jstate,istate)/obs_tot(jwtg,jstate)*sr_lambda(istate,jstate)
                 do iconf=1,nconf_n
                    do i=1,nparm
                       obs(jfjsi+i-1,istate)=obs(jfjsi+i-1,istate)&
                               +sr_o(i,iconf,istate)*sr_o(nparm+2,iconf,istate)*sr_o(nparm+2,iconf,jstate)*dum
                    enddo
                 enddo
              enddo
              do jstate=istate+1,nstates
                 !dum=obs_tot(jwtg+jstate,istate)/obs_tot(jwtg,jstate)
                 dum=obs_tot(jwtg+jstate,istate)/obs_tot(jwtg,jstate)*sr_lambda(istate,jstate)
                 do iconf=1,nconf_n
                    do i=1,nparm
                       obs(jfjsi+i-1,istate)=obs(jfjsi+i-1,istate)&
                               +sr_o(i,iconf,istate)*sr_o(nparm+2,iconf,istate)*sr_o(nparm+2,iconf,jstate)*dum
                    enddo
                 enddo
              enddo
              !call cpu_time(dend_time)
              !print *, "Old loop time (s): ", dend_time-dstart_time

              call MPI_REDUCE(obs(jfjsi,istate),obs_tot(jfjsi,istate),&
                      nparm,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ier)

       ! added to check with new loop
       !if(idtask.eq.0) print *, "state, first part of pen parm 1 tot: ", istate, obs_tot(jfjsi,istate)
           enddo
        endif

! end loop over main istate

        if(idtask.eq.0) then

           do istate=1,nstates
! normalize to current state via <psi(istate)^2/psig^2> (to elminate <psig|psig>)
              do i=2,n_obs
                 obs_tot(i,istate)=obs_tot(i,istate)/obs_tot(jwtg,istate)
              enddo

              do k=1,nparm
                 aux2=obs_tot(jfifj+k-1,istate)&
                         -obs_tot(jfj+k-1,istate)*obs_tot(jfj+k-1,istate)
                 s_diag(k,istate)=aux2*sr_adiag
                 s_ii_inv(k,istate)=aux2+s_diag(k,istate)

! gradient wrt parameters
                 h_sr(k,istate)=-2.0d0*(obs_tot(jefj+k-1,istate)&
                         -obs_tot(jfj+k-1,istate)*obs_tot(jelo,istate))
              enddo

              smax=0.0d0
              do k=1,nparm
                 if(s_ii_inv(k,istate).gt.smax) smax=s_ii_inv(k,istate)
              enddo
              write(6,'(''max S diagonal element '',t41,d9.2)') smax

              kk=0
              do k=1,nparm
                 if(s_ii_inv(k,istate)/smax.gt.eps_eigval) then
                    kk=kk+1
                    s_ii_inv(k,istate)=1.0d0/s_ii_inv(k,istate)
                 else
                    s_ii_inv(k,istate)=0.0d0
                 endif
              enddo
              write(6,'(''nparm, non-zero S diag'',t41,2i5)') nparm,kk

          ! Orthogonalization
              if(nstates.gt.1) then
! sum (<psi(jstate)|psi(istate)>/<psi(istate)|psi(istate)>)^2*<psi(istate)|psi(istate)>/<psi(jstate)|psi(jstate)>
                penalty=0.d0
                do jstate=1,nstates
                   if(jstate.ne.istate) then
                     !penalty=penalty+obs_tot(jwtg+jstate,istate)*obs_tot(jwtg+jstate,istate)&
                     !      *obs_tot(jwtg,istate)/obs_tot(jwtg,jstate)
                     penalty=penalty+obs_tot(jwtg+jstate,istate)*obs_tot(jwtg+jstate,istate)&
                           *obs_tot(jwtg,istate)/obs_tot(jwtg,jstate)*sr_lambda(istate,jstate)
                   endif
                enddo

                do k=1,nparm
                   !h_sr_penalty(k,istate)= -sr_lambda(istate,jstate)*2.d0*(obs_tot(jfjsi+k-1,istate)-penalty*obs_tot(jfj+k-1,istate))
                   h_sr_penalty(k,istate)=  -2.d0*(obs_tot(jfjsi+k-1,istate)-penalty*obs_tot(jfj+k-1,istate))
                enddo
                h_sr(1:nparm,istate)=h_sr(1:nparm,istate)+h_sr_penalty(1:nparm,istate)
                print *, "State", istate
                !print *, "lambda SR ortho", sr_lambda
                print *, "penalty ortho  ", penalty
              endif
           enddo
        endif

    end subroutine compute_gradients_sr_ortho

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
                tmp(i) = (sr_ho(i, l) - elocal(l, 1)*sr_o(i, l, 1))*sqrt(wtg(l, 1))
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
                        oloc(i) = oloc(i) + (sr_ho(i, l) - elocal(l, 1)*sr_o(i, l, 1)) &
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
