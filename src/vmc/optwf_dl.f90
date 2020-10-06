!------------------------------------------------------------------------------
!        Optimization routine using ML inspired optimizers
!------------------------------------------------------------------------------
!> @author
!> Claudia Filippi
!
! DESCRIPTION:
!> Opitmize the wave function parameters using the ADAM optimizer
!
! URL           : https://github.com/filippi-claudia/champ
!---------------------------------------------------------------------------

subroutine optwf_dl()

    use precision_kinds, only: dp
    use sr_mod, only: MPARM
    use optwf_contrl, only: ioptci, ioptjas, ioptorb, nparm
    use optwf_corsam, only: energy, energy_err, force
    use contrl, only: nblk
    use method_opt, only: method
    use optwf_contrl, only: idl_flag

    implicit real*8(a - h, o - z)

    character(20) dl_alg
    real(dp), dimension(:), allocatable :: deltap
    real(dp), dimension(:), allocatable :: dl_momentum
    real(dp), dimension(:), allocatable :: dl_EG_sq
    real(dp), dimension(:), allocatable :: dl_EG
    real(dp), dimension(:), allocatable :: parameters

    allocate (deltap(nparm))
    allocate (dl_momentum(nparm))  !< 'momentum' variables
    allocate (dl_EG_sq(nparm))     !< moving average of past squared gradients
    allocate (dl_EG(nparm))        !< moving average of past gradients
    allocate (parameters(nparm))   !< vector of wave function parameters

    if (method .ne. 'sr_n' .or. idl_flag .eq. 0) return

    write (6, '(''Started dl optimization'')')

    call set_nparms_tot

    if (nparm .gt. MPARM) call fatal_error('SR_OPTWF: nparmtot gt MPARM')

    call p2gtid('optwf:nopt_iter', nopt_iter, 6, 1)
    call p2gtid('optwf:nblk_max', nblk_max, nblk, 1)
    call p2gtfd('optwf:energy_tol', energy_tol, 1.d-3, 1)

    call p2gtfd('optwf:dparm_norm_min', dparm_norm_min, 1.0d0, 1)
    write (6, '(''Starting dparm_norm_min'',g12.4)') dparm_norm_min

    call p2gtfd('optwf:sr_tau', sr_tau, 0.02, 1)
    call p2gtfd('optwf:sr_adiag', sr_adiag, 0.01, 1)
    call p2gtfd('optwf:sr_eps', sr_eps, 0.001, 1)
    call p2gtfd('optwf:dl_mom', dl_mom, 0.0, 1)
    call p2gtid('optwf:idl_flag', idl_flag, 0, 1)
    call p2gtad('optwf:dl_alg', dl_alg, 'nag', 1)

    write (6, '(/,''SR adiag: '',f10.5)') sr_adiag
    write (6, '(''SR tau:   '',f10.5)') sr_tau
    write (6, '(''SR eps:   '',f10.5)') sr_eps
    write (6, '(''DL flag:   '',I10)') idl_flag

    inc_nblk = 0

    ! Initialize vectors to zero
    dl_momentum(:) = 0.d0
    dl_EG_sq(:) = 0.d0
    dl_EG(:) = 0.d0
    parameters(:) = 0.d0

    call save_nparms

    call fetch_parameters(parameters)

    ! do iteration
    do iter = 1, nopt_iter
        write (6, '(/,''DL Optimization iteration'',i5,'' of'',i5)') iter, nopt_iter

        call qmc

        write (6, '(/,''Completed sampling'')')

        call dl_more(iter, nparm, dl_momentum, dl_EG_sq, dl_EG, deltap, parameters)

        ! historically, we input -deltap in compute_parameters,
        ! so we multiply actual deltap by -1
        call dscal(nparm, -1.d0, deltap, 1)

        call test_solution_parm(nparm, deltap, dparm_norm, dparm_norm_min, sr_adiag, iflag)
        write (6, '(''Norm of parm variation '',g12.5)') dparm_norm
        if (iflag .ne. 0) then
            write (6, '(''Warning: dparm_norm>1'')')
            stop
        endif

        call compute_parameters(deltap, iflag, 1)
        call write_wf(1, iter)
        call save_wf

        if (iter .ge. 2) then
            denergy = energy(1) - energy_sav
            denergy_err = sqrt(energy_err(1)**2 + energy_err_sav**2)
            nblk = nblk*1.2
            nblk = min(nblk, nblk_max)
        endif
        write (6, '(''nblk = '',i6)') nblk

        energy_sav = energy(1)
        energy_err_sav = energy_err(1)
    enddo
    ! enddo iteration

    write (6, '(/,''Check last iteration'')')

    ioptjas = 0
    ioptorb = 0
    ioptci = 0

    call set_nparms_tot

    call qmc

    call write_wf(1, -1)

    deallocate (deltap)
    deallocate (dl_momentum)
    deallocate (dl_EG_sq)
    deallocate (dl_EG)
    deallocate (parameters)

    return
end

!=================================================
