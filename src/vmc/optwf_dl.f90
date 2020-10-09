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

module optwf_dl_mod

    use precision_kinds, only: dp
    implicit None

    type OptWFDLData
        integer :: nopt_iter, nblk_max
        real(dp) ::  energy_tol
        real(dp) :: dparm_norm_min
        real(dp) :: sr_tau, sr_adiag, sr_eps
        real(dp) :: dl_mom
        integer  :: idl_flag
        character(20) :: dl_alg

        real(dp), dimension(:), allocatable :: deltap
        real(dp), dimension(:), allocatable :: dl_momentum
        real(dp), dimension(:), allocatable :: dl_EG_sq
        real(dp), dimension(:), allocatable :: dl_EG
        real(dp), dimension(:), allocatable :: parameters
    end type OptWFDLData

    type(OptWFDLData)           :: opt_data
    private
    public :: optwf_dl
    save

contains

    subroutine optwf_dl()

        use precision_kinds, only: dp
        use sr_mod, only: MPARM
        use optwf_contrl, only: ioptci, ioptjas, ioptorb, nparm
        use optwf_corsam, only: energy, energy_err, force
        use contrl, only: nblk
        use method_opt, only: method

        implicit None

        integer :: iter, iflag
        real(dp) :: dparm_norm
        real(dp) :: denergy, energy_sav, denergy_err, energy_err_sav

        call read_input(opt_data)

        if (method .ne. 'sr_n' .or. opt_data%idl_flag .eq. 0) return

        write (6, '(''Started dl optimization'')')

        call set_nparms_tot

        if (nparm .gt. MPARM) call fatal_error('SR_OPTWF: nparmtot gt MPARM')

        call init_arrays(opt_data)

        call save_nparms()

        call fetch_parameters(opt_data%parameters)

        ! do iteration
        do iter = 1, opt_data%nopt_iter
            write (6, '(/,''DL Optimization iteration'',i5,'' of'',i5)') iter, opt_data%nopt_iter

            call qmc

            write (6, '(/,''Completed sampling'')')

            call optimization_step(iter, nparm, opt_data)

            ! historically, we input -deltap in compute_parameters,
            ! so we multiply actual deltap by -1
            call dscal(nparm, -1.d0, opt_data%deltap, 1)

            call test_solution_parm(nparm, opt_data%deltap, dparm_norm, opt_data%dparm_norm_min, opt_data%sr_adiag, iflag)
            write (6, '(''Norm of parm variation '',g12.5)') dparm_norm
            if (iflag .ne. 0) then
                write (6, '(''Warning: dparm_norm>1'')')
                stop
            endif

            call compute_parameters(opt_data%deltap, iflag, 1)
            call write_wf(1, iter)
            call save_wf

            if (iter .ge. 2) then
                denergy = energy(1) - energy_sav
                denergy_err = sqrt(energy_err(1)**2 + energy_err_sav**2)
                nblk = nblk*1.2
                nblk = min(nblk, opt_data%nblk_max)
            endif
            write (6, '(''nblk = '',i6)') nblk

            energy_sav = energy(1)
            energy_err_sav = energy_err(1)
        enddo
        ! enddo iteration

        write (6, '(/,''Check last iteration'')')

        ! Why ??!!
        ioptjas = 0
        ioptorb = 0
        ioptci = 0

        call set_nparms_tot

        call qmc

        call write_wf(1, -1)

        call deallocate_arrays(opt_data)

        return
    end subroutine optwf_dl

    subroutine init_arrays(options)
        !> Allocate and initialize to 0 all arrays
        use optwf_contrl, only: nparm
        type(OptWFDLData), INTENT(INOUT) :: options

        !> allocate
        allocate (options%deltap(nparm))
        allocate (options%dl_momentum(nparm))  !< 'momentum' variables
        allocate (options%dl_EG_sq(nparm))     !< moving average of past squared gradients
        allocate (options%dl_EG(nparm))        !< moving average of past gradients
        allocate (options%parameters(nparm))   !< vector of wave function parameters

        !> init
        options%dl_momentum(:) = 0.d0
        options%dl_EG_sq(:) = 0.d0
        options%dl_EG(:) = 0.d0
        options%parameters(:) = 0.d0

    end subroutine init_arrays

    subroutine deallocate_arrays(options)
        !> Deallocate arrays
        type(OptWFDLData), INTENT(INOUT) :: options

        deallocate (options%deltap)
        deallocate (options%dl_momentum)
        deallocate (options%dl_EG_sq)
        deallocate (options%dl_EG)
        deallocate (options%parameters)
    end subroutine deallocate_arrays

    subroutine read_input(options)
        !> Read the inputs
        use contrl, only: nblk
        type(OptWFDLData), INTENT(INOUT) :: options

        call p2gtid('optwf:nopt_iter', options%nopt_iter, 6, 1)
        call p2gtid('optwf:nblk_max', options%nblk_max, nblk, 1)
        call p2gtfd('optwf:energy_tol', options%energy_tol, 1.d-3, 1)

        call p2gtfd('optwf:dparm_norm_min', options%dparm_norm_min, 1.0d0, 1)
        write (6, '(''Starting dparm_norm_min'',g12.4)') options%dparm_norm_min

        call p2gtfd('optwf:sr_tau', options%sr_tau, 0.02, 1)
        call p2gtfd('optwf:sr_adiag', options%sr_adiag, 0.01, 1)
        call p2gtfd('optwf:sr_eps', options%sr_eps, 0.001, 1)
        call p2gtfd('optwf:dl_mom', options%dl_mom, 0.0, 1)
        call p2gtid('optwf:idl_flag', options%idl_flag, 0, 1)
        call p2gtad('optwf:dl_alg', options%dl_alg, 'nag', 1)

        write (6, '(/,''SR adiag: '',f10.5)') options%sr_adiag
        write (6, '(''SR tau:   '',f10.5)') options%sr_tau
        write (6, '(''SR eps:   '',f10.5)') options%sr_eps
        write (6, '(''DL flag:   '',I10)') options%idl_flag

    end subroutine read_input

    subroutine optimization_step(iter, nparm, opt)
        !> do 1 optimization step
        use mpi
        use precision_kinds, only: dp
        use mpiconf, only: idtask
        use optwf_sr_mod, only: sr_hs

        ! in/out variable
        integer, intent(in) :: iter, nparm
        type(OptWFDLData), intent(inout) :: opt
        integer :: ierr

        ! we only need h_sr = - grad_parm E
        call sr_hs(nparm, opt%sr_adiag)

        if (idtask .eq. 0) then
            call one_iter(iter, nparm, opt)
        endif

        call MPI_BCAST(opt%deltap, nparm, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)

        return
    end subroutine optimization_step

    subroutine one_iter(iter, nparm, opt)
        !> Individual routine to optimize the parameters
        use precision_kinds, only: dp
        use sr_mat_n, only: h_sr

        integer, intent(in) :: iter
        integer, intent(in) :: nparm
        type(OptWFDLData), intent(inout) :: opt

        integer :: i
        real(dp) :: dl_EG_corr, dl_EG_sq_corr, dl_momentum_prev, parm_old, dl_EG_old
        real(dp) :: v_corr
        real(dp) :: damp

        ! Damping parameter for Nesterov gradient descent
        damp = 10.d0
        select case (opt%dl_alg)
        case ('mom')
            do i = 1, nparm
                opt%dl_momentum(i) = opt%dl_mom*opt%dl_momentum(i) + opt%sr_tau*h_sr(i)
                opt%deltap(i) = opt%dl_momentum(i)
                opt%parameters(i) = opt%parameters(i) + opt%deltap(i)
            enddo
        case ('nag')
            do i = 1, nparm
                dl_momentum_prev = opt%dl_momentum(i)
                opt%dl_momentum(i) = opt%dl_mom*opt%dl_momentum(i) - opt%sr_tau*h_sr(i)
                opt%deltap(i) = -(opt%dl_mom*dl_momentum_prev + (1 + opt%dl_mom)*opt%dl_momentum(i))
                opt%parameters(i) = opt%parameters(i) + opt%deltap(i)
            enddo
        case ('rmsprop')
            ! Actually an altered version of rmsprop that uses nesterov momentum as well
            ! magic numbers: gamma = 0.9
            do i = 1, nparm
                opt%dl_EG_sq(i) = 0.9*opt%dl_EG_sq(i) + 0.1*(-h_sr(i))**2
                opt%parameters(i) = opt%parameters(i) + opt%deltap(i)
                parm_old = opt%parameters(i)
                dl_momentum_prev = opt%dl_momentum(i)
                opt%dl_momentum(i) = opt%parameters(i) + opt%sr_tau*h_sr(i)/sqrt(opt%dl_EG_sq(i) + 10.d0**(-8.d0))

                ! To avoid declaring more arrays, use dl_EG for \lambda, dl_EG_sq = gamma
                ! Better solution needed (custom types for each iterator a la Fortran 2003 or C++ classes?)
                dl_EG_old = opt%dl_EG(i)
                opt%dl_EG(i) = 0.5 + 0.5*sqrt(1 + 4*opt%dl_EG(i)**2)
                opt%dl_EG_sq(i) = (1 - dl_EG_old)/opt%dl_EG(i)
                dl_EG_sq_corr = opt%dl_EG_sq(i)*exp(-(iter - 1)/damp)
                opt%parameters(i) = (1 - v_corr)*opt%dl_momentum(i) + v_corr*dl_momentum_prev
                opt%deltap(i) = opt%parameters(i) - parm_old
            enddo
        case ('adam')
            ! Magic numbers: beta1 = 0.9, beta2 = 0.999
            do i = 1, nparm
                opt%dl_EG(i) = 0.9*opt%dl_EG(i) + 0.1*(-h_sr(i))
                opt%dl_EG_sq(i) = 0.999*opt%dl_EG_sq(i) + 0.001*(-h_sr(i))**2
                dl_EG_corr = opt%dl_EG(i)/(1 - 0.9**iter)
                dl_EG_sq_corr = opt%dl_EG_sq(i)/(1 - 0.999**iter)
                opt%deltap(i) = -opt%sr_tau*dl_EG_corr/(sqrt(dl_EG_sq_corr) + 10.d0**(-8.d0))
                opt%parameters(i) = opt%parameters(i) + opt%deltap(i)
            enddo
        case ('cnag')
            do i = 1, nparm
                parm_old = opt%parameters(i)
                dl_momentum_prev = opt%dl_momentum(i)
                opt%dl_momentum(i) = opt%parameters(i) + opt%sr_tau*h_sr(i)
                ! To avoid declaring more arrays, use dl_EG for \lambda, dl_EG_sq = gamma
                dl_EG_old = opt%dl_EG(i)
                opt%dl_EG(i) = 0.5 + 0.5*sqrt(1 + 4*opt%dl_EG(i)**2)
                opt%dl_EG_sq(i) = (1 - dl_EG_old)/opt%dl_EG(i)
                opt%parameters(i) = (1 - opt%dl_EG_sq(i))*opt%dl_momentum(i) + opt%dl_EG_sq(i)*dl_momentum_prev
                opt%deltap(i) = opt%parameters(i) - parm_old
            enddo
        end select

        return
    end subroutine one_iter

end module optwf_dl_mod
!=================================================
