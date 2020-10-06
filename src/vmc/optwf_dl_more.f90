subroutine dl_more(iter, nparm, dl_momentum, dl_EG_sq, dl_EG, deltap, parameters)

    use mpi
    use precision_kinds, only: dp
    use mpiconf, only: idtask

    implicit real*8(a - h, o - z)

    ! in/out variable
    integer, intent(in) :: iter, nparm
    real(dp), dimension(nparm), intent(in) :: dl_momentum
    real(dp), dimension(nparm), intent(in) :: dl_EG_sq
    real(dp), dimension(nparm), intent(in) :: dl_EG
    real(dp), dimension(nparm), intent(in) :: parameters
    real(dp), dimension(nparm), intent(inout) :: deltap

    character(20) dl_alg

    call p2gtfd('optwf:sr_adiag', sr_adiag, 0.01, 1)
    call p2gtfd('optwf:sr_tau', sr_tau, 0.02, 1)

    call p2gtfd('optwf:dl_mom', dl_mom, 0.0, 1)
    call p2gtid('optwf:idl_flag', idl_flag, 0, 1)
    call p2gtad('optwf:dl_alg', dl_alg, 'nag', 1)

    ! we only need h_sr = - grad_parm E
    call sr_hs(nparm, sr_adiag)

    if (idtask .eq. 0) then
        call dl_iter(iter, nparm, dl_alg, dl_mom, sr_tau, dl_momentum, dl_EG_sq, dl_EG, deltap, parameters)
    endif

    call MPI_BCAST(deltap, nparm, MPI_REAL8, 0, MPI_COMM_WORLD, ier)

    return
end

! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
