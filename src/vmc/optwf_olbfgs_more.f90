module olbfgs_more_mod
contains
      subroutine olbfgs_more(iter, nparm, deltap, parameters)

      use mpi
      use mpiconf, only: idtask
      use olbfgs,  only: olbfgs_iteration,update_hessian
      use optwf_control, only: sr_adiag,sr_tau
      use optwf_sr_mod, only: sr_hs
      use precision_kinds, only: dp
      use sr_mat_n, only: h_sr

      implicit none

      integer :: ier, iter, nparm

      real(dp), dimension(*) :: deltap
      real(dp), dimension(*) :: parameters
      character(len=20) dl_alg
      real(dp), dimension(1), allocatable :: parameters_old(:)
      real(dp), dimension(1), allocatable :: parms_lbfgs(:)

      allocate(parameters_old(nparm))
      allocate(parms_lbfgs(nparm))

      parms_lbfgs = parameters(1:nparm)
      parameters_old = parms_lbfgs


! we only need h_sr = - grad_parm E
      call sr_hs(nparm,sr_adiag)

      if(idtask.eq.0) then
! update stored Hessian approximation
        call update_hessian(parms_lbfgs, -h_sr(:,1))

! perform actual oLBFGS iteration
        call olbfgs_iteration(parms_lbfgs, -h_sr(:,1), sr_tau, iter)

        deltap(1:nparm) = parms_lbfgs - parameters_old
        parameters(1:nparm) = parameters(1:nparm) + deltap(1:nparm)
      end if

! Update parameter changes

      call MPI_BCAST(deltap,nparm,MPI_REAL8,0,MPI_COMM_WORLD,ier)

      deallocate(parameters_old)
      deallocate(parms_lbfgs)

      return
      end
end module
