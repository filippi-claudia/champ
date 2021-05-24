      subroutine olbfgs_more(iter, nparm, deltap, parameters)

      use olbfgs, only: update_hessian, olbfgs_iteration
      use sr_mat_n, only: h_sr
      use optwf_sr_mod, only: sr_hs
      use mpiconf, only: idtask
      use optwf_contrl, only: sr_tau , sr_adiag
      use mpi
      use precision_kinds, only: dp

      implicit none

      integer :: ier, iter, nparm

      real(dp), dimension(*) :: deltap
      real(dp), dimension(*) :: parameters
      character*20 dl_alg
      real(kind=8), dimension(1), allocatable :: parameters_old(:)
      real(kind=8), dimension(1), allocatable :: parms_lbfgs(:)

      allocate(parameters_old(nparm))
      allocate(parms_lbfgs(nparm))

      parms_lbfgs = parameters(1:nparm)


c we only need h_sr = - grad_parm E
      call sr_hs(nparm,sr_adiag)

      if(idtask.eq.0) then 
c update stored Hessian approximation
        call update_hessian(parms_lbfgs, -h_sr)

c perform actual oLBFGS iteration
        call olbfgs_iteration(parms_lbfgs, -h_sr, sr_tau, iter)

        deltap(1:nparm) = parms_lbfgs - parameters_old
      end if

c Update parameter changes

      call MPI_BCAST(deltap,nparm,MPI_REAL8,0,MPI_COMM_WORLD,ier)

      deallocate(parameters_old)
      deallocate(parms_lbfgs)

      return
      end
