      subroutine lbfgs_iter(iter, nparm, deltap, parameters, step_size)
      !use lbfgs_wrapper, only: lbfgs_iteration
      use olbfgs, only: olbfgs_iteration, update_hessian
      implicit real*8 (a-h,o-z)

      include 'vmc.h'
      include 'force.h'
      include 'mstates.h'
      include 'sr.h'

      common /sr_mat_n/ sr_o(MPARM,MCONF),sr_ho(MPARM,MCONF),obs(MOBS,MSTATES),s_diag(MPARM,MSTATES),s_ii_inv(MPARM),h_sr(MPARM),wtg(MCONF,MSTATES),elocal(MCONF,MSTATES),jfj,jefj,jhfj,nconf

      dimension deltap(*), parameters(*)

      real(kind=8), dimension(1), allocatable :: parameters_old(:)
      real(kind=8), dimension(1), allocatable :: parms_lbfgs(:)

      allocate(parameters_old(nparm))
      allocate(parms_lbfgs(nparm))

      parms_lbfgs = parameters(1:nparm)
      parameters_old = parms_lbfgs

      call update_hessian(parms_lbfgs, -h_sr)
      ! TODO make step size configurable
      print *, 'Doing o-lbfgs iteration'
      print *, step_size
      call olbfgs_iteration(parms_lbfgs, -h_sr, step_size, iter)

      deltap(1:nparm) = parms_lbfgs - parameters_old
      parameters(1:MPARM) = parameters(1:MPARM) + deltap(1:MPARM)

      deallocate(parms_lbfgs)
      deallocate(parameters_old)

      return
      end
