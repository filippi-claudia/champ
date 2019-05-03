      subroutine lbfgs_iter(iter, nparm, deltap, parameters, energy, diag, workspace, step_size)
      use lbfgs_wrapper, only: lbfgs_iteration
      implicit real*8 (a-h,o-z)

      include 'vmc.h'
      include 'force.h'
      include 'mstates.h'
      include 'sr.h'

      common /sr_mat_n/ sr_o(MPARM,MCONF),sr_ho(MPARM,MCONF),obs(MOBS,MSTATES),s_diag(MPARM,MSTATES),s_ii_inv(MPARM),h_sr(MPARM),wtg(MCONF,MSTATES),elocal(MCONF,MSTATES),jfj,jefj,jhfj,nconf

      dimension deltap(*), parameters(*), diag(*), workspace(*)

      real(kind=8), dimension(1), allocatable :: parameters_old(:)
      real(kind=8), dimension(1), allocatable :: parms_lbfgs(:)
      !dimension parameters_old(MPARM)
      !dimension parms_lbfgs(MPARM)

      allocate(parameters_old(nparm))
      allocate(parms_lbfgs(nparm))

      parms_lbfgs = parameters(1:nparm)
      parameters_old = parms_lbfgs

      ! TODO: make num_history configurable
      call lbfgs_iteration(energy, parms_lbfgs, -h_sr, nparm, 10, diag, workspace, step_size)

      deltap(1:nparm) = parms_lbfgs - parameters_old
      parameters(1:MPARM) = parameters(1:MPARM) + deltap(1:MPARM)

      deallocate(parms_lbfgs)
      deallocate(parameters_old)

      return
      end
