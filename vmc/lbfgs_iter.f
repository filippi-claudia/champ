      subroutine lbfgs_iter(iter, nparm, deltap, parameters, energy, diag, workspace)
      use lbfgs_wrapper, only: lbfgs_iteration
      implicit real*8 (a-h,o-z)

      include 'vmc.h'
      include 'force.h'
      include 'mstates.h'
      include 'sr.h'

      common /sr_mat_n/ sr_o(MPARM,MCONF),sr_ho(MPARM,MCONF),obs(MOBS,MSTATES),s_diag(MPARM,MSTATES),s_ii_inv(MPARM),h_sr(MPARM),wtg(MCONF,MSTATES),elocal(MCONF,MSTATES),jfj,jefj,jhfj,nconf
      ! need energy for LBFGS iteration
      ! TODO figure out a cleaner way to get eold (as a subroutine 
      ! parameter is best)

      dimension deltap(*), parameters(MPARM), parameters_old(MPARM), diag(MPARM), workspace(MPARM*11 + 10)

      ! save old parameters to compute deltap after lbfgs step
      parameters_old = parameters

      ! TODO: make num_history configurable
      call lbfgs_iteration(energy, parameters, -h_sr, nparm, 5, diag, workspace)

      do i=1,nparm
        deltap(i) = parameters(i) - parameters_old(i)
      end do

      return
      end
