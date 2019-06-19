      subroutine olbfgs_more(iter, nparm, deltap, parameters)
      use olbfgs, only: update_hessian, olbfgs_iteration

      implicit real*8 (a-h,o-z)
      character*20 dl_alg

      include 'vmc.h'
      include 'force.h'
      include 'mstates.h'
      include 'sr.h'

      common /sr_mat_n/ sr_o(MPARM,MCONF),sr_ho(MPARM,MCONF),obs(MOBS,MSTATES),s_diag(MPARM,MSTATES)
     &,s_ii_inv(MPARM),h_sr(MPARM),wtg(MCONF,MSTATES),elocal(MCONF,MSTATES),jfj,jefj,jhfj,nconf

      dimension deltap(*), parameters(*)

c we only need h_sr = - grad_parm E
      call sr_hs(nparm,sr_adiag)

      call olbfgs_iter(iter, nparm, deltap, parameters)

      return
      end
