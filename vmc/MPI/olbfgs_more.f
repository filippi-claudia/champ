      subroutine olbfgs_more(iter, nparm, deltap, parameters)

      implicit real*8 (a-h,o-z)
      character*20 dl_alg

      include 'vmc.h'
      include 'force.h'
      include 'mstates.h'
      include 'sr.h'

      include 'mpif.h'

      common /sr_mat_n/ sr_o(MPARM,MCONF),sr_ho(MPARM,MCONF),obs(MOBS,MSTATES),s_diag(MPARM,MSTATES)
     &,s_ii_inv(MPARM),h_sr(MPARM),wtg(MCONF,MSTATES),elocal(MCONF,MSTATES),jfj,jefj,jhfj,nconf

      dimension deltap(*), parameters(*)

c we only need h_sr = - grad_parm E
      call sr_hs(nparm,sr_adiag)

      if(idtask.eq.0) then 
        call olbfgs_iter(iter, nparm, deltap, parameters)
      end if

c Update parameter changes
      deltap(1:nparm) = parms_lbfgs - parameters_old
      parameters(1:nparm) = parameters(1:nparm) + deltap(1:nparm)

      call MPI_BCAST(deltap,nparm,MPI_REAL8,0,MPI_COMM_WORLD,ier)

      return
      end
