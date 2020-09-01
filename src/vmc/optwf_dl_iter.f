      subroutine dl_iter(iter,nparm,dl_alg,dl_mom,sr_tau,dl_momentum,dl_EG_sq,dl_EG,deltap,parameters)
      use sr_mat_n, only: elocal, h_sr, jefj, jfj, jhfj, nconf_n, obs, s_diag, s_ii_inv, sr_ho,
     &sr_o, wtg
      implicit real*8(a-h,o-z)


      include 'vmc.h'
      include 'force.h'
      include 'mstates.h'
      include 'sr.h'

      character*20 dl_alg
      real*8 dl_EG_corr, dl_EG_sq_corr, dl_momentum_prev, parm_old, dl_EG_old
      real*8 damp


      dimension deltap(*),dl_momentum(*),dl_EG(*),dl_EG_sq(*),parameters(*)
c Damping parameter for Nesterov gradient descent
      damp = 10.d0
      if(dl_alg.eq.'mom') then
        do i=1,nparm
          dl_momentum(i) = dl_mom*dl_momentum(i) + sr_tau*h_sr(i)
          deltap(i) = dl_momentum(i)
          parameters(i) = parameters(i) + deltap(i)
        enddo
       elseif(dl_alg.eq.'nag') then
        do i=1,nparm
          dl_momentum_prev = dl_momentum(i)
          dl_momentum(i) = dl_mom*dl_momentum(i) - sr_tau*h_sr(i)
          deltap(i) = -(dl_mom*dl_momentum_prev + (1 + dl_mom)*dl_momentum(i))
          parameters(i) = parameters(i) + deltap(i)
        enddo
      elseif(dl_alg.eq.'rmsprop') then
c Actually an altered version of rmsprop that uses nesterov momentum as well
c magic numbers: gamma = 0.9
        do i=1,nparm
            dl_EG_sq(i) = 0.9*dl_EG_sq(i) + 0.1*(-h_sr(i))**2
            parameters(i) = parameters(i) + deltap(i) 
            parm_old = parameters(i)
            dl_momentum_prev = dl_momentum(i)
            dl_momentum(i) = parameters(i) + sr_tau*h_sr(i)/sqrt(dl_EG_sq(i) + 10.d0**(-8.d0))
c To avoid declaring more arrays, use dl_EG for \lambda, dl_EG_sq = gamma
c Better solution needed (custom types for each iterator a la Fortran 2003 or C++ classes?)
            dl_EG_old = dl_EG(i)
            dl_EG(i) = 0.5 + 0.5*sqrt(1 + 4*dl_EG(i)**2)
            dl_EG_sq(i) = (1 - dl_EG_old)/dl_EG(i)
            dl_EG_sq_corr = dl_EG_sq(i)*exp(-(iter-1)/damp)
            parameters(i) = (1 - v_corr)*dl_momentum(i) +v_corr*dl_momentum_prev
            deltap(i) = parameters(i) - parm_old
        enddo
      elseif(dl_alg.eq.'adam') then
C Magic numbers: beta1 = 0.9, beta2 = 0.999
        do i=1,nparm
          dl_EG(i) = 0.9*dl_EG(i) + 0.1*(-h_sr(i))
          dl_EG_sq(i) = 0.999*dl_EG_sq(i) + 0.001*(-h_sr(i))**2
          dl_EG_corr = dl_EG(i)/(1-0.9**iter)
          dl_EG_sq_corr = dl_EG_sq(i)/(1-0.999**iter)
          deltap(i) = -sr_tau*dl_EG_corr/(sqrt(dl_EG_sq_corr)+10.d0**(-8.d0))
          parameters(i) = parameters(i) + deltap(i)
        enddo 
      elseif(dl_alg.eq.'cnag') then
        do i=1,nparm
          parm_old = parameters(i)
          dl_momentum_prev = dl_momentum(i)
          dl_momentum(i) = parameters(i) + sr_tau*h_sr(i)
c To avoid declaring more arrays, use dl_EG for \lambda, dl_EG_sq = gamma
          dl_EG_old = dl_EG(i)
          dl_EG(i) = 0.5 + 0.5*sqrt(1 + 4*dl_EG(i)**2)
          dl_EG_sq(i) = (1 - dl_EG_old)/dl_EG(i)
          parameters(i) = (1 - dl_EG_sq(i))*dl_momentum(i) + dl_EG_sq(i)*dl_momentum_prev
          deltap(i) = parameters(i) - parm_old
        enddo
      endif

      return
      end
