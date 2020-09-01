      subroutine store_diag_hs(nparm_p1,hii,sii)

      use csfs, only: ccsf, cxdet, iadet, ibdet, icxdet, ncsf, nstates
      use mpiconf, only: idtask, nproc
      use optwf_contrl, only: ioptci, ioptjas, ioptorb, nparm
      use optwf_func, only: ifunc_omega, omega, omega_hes
      use sa_weights, only: iweight, nweight, weights
      use sr_index, only: jelo, jelo2, jelohfj
      use sr_mat_n, only: elocal, h_sr, jefj, jfj, jhfj, nconf_n, obs, s_diag, s_ii_inv, sr_ho,
     &sr_o, wtg, obs_tot

      implicit real*8(a-h,o-z)

      include 'mpif.h'
      include 'sr.h'
      include 'vmc.h'
      include 'force.h'
      include 'mstates.h'
      include 'optorb.h'

      dimension obs_wtg(MSTATES),obs_wtg_tot(MSTATES)
      dimension hii(MPARM),sii(MPARM)

      write(6,*) 'nparm_p1,nparm',nparm_p1,nparm

      jwtg=1
      jelo=2
      n_obs=2
      jfj=n_obs+1
      n_obs=n_obs+nparm
      jefj=n_obs+1
      n_obs=n_obs+nparm
      jfifj=n_obs+1
      n_obs=n_obs+nparm

      jhfj=n_obs+1
      n_obs=n_obs+nparm
      jfhfj=n_obs+1
      n_obs=n_obs+nparm

      jelo2=n_obs+1
      n_obs=n_obs+1
      jelohfj=n_obs+1
      n_obs=n_obs+nparm

      ish=0
      if(ioptorb+ioptjas.gt.0) then
        ish=1
        hii( 1)=obs_tot(jelo, 1)
        sii( 1)=1
      endif
      do i=1,nparm
        hii( i+ish)=obs_tot(jfhfj-1+i,1)
        sii( i+ish)=obs_tot(jfifj-1+i,1)
      enddo
      
      return
      end
