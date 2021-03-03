      subroutine store_diag_hs(nparm_p1,hii,sii)

      use sr_mod, only: MPARM
      use optwf_contrl, only: ioptjas, ioptorb, nparm
      use mstates_mod, only: MSTATES
      use sr_index, only: jelo, jelo2, jelohfj
      use sr_mat_n, only: jefj, jfj, jhfj
      use sr_mat_n, only: obs_tot
      use mpi

      implicit real*8(a-h,o-z)

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
