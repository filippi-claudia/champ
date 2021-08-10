      subroutine store_diag_hs(nparm_p1,hii,sii)

      use sr_mod, only: MPARM
      use optwf_contrl, only: ioptjas, ioptorb, nparm

      use sr_index, only: jelo, jelo2, jelohfj
      use sr_mat_n, only: jefj, jfj, jhfj
      use sr_mat_n, only: obs_tot
      use mpi

      use precision_kinds, only: dp
      implicit none

      integer :: i, ish, jfhfj, jfifj, jwtg
      integer :: n_obs, nparm_p1

      real(dp), dimension(MPARM) :: hii
      real(dp), dimension(MPARM) :: sii


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
