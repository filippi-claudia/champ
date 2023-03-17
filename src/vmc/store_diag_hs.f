      module store_diag_hs_mod
      contains
      subroutine store_diag_hs(nparm_p1,hii,sii)

      use contrl_file, only: errunit,ounit
      use mpi
      use optwf_control, only: ioptjas,ioptorb,nparm
      use precision_kinds, only: dp
      use sr_mat_n, only: obs_tot
      use sr_mod,  only: mparm

      implicit none

      integer :: i, ish, jfhfj, jfifj, jwtg
      integer :: jelo,jelo2,jelohfj
      integer :: jefj,jfj,jhfj
      integer :: n_obs, nparm_p1

      real(dp), dimension(:) :: hii
      real(dp), dimension(:) :: sii
      real(dp) :: aux0, aux1, aux2


      write(ounit,*) 'nparm_p1,nparm',nparm_p1,nparm

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
        aux0=obs_tot(jfj-1+i,1)
        aux1=obs_tot(jefj-1+i,1)
        aux2=obs_tot(jhfj-1+i,1)
        hii(i+ish)=aux0*obs_tot(jelo,1)*obs_tot(jfj+i-1,1)
     &     -0.5d0*(aux0*(obs_tot(jefj+i-1,1)+obs_tot(jhfj+i-1,1))
     &                  +obs_tot(jfj+i-1,1)*(aux1+aux2))
        sii( i+ish)=obs_tot(jfifj-1+i,1) -
     &              obs_tot(jfj-1+i,1)**2
      enddo

 

      return
      end
      end module
