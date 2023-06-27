module nonloc_pot_mod
contains
      subroutine nonloc_pot(x,rshift,rvec_en,r_en,pe,vpsp_det,dvpsp_dj,t_vpsp,i_vpsp,ifr)
! Written by Claudia Filippi; modified by Cyrus Umrigar
! Calculates the local and nonlocal components of the pseudopotential
! Calculates non-local potential derivatives
! pe_en(loc) is computed in distances and pe_en(nonloc) here in nonloc_pot if nloc !=0 and iperiodic!=0.
      use contrl_per, only: iperiodic
      use nonloc_mod, only: nonloc
      use precision_kinds, only: dp
      use pseudo,  only: lpot,nloc,vps
      use pseudo_mod, only: MPS_QUAD
      use readps_gauss, only: getvps_gauss
      use system,  only: iwctype,ncent,ncent_tot,nelec
      use error,   only: fatal_error
      use vmc_mod, only: nbjx
      use optwf_parms, only: nparmj

      implicit none

      integer :: i, i1, i2, i_vpsp, ic
      integer :: ifr
      real(dp) :: pe
      real(dp), dimension(3, *) :: x
      real(dp), dimension(3, nelec, ncent_tot) :: rshift
      real(dp), dimension(3, nelec, ncent_tot) :: rvec_en
      real(dp), dimension(nelec, ncent_tot) :: r_en
      real(dp), dimension(2, nbjx) :: vpsp_det
      real(dp), dimension(nparmj, nbjx) :: dvpsp_dj
      real(dp), dimension(ncent_tot, MPS_QUAD, *) :: t_vpsp

      if(i_vpsp.gt.0)then
        i1=i_vpsp
        i2=i_vpsp
       else
        i1=1
        i2=nelec
      endif
      if(nloc.eq.4) then
         do i=i1,i2
            call getvps_gauss(rvec_en,r_en,i)
         enddo
      else
         call fatal_error('nonloc different to 4 is not supported')
      endif


! local component (highest angular momentum)
      if(iperiodic.eq.0) then
        do ic=1,ncent
          do i=i1,i2
            pe=pe+vps(i,ic,lpot(iwctype(ic)))
          enddo
        enddo
      endif
! non-local component (division by the Jastrow already in nonloc)
      call nonloc(x,rshift,rvec_en,r_en,vpsp_det,dvpsp_dj,t_vpsp,i_vpsp)

      return
      end
end module 
