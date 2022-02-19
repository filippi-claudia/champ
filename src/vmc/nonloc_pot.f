      module nonloc_pot_mod
      contains
      subroutine nonloc_pot(x,rshift,rvec_en,r_en,pe,vpsp_det,dvpsp_dj,t_vpsp,i_vpsp,ifr)
c Written by Claudia Filippi; modified by Cyrus Umrigar
c Calculates the local and nonlocal components of the pseudopotential
c Calculates non-local potential derivatives
c pe_en(loc) is computed in distances and pe_en(nonloc) here in nonloc_pot if nloc !=0 and iperiodic!=0.
      use pseudo_mod, only: MPS_QUAD
      use atom, only: iwctype, ncent, ncent_tot
      use const, only: nelec
      use contrl_per, only: iperiodic

      use pseudo, only: lpot, nloc, vps

      use precision_kinds, only: dp
      use readps_mod, only: getvps
      use readps_gauss, only: getvps_gauss
      use readps_tm_mod, only: getvps_tm
      use readps_champ_mod, only: getvps_champ
      use nonloc_mod, only: nonloc
      implicit none

      integer :: i, i1, i2, i_vpsp, ic
      integer :: ifr
      real(dp) :: pe
      real(dp), dimension(3, *) :: x
      real(dp), dimension(3, nelec, ncent_tot) :: rshift
      real(dp), dimension(3, nelec, ncent_tot) :: rvec_en
      real(dp), dimension(nelec, ncent_tot) :: r_en
      real(dp), dimension(*) :: vpsp_det
      real(dp), dimension(*) :: dvpsp_dj
      real(dp), dimension(ncent_tot, MPS_QUAD, *) :: t_vpsp



      if(i_vpsp.gt.0)then
        i1=i_vpsp
        i2=i_vpsp
       else
        i1=1
        i2=nelec
      endif
      do i=i1,i2
        if(nloc.eq.1) then
c         call getvps(r_en,i)
         elseif(nloc.eq.2.or.nloc.eq.3) then
          call getvps_tm(r_en,i)
         elseif(nloc.eq.4) then
          call getvps_gauss(rvec_en,r_en,i)
         elseif(nloc.eq.5) then
c         call getvps_champ(r_en,i)
        endif
      enddo
      
c local component (highest angular momentum)
      if(iperiodic.eq.0) then
        do ic=1,ncent
          do i=i1,i2
            pe=pe+vps(i,ic,lpot(iwctype(ic)))
          enddo
        enddo
      endif
      
c non-local component (division by the Jastrow already in nonloc)
      call nonloc(x,rshift,rvec_en,r_en,vpsp_det,dvpsp_dj,t_vpsp,i_vpsp)
      
      return
      end
      end module 
