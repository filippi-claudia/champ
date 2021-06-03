      subroutine nonloc_pot(x,rshift,rvec_en,r_en,pe,vpsp_det,dvpsp_dj,t_vpsp,i_vpsp,ifr)
c     Written by Claudia Filippi; modified by Cyrus Umrigar
c     Calculates the local and nonlocal components of the pseudopotential
c     Calculates non-local potential derivatives
c     pe_en(loc) is computed in distances and pe_en(nonloc) here in nonloc_pot if nloc !=0 and iperiodic!=0.
      use mstates_mod, only: MSTATES
      use optjas, only: MPARMJ
      use pseudo_mod, only: MPS_QUAD
      use vmc_mod, only: MELEC, MCENT
      use atom, only: iwctype, ncent
      use const, only: nelec
      use contrl_per, only: iperiodic
      use pseudo, only: lpot, nloc, vps

      implicit real*8(a-h,o-z)

      dimension x(3,*),rshift(3,MELEC,MCENT),rvec_en(3,MELEC,MCENT),r_en(MELEC,MCENT)
     &     ,dvpsp_dj(MPARMJ,MSTATES),t_vpsp(MCENT,MPS_QUAD,*)
      dimension vpsp_det(2,MSTATES)

      if(i_vpsp.gt.0)then
         i1=i_vpsp
         i2=i_vpsp
      else
         i1=1
         i2=nelec
      endif

      do i=i1,i2
         if(nloc.eq.1) then
            call getvps(r_en,i)
         elseif(nloc.eq.2.or.nloc.eq.3) then
            call getvps_tm(r_en,i)
         elseif(nloc.eq.4) then
            call getvps_gauss(rvec_en,r_en,i)
         elseif(nloc.eq.5) then
            call getvps_champ(r_en,i)
         endif
      enddo

c     local component (highest angular momentum)

      if(iperiodic.eq.0) then
         do ic=1,ncent
            do i=i1,i2
               pe=pe+vps(i,ic,lpot(iwctype(ic)))
            enddo
         enddo
      endif

c     non-local component (division by the Jastrow already in nonloc)

      call nonloc(x,rshift,rvec_en,r_en,vpsp_det,dvpsp_dj,t_vpsp,i_vpsp)

      end subroutine
