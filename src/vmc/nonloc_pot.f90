      module nonloc_pot_mod
      contains
      subroutine nonloc_pot(x,rvec_en,r_en,pe,vpsp_det,dvpsp_dj,t_vpsp,i_vpsp,ifr)
! Written by Claudia Filippi; modified by Cyrus Umrigar
! Calculates the local and nonlocal components of the pseudopotential
! Calculates non-local potential derivatives
! pe_en(loc) is computed in distances and pe_en(nonloc) here in nonloc_pot if nloc !=0 and iperiodic!=0.
      use contrl_per, only: iperiodic
      use contrldmc, only: icut_e
      use fragments, only: eloc_i, elocfrag, ifragelec, ifragcent, nfrag
      use error,   only: fatal_error
      use ewald, only: cos_n_sum, sin_n_sum
      use ewald_breakup, only: pot_en_ewald
      use m_force_analytic, only: iforce_analy
      use nonloc_mod, only: nonloc
      use optwf_parms, only: nparmj
      use precision_kinds, only: dp
      use pseudo,  only: lpot,nloc,vps
      use pseudo_mod, only: MPS_QUAD
      use da_pseudo, only: da_vps
      use readps_gauss, only: getvps_gauss, gauss_pot
      use system,  only: iwctype,ncent,ncent_tot,nelec, znuc
      use vmc_mod, only: nbjx

      implicit none

      integer :: i, i1, i2, i_vpsp, ic,k
      integer :: ifr, ict
      real(dp) :: pe, pe_en, r, ri, ri2,vcoule,dvpot, tmp
      real(dp), dimension(3, *) :: x
      real(dp), dimension(3, nelec, ncent_tot) :: rvec_en
      real(dp), dimension(nelec, ncent_tot) :: r_en
      real(dp), dimension(2, nbjx) :: vpsp_det
      real(dp), dimension(nparmj, nbjx) :: dvpsp_dj
      real(dp), dimension(ncent_tot, MPS_QUAD, *) :: t_vpsp
      real(dp), dimension(3,ncent) :: da_pe_en

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
            tmp = vps(i,ic,lpot(iwctype(ic)))
            pe=pe+tmp
            if (icut_e.lt.0) then
              eloc_i(i)=eloc_i(i)+tmp
            endif
            if (nfrag.gt.1) then
              elocfrag(ifragcent(ic)) = elocfrag(ifragcent(ic)) + 0.5d0 * tmp
              elocfrag(ifragelec(i)) = elocfrag(ifragelec(i)) + 0.5d0 * tmp
            endif
          enddo
        enddo

      else

! this add coulumb pe_en contribution from PBC/Ewald split Natoli-Ceperley algorithm

         call pot_en_ewald(x,pe_en,cos_n_sum(1,ifr),sin_n_sum(1,ifr))
         pe=pe+pe_en

!     Add and fix remaining local component from gaussian (BFD) pseudo
!     It can be done in two ways
!     this to be improved and updated later

         do ic=1,ncent
            ict=iwctype(ic)
            do i=i1,i2
               r=max(1.0d-10,r_en(i,ic))
               ri=1.d0/r
!     1- subtract Coulomb term obc in local component vps(:,:lpot(ict)
!     current simplest one /subtract local coulomb term which was added at the beginning at pe for periodic
!               pe=pe+vps(i,ic,lpot(ict))+(znuc(ict)*ri)

!     2- get local component without coulomb regarding local was added from pot_loc with pot_en_ewald
!     implies recompute gauss pot need to be a better way
               call gauss_pot(r,lpot(ict),ict,vcoule,dvpot)
               pe=pe+vcoule
               if (iforce_analy.gt.0) then
                 do k=1,3
                   da_vps(k,i,ic,lpot(ict))=-dvpot*rvec_en(k,i,ic)*ri
                 enddo
               endif
            enddo
         enddo
      endif


!     non-local component (division by the Jastrow already in nonloc)
      call nonloc(x,rvec_en,r_en,vpsp_det,dvpsp_dj,t_vpsp,i_vpsp)

      return
      end
end module
