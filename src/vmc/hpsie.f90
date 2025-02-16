module hpsie
contains
      subroutine psie(iel,coord,psid,psij,ipass,iflag)
! Written by Claudia Filippi by modifying hpsi

      use config,  only: xold
      use contrl_file, only: ounit
      use csfs, only: nstates
      use determinante_mod, only: determinante
      use determinante_psit_mod, only: determinante_psit
      use distance_mod, only: r_en, rvec_en
      use distances_mod, only: distances
      use error, only: fatal_error
      use estpsi, only: apsi, aref
      use jastrowe_mod, only: jastrowe
      use mstates_mod, only: MSTATES
      use multideterminante_mod, only: multideterminante
      use multiple_geo, only: iwf, iwftype
      use multislatern, only: detn
      use precision_kinds, only: dp
      use slater, only: kref
      use system, only: nelec
      use velocity_jastrow, only: vjn
      use vmc_mod, only: nwftypejas, stoo

      implicit none

      integer :: i,k
      integer :: iel, iflag, ipass, istate, icheck
      real(dp) :: apsi_now, aref_now, check_apsi, check_apsi_min, check_dref
      real(dp), dimension(3, nelec) :: coord
      real(dp), dimension(MSTATES) :: psid
      real(dp), dimension(nwftypejas) :: psij
      real(dp), dimension(nwftypejas) :: d2j

! Calculates wave function

      iwf=iwftype(1)

      call distances(iel,coord)

      call jastrowe(iel,coord,vjn,d2j,psij,iflag)

! compute all determinants

      call determinante(iel,coord,rvec_en,r_en,iflag)

      icheck=0
      do istate=1,nstates
        if(detn(kref,stoo(istate)).eq.0.d0) then
          psid(istate)=0.d0
          icheck=1
        endif
      enddo
      if(icheck.eq.1) return

      call multideterminante(iel)

! combine determinantal quantities to obtain trial wave function

      do istate=1,nstates
        call determinante_psit(iel,psid(istate),istate)
      enddo

      if(ipass.gt.2) then

        do istate=1,nstates
          check_apsi_min=1.d+99
          apsi_now=apsi(istate)/(ipass-1)
          check_apsi=abs(psid(istate))/apsi_now
          check_apsi_min=min(check_apsi,check_apsi_min)
          aref_now=aref(stoo(istate))/(ipass-1)
          check_dref=abs(detn(kref,stoo(istate)))/aref_now
        enddo

      endif

      return
      end
end module 
