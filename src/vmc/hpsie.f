      module hpsie
      contains
      subroutine psie(iel,coord,psid,psij,ipass,iflag)
c Written by Claudia Filippi by modifying hpsi

      use contrl_file, only: ounit
      use csfs,    only: nstates
      use determinante_mod, only: determinante
      use determinante_psit_mod, only: determinante_psit
      use distance_mod, only: r_en,rvec_en
      use distances_mod, only: distances
      use error,   only: fatal_error
      use estpsi,  only: apsi,aref
      use jastrow, only: ianalyt_lap
      use jastrowe_mod, only: jastrowe
      use mstates_mod, only: MSTATES
      use multideterminante_mod, only: multideterminante
      use multiple_geo, only: iwf,iwftype
      use multislatern, only: detn
      use precision_kinds, only: dp
      use slater,  only: kref
      use system,  only: nelec
      use velocity_jastrow, only: vjn

      implicit none

      integer :: iel, iflag, ipass, istate
      real(dp) :: apsi_now, aref_now, check_apsi, check_apsi_min, check_dref
      real(dp) :: d2j, psij
      real(dp), dimension(3, nelec) :: coord
      real(dp), dimension(MSTATES) :: psid

c Calculates wave function

      ! dimension coord(3,*),psid(*)

      iwf=iwftype(1)


      call distances(iel,coord)


      if(ianalyt_lap.eq.1) then
        call jastrowe(iel,coord,vjn,d2j,psij,iflag)
       else
        call fatal_error('HPSIE: numerical one-electron move not implemented')
      endif

c compute all determinants

      call determinante(iel,coord,rvec_en,r_en,iflag)

      if(detn(kref).eq.0.d0) then
        do istate=1,nstates
          psid(istate)=0.d0
        enddo
        return
      endif


      call multideterminante(iel)

c combine determinantal quantities to obtain trial wave function

      do istate=1,nstates
        call determinante_psit(iel,psid(istate),istate)
      enddo


      if(ipass.gt.2) then

        check_apsi_min=1.d+99
        do istate=1,nstates
          apsi_now=apsi(istate)/(ipass-1)
          check_apsi=abs(psid(istate))/apsi_now

         check_apsi_min=min(check_apsi,check_apsi_min)
        enddo

        aref_now=aref/(ipass-1)
        check_dref=abs(detn(kref))/aref_now

      endif

      return
      end
      end module 
