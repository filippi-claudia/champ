      module pot
      contains
      subroutine pot_nn(cent,znuc,iwctype,ncent,pecent)
c Written by Cyrus Umrigar
c get nuclear potential energy
      use contrl_per , only:  iperiodic, ibasis
      use contrl_per, only: iperiodic, ibasis
      use contrl_per, only: iperiodic, ibasis
      use contrl_per, only: iperiodic
      USE da_pseudo , only:  da_pecent, da_vps
      use da_pseudo, only: da_pecent, da_vps
      use da_pseudo, only: da_pecent, da_vps
      use da_pseudo, only: da_pecent, da_vps
      use system, only: ncent_tot, nctype_tot
      use m_force_analytic, only: iforce_analy
      use da_pseudo, only: da_pecent, da_vps
      use precision_kinds, only: dp
      use pw_ewald, only: pot_nn_ewald

      implicit none

      integer :: i, j, j1, k, ncent
      integer, dimension(ncent_tot) :: iwctype
      real(dp) :: pecent, r, r2, ri, ri2
      real(dp) :: term, zij
      real(dp), dimension(nctype_tot) :: znuc
      real(dp), dimension(3,ncent_tot) :: cent

!
      if(iperiodic.eq.0) then
        pecent=0
        do i=2,ncent
          j1=i-1
          do j=1,j1
            r2=0
            do k=1,3
              r2=r2+(cent(k,i)-cent(k,j))**2
            enddo
            r=dsqrt(r2)
            pecent=pecent+znuc(iwctype(i))*znuc(iwctype(j))/r
          enddo
        enddo
       else
        call pot_nn_ewald
      endif

      if(iforce_analy.eq.0) return

      if(iperiodic.eq.0) then
        do i=1,ncent
          do k=1,3
            da_pecent(k,i)=0.d0
          enddo
        enddo
        do i=2,ncent
          j1=i-1
          do j=1,j1
            zij=znuc(iwctype(i))*znuc(iwctype(j))
            r2=0
            do k=1,3
              r2=r2+(cent(k,i)-cent(k,j))**2
            enddo
            r=dsqrt(r2)
            ri=1/r
            ri2=ri*ri
            do k=1,3
              term=zij*(cent(k,i)-cent(k,j))*ri*ri2
              da_pecent(k,i)=da_pecent(k,i)-term
              da_pecent(k,j)=da_pecent(k,j)+term
            enddo
          enddo
        enddo
      endif

      return
      end
      end module
