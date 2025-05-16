module pot
contains
      subroutine pot_nn(cent,znuc,iwctype,ncent,pecent,cos_n_sum,sin_n_sum)
! Written by Cyrus Umrigar
! get nuclear potential energy
      use fragments, only: potnnfrag, ifragcent, nfrag
      use contrl_per, only: ibasis,iperiodic
      use contrldmc, only: icut_e
      use da_pseudo, only: da_pecent,da_vps
      use m_force_analytic, only: iforce_analy
      use precision_kinds, only: dp
      use ewald_breakup, only: pot_nn_ewald
      use system,  only: ncent_tot,nctype_tot

      implicit none

      integer :: i, j, j1, k, ncent
      integer, dimension(ncent_tot) :: iwctype
      real(dp) :: pecent, r, r2, ri, ri2
      real(dp) :: tmp, term, zij
      real(dp), dimension(nctype_tot)  :: znuc
      real(dp), dimension(3,ncent_tot) :: cent
      real(dp), dimension(*), optional :: cos_n_sum   ! required for pot_nn_ewald
      real(dp), dimension(*), optional :: sin_n_sum   ! required for pot_nn_ewald
      
      if (nfrag.gt.1) potnnfrag = 0.0d0

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
            tmp = znuc(iwctype(i))*znuc(iwctype(j))/r
            pecent=pecent+tmp
            if ( (nfrag.gt.1) ) then
              potnnfrag(ifragcent(i)) = potnnfrag(ifragcent(i)) + 0.5d0 * tmp
              potnnfrag(ifragcent(j)) = potnnfrag(ifragcent(j)) + 0.5d0 * tmp
            endif
          enddo
        enddo
       else
        call pot_nn_ewald(cent,znuc,ncent,pecent,cos_n_sum,sin_n_sum)
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
