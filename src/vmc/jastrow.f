      module jastrow_mod
      contains
      subroutine jastrow_factor(x,v,d2,value,ifr)
c Written by Cyrus Umrigar

      use deriv_jastrow1_mod, only: deriv_jastrow1
      use deriv_jastrow4_mod, only: deriv_jastrow4
      use jastrow1_mod, only: jastrow_factor1
      use jastrow4_mod, only: jastrow_factor4
      use optwf_control, only: ioptjas
      use precision_kinds, only: dp
      use system,  only: nelec
      use jastrow, only: ijas
      use contrl_file, only: ounit
      implicit none

      integer :: i, ifr
      real(dp) :: d2, value
      real(dp), dimension(3, *) :: x
      real(dp), dimension(3, *) :: v

      do i=1,nelec
        v(1,i)=0
        v(2,i)=0
        v(3,i)=0
      enddo
      d2=0

c       write(ounit,*),"d2o before", d2
c       write(ounit,*),"fsumo before", value
      
      if(ifr.gt.1.or.ioptjas.eq.0) then
        

        if(ijas.eq.1) then
          call jastrow_factor1(x,v,d2,value)
         else
          call jastrow_factor4(x,v,d2,value)
       endif
       
      else
        if(ijas.eq.1) then
          call deriv_jastrow1(x,v,d2,value)
         else
          call deriv_jastrow4(x,v,d2,value)
        endif
      endif


c      write(ounit,*),"d2o after", d2
c      write(ounit,*),"fsumo after", value
      
      return
      end
      end module
