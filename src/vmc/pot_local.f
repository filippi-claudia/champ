      module pot_local_mod
      contains
      subroutine pot_local(x, pe)
      use contrl_file, only: ounit
      use contrl_per, only: iperiodic
      use control, only: ipr
      use distance_mod, only: r_ee,r_en
      use multiple_geo, only: pecent
      use precision_kinds, only: dp
      use pseudo,  only: nloc
      use pw_ewald, only: pot_ee_ewald,pot_en_ewald
      use system,  only: iwctype,ncent,nelec,nghostcent,znuc
      implicit none

      integer :: i, ic, ij, j
      real(dp) :: pe, pe_ee, pe_en
      real(dp), dimension(3,*) :: x

c  pe from nucleus-nucleus repulsion
      pe=pecent
      pe_ee=0.d0
      pe_en=0.d0
      if(iperiodic.eq.0) then

         if(nloc.eq.0) then
            do i=1,nelec
               do ic=1,ncent
                  pe=pe-znuc(iwctype(ic))/r_en(i,ic)
               enddo
            enddo
         endif
         
         ij=0
         do i=2,nelec
            do j=1,i-1
               ij=ij+1
               pe=pe+1/r_ee(ij)
            enddo
         enddo
         
      else

         
         call pot_ee_ewald(x,pe_ee)
         pe=pe+pe_en+pe_ee

      endif
      if(ipr.ge.3) write(ounit,'(''pe,pe_en(loc),pe_ee'',9f9.5)') pe,pe_en,pe_ee
      return
      end
      end module
