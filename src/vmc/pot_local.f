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

c         do i=1,nelec
c            do ic=1,ncent+nghostcent
c               if(nloc.eq.0.and.ic.le.ncent) pe=pe-znuc(iwctype(ic))/r_en(i,ic)
c            enddo
c        enddo
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
c         write(ounit,'(''00 pe'',20f12.4)') pe
c         call pot_en_ewald(x,pe_en)
         call pot_ee_ewald(x,pe_ee)
c     write(ounit,'(''02 pe_ee'',20d12.4)') pe_ee
         pe=pe+pe_en+pe_ee
      endif
c     write(ounit,'(''03 pe,pe_en(loc),pe_ee'',20d12.4)') pe,pe_en,pe_ee
      if(ipr.ge.3) write(ounit,'(''pe,pe_en(loc),pe_ee'',9f9.5)') pe,pe_en,pe_ee
      return
      end
      end module
