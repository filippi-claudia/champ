      module pot_local_mod
      contains
      subroutine pot_local(x, pe)
      use fragments, only: eloc_i, elocfrag, ifragelec, ifragcent, nfrag
      use contrl_file, only: ounit
      use contrl_per, only: iperiodic
      use contrldmc, only: icut_e
      use control, only: ipr
      use distance_mod, only: r_ee,r_en
      use multiple_geo, only: pecent
      use precision_kinds, only: dp
      use pseudo,  only: nloc
      use ewald_breakup, only: pot_ee_ewald,pot_en_ewald
      use system,  only: iwctype,ncent,nelec,nghostcent,znuc
      implicit none

      integer :: i, ic, ij, j
      real(dp) :: pe, pe_ee, pe_en
      real(dp) :: tmp
      real(dp), dimension(3,*) :: x
      

!  pe from nucleus-nucleus repulsion
      pe=pecent
      pe_ee=0.d0
      pe_en=0.d0

      if (icut_e.lt.0) then
         do i=1,nelec
            eloc_i(i) = eloc_i(i) + pecent/nelec
         enddo
      endif

      if(iperiodic.eq.0) then

         if(nloc.eq.0) then
            do i=1,nelec
               do ic=1,ncent
                  tmp=-znuc(iwctype(ic))/r_en(i,ic)
                  pe=pe+tmp
                  if (icut_e.lt.0) then
                     eloc_i(i)=eloc_i(i)+tmp
                  endif
                  if (nfrag.gt.1) then
                     elocfrag(ifragelec(i)) = elocfrag(ifragelec(i)) + 0.5d0*tmp
                     elocfrag(ifragcent(ic)) = elocfrag(ifragcent(ic)) + 0.5d0*tmp
                  endif
               enddo
            enddo
         endif

         ij=0
         do i=2,nelec
            do j=1,i-1
               ij=ij+1
               tmp = 1/r_ee(ij)
               pe=pe+tmp
               
               if (icut_e.lt.0) then
                  eloc_i(i) = eloc_i(i) + 0.5d0 * tmp
                  eloc_i(j) = eloc_i(j) + 0.5d0 * tmp
               endif
               if (nfrag.gt.1) then
                  elocfrag(ifragelec(i)) = elocfrag(ifragelec(i)) + 0.5d0 * tmp
                  elocfrag(ifragelec(j)) = elocfrag(ifragelec(j)) + 0.5d0 * tmp
               endif
            enddo
         enddo

      else

         call pot_ee_ewald(x,pe_ee)

         pe=pe+pe_ee
      endif
      if(ipr.ge.3) write(ounit,'(''pe,pe_en(loc),pe_ee'',9f9.5)') pe,pe_en,pe_ee
      return
      end
end module
