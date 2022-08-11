      module efield_f_mod
      contains
      subroutine efield_extpot_ene(coord,nelec,efield_pot)

      use efield,  only: iscreen,ncharges
      use efield_blk, only: ascreen,bscreen,qcharge,xcharge,ycharge
      use efield_blk, only: zcharge
      use precision_kinds, only: dp


      implicit none

      integer :: i, j, nelec
      real(dp) :: atmp, efield_pot, rij, rij2
      real(dp), dimension(3, *) :: coord




      atmp=0
      efield_pot=0
      do i=1,ncharges
         do j=1,nelec
           rij2=(xcharge(i)-coord(1,j))**2+(ycharge(i)-coord(2,j))**2+(zcharge(i)-coord(3,j))**2
           rij=sqrt(rij2)
           if(iscreen.gt.0) atmp=ascreen(i)*dexp(-bscreen(i)*rij**2)
           efield_pot=efield_pot-qcharge(i)*(1.d0-atmp)/rij
         enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine efield_compute_extint

      use contrl_file, only: ounit
      use efield,  only: iscreen,ncharges
      use efield_blk, only: ascreen,bscreen,qcharge,xcharge,ycharge
      use efield_blk, only: zcharge
      use precision_kinds, only: dp
      use system,  only: cent,iwctype,ncent,znuc

      implicit none

      integer :: i, j
      real(dp) :: atmp, efield_extext, efield_nucext, rij, rij2




      efield_nucext=0
      atmp=0
      do i=1,ncharges
c         write(ounit,*) 'HELLO',i,ascreen(i),bscreen(i)
         do j=1,ncent
           rij2=(xcharge(i)-cent(1,j))**2+(ycharge(i)-cent(2,j))**2+(zcharge(i)-cent(3,j))**2
           rij=sqrt(rij2)
           if(iscreen.gt.0) atmp=ascreen(i)*dexp(-bscreen(i)*rij**2)
           efield_nucext=efield_nucext+znuc(iwctype(j))*qcharge(i)*(1.d0-atmp)/rij
         enddo
      enddo

      efield_extext=0
      do i=1,ncharges
         do j=i+1,ncharges
           rij2=(xcharge(i)-xcharge(j))**2+(ycharge(i)-ycharge(j))**2+(zcharge(i)-zcharge(j))**2
           rij=sqrt(rij2)
           efield_extext=efield_extext+qcharge(i)*qcharge(j)/rij
         enddo
      enddo

      write(ounit,'(''Number of charges: '',i5)') ncharges
      write(ounit,'(''Number of centers: '',i5)') ncent

c      write(ounit,*) 'znuc=',znuc(1:8)

      write(ounit,'(''Nuclear-external interaction: '',f15.9)') efield_nucext
      write(ounit,'(''External-external interaction:'',f15.9)') efield_extext
      return
      end
c-----------------------------------------------------------------------
      end module
