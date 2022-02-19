module qmmm_pot
  integer, parameter :: dbl = kind(1.0d0)
contains
!*********************************************************************
        subroutine qmmm_extpot_read
!*********************************************************************
     
        use qmmm_extpot

        implicit none
        integer :: i,ix,iy,iz
        double precision :: tmp(3)

        filename_cube= 'potential.cube'
        open (12,file=filename_cube,form='formatted',status='old')
        read (12,'(a80)') title_cube(1)
        read (12,'(a80)') title_cube(2)
        read (12,*) n_atoms,x0(1),x0(2),x0(3)
        read (12,*) n_x,delta(1),tmp(2),tmp(3)
        read (12,*) n_y,tmp(1),delta(2),tmp(3)
        read (12,*) n_z,tmp(1),tmp(2),delta(3)
        allocate(x_atom(n_atoms,3), source=0.0_dbl)
        allocate(id_atom(n_atoms), source=0)
        allocate(chrg_atom(n_atoms), source=0.0_dbl)
        allocate(pot(n_x,n_y,n_z), source=0.0_dbl)
        do i=1,n_atoms
          read(12,*) id_atom(i),chrg_atom(i),x_atom(i,:)
        enddo
        do ix=1,n_x
          do iy=1,n_y
!           read(12,'(6e13.5)') pot(ix,iy,:)
!           read(12,'(90e13.5)') pot(ix,iy,:)
            read(12,*) (pot(ix,iy,iz),iz=1,n_z)
          enddo
        enddo

        allocate(xdata(n_x),ydata(n_y),zdata(n_z), source=0.0_dbl)
         do i = 1, n_x
            xdata(i) = x0(1)+dble(i-1)*delta(1)
!           write (6,*) xdata(i)
         end do
         do i = 1, n_y
            ydata(i) = x0(2)+dble(i-1)*delta(2)
!           write (6,*) ydata(i)
         end do
         do i = 1, n_z
            zdata(i) = x0(3)+dble(i-1)*delta(3)
!           write (6,*) zdata(i)
         end do

!       do i=1,n_x
!        iy = (n_y-1)/2+1
!        iz = (n_z-1)/2+1
!        write(111,'(4f12.6)') xdata(i),ydata(iy),zdata(iz), &
!    &                  pot(i,iy,iz)
!       enddo

       write (6,*) "Reading cube file ...",filename_cube
       write (6,'(a15,3i7)') 'Dimensions  :' ,n_x,n_y,n_z
       write (6,'(a10,3f12.6)') 'Mesh  : ',delta

!      call qmmm_writecube("read.cube",title_cube,n_atoms,x0,n_x,n_y,n_z &
!   &         ,delta,x_atom,id_atom,chrg_atom,pot)

       call qmmm_interpolate

        return
        end

!*********************************************************************
        subroutine qmmm_interpolate
!*********************************************************************
! set up the interpolation parameters of the splines

        use qmmm_extpot
        use qmmm_splines
        use qmmm_bspline

        implicit none
      
        kxord=5
        kyord=5
        kzord=5

        nxknot=n_x+kxord
        nyknot=n_y+kyord
        nzknot=n_z+kzord

        allocate(bscoef(n_x,n_y,n_z), source=0.0_dbl)
        allocate(xknot(nxknot),yknot(nyknot),zknot(nzknot), source=0.0_dbl)

!        generate knots

        call dbsnak (n_x, xdata, kxord, xknot)
        call dbsnak (n_y, ydata, kyord, yknot)
        call dbsnak (n_z, zdata, kzord, zknot)

!        interpolate
        write (6,*) 
        write (6,'(a50,3i4)') &
     &  'Interpolating with 3D spline of order (x,y,z) : '  &
     &  ,kxord,kyord,kzord

        call dbs3in (n_x, xdata, n_y, ydata, n_z, zdata, pot,  &
     &            n_x, n_y, kxord, kyord, kzord, xknot, yknot, zknot,   &
     &            bscoef)
  
      write (6,*) 'Interpolation done.'
      nxcoef = n_x
      nycoef = n_y
      nzcoef = n_z

      return
      end

!*********************************************************************
      subroutine qmmm_evaluate(nxvec,nyvec,nzvec,xvec,yvec,zvec,value)
!*********************************************************************
! evaluate the function on the points

        use qmmm_extpot
        use qmmm_splines
        use qmmm_bspline

        implicit none
        integer :: nxvec,nyvec,nzvec
        integer :: i,j,k
        double precision ::  f
        double precision :: xvec(nxvec),yvec(nyvec),zvec(nzvec)
        double precision :: value(nxvec,nyvec,nzvec)


      call dbs3gd (0, 0, 0, nxvec, xvec, nyvec, yvec, nzvec, zvec,      &
     &            kxord, kyord, kzord, xknot, yknot, zknot, nxcoef,     &
     &            nycoef, nzcoef, bscoef, value, nxvec, nyvec)


      do i = 1, nxvec
         do j = 1, nyvec
            do k = 1, nzvec
!              write (6,'(4f13.4, f13.6)') xvec(i), yvec(j),            &
!    &                                       zvec(k), value(i,j,k),     &
!    &                                       f(xvec(i),yvec(j),zvec(k)) &
!    &                                        - value(i,j,k)
            end do
         end do
      end do


      return
      end

!*********************************************************************
      subroutine qmmm_extpot_ene(coord,nelec,ext_ene)
!*********************************************************************

       use qmmm_extpot
       use qmmm_splines
       use qmmm_vector
       use qmmm_bspline

       implicit none
 
       integer :: nelec
       double precision :: coord(3,nelec), ext_ene
       integer :: n
       double precision :: point(3),valuep(1,1,1),x1(3)
     
       x1(1)=x0(1)+(n_x-1)*delta(1)
       x1(2)=x0(2)+(n_y-1)*delta(2)
       x1(3)=x0(3)+(n_z-1)*delta(3)
       ext_ene=0.d0 
       do n=1, nelec
         point(:)=coord(:,n)
! check if the point is out of the potential box
         if(point(1).gt.x0(1) .and. point(1).lt.x1(1) .and. &
&           point(2).gt.x0(2) .and. point(2).lt.x1(2) .and. &
&           point(3).gt.x0(3) .and. point(3).lt.x1(3) ) then
           call qmmm_evaluate(1,1,1,point(1),point(2),point(3),valuep)
!          write (112,'(4f12.6)') point(:),valuep
           ext_ene = ext_ene + valuep(1,1,1)
         else
           ext_ene = 0.d0
           nout = nout + 1
         endif 
       enddo 
       ncount = ncount + 1
       ave = ave + ext_ene
       ave2 = ave2 + ext_ene*ext_ene

       return
       end

!*********************************************************************
        subroutine qmmm_extpot_final(nelec)
!*********************************************************************
        
        use qmmm_extpot

        implicit none
      
        integer :: nelec
 
        ave = ave /ncount
        err = dsqrt((ave2/ncount - ave**2)/ncount)
        
        write (6,'(''QMMM - 1proc Total number of N-el evaluations: '',i6)') ncount
        write (6,'(''QMMM - 1proc Rate of out of box electrons: '',d12.6)') nout/dble(ncount*nelec)
        write (6,'(''QMMM - 1proc Average QMC/MM energy: '',f12.6,''('',f12.6,'')'')') ave,err

        return
        end

!*********************************************************************
        subroutine qmmm_test_points
!*********************************************************************
! set up a grid of points where to evaluate the function

        use qmmm_extpot
        use qmmm_splines
        use qmmm_vector
        use qmmm_bspline

        implicit none
        integer :: i
      
        double precision ::  f
 
        nxvec=83
        nyvec=83
        nzvec=83

        allocate(value(nxvec,nyvec,nzvec), source=0.0_dbl)
        allocate(xvec(nxvec),yvec(nyvec),zvec(nzvec), source=0.0_dbl)


!     write (6,99999)
 
        deltavec(1)=delta(1)/dble(nxvec-1)*dble(n_x-1)
        deltavec(2)=delta(2)/dble(nyvec-1)*dble(n_y-1)
        deltavec(3)=delta(3)/dble(nzvec-1)*dble(n_z-1)

         do i = 1, nxvec
            xvec(i) = x0(1)+dble(i-1)*deltavec(1)
!           write (6,*) xvec(i)
         end do
         do i = 1, nyvec
            yvec(i) = x0(2)+dble(i-1)*deltavec(2)
!           write (6,*) yvec(i)
         end do
         do i = 1, nzvec
            zvec(i) = x0(3)+dble(i-1)*deltavec(3)
!           write (6,*) zvec(i)
         end do

       write (6,*)
       write (6,'(a35,3i7)') 'Dimensions of the evaluating box :' &
     &  ,nxvec,nyvec,nzvec
       write (6,'(a30,3f12.6)') 'Mesh of the evaluating box : ',deltavec
       write (6,'(a30,3f12.6)') 'Box dimensions : ' &
     &   ,(nxvec-1)*deltavec(1), (nyvec-1)*deltavec(2),(nzvec-1)*deltavec(3)

!     99999 format (10x, 'x', 11x, 'y', 10x, 'z', 10x, 's(x,y,z)', 7x,        &
!    &       'error')


      return
      end

end module
