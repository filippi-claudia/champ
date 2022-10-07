      module grid3d
      use error, only: fatal_error
      contains
c 3d grid module
c Written by Anthony Scemama
c
c Example of input:
c
c &3dgrid nstepx 30 nstepy 30 nstepz 30
c &3dgrid stepx .1 stepy .1 stepz .1
c
c The following are computed if they are not present in the input.
c &3dgrid x0 0. y0 0. z0 0.
c &3dgrid xn 1. yn 1. zn 1.

c----------------------------------------------------------------------

      subroutine setup_grid()

      use grid_mod, only: MXNSTEP
      use grid_mod, only: IUNDEFINED, UNDEFINED, SHIFT
      use grid_mod, only: grid3d, cart_from_int
      use system, only: cent, ncent
!      use contrl, only: irstar
      use control_vmc, only: vmc_irstar
      use grid3d_param, only: endpt, nstep3d, origin, step3d

      implicit none

      integer :: i, input_ok, j, k


c     Test if the input is consistent

       input_ok=1
       do i=1,3
        if ( (nstep3d(i).eq.IUNDEFINED).and.(step3d(i).eq.UNDEFINED) ) then
          write (6,*) 'nstep3d(',i,') and step3d(',i,') are undefined'
          input_ok=0
        endif
        if ( origin(i).gt.endpt(i) ) then
         write (6,*)
     >  'The origin coordinates have to be smaller than the end points'
         input_ok=0
        endif
       enddo

       if (input_ok.ne.1)
     >   call fatal_error('3D Grid input inconsistent')

c      Origin and end of the grid. If not in input, compute it from the atomic
c      coordinates.

       do i=1, 3

        if ( origin(i).eq.UNDEFINED ) then
          origin(i) = cent(i,1)
          do j=2, ncent
            if ( cent(i,j).lt.origin(i) ) origin(i) = cent(i,j)
          enddo
          origin(i) = origin(i) - SHIFT
        endif

        if ( endpt(i).eq.UNDEFINED ) then
          endpt(i) = cent(i,1)
          do j=2, ncent
            if ( cent(i,j).gt.endpt(i) ) endpt(i) = cent(i,j)
          enddo
          endpt(i) = endpt(i) + SHIFT
        endif

       enddo

c      If the step is undefined, use the value in nstep3d to compute it.
c      Else, compute the value of nstep3d

       do i=1, 3

        if ( step3d(i).eq.UNDEFINED ) then
         step3d(i) = ( endpt(i) - origin(i) ) / (nstep3d(i)-1)
        else
         nstep3d(i) = int( (endpt(i) - origin(i))/step3d(i) )+1
        endif

        if ( nstep3d(i).gt.MXNSTEP ) then
         write (6,*) 'Warning: using nstep3d(',i,') = ',MXNSTEP
         nstep3d(i) = MXNSTEP
         endpt(i) = MXNSTEP*step3d(i) + origin(i)
        endif

       enddo

c      Prepare the integer->cartesian array

       do i=1,3
        do j=1, nstep3d(i)
         cart_from_int(j,i) = origin(i) + (j-1)*step3d(i)
        enddo
       enddo

c      Update the end point

       do i=1,3
         endpt(i) = cart_from_int(nstep3d(i),i)
       enddo

       do i=1,nstep3d(1)
        do j=1,nstep3d(2)
         do k=1,nstep3d(3)
          grid3d(i,j,k)=0.d0
         enddo
        enddo
       enddo

CACTIVATE
c     endif

c     Print the parameters to the output file

      write(45,*)
      write(45,*) '3D grid parameters'
      write(45,*) '------------------'
      write(45,*)
      write(45,*) 'Origin and end points'
      write(45,'(3(F10.6, 3X))') ( origin(i), i=1,3 )
      write(45,'(3(F10.6, 3X))') ( endpt (i), i=1,3 )
      write(45,*)
      write(45,*) 'Number of steps'
      write(45,'(3(I5, 3X))') ( nstep3d(i), i=1,3 )
      write(45,*)
      write(45,*) 'Step sizes'
      write(45,'(3(F10.6, 3X))') ( step3d (i), i=1,3 )
      write(45,*)

      end ! subroutine setup_grid

c----------------------------------------------------------------------

      function int_from_cart(value,iaxis)
      use grid_mod, only: IUNDEFINED
      use grid3d_param, only: endpt, origin, step3d
      use precision_kinds, only: dp
      implicit none

      integer :: iaxis
      real(dp) :: value
      integer :: int_from_cart


      if (( value.lt.origin(iaxis) ).or.
     &    ( value.ge.endpt (iaxis) )) then
        int_from_cart = IUNDEFINED
      else
        int_from_cart = int(( value-origin(iaxis) )/step3d(iaxis) +1.0)
      endif

      end ! subroutine int_from_cart

c----------------------------------------------------------------------

      subroutine write_cube(cube_file)

      use grid_mod, only: grid3d
      use system, only: znuc, cent, iwctype, ncent
      use system, only: nghostcent
      use grid3d_param, only: nstep3d, origin, step3d

      implicit none

      integer :: i, iu2, ix, iy, iz
      integer :: l

      character*(*) cube_file

c    Header
      iu2=13
      open (iu2,file=cube_file)
      write (iu2,*) ' Cube File'
      write (iu2,*) 'test'
 10   format (2X,I3,3(2X,F10.6))
 11   format (2X,I3,4(2X,F10.6))
      write (iu2,10) ncent, (origin(i), i=1,3)
      write (iu2,10) nstep3d(1), step3d(1), 0.d0, 0.d0
      write (iu2,10) nstep3d(2), 0.d0, step3d(2), 0.d0
      write (iu2,10) nstep3d(3), 0.d0, 0.d0, step3d(3)
c    Molecule
      do i=1, ncent+nghostcent
        write (iu2,11) int(znuc(iwctype(i))), znuc(iwctype(i)), (cent(l,i),l=1,3)
      enddo
c    Values
 20   format (6(E13.5))
      do ix=1,nstep3d(1)
       do iy=1,nstep3d(2)
         write (iu2,20) (grid3d(ix,iy,iz), iz=1,nstep3d(3))
       enddo
      enddo
      close(iu2)

      end

      end module
