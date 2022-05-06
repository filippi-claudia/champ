subroutine init_md_vel(md_T)

     use precision_kinds, only: dp
     use atom, only: ncent
     use md_mass,only: amass
     use md_var,only:vel

  implicit none

     logical :: file_vel
     integer :: i, j
     real(kind=dp) :: rand, e_kin, md_T

    vel = 0.0000

    INQUIRE(FILE="input_vel", EXIST=file_vel)
    if ( file_vel) then
      do i = 1, ncent
        open(9 , FILE = 'input_vel',form='formatted',status='old')
        read(9,*) vel(1,i), vel(2,i), vel(3,i)
      enddo
      write(6,*) "READING INPUT VELOCITIES"
      write(6,*) vel
    else
      vel = 0.0000
      call init_random_seed()
      do i = 1, ncent
        do j = 1,3
          call gasdev_s(rand)
          vel(j,i) = rand*sqrt(3.1668105e-06*md_T/amass(i))
        enddo
      enddo
    endif

    if( md_T /= 0.0d0) call remove_cm()

    if ( .not. file_vel .and. md_T /= 0) then
      vel = vel*sqrt(md_T/(dot_product(amass,sum(vel**2, dim = 1))/(3*3.1668105e-06*ncent)))
    endif

    e_kin =  0.5*dot_product(amass,sum(vel**2, dim = 1))

!    write(6,*) "STARTING VEL", vel
    write(6,*) "TEMPERATURE", (2.0d0*e_kin/(3*3.1668105e-06*ncent))     
!    write(6,*) "EKIN",e_kin 

   return
end


     subroutine remove_cm

     use precision_kinds, only: dp
     use atom, only: ncent
     use md_mass,only: amass
     use md_var,only: vel, pos

     implicit none
     integer i,j
     real(kind=dp):: etrans,erot
     real(kind=dp):: weigh,totmass,eps
     real(kind=dp):: xx,yy,zz,xy,xz,yz
     real(kind=dp):: xtot,ytot,ztot
     real(kind=dp):: xdel,ydel,zdel
     real(kind=dp):: mang(3),vang(3)
     real(kind=dp):: vtot(3),tensor(3,3)

     real(kind=dp), dimension(:),   allocatable :: x, y, z
     real(kind=dp), dimension(:,:), allocatable :: v
     
     allocate(x(ncent), y(ncent),z(ncent))
     allocate(v(3,ncent))

     x(:) = pos(1,:) 
     y(:) = pos(2,:) 
     z(:) = pos(3,:) 
     v = vel

     totmass = 0.0d0
     do j = 1, 3
        vtot(j) = 0.0d0
     end do

!    compute linear velocity of the system center of mass
      do i = 1, ncent
         weigh = amass(i)
         totmass = totmass + weigh
         do j = 1, 3
            vtot(j) = vtot(j) + v(j,i)*weigh
         end do
      end do

!     compute translational kinetic energy of overall system
      etrans = 0.0d0
      do j = 1, 3
         vtot(j) = vtot(j) / totmass
         etrans = etrans + vtot(j)**2
      end do
      etrans = 0.5d0 * etrans * totmass 

!     find the center of mass coordinates of the overall system
      xtot = 0.0d0
      ytot = 0.0d0
      ztot = 0.0d0
         do i = 1, ncent
            weigh = amass(i)
            xtot = xtot + x(i)*weigh
            ytot = ytot + y(i)*weigh
            ztot = ztot + z(i)*weigh
         end do
      xtot = xtot / totmass
      ytot = ytot / totmass
      ztot = ztot / totmass

!     compute the angular momentum of the overall system

      do j = 1, 3
         mang(j) = 0.0d0
      end do
      do i = 1, ncent
         weigh = amass(i)
         mang(1) = mang(1) + (y(i)*v(3,i)-z(i)*v(2,i))*weigh
         mang(2) = mang(2) + (z(i)*v(1,i)-x(i)*v(3,i))*weigh
         mang(3) = mang(3) + (x(i)*v(2,i)-y(i)*v(1,i))*weigh
      end do
      mang(1) = mang(1) - (ytot*vtot(3)-ztot*vtot(2))*totmass
      mang(2) = mang(2) - (ztot*vtot(1)-xtot*vtot(3))*totmass
      mang(3) = mang(3) - (xtot*vtot(2)-ytot*vtot(1))*totmass

!     calculate the moment of inertia tensor

      xx = 0.0d0
      xy = 0.0d0
      xz = 0.0d0
      yy = 0.0d0
      yz = 0.0d0
      zz = 0.0d0
      do i = 1, ncent
         weigh = amass(i)
         xdel = x(i) - xtot
         ydel = y(i) - ytot
         zdel = z(i) - ztot
         xx = xx + xdel*xdel*weigh
         xy = xy + xdel*ydel*weigh
         xz = xz + xdel*zdel*weigh
         yy = yy + ydel*ydel*weigh
         yz = yz + ydel*zdel*weigh
         zz = zz + zdel*zdel*weigh
      end do
      tensor(1,1) = yy + zz
      tensor(2,1) = -xy
      tensor(3,1) = -xz
      tensor(1,2) = -xy
      tensor(2,2) = xx + zz
      tensor(3,2) = -yz
      tensor(1,3) = -xz
      tensor(2,3) = -yz
      tensor(3,3) = xx + yy

!     fix to avoid singularity for one- or two-body systems

       if (ncent .le. 2) then
          eps = 0.000001d0
          tensor(1,1) = tensor(1,1) + eps
          tensor(2,2) = tensor(2,2) + eps
          tensor(3,3) = tensor(3,3) + eps
       end if

!     diagonalize the moment of inertia tensor
       call invert (3,tensor)

!     compute angular velocity and rotational kinetic energy
      erot = 0.0d0
      do i = 1, 3
         vang(i) = 0.0d0
         do j = 1, 3
            vang(i) = vang(i) + tensor(i,j)*mang(j)
         end do
         erot = erot + vang(i)*mang(i)
      end do
      erot = 0.5d0 * erot 

!     eliminate any translation of the overall system

      do i = 1, ncent
         do j = 1, 3
            v(j,i) = v(j,i) - vtot(j)
         end do
      end do

!     print the translational velocity of the overall system

!      if (debug) then
!         write (iout,10)  (vtot(i),i=1,3),etrans
!   10    format (' System Linear Velocity :  ',3d12.2,
!     &           /,' Translational Kinetic Energy :',10x,f12.4,
!     &              ' Kcal/mole')
!      end if

!     eliminate any rotation about the system center of mass

        do i = 1, ncent
           xdel = x(i) - xtot
           ydel = y(i) - ytot
           zdel = z(i) - ztot
           v(1,i) = v(1,i) - vang(2)*zdel + vang(3)*ydel
           v(2,i) = v(2,i) - vang(3)*xdel + vang(1)*zdel
           v(3,i) = v(3,i) - vang(1)*ydel + vang(2)*xdel
       end do

!     print the angular velocity of the overall system

!         if (debug) then
!            write (iout,20)  (vang(i),i=1,3),erot
!   20       format (' System Angular Velocity : ',3d12.2,
!     &              /,' Rotational Kinetic Energy :',13x,f12.4,
!     &                 ' Kcal/mole')
!         end if

     
     write(6,*) "center of mass removed"
!     write(6,*) "vel in", vel

     pos(1,:) = x(:) 
     pos(2,:) = y(:) 
     pos(3,:) = z(:) 
     vel = v

!     write(6,*) "vel out", vel
     deallocate(x, y,z)
     deallocate(v)

      return
      end


SUBROUTINE init_random_seed()
  INTEGER :: i, n, clock
  INTEGER, DIMENSION(:), ALLOCATABLE :: seed

  CALL RANDOM_SEED(size = n)
  ALLOCATE(seed(n))

  CALL SYSTEM_CLOCK(COUNT=clock)

  seed = clock + 37 * (/ (i - 1, i = 1, n) /)
  CALL RANDOM_SEED(PUT = seed)

  DEALLOCATE(seed)
END SUBROUTINE
