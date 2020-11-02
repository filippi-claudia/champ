subroutine md

  implicit none

     integer, parameter :: dp = selected_real_kind(15, 307)
     integer :: i, j,md_step, md_tot_step, ncent
     real(kind=dp) :: md_dt, e_kin, bond_dist_sq
     real(kind=dp), dimension(3) :: vel_cm, pos_cm
     real(kind=dp), dimension(:), allocatable:: mass
     real(kind=dp), dimension(:,:), allocatable:: acc, acc_new, vel
     real(kind=dp), dimension(:,:), allocatable:: pos, forces_ave, dumm

  !get variables from input     
  call p2gtfd('mdyn:md_dt', md_dt, 0.5,1)         !molecular dynamcis time step
  call p2gtid('mdyn:md_step',md_tot_step,100,1)   !number of molecular dynamics iterations
  call p2gti('atoms:natom',ncent,1)               !number of atoms

  write(6,*) "TIME STEP", md_dt
  write(6,*) "TOTAL # of MD STEPS", md_tot_step, "# of ATOMS", ncent

  allocate(mass(ncent)) 
  allocate(acc(3,ncent), acc_new(3,ncent), vel(3,ncent))
  allocate(pos(3,ncent), forces_ave(3,ncent), dumm(3,ncent))

  open(13, FILE = 'geometry_md.xyz',form='formatted',status='unknown')
!  open(15, FILE = 'write_forces',form='formatted',status='unknown')
  open(99, FILE = 'energy_md',form='formatted',status='unknown')
 
  vel = 0.0000
  call random_number(vel)
  vel = (vel-0.5)*0.0008
!  write(6,*) "VELO", vel
!  vel(3,1) = 0.00035
!  vel(3,2) = -0.00035

  call init_mass(mass)
  write(6,*) "READ MASS", mass
  e_kin = 0.0


  call optwf_sr

  ! subroutine in md_more.f to get positions and forces
  call init_pos_md(pos, forces_ave) 

  ! computes accelerations
  do i = 1, ncent
    acc(:,i) = -forces_ave(:,i)/mass(i)
  enddo

  ! subtract center of mass velocity and positions
  do i =1,3
  vel_cm(i) = dot_product(mass,vel(i,:))/sum(mass)
  vel(i,:) = vel(i,:) - vel_cm(i)
  pos_cm(i) = dot_product(mass,pos(i,:))/sum(mass)
  pos(i,:) = pos(i,:) - pos_cm(i)
  enddo


  write(6,*) "STARTING VEL", vel
  write(6,*) "EKIN",  0.5*dot_product(mass,sum(vel**2, dim = 1))     
  write(6,*) "TEMPERATURE", dot_product(mass,sum(vel**2, dim = 1))/(3*3.1668105e-06*ncent)     
  
  do md_step = 1, md_tot_step

    ! compute new positions
    pos = pos + vel * md_dt + 0.5 * md_dt**2 * acc 
 
    call broadcas_md_pos(pos)
    call get_forces_md(forces_ave)

    ! computes  new accelerations
    do i = 1, ncent
      acc_new(:,i) = -forces_ave(:,i)/mass(i)
    enddo

    ! computes velocities
    vel = vel + 0.5 * (acc + acc_new) * md_dt

    ! sum is summing on x,y,z 
    e_kin = 0.5*dot_product(mass,sum(vel**2, dim = 1))     
!    write(6,*) "ekin", e_kin
    acc = acc_new

    ! write time step and current total energy and write geometry for visualizaiton. subroutines in md_more.f90
    call write_md(md_step,md_dt,e_kin)
    call write_md_geo(md_step,md_dt, pos)


   !temporary checks
    write (6,*) 'c2_bondlenght', sqrt(sum((pos(:,1)-pos(:,2))**2))
    write (6,*) 'velocities', vel
    write (6,*) 'accelerations_new', acc_new

  enddo
 
  deallocate(mass,acc,acc_new, vel)
  deallocate(pos, forces_ave, dumm)

  write(6,*) "END OF MD"
 
  close(13) 
  close(99)
!  close(15) 
  
 return
end
