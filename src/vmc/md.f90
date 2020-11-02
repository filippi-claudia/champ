subroutine md
!uses the velocity verlet scheme as in the wikipedia page https://en.wikipedia.org/wiki/Verlet_integration#Velocity_Verlet 
  implicit none

     logical :: file_pos
     integer, parameter :: dp = selected_real_kind(15, 307)
     integer :: i, j,md_step, md_tot_step, ncent,imove
     real(kind=dp) :: rand, md_dt, e_kin, md_T
     real(kind=dp), dimension(3) :: vel_cm, pos_cm
     real(kind=dp), dimension(:), allocatable:: mass
     real(kind=dp), dimension(:,:), allocatable:: acc, acc_new, vel
     real(kind=dp), dimension(:,:), allocatable:: pos, forces_ave, sigma

  call p2gti('atoms:natom',ncent,1)               !gets number of atoms
  call p2gtid('mdyn:imove', imove,0,1)

  allocate(pos(3,ncent), forces_ave(3,ncent), sigma(3,ncent))
 
  !check if there is an external file to read the positions from
  INQUIRE(FILE="qmcoords", EXIST=file_pos)
  
  if (imove == 0) then
    if (.not. file_pos) stop 'ERROR: THERE IS NOT TINKER COORD FILE!'

    open(75, FILE = 'write_forces',form='formatted',status='replace')
    write(6,*) "COUPLING WITH TINKER"

    do i = 1, ncent
      open(9 , FILE = 'qmcoords',form='formatted',status='old')
      read(9,*) pos(1,i), pos(2,i), pos(3,i)
    enddo
    write(6,*) "postions", pos

    call broadcas_md_pos(pos)
    call get_forces_md(forces_ave,sigma,imove)
    close(75)

  elseif(imove == 2) then

    write(6,*) "USING STEEPEST DECENT TO MOVE ATOMS"

    call p2gtid('mdyn:md_step',md_tot_step,100,1)               !number of iterations
    write(6,*) "TOTAL # of STEPS", md_tot_step, "# of ATOMS", ncent
    call init_pos_md(pos, forces_ave) 
    
    do i = i, md_tot_step
      call broadcas_md_pos(pos)
      call get_forces_md(forces_ave,sigma,imove)
    enddo
  
  else

    write(6,*) "INITIALIZING MOLECULAR DYNAMICS"

    call p2gtid('mdyn:md_step',md_tot_step,100,1)   !number of iterations
    call p2gtfd('mdyn:md_dt', md_dt, 0.5,1)         !molecular dynamcis time step
    call p2gtfd('mdyn:md_T', md_T,0, 300)           !Temperature in Kelvin

    write(6,*) "TIME STEP", md_dt, "INPUT TEMPERATURE", md_T
    write(6,*) "TOTAL # of STEPS", md_tot_step, "# of ATOMS", ncent

    allocate(mass(ncent)) 
    allocate(acc(3,ncent), acc_new(3,ncent), vel(3,ncent))

    open(13, FILE = 'geometry_md.xyz',form='formatted',status='unknown')
    open(99, FILE = 'energy_md',form='formatted',status='unknown')

    call init_mass(mass)
    write(6,*) "READ MASS", mass
    vel = 0.0000
    call init_random_seed()
    do i = 1, ncent
      do j = 1,3
        call gasdev_s(rand)
        vel(j,i) = rand*sqrt(3.1668105e-06*md_T/mass(i)) 
      enddo
    enddo 
    write(6,*) "VEL GAUSS", vel
 
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
 
    e_kin = 0.0
    write(6,*) "STARTING VEL", vel
    write(6,*) "TEMPERATURE", dot_product(mass,sum(vel**2, dim = 1))/(3*3.1668105e-06*ncent)     

    vel = vel*sqrt(md_T/(dot_product(mass,sum(vel**2, dim = 1))/(3*3.1668105e-06*ncent)))
    write(6,*) "vel", vel  
    write(6,*) "RESETTED TEMPERATURE", dot_product(mass,sum(vel**2, dim = 1))/(3*3.1668105e-06*ncent)
    write(6,*) "EKIN",  0.5*dot_product(mass,sum(vel**2, dim = 1))     
 
    do md_step = 1, md_tot_step
 
      ! compute new positions
      pos = pos + vel * md_dt + 0.5 * md_dt**2 * acc 
   
      call broadcas_md_pos(pos)

      call get_forces_md(forces_ave,sigma,imove)
 
       ! computes  new accelerations
       do i = 1, ncent
         acc_new(:,i) = -forces_ave(:,i)/mass(i)
       enddo

       vel = vel + 0.5 * (acc + acc_new) * md_dt
       e_kin = 0.5*dot_product(mass,sum(vel**2, dim = 1))     
       acc = acc_new
      ! writes molecular dynamics infos and geometry
      write(6,*) "MD STEP", md_step
      call write_md(md_step,md_dt,e_kin)
      call write_md_geo(md_step,md_dt, pos)
    enddo

    deallocate(mass,acc,acc_new, vel)
    write(6,*) "END OF MD"
    close(99)
    close(13) 
  endif

  deallocate(pos, forces_ave, sigma)

 return
end
