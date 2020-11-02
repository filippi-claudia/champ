!LANGEVINES DYNAMICS SCHEME OF SORELLA https://doi.org/10.1063/1.4901430
!COULD ALSO TRY FOLLOWING THE SCHEME OF GROMACS 
!http://manual.gromacs.org/current/reference-manual/algorithms/stochastic-dynamics.html
!
subroutine langev
  implicit none

     integer, parameter :: dp = selected_real_kind(15, 307)
     integer :: i, j
     integer  :: md_step, md_tot_step, ncent
     real(kind=dp) :: rand, md_dt, e_kin, bond_dist_sq, md_T
     real(kind=dp), dimension(3) :: vel_cm, pos_cm
     real(kind=dp), dimension(:), allocatable:: mass
     real(kind=dp), dimension(:,:), allocatable:: gamman ,alpha,beta
     real(kind=dp), dimension(:,:), allocatable:: acc, vel
     real(kind=dp), dimension(:,:), allocatable:: pos, forces_ave, sigma, dumm

  !get variables from input     
  call p2gtfd('mdyn:md_dt', md_dt, 0.5,1)         !molecular dynamcis time step
  call p2gtfd('mdyn:md_T', md_T,0, 300)          !molecular dynamcis time step
  call p2gtid('mdyn:md_step',md_tot_step,100,1)   !number of molecular dynamics iterations
  call p2gti('atoms:natom',ncent,1)               !number of atoms

  write(6,*) "TIME STEP", md_dt, "INPUT TEMPERATURE", md_T
  write(6,*) "TOTAL # of MD STEPS", md_tot_step, "# of ATOMS", ncent

  allocate(mass(ncent), alpha(3,ncent), beta(3,ncent)) 
  allocate(acc(3,ncent), vel(3,ncent), gamman(3,ncent))
  allocate(pos(3,ncent), forces_ave(3,ncent), sigma(3,ncent), dumm(3,ncent))

  open(13, FILE = 'geometry_md.xyz',form='formatted',status='unknown')
  open(99, FILE = 'energy_md',form='formatted',status='unknown')
  open(15) 

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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! uniformly dsitrubute velocities
!  call random_number(vel)
!  do i = 1, ncent
!    vel(:,i) = (vel(:,i)-0.5)*2*sqrt(3.)*sqrt(3.1668105e-06*md_T/mass(i))     !to be checked!   
!  enddo
!  write(6,*) "VEL UNIFORM", vel
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  call optwf_sr
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
  write(6,*) "EKIN",  0.5*dot_product(mass,sum(vel**2, dim = 1))     
  write(6,*) "TEMPERATURE", dot_product(mass,sum(vel**2, dim = 1))/(3*3.1668105e-06*ncent)     
  write(6,*) " STARTING LANGEVIN DYNAMICS" 

  vel = vel*sqrt(md_T/(dot_product(mass,sum(vel**2, dim = 1))/(3*3.1668105e-06*ncent)))
  write(6,*) "vel", vel  
  write(6,*) "RESETTED TEMPERATURE", dot_product(mass,sum(vel**2, dim = 1))/(3*3.1668105e-06*ncent)     
  
  do md_step = 1, md_tot_step

    call get_forces_md(forces_ave,sigma)
    
    !set langevin parameters
    gamman = md_dt*sigma**2/(2*3.1668105e-06*md_T)                       !0.168413 assuming sigma = 0.004
    write(6,*) "gammn", gamman
    
    do i = 1, ncent
      alpha(:,i) = (1.0/gamman(:,i))*(1.0 - exp(-gamman(:,i)*md_dt/mass(i)))      !capital gamma in Sorella
      beta(:,i)  = exp(-gamman(:,i)*md_dt/mass(i))
    enddo 
    write(6,*) "alpha", alpha
    write(6,*) "beta", beta

!    do i = 1, ncent
    ! computes velocities
!    vel(:,i) = beta(i) * vel(:,i) + alpha(i) *(-forces_ave(:,i))
    vel = beta * vel + alpha *(-forces_ave)
!    enddo
    ! compute new positions
    pos = pos + vel * md_dt 
    call broadcas_md_pos(pos)

    ! sum is summing on x,y,z 
    e_kin = 0.5*dot_product(mass,sum(vel**2, dim = 1))     

    ! write time step and current total energy and write geometry for visualizaiton. subroutines in md_more.f90
    call write_md(md_step,md_dt,e_kin)
    call write_md_geo(md_step,md_dt, pos)


   !temporary checks
    write (6,*) 'c2_bondlenght', sqrt(sum((pos(:,1)-pos(:,2))**2))
    write (6,*) 'velocities', vel
    write (6,*) 'accelerations', acc

  enddo
 
  deallocate(mass,acc, vel, alpha, beta)
  deallocate(pos, forces_ave, dumm)

  write(6,*) "END OF MD"
 
  close(13) 
  close(99)
  close(15)  
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
