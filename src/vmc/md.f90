subroutine md
!uses the velocity verlet scheme as in the wikipedia page https://en.wikipedia.org/wiki/Verlet_integration#Velocity_Verlet 

     use atom, only: ncent
     use md_mass,only: amass
     use md_var
     use md_fit, only: nfit, allocate_mdfit, deallocate_mdfit
     use precision_kinds, only: dp

  implicit none

     logical :: file_pos
     integer :: i, j, md_step, md_tot_step
     integer :: mdfit, forces_kick,nopt
     integer :: imove, md_cm
     real(kind=dp) :: e_kin, md_T
     character*20 :: md_method

  call p2gtid('mdyn:imove', imove,0,1)
  call p2gtid('mdyn:md_cm', md_cm,0,1)

  
  if (imove == 0) then
    !check if there is an external file to read the positions from
    INQUIRE(FILE="qmcoords", EXIST=file_pos)
    if (.not. file_pos) stop 'ERROR: THERE IS NOT TINKER COORD FILE!'

    open(75, FILE = 'write_forces',form='formatted',status='replace')
    write(6,*) "COUPLING WITH TINKER"

    do i = 1, ncent
      open(9 , FILE = 'qmcoords',form='formatted',status='old')
      read(9,*) pos(1,i), pos(2,i), pos(3,i)
    enddo
    write(6,*) "postions", pos

    call broadcas_md_pos()
    call get_forces_md(forces_ave,sigma,1,imove)
    close(75)

  elseif(imove == 2) then

    write(6,*) "USING STEEPEST DECENT TO MOVE ATOMS"

    call p2gtid('mdyn:md_step',md_tot_step,100,1)               !number of iterations
    write(6,*) "TOTAL # of STEPS", md_tot_step, "# of ATOMS", ncent
    call init_pos_md() 
    
    do i = i, md_tot_step
      call broadcas_md_pos()
      call get_forces_md(imove,1)
    enddo
  
  else
    
    write(6,*) "INITIALIZING MOLECULAR DYNAMICS"

    call p2gtid('mdyn:md_step',md_tot_step,100,1)   !number of iterations
    call p2gtfd('mdyn:md_dt', md_dt, 0.5,1)         !molecular dynamcis time step
    call p2gtfd('mdyn:md_T', md_T,0, 300)           !Temperature in Kelvin
    call p2gtid('mdyn:nfit', nfit,0,10)             !how many points in the fit
    call p2gtid('mdyn:nopt', nopt,0,20)            ! Optimize after nopt
    call p2gtad('mdyn:md_method', md_method, 'verlet',1)

    write(6,*) "TIME STEP", md_dt, "INPUT TEMPERATURE", md_T
    write(6,*) "TOTAL # of STEPS", md_tot_step, "# of ATOMS", ncent

    open(13, FILE = 'geometry_md.xyz',form='formatted',status='unknown')
    open(75, FILE = 'forces_md.xyz',form='formatted',status='unknown')
    open(15, FILE = 'velocities_md.xyz',form='formatted',status='unknown')
    open(99, FILE = 'energy.md',form='formatted', status='unknown')

    md_tau = 300                       !Tau for kick in the forces
    write(6,*) "TAU", md_tau

    call init_mass()
    write(6,*) "MASS", amass

    call init_pos_md() 
    call broadcas_md_pos()

    call init_md_vel(md_T)
 
    e_kin = 0.5*dot_product(amass,sum(vel**2, dim = 1))     

    call write_md_vel(md_step,md_dt)
    call write_md(0,md_dt,e_kin)
    call write_md_geo(md_step,md_dt)
    call write_md_force(md_step,md_dt)

    if(md_method == 'verlet') then
       write(6,*) "VERLET INTEGRATION SCHEME"
    elseif(md_method == 'verlet_fit') then
       call allocate_mdfit
       write(6,*) "VERLET+FIT INTEGRATION SCHEME"
    elseif(md_method == 'langevin') then
       write(6,*) "LANGEVIN INTEGRATION SCHEME"
    elseif(md_method == 'mix_verlang') then
       call allocate_mdfit
       write(6,*) "MIX-VERLANG INTEGRATION SCHEME"
    elseif(md_method == 'mix_memory') then
       call allocate_mdfit
       write(6,*) "MIX-MEMORY INTEGRATION SCHEME"
    else
      stop 'ERROR: NO INTEGRATION SCHEME SELECTED!' 
    endif
    ! computes accelerations
    do i = 1, ncent
      acc(:,i) = -forces_ave(:,i)/amass(i)
    enddo
 
    do md_step = 1, md_tot_step
 
       if(md_method=='verlet' ) then       
         call verlet(md_step, nopt)
       elseif(md_method == 'verlet_fit') then
         call verlet_fit(md_step, nopt)
       elseif(md_method == 'langevin') then
         call langevin(md_step, nopt, md_T)
       elseif(md_method == 'mix_verlang') then
         call mix_verlang(md_step, nopt, md_T)
       elseif(md_method == 'mix_memory') then
         call mix_memory(md_step, nopt, md_T)
       endif
       
       if(md_cm ==1) call remove_cm()
       
       e_kin = 0.5*dot_product(amass,sum(vel**2, dim = 1))     
       
       write(6,*) "END OF MD STEP", md_step
       
       call write_md_vel(md_step,md_dt)
       call write_md(md_step,md_dt,e_kin)
       call write_md_geo(md_step,md_dt)
       call write_md_force(md_step,md_dt)
    enddo
    
    write(6,*) "END OF MD"

    close(99)
    close(15)
    close(75)
    close(13)

    deallocate(sigma_md) 
    if(mdfit==1) call deallocate_mdfit
  endif
         
  deallocate(sigma)

 return
end


