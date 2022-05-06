Subroutine verlet(md_step, nopt)
     use atom, only: ncent
     use md_mass,only: amass
     use md_var
     use precision_kinds, only: dp

  implicit none

     logical :: md_kick
     integer :: i, j, md_step
     integer :: forces_kick, nopt, md_opt

    
!     if (md_method .ne. 'verlet') return

     pos = pos + vel * md_dt + 0.5 * md_dt**2 * acc 
     call broadcas_md_pos()

     md_opt = 0
     if(mod(md_step, nopt) == 0.0d0) md_opt =1
     call get_forces_md(md_opt, 1)                 !imove is 1 becouse we are doing dynamics not steepest descent

     forces_kick = 1
     if (md_step .le. 2) then
       sigma_md = sigma
     else
       do while (forces_kick .ge. 1)
         md_kick = any(sum(sigma,dim=1) .ge. sum((2.5*sigma_md),dim=1))
         write(6,*) "MD KICK", md_kick
         if(md_kick) then
            call get_forces_md(0, 1)
         else
            forces_kick = 0
            sigma_md = (1-md_dt/md_tau)*sigma_md + md_dt/md_tau*sigma
         endif
       enddo
     endif
   
    ! computes  new accelerations
    do i = 1, ncent
      acc_new(:,i) = -forces_ave(:,i)/amass(i)
    enddo
  
    vel = vel + 0.5 * (acc + acc_new) * md_dt

    acc = acc_new

    return
    
end subroutine

subroutine verlet_fit(md_step, nopt)
     use atom, only: ncent
     use md_mass,only: amass
     use md_var
     use md_fit, only: nfit
     use precision_kinds, only: dp

  implicit none

     logical :: md_kick
     integer :: i, j, md_step
     integer :: forces_kick, nopt, md_opt
    
!     if (md_method .ne. 'verlet_fit') return

     pos = pos + vel * md_dt + 0.5 * md_dt**2 * acc 
     call broadcas_md_pos()

     md_opt = 0
     if(mod(md_step, nopt) == 0.0d0) md_opt =1
     call get_forces_md(md_opt, 1)                 !imove is 1 becouse we are doing dynamics not steepest descent

     forces_kick = 1
     if (md_step .le. 2) then
       sigma_md = sigma
     else
       do while (forces_kick .ge. 1)
         md_kick = any(sum(sigma,dim=1) .ge. sum((2.5*sigma_md),dim=1))
         write(6,*) "MD KICK", md_kick
         if(md_kick) then
            call get_forces_md(0, 1)
         else
            forces_kick = 0
            sigma_md = (1-md_dt/md_tau)*sigma_md + md_dt/md_tau*sigma
         endif
       enddo
     endif
   
    ! computes  new accelerations
    do i = 1, ncent
      acc_new(:,i) = -forces_ave(:,i)/amass(i)
    enddo
  
    call fit_force(md_step,md_dt)

    vel = vel + 0.5 * (acc + acc_new) * md_dt

    acc = acc_new

    return
    
end subroutine
