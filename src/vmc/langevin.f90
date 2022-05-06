Subroutine langevin(md_step, nopt, md_T)
     use atom, only: ncent
     use md_mass,only: amass
     use md_var
     use precision_kinds, only: dp

  implicit none

     logical :: md_kick
     integer :: i, j, md_step
     integer :: forces_kick, nopt, md_opt
     real(kind=dp), dimension(:,:), allocatable:: gamman ,alpha,beta
     real(kind=dp):: md_T
    
!     if (md_method .ne. 'langevin') return

     allocate(alpha(3,ncent), beta(3,ncent), gamman(3,ncent))

     call broadcas_md_pos()
     md_opt = 0
     if(mod(md_step, nopt) == 0.0d0) md_opt =1
     call get_forces_md(md_opt, 1)                 !imove is 1 becouse we are doing dynamics not steepest descent

     !check if there are kicks in the forces
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


    !set langevin parameters 
    gamman = md_dt*sigma_md**2/(2*3.1668105e-06*md_T)
    do i = 1, ncent
      alpha(:,i) = (1.0/gamman(:,i))*(1.0 - exp(-gamman(:,i)*md_dt/amass(i)))
      beta(:,i)  = exp(-gamman(:,i)*md_dt/amass(i))
    enddo

    vel = beta * vel + alpha *(-forces_ave)
    ! compute new positions
    pos = pos + vel * md_dt
    call broadcas_md_pos()

    deallocate(gamman, alpha, beta)
    return
end subroutine
