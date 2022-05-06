Subroutine mix_verlang(md_step, nopt, md_T)
  use atom, only: ncent
  use md_mass,only: amass
  use md_var
  use md_fit, only: nfit, wnoi
  use precision_kinds, only: dp

  implicit none

    logical :: md_kick
    integer :: k, j, md_step
    integer :: forces_kick, nopt, md_opt
    real(dp) :: dt_2, Treal,Tgam,rand
    real(kind=dp):: md_T,kb
    real(kind=dp), dimension(:,:), allocatable:: ynoi
    real(kind=dp), dimension(:,:), allocatable:: Gau_y, Gau_w
    real(kind=dp), dimension(:,:), allocatable:: gamman ,alpha,beta
    real(kind=dp), dimension(:,:), allocatable:: clan,eta
    
    kb = 3.1668105e-06
    dt_2 = 0.5d0 * md_dt
    allocate (ynoi(3,ncent))
    allocate(alpha(3,ncent), beta(3,ncent), gamman(3,ncent))
    allocate (clan(3,ncent), eta(3,ncent))
    allocate (Gau_y(3,ncent),Gau_w(3,ncent))

    Tgam  = md_T
    Treal = md_T 

    if (md_step .eq. 1) sigma_md = sigma
! Initialize Variables
    do k = 1, ncent
      do j = 1,3
        gamman(j,k) = md_dt*sigma_md(j,k)**2/2.0d0/kb/md_T/amass(k)
        eta(j,k)   = exp(-gamman(j,k)*md_dt)
        alpha(j,k) = (1.0d0 - eta(j,k))/gamman(j,k)
        beta(j,k)  = (1.0d0 - exp(-2.0d0*gamman(j,k)*md_dt))/gamman(j,k)
        clan(j,k)  = sqrt(2.0d0*kb*Treal/gamman(j,k)/amass(k))
        call gasdev_s(rand)
        Gau_y(j,k) = rand
        call gasdev_s(rand)
        Gau_w(j,k) = rand
      enddo
    enddo
    if (md_step == 1) write (6,*) "xi/m dt", gamman*md_dt

    if(md_step .le. nfit) then
      do k = 1, ncent
         do j = 1, 3
           vel(j,k) = vel(j,k) + acc(j,k)*dt_2
           wnoi(j,k) = 0.0d0
           pos(j,k) = pos(j,k) + vel(j,k)*md_dt
         enddo
      enddo
    else
    do k = 1, ncent
       do j = 1, 3
        vel(j,k) = vel(j,k) + acc(j,k) * dt_2

        ynoi(j,k) = wnoi(j,k)*alpha(j,k) + clan(j,k)*sqrt(md_dt-2.0d0*alpha(j,k)+0.5d0*beta(j,k)) &
             + clan(j,k)*sqrt(md_dt-2.0d0*alpha(j,k)+0.5d0*beta(j,k))*Gau_y(j,k)

        wnoi(j,k) = wnoi(j,k)*eta(j,k) + sqrt(0.5d0)*gamman(j,k)*clan(j,k)*sqrt(beta(j,k))  &
               + sqrt(0.5d0)*gamman(j,k)*clan(j,k)*sqrt(beta(j,k))*Gau_w(j,k)
        enddo
      enddo
      pos = pos + vel*md_dt + ynoi
    endif
    call broadcas_md_pos()
   
    !compute forces atnew positions
    md_opt = 0
    if(mod(md_step, nopt) == 0.0d0) md_opt =1
    call get_forces_md(md_opt, 1)                 !imove is 1 becouse we are not doing steepest descent
   
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
    do k = 1, ncent
      acc_new(:,k) = -forces_ave(:,k)/amass(k)
    enddo

    call fit_force(md_step,md_dt)

    do k = 1, ncent
       do j = 1, 3
          vel(j,k) = vel(j,k) + acc_new(j,k) * dt_2   
       end do
    end do

    acc = acc_new

!   perform deallocation of some local arrays

    deallocate (gamman,alpha, beta)
    deallocate (clan,eta)
    deallocate (Gau_y,Gau_w)

    return
    
end subroutine

