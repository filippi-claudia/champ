Subroutine mix_memory(md_step, nopt, md_T)
  use atom, only: iwctype, nctype, ncent
  use md_mass,only: ntype,amass
  use md_var
  use md_fit, only: nfit, wnoi
  use precision_kinds, only: dp

  implicit none

    logical :: md_kick
    integer :: k, j, md_step
    integer :: forces_kick, nopt, md_opt
    real(dp):: dt_2, rand,temp
    real(dp):: md_T, kb, ekinoi
    real(kind=dp), dimension(:), allocatable:: tau, tf, eta
    real(kind=dp), dimension(3):: tau_md
    real(kind=dp), dimension(:,:), allocatable:: ynoi
    real(kind=dp), dimension(:,:), allocatable:: Gau_y, Gau_w
    real(kind=dp), dimension(:,:), allocatable:: omeg, opom, omom
    real(kind=dp), dimension(:,:), allocatable:: sinf, cosf
    real(kind=dp), dimension(:,:), allocatable:: somf, somg
    real(kind=dp), dimension(:,:), allocatable:: syf, syg

    kb = 3.1668105e-06
    dt_2 = 0.5d0 * md_dt

    allocate (ynoi(3,ncent))
    allocate (Gau_y(3,ncent),Gau_w(3,ncent))
    allocate (tau(ncent), tf(ncent), eta(ncent))
    allocate (omeg(3,ncent), opom(3,ncent), omom(3,ncent))
    allocate (sinf(3,ncent), cosf(3,ncent))
    allocate (somf(3,ncent), somg(3,ncent))
    allocate (syf(3,ncent), syg(3,ncent))


    if (md_step .eq. 1) sigma_md = sigma


    tau_md(1) = 0.83d0
    tau_md(2) = 1.0d0
    tau_md(3) = 0.67d0

    do k = 1, ncent
      do j= 1,nctype
         if(iwctype(k).eq.ntype(j)) then
         tau(k) = tau_md(j) *md_dt
         endif
      enddo
      tf(k)  = md_dt/tau(k)
      eta(k) = exp(-tf(k))
      do j = 1,3
        call gasdev_s(rand)
        Gau_y(j,k) = rand
        call gasdev_s(rand)
        Gau_w(j,k) = rand
        temp = (sigma_md(j,k)*tau(k))**2/kb/md_T/amass(k)
      enddo
    enddo

    if (temp .le. 1.0d0 ) then 

    do k = 1, ncent
      do j = 1,3
        omeg(j,k) = sqrt(1.d0 -((sigma_md(j,k)*tau(k))**2/kb/md_T/amass(k)))
        opom(j,k)  = 1.d0 + omeg(j,k)**2
        omom(j,k)  = 1.d0 - omeg(j,k)**2         

        somf(j,k) = eta(k)*sigma_md(j,k)*tau(k)/amass(k)/omeg(j,k)/sqrt(2.d0)
        somg(j,k) = sqrt((omeg(j,k)**2-1.d0-2.d0*exp(2.d0*tf(k))*omeg(j,k)**2  &
                  + opom(j,k)*cosh(2.d0*omeg(j,k)*tf(k))           &
                  +  2.d0*omeg(j,k)*sinh(2.d0*omeg(j,k)*tf(k)))/(omeg(j,k)**2-1.d0))

        cosf(j,k) = cosh(omeg(j,k)*tf(k))
        sinf(j,k) = sinh(omeg(j,k)*tf(k))

        syf(j,k) = sigma_md(j,k)/amass(k)*sqrt((tau(k)/omom(j,k))**3)
        syg(j,k) = sqrt(4.d0*md_dt*omom(j,k) - 2.d0*tau(k)*(6.d0-omom(j,k)) -tau(k)*(eta(k)*(2.d0*cosf(j,k) &
             + (1.d0+omeg(j,k)**2)/omeg(j,k)*sinf(j,k)))**2 +eta(k)*2.d0*tau(k)*((8.d0-omom(j,k))*cosf(j,k) &
             + (8.d0-5.d0*omom(j,k))/omeg(j,k)*sinf(j,k)))
        if (md_step == 1)  write (6,*) 'Omega', omeg(j,k)
      enddo
    enddo
    else
    do k = 1, ncent
      do j = 1,3
        omeg(j,k) = sqrt(((sigma_md(j,k)*tau(k))**2/kb/md_T/amass(k))-1.d0 )
        opom(j,k)  = 1.d0 + omeg(j,k)**2
        omom(j,k)  = 1.d0 - omeg(j,k)**2

        somf(j,k) =eta(k)*sigma_md(j,k)*tau(k)/amass(k)/omeg(j,k)/sqrt(2.d0*opom(j,k))
        somg(j,k) = sqrt(-opom(j,k) + 2.d0*exp(2.d0*tf(k))*omeg(j,k)**2 &
                  + omom(j,k)*cos(2.d0*omeg(j,k)*tf(k))- 2.d0*omeg(j,k)*sin(2.d0*omeg(j,k)*tf(k)))

        cosf(j,k) = cos(omeg(j,k)*tf(k))
        sinf(j,k) = sin(omeg(j,k)*tf(k))

        syf(j,k) = sigma_md(j,k)/amass(k)*sqrt((tau(k)/opom(j,k))**3)
        syg(j,k) = sqrt(4.d0*md_dt*opom(j,k) - 2.d0*tau(k)*(6.d0-opom(j,k))           &
                 -tau(k)*(eta(k)*(2.d0*cosf(j,k) + sinf(j,k)*omom(j,k)/omeg(j,k)))**2 &
                 +  eta(k)*2.d0*tau(k)*((8.d0-opom(j,k))*cosf(j,k) + (5.d0*omom(j,k)-2.d0)/omeg(j,k)*sinf(j,k)))
        if (md_step == 1)  write (6,*) 'Omega', omeg(j,k)
       enddo
     enddo
     endif


      ynoi = 0.0d0
      if(md_step .le. nfit) then
        do k = 1, ncent
           do j = 1, 3
             vel(j,k) = vel(j,k) + acc(j,k)*dt_2
             wnoi(j,k) = 0.0d0
             pos(j,k) = pos(j,k) + vel(j,k)*md_dt
           enddo
        enddo
        ekinoi = 0.0d0
      else
        ekinoi = 0.0d0
        do k = 1, ncent
          do j = 1, 3
           vel(j,k) = vel(j,k) + acc(j,k) * dt_2
           if (temp .le. 1.0d0 ) then
             ynoi(j,k) = tau(k)/omom(j,k)*(2.d0-eta(k)*(2.d0*cosf(j,k)+opom(j,k)*sinf(j,k)/omeg(j,k)))
             ynoi(j,K) = ynoi(j,k)*wnoi(j,k)+ syf(j,k)*syg(j,k)*Gau_y(j,k)
             wnoi(j,k) = wnoi(j,k)*eta(k)*(cosf(j,k)+sinf(j,k)/omeg(j,k)) + somf(j,k)*somg(j,k)*Gau_w(j,k)
           else
             ynoi(j,k) = tau(k)/opom(j,k)*(2.d0-eta(k)*(2.d0*cosf(j,k) + omom(j,k)*sinf(j,k)/omeg(j,k)))
             ynoi(j,K) = ynoi(j,k)*wnoi(j,k)+ syf(j,k)*syg(j,k)*Gau_y(j,k)
             wnoi(j,k) = wnoi(j,k)*eta(k)*(cosf(j,k)+sinf(j,k)/omeg(j,k)) + somf(j,k)*somg(j,k)*Gau_w(j,k)
           endif
           ekinoi = ekinoi + 0.50d0*amass(k)*wnoi(j,k)*wnoi(j,k)
          enddo
        enddo
        pos = pos + vel*md_dt + ynoi
      endif

! Initialize Variables
    write(6,*) "HERE  1"
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


    deallocate (ynoi)
    deallocate (Gau_y,Gau_w)
    deallocate (tf, eta)
    deallocate (omeg, opom, omom)
    deallocate (sinf, cosf)
    deallocate (somf, somg)
    deallocate (syf, syg)
    return
    
end subroutine
