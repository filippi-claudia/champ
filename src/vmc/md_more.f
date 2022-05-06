ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      subroutine init_mass(amass)
      subroutine init_mass()

      use vmc_mod
      use md_mass 
      use atom

       do k= 1,ncent
         amass(k) = 0
        do j= 1,nctype
         if(iwctype(k).eq.ntype(j)) then
c           write(6,*) "IWCTYPE, NTYPE", ntype(j)
c           write(6,*) "ASYMB(J) = ", asymb(j)
           call get_masses(asymb(j), amass(k))
         endif
        enddo
       enddo 
       
      return
      end


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c To have the updated positions and gradients of the enrgy withouth loading the common blocks in md.f90
      subroutine init_pos_md()

      use vmc_mod
      use force_mod
      use mstates_mod
      use sr_mod
      use md_mass
      use md_var,only:pos,forces_ave,sigma

      use atom
      use inputflags, only: node_cutoff, eps_node_cutoff
      use force_analy, only: iforce_analy
      use force_fin, only: da_energy_ave, da_energy_err
      use csfs, only : nstates 
      use optwf_contrl, only : ioptjas,ioptorb,ioptci,nparm
      use optwf_sr_ortho_mod, only: optwf_sr_ortho, optwf_sr_ortho_nogeo
      use config, only: anormo

      implicit real*8 (a-h,o-z)

      anormo=1.0d0
      write(6,*) "NSTATES IN INIT", nstates
 
     call set_nparms_tot
     call save_nparms

      ioptjas_sav=ioptjas; ioptorb_sav=ioptorb; ioptci_sav=ioptci
      iforce_analy=0; nstates_sav=nstates  
      node_cutoff_save = node_cutoff
 
      if (nstates.ge.2) then
        node_cutoff = 0
        call optwf_sr_ortho
        iguiding_sav=iguiding;      
        iroot_geo = 2
        write(6,*) "SELECTING ROOT", iroot_geo
        call select_wf_root(iroot_geo)
        nstates=1; iguiding=0
      else
        call optwf_sr_ortho
      endif

      node_cutoff = node_cutoff_save
      iforce_analy=1
      ioptjas=0; ioptorb=0; ioptci=0 
      call set_nparms

      call vmc
      call compute_position_bcast

      pos(:,:) = cent(:, :ncent)
      forces_ave(:,:) = da_energy_ave(:,:ncent,1)
      sigma(:,:) = da_energy_err(:,:ncent,1)

      if (nstates_sav .ge. 2) then 
         nstates=nstates_sav
         iguiding=iguiding_sav
         node_cutoff = 0
      endif

      ioptjas=ioptjas_sav; ioptorb=ioptorb_sav; ioptci=ioptci_sav
      call set_nparms

      write(6,*) "AFTER NSTATES", nstates
      if (nstates .ge. 2)  call restore_wf(1)

      return
      end


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Write md velocities

      subroutine write_md_vel(md_step,md_dt)

      use vmc_mod
      use atom, only: iwctype,ncent,nctype
      use md_mass
      use md_var, only: vel
      use mpiconf

c      implicit real*8 (a-h,o-z)
      real*8 md_dt


      if(wid) then
       write(15,*) ncent                         !number of atoms
       write(15,*) "time step", md_dt*md_step    !random comment line
       do i= 1, ncent
         do j= 1,nctype
          if(iwctype(i)== ntype(j)) then
          write(15,*) asymb(j), vel(:,i)       !write velocities in atomic units
         endif
        enddo
       enddo
      endif
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Write md forces

      subroutine write_md_force(md_step,md_dt)

      use vmc_mod
      use atom, only: iwctype,ncent,nctype
      use md_mass
      use md_var, only: forces_ave
      use mpiconf

c      implicit real*8 (a-h,o-z)
      real*8 md_dt


      if(wid) then
       write(75,*) ncent                         !number of atoms
       write(75,*) "time step", md_dt*md_step    !random comment line
       do i= 1, ncent
         do j= 1,nctype
          if(iwctype(i)== ntype(j)) then
          write(75,*) asymb(j), forces_ave(:,i)       !write Forces atomic units
         endif
        enddo
       enddo
      endif
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c To broadcast the new positions MD positions to get new gradients
      subroutine broadcas_md_pos()

      use vmc_mod
      use mpiconf
      use mpi
      use md_var, only:pos
      use atom, only: cent,ncent

c      implicit real*8 (a-h,o-z)
c      dimension positions(3,ncent)

       cent(:, :ncent) =  pos(:,:)
c       write(6,*) "CENT", cent
       call MPI_BCAST(cent,3*ncent,MPI_REAL8,0,MPI_COMM_WORLD,ier)
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Write md energy and time step
      subroutine write_md(md_step, md_dt, e_kin)

      use vmc_mod
      use force_mod
      use mstates_mod            

      use atom, only : ncent
      use sa_check 
      use mpiconf

      implicit real*8 (a-h,o-z)
      real*8 md_dt

      if(wid) then
        write (99,*) 'ekin', e_kin, 'epot', energy_all
        write (99,*) 'T', 2*e_kin/(3*3.1668105e-06*ncent)
        write (99,*) 'total ene', e_kin + energy_all, 'potential error', energy_err_all
      endif

      return
      end


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Write md geometries

      subroutine write_md_geo(md_step,md_dt)

      use vmc_mod

      use atom, only: iwctype, nctype, ncent
      use md_mass
      use mpiconf
      use md_var,only: pos

c      implicit real*8 (a-h,o-z)
      real*8 md_dt
c      dimension positions(3,ncent)

      if(wid) then
       write(13,*) ncent                         !number of atoms
       write(13,*) "time step", md_dt*md_step    !random comment line
       do i= 1, ncent
         do j= 1,nctype
          if(iwctype(i)== ntype(j)) then
          write(13,*) asymb(j), pos(:,i)*0.529177       !write positions in angstrom for vmd
         endif
        enddo
       enddo
      endif
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine get_forces_md( md_opt, imove)

      use vmc_mod
      use force_mod
      use mstates_mod            
      use md_mass
      use md_var, only: forces_ave,sigma

      use inputflags, only: node_cutoff, eps_node_cutoff
      use atom, only : ncent
      use force_fin, only: da_energy_ave, da_energy_err
      use force_analy, only: iforce_analy
      use sa_check 
      use mpiconf
      use csfs, only : nstates 
      use optwf_contrl, only : ioptjas,ioptorb,ioptci,nparm
      use optwf_sr_ortho_mod, only: optwf_sr_ortho, optwf_sr_ortho_nogeo
      use config, only: anormo

      implicit real*8 (a-h,o-z)

      dimension force_sig_write(6,ncent)

      ioptjas_sav=ioptjas; ioptorb_sav=ioptorb; ioptci_sav=ioptci
      iforce_analy=0; nstates_sav=nstates
      node_cutoff_save = node_cutoff

      anormo=1.0d0

      call set_nparms

      write(6,*) "START MDMORE, NSTATES =", nstates
      write(6,*) "md_opt", md_opt

      if (nstates.ge.2) then
        node_cutoff = 0
        iguiding_sav=iguiding        
        if(md_opt == 1) then
           call optwf_sr_ortho
        else
           call vmc
        endif
        iroot_geo = 2
        write(6,*) "SELECTING ROOT", iroot_geo
        call select_wf_root(iroot_geo)
        nstates=1; iguiding=0
      else
        if(md_opt == 1) call optwf_sr_ortho
      endif

      node_cutoff = node_cutoff_save
      iforce_analy=1
      ioptjas=0; ioptorb=0; ioptci=0 
      call set_nparms

      call vmc
      call compute_position_bcast

c some if it's tinker or md we just take forces and sigma, if it's steepest descent it moves
      if (imove == 2)then
        call compute_positions      
      else
         forces_ave(:,:) = da_energy_ave( :,:ncent, 1)
         sigma(:,:) = da_energy_err(:,:ncent,1)
         if(imove==0) then
           force_sig_write(:3,:)= -forces_ave(:,:)
           force_sig_write(4:,:)= sigma(:,:)
           if(wid) write(75,'(1p6e14.5)') (force_sig_write) 
         endif
      endif

      if (nstates_sav .ge. 2)then 
         nstates=nstates_sav
         iguiding=iguiding_sav
         node_cutoff = 0
      endif

      ioptjas=ioptjas_sav; ioptorb=ioptorb_sav; ioptci=ioptci_sav
      call set_nparms

      if (nstates .ge. 2) call restore_wf(1)
      write(6,*) "END OF MDMORE WITH #STATES:", nstates

      return
      end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
