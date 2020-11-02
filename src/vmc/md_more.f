ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine init_mass(amass)

      implicit real*8 (a-h,o-z)
      include 'vmc.h'
      character*5 asymb
      common /mass_type/ ntype(MCTYPE)
      common/mass_symb/ asymb(MCTYPE)
      common /atom/ znuc(MCTYPE),cent(3,MCENT),pecent, iwctype(MCENT),nctype,ncent
 
      dimension amass(ncent)

       do k= 1,ncent
        atom_mass= 0
        do j= 1,nctype
         if(iwctype(k).eq.ntype(j)) then
           call get_masses(asymb(j), atom_mass)
           amass(k)= atom_mass      
         endif
        enddo
       enddo 
       
      return
      end


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c To have the updated positions and gradients of the enrgy withouth loading the common blocks in md.f90
      subroutine init_pos_md(positions,forces_average)

      implicit real*8 (a-h,o-z)
      include 'vmc.h'
      include 'force.h'
      include 'mstates.h'
      include 'sr.h'

      common /atom/ znuc(MCTYPE),cent(3,MCENT),pecent,iwctype(MCENT),nctype,ncent
      common /force_fin/ da_energy_ave(3,MCENT),da_energy_err(3,MCENT)
      dimension positions(3,ncent)
      dimension forces_average(3,ncent)

      positions(:,:) = cent(:, :ncent)
      forces_average(:,:) = da_energy_ave(:,:ncent)
      return
      end


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c To broadcast the new positions MD positions to get new gradients
      subroutine broadcas_md_pos(positions)
      implicit real*8 (a-h,o-z)
      include 'vmc.h'
      include 'mpif.h'

      common /mpiconf/ idtask,nproc,wid
      common /atom/ znuc(MCTYPE),cent(3,MCENT),pecent,iwctype(MCENT),nctype,ncent
      dimension positions(3,ncent)

       cent(:, :ncent) =  positions(:,:)
       call MPI_BCAST(cent,3*ncent,MPI_REAL8,0,MPI_COMM_WORLD,ier)
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Write md energy and time step
      subroutine write_md(md_step, md_dt, e_kin)

      implicit real*8 (a-h,o-z)
      include 'vmc.h'
      include 'mpif.h'
      include 'force.h'
      include 'mstates.h'            

      common /atom/ znuc(MCTYPE),cent(3,MCENT),pecent, iwctype(MCENT),nctype,ncent
c      common /optwf_corsam/ add_diag(MFORCE),energy(MFORCE),energy_err(MFORCE),
c     & force(MFORCE),force_err(MFORCE),sigma
      common /sa_check/ energy_all(MSTATES),energy_err_all(MSTATES) 
      real*8 md_dt
      common /mpiconf/ idtask,nproc,wid
      logical wid

      dimension amass(ncent)

      if(wid) then
        write (99,*) 'ekin', e_kin, 'epot', energy_all
        write(6,*) "epot", energy_all
        write (99,*)  'T', 2*e_kin/(3*3.1668105e-06*ncent)
        write (99,*) 'total ene', e_kin + energy_all, 'potential error', energy_err_all(1)
        write (99,*) 'Ene per atom', (e_kin + energy_all(1))/ncent
      endif
      write (99,*) 'End of molueclar dynamics step'

      return
      end


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Write md geometries

      subroutine write_md_geo(md_step,md_dt, positions)

      implicit real*8 (a-h,o-z)
      real*8 md_dt
      include 'mpif.h'
      include 'vmc.h'
      character*5 asymb

      common /atom/ znuc(MCTYPE),cent(3,MCENT),pecent, iwctype(MCENT),nctype,ncent
      common /mass_type/ ntype(MCTYPE)
      common/mass_symb/ asymb(MCTYPE)
      common /mpiconf/ idtask,nproc,wid
      logical wid
      dimension positions(3,ncent)

      if(wid) then
       write(13,*) ncent                         !number of atoms
       write(13,*) "time step", md_dt*md_step    !random comment line
       do i= 1, ncent
         do j= 1,nctype
          if(iwctype(i)== ntype(j)) then
          write(13,*) asymb(j), positions(:,i)*0.529177       !write positions in angstrom for vmd
         endif
        enddo
       enddo
      endif
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine get_forces_md(forces_average, sigma, imove)

      implicit real*8 (a-h,o-z)
      include 'vmc.h'
      include 'mpif.h'
      include 'force.h'
      include 'mstates.h'            

      common /atom/ znuc(MCTYPE),cent(3,MCENT),pecent,iwctype(MCENT),nctype,ncent
      common /force_fin/ da_energy_ave(3,MCENT),da_energy_err(3,MCENT)
      common /force_analy/ iforce_analy
      common /mpiconf/ idtask,nproc,wid
      common /csfs/ ccsf(MDET,MSTATES,MWF),cxdet(MDET*MDETCSFX)
     &,icxdet(MDET*MDETCSFX),iadet(MDET),ibdet(MDET),ncsf,nstates
      common /optwf_contrl/ ioptjas,ioptorb,ioptci,nparm
      dimension forces_average(3,ncent), sigma(3,ncent)
      dimension force_sig_write(6,ncent)
      logical wid


c     iforce_analy is set to zero in the SR

      ioptjas_sav=ioptjas; ioptorb_sav=ioptorb; ioptci_sav=ioptci
      write(6,*) "START MDMORE", nstates
       
      ifroce_analy=0
      if (nstates.ge.2) then
        call optwf_mix_nogeo
        iroot_geo = 2
        write(6,*) "SELECTING ROOT", iroot_geo
        call p2gtid('optgeo:iroot_geo',iroot_geo,0,0)
        call select_ci_root(iroot_geo)
        iguiding_sav=iguiding; nstates_sav=nstates        
        nstates=1; iguiding=0
      else
        call optwf_sr_nogeo
      endif

      call p2gtid('vmc:node_cutoff',node_cutoff,0,1)
      write(6,*) "NODE CUTOFF", node_cutoff

      iforce_analy=1
      ioptjas=0; ioptorb=0; ioptci=0 
      call set_nparms

      call vmc

c some if it's tinker or md we just take forces and sigma, if it's steepest descent it moves
      if (imove == 2)then
        call compute_positions      
      else
         forces_average(:,:) = da_energy_ave(:,:ncent)
         sigma(:,:) = da_energy_err(:,:ncent)
         if(imove==0) then
           force_sig_write(:3,:)= -forces_average(:,:)
           force_sig_write(4:,:)= sigma(:,:)
           if(wid) write(75,'(1p6e14.5)') (force_sig_write) 
         endif
      endif

      if (nstates_sav .ge. 2)then 
         nstates=nstates_sav
         iguiding=iguiding_sav
      endif

      ioptjas=ioptjas_sav; ioptorb=ioptorb_sav; ioptci=ioptci_sav
      call set_nparms
      call restore_wf(1)
      write(6,*) "MDMORE", nstates

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
