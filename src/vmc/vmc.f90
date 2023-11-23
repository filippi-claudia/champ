module vmc_f_mod
contains
      subroutine vmc
! Written by Cyrus Umrigar and Claudia Filippi

! Program to do variational Monte Carlo calculations
! on atoms and molecules.
! Various types of Metropolis moves can be done, including a few
! versions of directed Metropolis in spherical polar coordinates.
! Also, one or all electrons can be moved at once.
! Currently this program contains
! 1s, 2s, 2p, 3s, 3p, 3d, 4s,  and 4p  Slater basis states.
! and sa, pa, da asymptotic functions

      use acuest_mod, only: acuest,zerest
      use coefs,   only: nbasis
      use config,  only: eold,psido,psijo,xold
      use contrl_file, only: ounit
      use control_vmc, only: vmc_idump,vmc_irstar,vmc_nblk,vmc_nblkeq
      use control_vmc, only: vmc_nconf,vmc_nconf_new,vmc_nstep
      use dumper_mod, only: dumper,startr
      use error,   only: fatal_error
      use finwrt_mod, only: finwrt
      use mc_configs, only: mc_configs_start,mc_configs_write
      use metrop_mov1_slat, only: metrop6
      use mpitimer, only: elapsed_time
      use multiple_geo, only: iwftype,nforce,nwftype
      use precision_kinds, only: dp
      use pseudo,  only: nloc
      use rotqua_mod, only: rotqua
      use slater,  only: coef
      use strech_mod, only: setup_force
      use system,  only: nelec
!      use contrl, only: idump, irstar, nconf, nblk, nblkeq, nconf_new, nstep

      implicit none

      integer :: i, ii, j, jj, l
      integer :: ngfmc
      real(dp) ::err




      character(len=25) fmt

! common block variables:

!   /const/
!        nelec  = number of electrons
!        pi     = 3.14159...
!        hb     = hbar**2/(2m)
!        delta  = side of box in which metropolis steps are made
!        deltai = 1/delta
!        fbias  = force bias parameter
!   /contrl/
!        nstep  = number of metropolis steps/block
!        nblk   = number of blocks od nstep steps after the
!                equilibrium steps
!        nblkeq = number of equilibrium blocks
!        nconf  = target number of mc configurations (dmc only)
!        nconf_new = number of mc configurations generated for optim and dmc
!        idump  =  1 dump out stuff for a restart
!        irstar =  1 pick up stuff for a restart
!   /config/
!        xold   = current position of the electrons
!        xnew   = new position after a trial move
!        vold   = grad(psi)/psi at current position
!        vnew   = same after trial move
!        psi2o  = psi**2 at current position
!        psi2n  = same after trial move
!        eold   = local energy at current position
!        enew   = same after trial move
!        ekino   = local kinetic energy at current position
!        psido  = determinantal part of wave function
!        psijo  = log(Jastrow)
!   /coefs/
!        coef   = read in coefficients of the basis functions
!                 to get the molecular orbitals used in determinant
!        nbasis = number of basis functions read in
!   /dets/
!        cdet   = coefficients of the determinants
!        ndet   = number of determinants of molecular orbitals
!                 used
!        nup    = number of up spin electrons
!        ndn    = number of down spin electrons

      if(nforce.gt.1) then
! force parameters
        call setup_force
       else
        nwftype=1
        iwftype(1)=1
      endif

! initialize the walker configuration
      call mc_configs_start
      if (vmc_nconf_new.eq.0) then
        ngfmc=2*vmc_nstep*vmc_nblk
       else
        ngfmc=(vmc_nstep*vmc_nblk+vmc_nconf_new-1)/vmc_nconf_new
      endif
      call elapsed_time("VMC : initial walker configuration : ")

! zero out estimators and averages

      if (vmc_irstar.ne.1) then
            call zerest
            call elapsed_time("VMC : zero out estimators and averages : ")
      endif


! check if restart flag is on. If so then read input from
! dumped data to restart

      if (vmc_irstar.eq.1) then
        open(10,err=401,form='unformatted',file='restart_vmc')
        goto 402
      401   call fatal_error('VMC: restart_vmc empty, not able to restart')
      402   rewind 10
        call startr
        close(10)
        call elapsed_time("VMC : reading restart files : ")
      endif

! if there are equilibrium steps to take, do them here
! skip equilibrium steps if restart run
! imetro = 6 spherical-polar with slater T
      if (vmc_nblkeq.ge.1.and.vmc_irstar.ne.1) then
        l=0
        do i=1,vmc_nblkeq
          do j=1,vmc_nstep
            l=l+1
            if (nloc.gt.0) call rotqua
            call metrop6(l,0)
          enddo

         call acuest
        enddo

!       Equilibration steps done. Zero out estimators again.

        call elapsed_time("VMC : equilibrium CP : ")

        call zerest
      endif
! now do averaging steps

      l=0
      do i=1,vmc_nblk
        do j=1,vmc_nstep
        l=l+1
!   write(ounit, *) i, nblk, j, nstep
        if (nloc.gt.0) call rotqua
        call metrop6(l,1)
! write out configuration for optimization/dmc/gfmc here
        if (mod(l,ngfmc).eq.0 .or. ngfmc.eq.1) then
          if(3*nelec.lt.100) then
           write(fmt,'(a1,i2,a21)')'(',3*nelec,'f13.8,i3,d12.4,f12.5)'
          else
           write(fmt,'(a1,i3,a21)')'(',3*nelec,'f13.8,i3,d12.4,f12.5)'
          endif
          write(7,fmt) ((xold(ii,jj),ii=1,3),jj=1,nelec), &
          int(sign(1.d0,psido(1))),log(dabs(psido(1)))+psijo,eold(1,1)
        endif
        enddo

      call acuest
      enddo

      call elapsed_time("VMC : all CP : ")

! write out last configuration to mc_configs_start
! call fin_reduce to write out additional files for efpci, embedding etc.
! collected over all the run and to reduce cum1 in mpi version
      call mc_configs_write
      call elapsed_time("writing configs to a file : ")

! print out final results
      call finwrt
      call elapsed_time("writing final results : ")

! if dump flag is on then dump out data for a restart
      if (vmc_idump.eq.1) then
        open(10,form='unformatted',file='restart_vmc')
        rewind 10
        call dumper
        close(10)
        call elapsed_time("dumping restart files : ")
      endif
      if(vmc_nconf_new.ne.0) close(7)

      return
      end
end module
