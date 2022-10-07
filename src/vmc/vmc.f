      module vmc_f_mod
      contains
      subroutine vmc
c Written by Cyrus Umrigar and Claudia Filippi

c Program to do variational Monte Carlo calculations
c on atoms and molecules.
c Various types of Metropolis moves can be done, including a few
c versions of directed Metropolis in spherical polar coordinates.
c Also, one or all electrons can be moved at once.
c Currently this program contains
c 1s, 2s, 2p, 3s, 3p, 3d, 4s,  and 4p  Slater basis states.
c and sa, pa, da asymptotic functions

      use system, only: nelec
      use config, only: eold
      use config, only: psido, psijo
      use config, only: xold
      use forcepar, only: nforce
      use wfsec, only: iwftype, nwftype
      use coefs, only: coef, nbasis
!      use contrl, only: idump, irstar, nconf, nblk, nblkeq, nconf_new, nstep
      use control_vmc, only: vmc_idump, vmc_irstar, vmc_nconf, vmc_nblk
      use control_vmc, only: vmc_nblkeq, vmc_nconf_new, vmc_nstep
      use pseudo, only: nloc
      use mpitimer,    only: elapsed_time
      use contrl_file,    only: ounit

      use precision_kinds, only: dp
      use strech_mod, only: setup_force
      use dumper_mod, only: dumper, startr
      use mc_configs, only: mc_configs_start, mc_configs_write
      use acuest_mod, only: zerest
      use error, only: fatal_error
      use finwrt_mod, only: finwrt
      use rotqua_mod, only: rotqua
      use acuest_mod, only: acuest
      use metrop_mov1_slat, only: metrop6
      implicit none

      integer :: i, ii, j, jj, l
      integer :: ngfmc
      real(dp) ::err




      character*25 fmt

c common block variables:

c   /const/
c        nelec  = number of electrons
c        pi     = 3.14159...
c        hb     = hbar**2/(2m)
c        delta  = side of box in which metropolis steps are made
c        deltai = 1/delta
c        fbias  = force bias parameter
c   /contrl/
c        nstep  = number of metropolis steps/block
c        nblk   = number of blocks od nstep steps after the
c                equilibrium steps
c        nblkeq = number of equilibrium blocks
c        nconf  = target number of mc configurations (dmc only)
c        nconf_new = number of mc configurations generated for optim and dmc
c        idump  =  1 dump out stuff for a restart
c        irstar =  1 pick up stuff for a restart
c   /config/
c        xold   = current position of the electrons
c        xnew   = new position after a trial move
c        vold   = grad(psi)/psi at current position
c        vnew   = same after trial move
c        psi2o  = psi**2 at current position
c        psi2n  = same after trial move
c        eold   = local energy at current position
c        enew   = same after trial move
c        peo    = local potential at current position
c        pen    = same after trial move
c        tjfo   = Jackson Feenberg kinetic energy at current position
c        tjfn   = same after trial move
c        psido  = determinantal part of wave function
c        psijo  = log(Jastrow)
c   /coefs/
c        coef   = read in coefficients of the basis functions
c                 to get the molecular orbitals used in determinant
c        nbasis = number of basis functions read in
c   /dets/
c        cdet   = coefficients of the determinants
c        ndet   = number of determinants of molecular orbitals
c                 used
c        nup    = number of up spin electrons
c        ndn    = number of down spin electrons

      if(nforce.gt.1) then
c force parameters
        call setup_force
       else
        nwftype=1
        iwftype(1)=1
      endif

c initialize the walker configuration
      call mc_configs_start
      if (vmc_nconf_new.eq.0) then
        ngfmc=2*vmc_nstep*vmc_nblk
       else
        ngfmc=(vmc_nstep*vmc_nblk+vmc_nconf_new-1)/vmc_nconf_new
      endif
      call elapsed_time("VMC : initial walker configuration : ")

c zero out estimators and averages

      if (vmc_irstar.ne.1) then
            call zerest
            call elapsed_time("VMC : zero out estimators and averages : ")
      endif


c check if restart flag is on. If so then read input from
c dumped data to restart

      if (vmc_irstar.eq.1) then
        open(10,err=401,form='unformatted',file='restart_vmc')
        goto 402
  401   call fatal_error('VMC: restart_vmc empty, not able to restart')
  402   rewind 10
        call startr
        close(10)
        call elapsed_time("VMC : reading restart files : ")
      endif

c if there are equilibrium steps to take, do them here
c skip equilibrium steps if restart run
c imetro = 6 spherical-polar with slater T
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

c       Equilibration steps done. Zero out estimators again.

        call elapsed_time("VMC : equilibrium CP : ")

        call zerest
      endif
c now do averaging steps

      l=0
      do i=1,vmc_nblk
        do j=1,vmc_nstep
        l=l+1
      !   write(ounit, *) i, nblk, j, nstep
        if (nloc.gt.0) call rotqua
        call metrop6(l,1)
c write out configuration for optimization/dmc/gfmc here
        if (mod(l,ngfmc).eq.0 .or. ngfmc.eq.1) then
          if(3*nelec.lt.100) then
           write(fmt,'(a1,i2,a21)')'(',3*nelec,'f13.8,i3,d12.4,f12.5)'
          else
           write(fmt,'(a1,i3,a21)')'(',3*nelec,'f13.8,i3,d12.4,f12.5)'
          endif
          write(7,fmt) ((xold(ii,jj),ii=1,3),jj=1,nelec),
     &    int(sign(1.d0,psido(1))),log(dabs(psido(1)))+psijo,eold(1,1)
        endif
        enddo

      call acuest
      enddo

      call elapsed_time("VMC : all CP : ")

c write out last configuration to mc_configs_start
c call fin_reduce to write out additional files for efpci, embedding etc.
c collected over all the run and to reduce cum1 in mpi version
      call mc_configs_write
      call elapsed_time("writing configs to a file : ")

c print out final results
      call finwrt
      call elapsed_time("writing final results : ")

c if dump flag is on then dump out data for a restart
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
