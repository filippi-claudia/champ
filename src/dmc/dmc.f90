module dmc_f_mod
contains
      subroutine dmc
! Written by Cyrus Umrigar with major contributions by Claudia Filippi.
! Uses the diffusion Monte Carlo algorithm described in:
! 1) A Diffusion Monte Carlo Algorithm with Very Small Time-Step Errors,
!    C.J. Umrigar, M.P. Nightingale and K.J. Runge, J. Chem. Phys., 99, 2865 (1993).

      use acues1_mod, only: acues1
      use acues1_reduce_mod, only: acues1_reduce
      use acuest_mod, only: acuest
      use averages, only: average,average_write,init_averages_index
      use constants, only: pi
      use contrl_file, only: ounit
      use contrldmc, only: idmc
      use control_dmc, only: dmc_idump,dmc_irstar,dmc_nblk,dmc_nblkeq
      use control_dmc, only: dmc_nconf,dmc_nstep
      use dmc_ps_mov1, only: dmc_ps
      use dumper_mod, only: dumper
      use error,   only: fatal_error
      use estcum,  only: ipass
      use finwrt_mod, only: finwrt
      use init_mod, only: init
      use mc_configs_mod, only: mc_configs,mc_configs_write
      use mpitimer, only: elapsed_time
      use multiple_geo, only: iwftype,nforce,nwftype,nwprod
      use precision_kinds, only: dp
      use pseudo,  only: nloc
      use rotqua_mod, only: rotqua
      use strech_mod, only: setup_force
      use zerest_mod, only: zerest
!      use contrl, only: idump, irstar, nblk, nblkeq, nconf, nstep


      implicit none

      integer :: i, j, irun, lpass
      real(dp), parameter :: one = 1.d0
      real(dp), parameter :: four = 4.d0


! variables:
!        nelec  = number of electrons
!        pi     = 3.14159...
!        hb     = hbar**2/(2m)
!        delta  = side of box in which metropolis steps are made
!        deltai = 1/delta
!        fbias  = force bias parameter
!        nstep  = number of metropolis steps/block
!        nblk   = number of blocks od nstep steps after the
!                 equilibrium steps
!        nblkeq = number of equilibrium blocks
!        nconf  = initial and target number of dmc configurations
!        idump  =  1 dump out stuff for a restart
!        irstar =  1 pick up stuff for a restart
!        tau    = time-step
!        nfprod = number of f's used for removing finite popul. bias
! Control variables are:
! idmc         < 0     VMC
!              > 0     DMC
! abs(idmc)    = 1 **  simple kernel using dmc.brock.f
!              = 2     good kernel using dmc_good or dmc_good_inhom
! ipq         <= 0 *   do not use expected averages
!             >= 1     use expected averages (mostly in all-electron move algorithm)
! itau_eff    <=-1 *   always use tau in branching (not implemented)
!              = 0     use 0 / tau for acc /nacc moves in branching
!             >= 1     use tau_eff (calcul in equilibration runs) for all moves
! iacc_rej    <=-1 **  accept all moves (except possibly node crossings)
!              = 0 **  use weights rather than accept/reject
!             >= 1     use accept/reject
! icross      <=-1 **  kill walkers that cross nodes (not implemented)
!              = 0     reject walkers that cross nodes
!             >= 1     allow walkers to cross nodes
!                      (OK since crossing prob. goes as tau^(3/2))
! icuspg      <= 0     approximate cusp in Green function
!             >= 1     impose correct cusp condition on Green function
! icut_br     <= 0     do not limit branching
!             >= 1 *   use smooth formulae to limit branching to (1/2,2)
!                      (bad because it makes energies depend on E_trial)
! icut_e      <= 0     do not limit energy
!             >= 1 *   use smooth formulae to limit energy (not implemented)
! *  => bad option, modest deterioration in efficiency or time-step error
! ** => very bad option, big deterioration in efficiency or time-step error
! So, idmc=6,66 correspond to the foll. two:
! 2 1 1 1 0 0 0 0 0  idmc,ipq,itau_eff,iacc_rej,icross,icuspg,idiv_v,icut_br,icut_e
! 2 1 0 1 1 0 0 0 0  idmc,ipq,itau_eff,iacc_rej,icross,icuspg,idiv_v,icut_br,icut_e
! Another reasonable choice is:
! 2 1 0 1 1 1 1 0 0  idmc,ipq,itau_eff,iacc_rej,icross,icuspg,idiv_v,icut_br,icut_e
!        xold_dmc   = current position of the electrons
!        xnew   = new position after a trial move
!        vold_dmc   = grad(psi)/psi at current position
!        vnew   = same after trial move
!        psi2o  = psi**2 at current position
!        psi2n  = same after trial move
!        eold   = local energy at current position
!        enew   = same after trial move
!        peo_dmc    = local potential at current position
!        pen    = same after trial move
!        psio   = psi
!        coef   = read in coefficients of the basis functions
!                 to get the molecular orbitals used in determinant
!        nbasis = number of basis functions read in
!        cdet   = coefficients of the determinants
!        ndet   = number of determinants of molecular orbitals
!                 used
!        nup    = number of up spin electrons
!        ndn    = number of down spin electrons
!        Jastrow function is dexp(cjas1*rij/(1+cjas2*rij))
!   Other variables main program
!        title  = title of run
!        date   = date of run
!        eunit  = energy units
!        sitsca = scaling factor to set up initial configuration of
!                 electrons on sites
!        nsite  = number of electrons to put on each site initially
!        isite  = flag if 1 then take initial configuration from
!                 sites routine
!        cjasa  = simple Jastrow a (should be .5 or .25)
!        cjasb  = simple Jastrow b
!        ijas   = form of wavefunction
!        isc    = form of scaled variables
!        nspin12= If (11) Parallel-spin a's = 1/2 antipar-spin a's
!                         Parallel-spin b's = antipar-spin b's
!                    (12) a's and b's for par and antipar are indep.
!                 The above applies to good psi.

      open(unit=8,form='formatted',file='tape8')
      rewind 8

!     call flush(6)

      if(nforce.gt.1) then
        call setup_force
       else
        ! debug line. ravindra
        if (.not. allocated(iwftype)) allocate (iwftype(nforce), source=0)
        nwprod=1
        nwftype=1
        iwftype(1)=1
      endif

! read walker configurations
      call mc_configs

! get initial value of cpu time

      call elapsed_time("DMC : Reading initial walker configuration : ")



! initialize sums and averages
      call init_averages_index
      if(dmc_irstar.ne.1) call init
      if(dmc_irstar.ne.1) call average(0)

! forces implemented only for certain dmc control options
      if(nforce.gt.1) write(ounit,'(''Possible Warning: force implemented for certain dmc control options'')')

!     call flush(6)
! loops for dmc calculation
      do i=1,dmc_nblk+2*dmc_nblkeq
        if((i.eq.dmc_nblkeq+1.or.i.eq.2*dmc_nblkeq+1).and.dmc_irstar.ne.1) then

          call elapsed_time("DMC : equilibrium CP : ")

          call zerest
          call average(0)
          call elapsed_time("DMC : zero out estimators and averages : ")
        endif
        do j=1,dmc_nstep
          ipass=ipass+1
          if (nloc.gt.0) call rotqua
          if(iabs(idmc).eq.1) then
!           call dmc_brock
           elseif(iabs(idmc).eq.2) then
            if (nloc.gt.0) then
              call dmc_ps(lpass,irun)
             else
!             call dmc_good
            endif
           else
            call fatal_error('DMC: iabs(idmc) must be 1 or 2')
          endif
          call mc_configs_write(i,ipass)

          call acues1
        enddo
        call average(2)
        call average_write
        call acuest
      enddo

      call acues1_reduce

      call finwrt
      call elapsed_time("DMC : all CP : ")

      if (dmc_idump.eq.1) call dumper
      close (unit=9)
      if (dmc_nconf.ne.0) close (unit=7)
      call elapsed_time("dumping restart files : ")

      end
end module
