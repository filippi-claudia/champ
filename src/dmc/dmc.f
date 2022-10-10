      module dmc_f_mod
      contains
      subroutine dmc
c Written by Cyrus Umrigar with major contributions by Claudia Filippi.
c Uses the diffusion Monte Carlo algorithm described in:
c 1) A Diffusion Monte Carlo Algorithm with Very Small Time-Step Errors,
c    C.J. Umrigar, M.P. Nightingale and K.J. Runge, J. Chem. Phys., 99, 2865 (1993).

      use precision_kinds, only: dp
      use const, only: pi
      use forcepar, only: nforce
      use contrldmc, only: idmc
      use estcum, only: ipass
      use force_dmc, only: nwprod
      use pseudo, only: nloc
      use wfsec, only: iwftype, nwftype
!      use contrl, only: idump, irstar, nblk, nblkeq, nconf, nstep
      use control_dmc, only: dmc_idump, dmc_irstar, dmc_nblk, dmc_nblkeq
      use control_dmc, only: dmc_nconf, dmc_nstep
      use mpitimer,    only: elapsed_time
      use contrl_file,    only: ounit

      use strech_mod,     only: setup_force
      use dumper_mod,     only: dumper
      use mc_configs_mod, only: mc_configs, mc_configs_write
      use averages,       only: init_averages_index, average, average_write
      use init_mod,       only: init
      use zerest_mod,     only: zerest
      use rotqua_mod,     only: rotqua
      use error,          only: fatal_error
      use dmc_ps_mov1,    only: dmc_ps
      use acues1_mod,     only: acues1
      use acuest_mod,     only: acuest
      use acues1_reduce_mod,only: acues1_reduce
      use finwrt_mod,     only: finwrt
      implicit none

      integer :: i, j
      real(dp), parameter :: one = 1.d0
      real(dp), parameter :: four = 4.d0


c variables:
c        nelec  = number of electrons
c        pi     = 3.14159...
c        hb     = hbar**2/(2m)
c        delta  = side of box in which metropolis steps are made
c        deltai = 1/delta
c        fbias  = force bias parameter
c        nstep  = number of metropolis steps/block
c        nblk   = number of blocks od nstep steps after the
c                 equilibrium steps
c        nblkeq = number of equilibrium blocks
c        nconf  = initial and target number of dmc configurations
c        idump  =  1 dump out stuff for a restart
c        irstar =  1 pick up stuff for a restart
c        tau    = time-step
c        nfprod = number of f's used for removing finite popul. bias
c Control variables are:
c idmc         < 0     VMC
c              > 0     DMC
c abs(idmc)    = 1 **  simple kernel using dmc.brock.f
c              = 2     good kernel using dmc_good or dmc_good_inhom
c ipq         <= 0 *   do not use expected averages
c             >= 1     use expected averages (mostly in all-electron move algorithm)
c itau_eff    <=-1 *   always use tau in branching (not implemented)
c              = 0     use 0 / tau for acc /nacc moves in branching
c             >= 1     use tau_eff (calcul in equilibration runs) for all moves
c iacc_rej    <=-1 **  accept all moves (except possibly node crossings)
c              = 0 **  use weights rather than accept/reject
c             >= 1     use accept/reject
c icross      <=-1 **  kill walkers that cross nodes (not implemented)
c              = 0     reject walkers that cross nodes
c             >= 1     allow walkers to cross nodes
c                      (OK since crossing prob. goes as tau^(3/2))
c icuspg      <= 0     approximate cusp in Green function
c             >= 1     impose correct cusp condition on Green function
c icut_br     <= 0     do not limit branching
c             >= 1 *   use smooth formulae to limit branching to (1/2,2)
c                      (bad because it makes energies depend on E_trial)
c icut_e      <= 0     do not limit energy
c             >= 1 *   use smooth formulae to limit energy (not implemented)
c *  => bad option, modest deterioration in efficiency or time-step error
c ** => very bad option, big deterioration in efficiency or time-step error
c So, idmc=6,66 correspond to the foll. two:
c 2 1 1 1 0 0 0 0 0  idmc,ipq,itau_eff,iacc_rej,icross,icuspg,idiv_v,icut_br,icut_e
c 2 1 0 1 1 0 0 0 0  idmc,ipq,itau_eff,iacc_rej,icross,icuspg,idiv_v,icut_br,icut_e
c Another reasonable choice is:
c 2 1 0 1 1 1 1 0 0  idmc,ipq,itau_eff,iacc_rej,icross,icuspg,idiv_v,icut_br,icut_e
c        xold_dmc   = current position of the electrons
c        xnew   = new position after a trial move
c        vold_dmc   = grad(psi)/psi at current position
c        vnew   = same after trial move
c        psi2o  = psi**2 at current position
c        psi2n  = same after trial move
c        eold   = local energy at current position
c        enew   = same after trial move
c        peo_dmc    = local potential at current position
c        pen    = same after trial move
c        tjfo   = Jackson Feenberg kinetic energy at current position
c        tjfn   = same after trial move
c        psio   = psi
c        coef   = read in coefficients of the basis functions
c                 to get the molecular orbitals used in determinant
c        nbasis = number of basis functions read in
c        cdet   = coefficients of the determinants
c        ndet   = number of determinants of molecular orbitals
c                 used
c        nup    = number of up spin electrons
c        ndn    = number of down spin electrons
c        Jastrow function is dexp(cjas1*rij/(1+cjas2*rij))
c   Other variables main program
c        title  = title of run
c        date   = date of run
c        eunit  = energy units
c        sitsca = scaling factor to set up initial configuration of
c                 electrons on sites
c        nsite  = number of electrons to put on each site initially
c        isite  = flag if 1 then take initial configuration from
c                 sites routine
c        cjasa  = simple Jastrow a (should be .5 or .25)
c        cjasb  = simple Jastrow b
c        ijas   = form of wavefunction
c        isc    = form of scaled variables
c        nspin12= If (11) Parallel-spin a's = 1/2 antipar-spin a's
c                         Parallel-spin b's = antipar-spin b's
c                    (12) a's and b's for par and antipar are indep.
c                 The above applies to good psi.

      pi=four*datan(one)

      open(unit=8,form='formatted',file='tape8')
      rewind 8

c     call flush(6)

      if(nforce.gt.1) then
        call setup_force
       else
        ! debug line. ravindra
        if (.not. allocated(iwftype)) allocate (iwftype(nforce), source=0)
        nwprod=1
        nwftype=1
        iwftype(1)=1
      endif

c read walker configurations
      call mc_configs

c get initial value of cpu time

      call elapsed_time("DMC : Reading initial walker configuration : ")



c initialize sums and averages
      call init_averages_index
      if(dmc_irstar.ne.1) call init
      if(dmc_irstar.ne.1) call average(0)

c forces implemented only for certain dmc control options
      if(nforce.gt.1) write(ounit,'(''Possible Warning: force implemented for certain dmc control options'')')

c     call flush(6)
c loops for dmc calculation
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
c           call dmc_brock
           elseif(iabs(idmc).eq.2) then
            if (nloc.gt.0) then
              call dmc_ps
             else
c             call dmc_good
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
