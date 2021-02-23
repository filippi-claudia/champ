      subroutine dmc
c Written by Cyrus Umrigar with major contributions by Claudia Filippi.
c Uses the diffusion Monte Carlo algorithm described in:
c 1) A Diffusion Monte Carlo Algorithm with Very Small Time-Step Errors,
c    C.J. Umrigar, M.P. Nightingale and K.J. Runge, J. Chem. Phys., 99, 2865 (1993).

      use vmc_mod, only: MELEC, MORB, MBASIS, MDET, MCENT, MCTYPE, MCTYP3X,
     &NSPLIN, nrad, MORDJ, MORDJ1, MMAT_DIM, MMAT_DIM2, MMAT_DIM20,
     &radmax, delri, NEQSX, MTERMS, MCENT3, NCOEF, MEXCIT
      use dmc_mod, only: MWALK, MFPROD, MFPRD1, MPATH
      use basis, only: zex, betaq, n1s, n2s, n2p, n3s, n3p, n3dzr, n3dx2, n3dxy, n3dxz, n3dyz,
     & n4s, n4p, n4fxxx, n4fyyy, n4fzzz, n4fxxy, n4fxxz, n4fyyx, n4fyyz,
     & n4fzzx, n4fzzy, n4fxyz, nsa, npa, ndzra, ndz2a, ndxya, ndxza, ndyza, ndx2a
      use const, only: delta, deltai, etrial, fbias, hb, imetro, ipr, nelec, pi
      use forcest, only: fgcm2, fgcum
      use forcepar, only: deltot, istrech, nforce
      use contrl_per, only: ibasis, iperiodic
      use contrldmc, only: iacc_rej, icross, icuspg, icut_br, icut_e, idiv_v, idmc, ipq,
     &itau_eff, nfprod, rttau, tau, taueff, tautot
      use atom, only: cent, iwctype, ncent, nctype, pecent, znuc
      use iterat, only: iblk, ipass
      use config, only: d2o, peo_dmc, psido_dmc, psijo_dmc, vold_dmc, xold_dmc
      use force_dmc, only: itausec, nwprod
      use force_mod, only: MFORCE, MFORCE_WT_PRD, MWF
      use ewald_mod, only: NCOEFX, NPX, IVOL_RATIO, IBIG_RATIO, NSYM,
     &NGNORMX, NGVECX, NG1DX, NGNORM_SIMX, NGVEC_SIMX, NGNORM_BIGX,
     &NGVEC_BIGX, NGNORM_SIM_BIGX, NGVEC_SIM_BIGX
      use numbas_mod, only: MRWF, MRWF_PTS
      use numbas, only: arg, d2rwf, igrid, iwrwf, nr, nrbas, numr, r0, rwf
      use pseudo_mod, only: MPS_L, MPS_QUAD, MPS_GRID, MGAUSS
      use qua, only: nquad, wq, xq, xq0, yq, yq0, zq, zq0
      use branch, only: eest, eigv, eold, ff, fprod, nwalk, pwt, wdsumo, wgdsumo, wt, wtgen,
     &wthist
      use bparm, only: nocuspb, nspin2b
      use contr2, only: i3body, ianalyt_lap, iaver, icusp, icusp2, ifock, ijas, irewgt,
     &isc, istrch
      use pseudo, only: lpot, nloc, vps, vpso
      use dorb_m, only: iworbd
      use pars, only: Z, a00, a20, a21, c0000, c1110, c2000, eps_fock, xm1, xm12, xm2, xma,
     &xms
      use rlobxy, only: rlobx, rloby, rloby2
      use rnyucm, only: ll, mm
      use coefs, only: coef, nbasis, norb
      use wfsec, only: iwf, iwftype, nwftype
      use jaspar, only: is, nspin1, nspin2, sspin, sspinn
      use jaspar1, only: cjas1, cjas2
      use jaspar2, only: a1, a2
      use jaspar3, only: a, b, c, fck, nord, scalek
      use jaspar4, only: a4, norda, nordb, nordc

      implicit real*8(a-h,o-z)

      parameter (one=1.d0,four=4.d0)

      common /contrl/ nstep,nblk,nblkeq,nconf,nconf_new,isite,idump,irstar

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
        nwprod=1
        nwftype=1
        iwftype(1)=1
      endif

c read walker configurations
      call mc_configs

c get initial value of cpu time
  350 call my_second(0,'begin ')

c initialize sums and averages
      call init_averages_index
      if(irstar.ne.1) call init
      if(irstar.ne.1) call average(0)

c forces implemented only for certain dmc control options
      if(nforce.gt.1) write(6,'(''Possible Warning: force implemented for certain dmc control options'')')

c     call flush(6)
c loops for dmc calculation
      do 360 i=1,nblk+2*nblkeq
        if((i.eq.nblkeq+1.or.i.eq.2*nblkeq+1).and.irstar.ne.1) then
          call my_second(2,'equilb')
          call zerest
          call average(0)
        endif
        do 355 j=1,nstep
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

  355     call acues1
        call average(2)
        call average_write
  360   call acuest

      call acues1_reduce

      call finwrt
      call my_second(2,'all   ')

      if (idump.eq.1) call dumper
      close (unit=9)
      if (nconf.ne.0) close (unit=7)

      end
