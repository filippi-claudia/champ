      subroutine read_input
c Written by Friedemann Schautz
      
      use contr3, only: mode
      implicit real*8(a-h,o-z)

      
c Initialize flags
      call flaginit
c Initialize input parser
      
      call p2init
      
c Parse input (standard input)
      call p2go(5,0)

c Compute the size of some matrices we need
      call preprocess_input()
      call compute_mat_size()

c Allocate memory of all arrays      
      call allocate_all_arrays()

c Transfer from lists to fortran variables, print out, check,
c and read in everything which is still in the old format
      call process_input

      if(index(mode,'mc').ne.0 ) call p2vin('input.log',1)
      
      return
      end subroutine read_input

      subroutine preprocess_input()
        !> read some parts of the input that 
        !> are needed for the dynamic allocation
        use elec, only: ndn, nup
        use const, only: nelec
        use ghostatom, only: newghostype, nghostcent
        use atom, only: nctype, ncent
        use contrl, only: nstep, nblk, nblk_max
        use contr3, only: mode
        use wfsec, only: nwftype
        use forcepar, only: nforce
        use force_mod, only: MFORCE
        use method_opt, only: method

        implicit none

        !> electrons
        call p2gti('electrons:nelec',nelec,1)
        call p2gti('electrons:nup',nup,1)
        ndn=nelec-nup
  
        !> atoms
        call p2gti('atoms:nctype',nctype,1)
        call p2gti('atoms:natom',ncent,1)
        call p2gtid('atoms:addghostype',newghostype,0,1)
        call p2gtid('atoms:nghostcent',nghostcent,0,1)

        !> force 
        call p2gtid('general:nforce',nforce,1,1)
        MFORCE = nforce

        !> wftype
        call p2gtid('general:nwftype',nwftype,1,1)

        !> sampling
        if(index(mode,'vmc').eq.0) call p2gtad('dmc:mode_dmc',mode,'dmc_one_mpi1',1)
        if(index(mode,'mc').ne.0 ) then

          if(index(mode,'vmc').ne.0) then
            call p2gti('blocking_vmc:nstep',nstep,1)
            call p2gti('blocking_vmc:nblk',nblk,1)
            call p2gtid('optwf:nblk_max',nblk_max,nblk,1)
          endif
          if(index(mode,'dmc').ne.0) then
            call p2gti('blocking_dmc:nstep',nstep,1)
            call p2gti('blocking_dmc:nblk',nblk,1)
            call p2gtid('optwf:nblk_max',nblk_max,nblk,1)
          endif
        endif

        !> optimization method
        call p2gtad('optwf:method', method, 'linear', 1)

      end subroutine preprocess_input

      subroutine compute_mat_size()
        !> compute various size that are derived from the input
        ! use vmc_mod, only: MMAT_DIM, MMAT_DIM2, MCTYP3X, MCENT3
        ! use const, only: nelec
        ! use atom, only: nctype_tot, ncent_tot
        use sr_mod, only: MPARM, MOBS, MCONF
        use contrl, only: nstep, nblk_max

        use vmc_mod, only: set_vmc_size
        use optci, only: set_optci_size
        use optorb_mod, only: set_optorb_size
        use gradhess_all, only: set_gradhess_all_size
        ! use sr_mod, only: set_sr_size

        implicit none 

        ! leads to circular dependecy of put in sr_mod ..
        MOBS = 10 + 6*MPARM
        MCONF = nstep * nblk_max

        call set_vmc_size
        call set_optci_size
        call set_optorb_size
        call set_gradhess_all_size
        ! call set_sr_size

      end subroutine

      subroutine allocate_all_arrays()
            !> massive dynamic allocation
            call allocate_m_common
            call allocate_m_basis
            call allocate_m_control
            call allocate_m_deriv
            call allocate_m_efield
            call allocate_m_estimators
            call allocate_m_ewald
            call allocate_m_force
            call allocate_m_gradhess
            call allocate_m_grdnt
            call allocate_m_grid
            call allocate_m_jastrow
            call allocate_m_mixderiv
            call allocate_m_mmpol
            call allocate_m_mstates
            call allocate_m_optci
            call allocate_m_optorb
            call allocate_m_optwf
            call allocate_m_pcm
            call allocate_m_prop
            call allocate_m_pseudo
            call allocate_m_sampling
            call allocate_m_sr
            call allocate_m_state_avrg
      end subroutine allocate_all_arrays

c-----------------------------------------------------------------------
      subroutine process_input
c Written by Cyrus Umrigar, Claudia Filippi, Friedemann Schautz,
c and Anthony Scemema
      use sr_mod, only: MCONF, MVEC
      use pseudo_mod, only: MPS_QUAD
      use properties, only: MAXPROP
      use optorb_mod, only: MXORBOP, MXREDUCED
      use optci, only: MXCITERM
      use mmpol_mod, only: mmpolfile_sites, mmpolfile_chmm
      use force_mod, only: MFORCE, MWF
      use vmc_mod, only: MELEC, MORB, MBASIS, MCENT, MCTYPE, MCTYP3X
      use atom, only: znuc, cent, pecent, iwctype, nctype, ncent, ncent_tot, nctype_tot
      use jaspar, only: nspin1, nspin2, is
      use ghostatom, only: newghostype, nghostcent
      use const, only: pi, hb, etrial, delta, deltai, fbias, nelec, imetro, ipr
      use jaspar1, only: cjas1, cjas2
      use general, only: pooldir, pp_id, bas_id
      use general, only: filenames_bas_num, wforce 
      use csfs, only: cxdet, ncsf, nstates
      use dets, only: cdet, ndet
      use elec, only: ndn, nup
      use forcepar, only: nforce
      use grdntspar, only: igrdtype, ngradnts
      use header, only: title
      use jaspar2, only: a1, a2
      use jaspar3, only: a, b, c, nord, scalek
      use jaspar4, only: a4, norda, nordb, nordc
      use jaspar6, only: asymp_jasa, asymp_jasb, asymp_r, c1_jas6, c1_jas6i, c2_jas6
      use jaspar6, only: cutjas, cutjasi
      use numbas, only: numr
      use numbas1, only: nbastyp
      use numbas2, only: ibas0, ibas1
      use optwf_contrl, only: ioptci, ioptjas, ioptorb, ioptwf
      use optwf_contrl, only: idl_flag, ilbfgs_flag, ilbfgs_m, dl_mom, dl_alg
      use optwf_contrl, only: ibeta, ratio_j, iapprox, ncore
      use optwf_contrl, only: iuse_orbeigv
      use optwf_parms, only: nparmj
      use optwf_sr_mod, only: i_sr_rescale, izvzb
      use pars, only: Z, a20, a21
      use rlobxy, only: rlobx
      use sa_weights, only: iweight, nweight, weights
      use wfsec, only: nwftype
      use zmatrix, only: izmatrix
      use bparm, only: nocuspb, nspin2b
      use casula, only: i_vpsp, icasula
      use coefs, only: coef, nbasis, norb
      use const2, only: deltar, deltat
      use contr2, only: ianalyt_lap, ijas
      use contr2, only: isc
      use contr3, only: mode
      use contrldmc, only: iacc_rej, icross, icuspg, icut_br, icut_e, idiv_v, idmc, ipq
      use contrldmc, only: itau_eff, nfprod, rttau, tau
      use contrl, only: idump, irstar, isite, nconf, nblk, nblkeq, nconf_new, nstep
      use contrl, only: icharged_atom, nblk_max, nblk_ci
      use dorb_m, only: iworbd
      use contrl_per, only: iperiodic, ibasis
      use force_analy, only: iforce_analy, iuse_zmat, alfgeo
      use force_dmc, only: itausec, nwprod
      use pseudo, only: nloc
      use optorb_cblock, only: idump_blockav
      use gradjerrb, only: ngrad_jas_blocks
      use qua, only: nquad, wq, xq, yq, zq
      use mmpol_cntrl, only: ich_mmpol, immpol, immpolprt, isites_mmpol
      use mmpol_parms, only: chmm
      use mmpol_fdc, only: a_cutoff, rcolm
      use grid3dflag, only: i3ddensity, i3dgrid, i3dlagorb, i3dsplorb
      use grid_mod, only: UNDEFINED, IUNDEFINED
      use efield, only: iefield, ncharges
      use mstates_ctrl, only: iefficiency, iguiding, nstates_psig
      use mstates3, only: iweight_g, weights_g
      use ci000, only: iciprt, nciprim, nciterm
      use pcm_cntrl, only: ichpol, ipcm, ipcmprt, isurf
      use pcm_unit, only: pcmfile_cavity, pcmfile_chs, pcmfile_chv
      use pcm_parms, only: eps_solv, iscov
      use pcm_parms, only: ncopcm, nscv, nvopcm
      use prp000, only: iprop, ipropprt, nprop
      use pcm_fdc, only: qfree, rcolv
      use pcm_grid3d_contrl, only: ipcm_3dgrid
      use pcm_grid3d_param, only: ipcm_nstep3d, pcm_step3d, pcm_origin, pcm_endpt
      use pcm_3dgrid, only: PCM_SHIFT
      use prp003, only: cc_nuc
      use method_opt, only: method
      use optorb_cblock, only: nefp_blocks, isample_cmat, iorbsample
      use orbval, only: ddorb, dorb, nadorb, ndetorb, orb
      use array_resize_utils, only: resize_tensor
      use grid3d_param, only: endpt, nstep3d, origin, step3d
      use inputflags, only: node_cutoff, eps_node_cutoff, iqmmm, scalecoef
      use optwf_contrl, only: energy_tol, dparm_norm_min, nopt_iter, micro_iter_sr
      use optwf_contrl, only: nvec, nvecx, alin_adiag, alin_eps, lin_jdav, multiple_adiag
      use optwf_contrl, only: ilastvmc, iroot_geo
      use optwf_contrl, only: sr_tau , sr_adiag, sr_eps 
      use optwf_func, only: ifunc_omega, omega0, n_omegaf, n_omegat
      use optwf_corsam, only: add_diag

      implicit real*8(a-h,o-z)

      parameter (zero=0.d0,one=1.d0,two=2.d0,four=4.d0)

c      include 'dmc.h' now emty

      character*20 fmt
      character*32 keyname
      character*10 eunit
      character*16 cseed
      dimension irn(4),cent_tmp(3),anorm(nbasis)

c Inputs:
c  title      title
c  irn        random number seeds (four 4-digit integers)
c  ijas       form of Jastrow. (between 1 and 6, mostly we use 4)
c  isc        form of scaling function for ri,rj,rij in Jastrow (between 1 and 7, mostly use 2,4,6,7,16,17)
c  iperiodic  0  finite system
c             >0 periodic system
c  ibasis     form of basis
c  hb         hbar=0.5 for Hartree units
c  etrial     guess for energy
c  eunit      'Hartree'
c  nstep      number of steps per block
c  nblk       number of blocks
c  nblkeq     number of equilibration blocks
c  nconf      target number of MC configurations in dmc
c  nconf_new  number of new MC configs. saved per processor.
c  idump      dump restart file
c  irstar     restart from restart file
c  isite      call sites to generate starting MC config. in vmc
c  ipr        print level
c  imetro     form of Metropolis (6 is most efficient choice for most systems)
c             1 simple algorithm with force-bias
c             6 accelerated Metropolis algorithm from Cyrus' 1993 PRL
c  delta      step-size for simple algorithm
c  deltar     radial step-size for accelerated algorithm
c  deltat     angular step-size for accelerated algorithm
c  fbias      force-bias.  (Use 1 always).
c  idmc       form of dmc algorithm
c             1  simple dmc
c             2  improved dmc from Umrigar, Nightingale, Runge 1993 JCP
c  ipq,itau_eff,iacc_rej,icross,icuspg,idiv_v
c  nfprod     number of products to undo for estimating population control bias in dmc
c  tau        time-step in dmc
c  nloc       nonlocal pseudopotential
c             0  local
c             1  in Fahy format
c             2  in Troullier-Martins format (unformatted)
c             3  in Troullier-Martins format (formatted)
c  nquad      number of angular quadrature points for nonlocal psp.
c  nelec      number of electrons
c  nup        number of up-spin electrons
c  nctype     number of atom/center types
c  ncent      number of atoms/centers
c  iwctype    specify atom-type for each atom
c  znuc       nuclear charge
c  cent       atom positions
c  ndet       number of determinants in wavefunction
c  nbasis     number of basis functions
c  norb       number of orbitals
c  cdet       coefficients of determinants
c  iworbd     which orbitals enter in which determinants
c  ianalyt_lap analytic laplacian or not
c  ijas     form of Jastrow. (between 1 and 6, mostly we use 4)
c  isc      form of scaling function for ri,rj,rij in Jastrow (between 1 and 7, mostly use 2,4,6,7)
c           2  [1-exp(scalek*r)]/scalek
c           3  [1-exp{-scalek*r-(scalek*r)**2/2}]/scalek
c           4  r/(1+scalek*r)
c           5  r/{1+(scalek*r)**2}**.5
c           6  Short-range version of 2 (range given bu cutjas)
c           7  Short-range version of 4 (range given bu cutjas)

c  nspin2   1,2,3,-1,-2 -> nspin2b=abs(nspin2)
c  nspin2   > 0  nspin2 sets of a, c parms, nspin2b sets of b parms
c              nocuspb=0  parallel e-e cusp conditions satisfied (b=1/2,1/4)
c  nspin2   < 0  -> nspin2=1
c                nspin2=1 sets of a and c parms, nspin2b sets of b parms
c                -1 nocuspb=1 parallel e-e cusp conditions not satisfied (1/2,1/2)
c                -2 nocuspb=0 parallel e-e cusp conditions satisfied (1/2,1/4)
c  nord     order of the polynmial
c  norda    order of the e-n polynmial in Jastrow4
c  nordb    order of the e-e polynmial in Jastrow4
c  nordc    order of the e-e-n polynmial in Jastrow4
c  cjas1    simple jastrow1 (0.5 to satisfy cusps, parallel-spins automatically take half this value)
c  cjas2    simple jastrow1 parameter
c  scalek   scale factor for Jastrow
c  a1,a2    Jastrow parameters for Jastrow2
c  a,b,c    Jastrow parameters for Jastrow3
c  a4,b,c   Jastrow parameters for Jastrow4,5,6
c  cutjas   cutoff for Jastrow4,5,6 if cutjas=6,7
c  rlobx(y) Lobachevsky parameters for Fock expansion

      pi=four*datan(one)

c Default dmc_one_mpi1 separate populations; dmc_one_mpi2 global population
      if(index(mode,'vmc').eq.0) call p2gtad('dmc:mode_dmc',mode,'dmc_one_mpi1',1)

c Check that necessary blocks are in input file
      call flagcheck

c Get weights for multiple states calculations
      keyname='weights:'
      call get_weights(keyname,weights,iweight,nweight)

c General section
      call p2gtad('general:title',title,' ',1)
      call p2gtad('general:pool',pooldir,'.',1)
      call p2gtad('general:pseudopot',pp_id,'none',1)
      call p2gtad('general:basis',bas_id,'none',1)

      call stripquotes(title)
      call stripquotes(pooldir)
      call stripquotes(bas_id)
      call stripquotes(pp_id)

      if(mode.eq.'vmc')
     & write(6,'(''Variational MC'',a40)') title
      if(mode.eq.'vmc_one_mpi')
     & write(6,'(''Variational MC one-electron move mpi'',a40)') title
      if(mode.eq.'dmc')
     & write(6,'(''Diffusion MC'',a40)') title
      if(mode.eq.'dmc_one_mpi1')
     & write(6,'(''Diffusion MC 1-electron move, mpi no global pop'',a40)') title
      if(mode.eq.'dmc_one_mpi2')
     & write(6,'(''Diffusion MC 1-electron move, mpi global pop comm'',a40)') title

      call p2gtid('general:iperiodic',iperiodic,0,1)
      call p2gtid('general:ibasis',ibasis,1,1)

      if(index(mode,'vmc').ne.0 .and. iperiodic.gt.0) 
     & call fatal_error('INPUT: VMC for periodic system -> run dmc/dmc.mov1 with idmc < 0')

      if(index(mode,'mc').ne.0 ) then
        call p2gta('general:seed',cseed,1)
        read(cseed,'(4i4)') irn
        write(6,'(/,''random number seeds'',t30,4i4)') irn
        call setrn(irn)

        call p2gtid('general:ipr',ipr,-1,1)
      endif

      call p2gtad('general:unit',eunit,'Hartrees',1)
      call p2gtfd('general:mass',hb,0.5d0,1)
      write(6,'(''all energies are in'',t30,a10)') eunit
      write(6,'(''hbar**2/(2.*m) ='',t30,f10.5)') hb

      call p2gtid('general:nforce',nforce,1,1)
      write(6,'(/,''number of geometries ='',t30,i10)') nforce
      
      ! if(nforce.gt.MFORCE) call fatal_error('INPUT: nforce > MFORCE')
      call p2gtid('general:nwftype',nwftype,1,1)
      write(6,'(/,''number of wave functions='',t30,i10)') nwftype
      if(nwftype.gt.nforce) call fatal_error('INPUT: nwftype gt nforce')
      !if(nwftype.gt.MWF) call fatal_error('INPUT: nwftype gt MWF')

c Electron section
      call p2gti('electrons:nelec',nelec,1)
      write(6,'(''number of electrons ='',t30,i10)') nelec
      ! if(nelec.gt.MELEC) call fatal_error('INPUT: nelec exceeds MELEC')

      call p2gti('electrons:nup',nup,1)
      if(nup.gt.nelec/2) call fatal_error('INPUT: nup exceeds nelec/2')
      ndn=nelec-nup
      write(6,'(''number of up,dn electrons ='',t30,2i5)') nup,ndn

c VMC parameters
      if(index(mode,'vmc').ne.0) then

        call p2gtid('vmc:imetro',imetro,6,1)
        write(6,'(''version of Metropolis ='',t30,i10)') imetro
        call p2gtid('vmc:node_cutoff',node_cutoff,0,1)
        call p2gtfd('vmc:enode_cutoff',eps_node_cutoff,1.d-7,1)

        if(imetro.eq.1) then
          call p2gtf('vmc:delta',delta,1)
          deltai=one/delta

          write(6,'(''step size ='',t30,f10.5)') delta
         else
          call p2gtf('vmc:deltar',deltar,1)
          call p2gtf('vmc:deltat',deltat,1)
          if(deltar.lt.one) then
            write(6,*) '**Warning value of deltar reset to 2.'
            deltar=two
          endif
          if(deltat.lt.zero .or. deltat.gt.two) then
            write(6,*) '**Warning value of deltat reset to 2.'
            deltat=two
          endif
          write(6,'(''radial step multiplier ='',t30,f10.5)') deltar
          write(6,'(''cos(theta) step size ='',t30,f10.5)') deltat
        endif
        call p2gtfd('vmc:fbias',fbias,1.0d0,1)
c Truncate fbias so that fbias and the sampled quantity are never negative
        fbias=dmin1(two,dmax1(zero,fbias))

        write(6,'(''force bias ='',t30,f10.5)') fbias
      endif

      if(index(mode,'mov1').ne.0 .and. imetro.eq.1) call fatal_error('INPUT: metrop_mov1 not updated')

c DMC parameters
      if(index(mode,'dmc').ne.0) then

        call p2gtid('dmc:idmc',idmc,2,1)
        call p2gtid('dmc:ipq',ipq,1,1)
        call p2gtid('dmc:itau_eff',itau_eff,1,1)
        call p2gtid('dmc:iacc_rej',iacc_rej,1,1)
        call p2gtid('dmc:icross',icross,1,1)
        call p2gtid('dmc:icuspg',icuspg,0,1)
        call p2gtid('dmc:idiv_v',idiv_v,0,1)
        call p2gtid('dmc:icut_br',icut_br,0,1)
        call p2gtid('dmc:icut_e',icut_e,0,1)

        call p2gtid('dmc:node_cutoff',node_cutoff,0,1)
        call p2gtfd('dmc:enode_cutoff',eps_node_cutoff,1.d-7,1)

        call p2gti('dmc:nfprod',nfprod,1)
        call p2gtf('dmc:tau',tau,1)
        rttau=dsqrt(tau)

        call p2gtf('dmc:etrial',etrial,1)
        call p2gtid('forces:nwprod',nwprod,200,1)
        call p2gtid('forces:itausec',itausec,1,1)

        write(6,'(''version of DMC ='',t30,i10)') idmc
        write(6,'(''nfprod,tau'',t30,i10,f10.5)') nfprod,tau
        write(6,'(''ipq,itau_eff,iacc_rej,icross,icuspg,idiv_v,icut_br,icut_e='',9i4)')
     &  idmc,ipq,itau_eff,iacc_rej,icross,icuspg,idiv_v,icut_br,icut_e
        write(6,'(''etrial'',t30,f10.6)') etrial

        call p2gtid('dmc:icasula',icasula,0,1)
        write(6,'(''casula'',i4)') icasula

        write(6,'(''node_cutoff'',i4,f12.6)') node_cutoff
        if(node_cutoff.gt.0) write(6,'(''enode_cutoff'',f12.6)') enode_cutoff

        if(idmc.ne.2) call fatal_error('INPUT: only idmc=2 supported')
        call p2gtid('pseudo:nloc',nloc,0,1)
        if(nloc.eq.0) call fatal_error('INPUT: no all-electron DMC calculations supported')

       else
        icasula=0

      endif
c Inizialized to zero for call to hpsi in vmc or dmc with no casula or/and in acuest
      i_vpsp=0

c Parameters for blocking/start/dump
      if(index(mode,'mc').ne.0 ) then

        if(index(mode,'vmc').ne.0) then
          call p2gti('blocking_vmc:nstep',nstep,1)
          call p2gti('blocking_vmc:nblk',nblk,1)
          call p2gti('blocking_vmc:nconf_new',nconf_new,1)
          call p2gtid('blocking_vmc:nblkeq',nblkeq,2,1)
        endif
        if(index(mode,'dmc').ne.0) then
          call p2gti('blocking_dmc:nstep',nstep,1)
          call p2gti('blocking_dmc:nblk',nblk,1)
          call p2gti('blocking_dmc:nconf',nconf,1)
          call p2gtid('blocking_dmc:nblkeq',nblkeq,2,1)
          call p2gtid('blocking_dmc:nconf_new',nconf_new,0,1)
        endif
        call p2gtid('startend:idump',idump,1,1)
        call p2gtid('startend:irstar',irstar,0,1)
        call p2gtid('startend:isite',isite,1,1)
        if (isite.eq.1) call p2gtid('startend:icharged_atom',icharged_atom,0,1) 

c Make sure that the printout is not huge
        if(nstep*(nblk+2*nblkeq).gt.104000) ipr=-1
        if(irstar.eq.1) nblkeq=0

        write(6,'(''no. of steps/block ='',t30,i10)') nstep
        write(6,'(''no. of blocks after eq.='',t30,i10)') nblk
        write(6,'(''no. of blocks before eq. ='',t30,i10)') nblkeq
        if(index(mode,'vmc').ne.0) then
          write(6,'(''no. configurations saved ='',t30,i10)') nconf_new
         elseif(index(mode,'dmc').ne.0) then
          write(6,'(''target walker population ='',t30,i10)') nconf
          if(nconf.le.0) call fatal_error('INPUT: target population <= 0')
          if(nconf.gt.MWALK) call fatal_error('INPUT: target population > MWALK')
          write(6,'(''no. configurations saved ='',t31,i10)') nconf_new
        endif

       endif
      write(6,*)

c Analytical forces flags (vmc only)
      if(index(mode,'vmc').ne.0) then
        call p2gtid('optgeo:iforce_analy',iforce_analy,0,0)
        call p2gtid('optgeo:iuse_zmat',iuse_zmat,0,0)
        call p2gtid('optgeo:izvzb',izvzb,0,1)
        if(iforce_analy.gt.0) then
          if(nordc.gt.0) call fatal_error('READ_INPUT: Nuclear analytic forces not implemented for 3-body J')
          call p2gtfd('optgeo:alfgeo',alfgeo,1.d0,0)
          write(6,'(''Geometry optimization with analytic gradients'')')
          if(iuse_zmat.gt.0) write(6,'(''use internal coordinates'')')
          write(6,'(''starting alfgeo ='',f10.4)') alfgeo
        endif
      endif

c Optimization flags (vmc/dmc only)
      if(index(mode,'vmc').ne.0.or.index(mode,'dmc').ne.0) then
        call p2gtid('optwf:ioptwf', ioptwf, 0, 1)
        call p2gtad('optwf:method', method, 'linear', 1)
        call p2gtid('optwf:idl_flag', idl_flag, 0, 1)
        call p2gtid('optwf:ilbfgs_flag', ilbfgs_flag, 0, 1)
        call p2gtid('optwf:ilbfgs_m',ilbfgs_m,5,1)
        call p2gtid('optwf:sr_rescale',i_sr_rescale,0,1)

        call p2gtid('optwf:ibeta',ibeta,-1,1)
        call p2gtfd('optwf:ratio',ratio,ratio_j,1)
        call p2gtid('optwf:approx',iapprox,0,1)
        call p2gtid('optwf:ncore',ncore,0,1)
        call p2gtid('optwf:iuse_orbeigv',iuse_orbeigv,0,1)

        call p2gtid('optwf:ioptjas',ioptjas,0,1)
CVARDOC flag: Jastrow derivatives will be sampled
        call p2gtid('optwf:ioptorb',ioptorb,0,1)
CVARDOC flag: ORB-PT derivatives will be sampled
        call p2gtid('optwf:ioptci',ioptci,0,1)
CVARDOC flag: CI derivatives will be sampled
        call p2gtid('optwf:idl_flag',idl_flag,0,1)
CVARDOC flag: Deep learning optimization algorithm wil be used
        call p2gtid('optwf:ilbfgs_flag',ilbfgs_flag,0,1)
CVARDOC flag: oLBFGS optimization algorithm wil be used

        call p2gtid('optwf:ioptwf',ioptwf,0,1)
        if(ioptwf.gt.0) then
          write(6,'(''Perform wave function optimization in vmc/dmc'')')
         elseif(ioptjas.eq.1.or.ioptorb.eq.1.or.ioptci.eq.1) then
          write(6,'(''Only sample derivatives of wave function for external use'')')
        endif

        if(ioptwf.gt.0.or.ioptjas+ioptorb+ioptci.ne.0) then

        call p2gtad('optwf:method',method,'linear',1)
        call p2gtid('optwf:nblk_max',nblk_max,nblk,1)

! lin_d, sr_n, mix_n and linear shared flags: 
        if((method.eq.'lin_d').or.(method.eq.'sr_n').or.(method.eq.'mix_n')
     #       .or.(method.eq.'linear')) then
           call p2gtfd('optwf:energy_tol', energy_tol, 1.d-3, 1)
           call p2gtfd('optwf:dparm_norm_min', dparm_norm_min, 1.0d0, 1)
           call p2gtfd('optwf:add_diag',add_diag(1),1.d-6,1)
           call p2gtid('optwf:nopt_iter', nopt_iter, 6, 1)
           call p2gtid('optwf:micro_iter_sr', micro_iter_sr, 1, 1)
        end if

! lin_d and sr_n shared flags: 
        if((method.eq.'lin_d').or.(method.eq.'sr_n')) then
           call p2gtid('optwf:func_omega', ifunc_omega, 0, 1)
           if(ifunc_omega.gt.0) then
             call p2gtfd('optwf:omega', omega0, 0.d0, 1)
             call p2gtid('optwf:n_omegaf', n_omegaf, nopt_iter, 1)
             call p2gtid('optwf:n_omegat', n_omegat, 0, 1)
           end if
        end if

! lin_d and mix_n shared flags: 
        if ((method.eq.'lin_d').or.(method.eq.'mix_n').or.(method.eq.'linear')) then
          call p2gtid('optwf:lin_nvec', nvec, 5, 1)
          call p2gtid('optwf:lin_nvecx', nvecx, MVEC, 1)
          call p2gtfd('optwf:lin_adiag', alin_adiag, 0.01, 1)
          call p2gtfd('optwf:lin_eps', alin_eps, 0.001, 1)
          call p2gtid('optwf:lin_jdav',lin_jdav,0,1)
          call p2gtid('optwf:multiple_adiag',multiple_adiag,0,1)
        end if

! sr_n and mix_n shared flags: 
        if ((method.eq.'sr_n').or.(method.eq.'mix_n')) then
          call p2gtfd('optwf:sr_tau', sr_tau, 0.02, 1)
          call p2gtfd('optwf:sr_adiag', sr_adiag, 0.01, 1)
          call p2gtfd('optwf:sr_eps', sr_eps, 0.001, 1)
        end if
! mix_n and linear flags:
        if ((method.eq.'mix_n').or.(method.eq.'linear')) then
           if(iforce_analy.gt.0) call p2gtid('optgeo:iroot_geo',iroot_geo,0,0)
           call p2gtid('optwf:nblk_ci',nblk_ci,nblk,1)
           call p2gtid('optwf:ilastvmc',ilastvmc,1,1)
        end if
! dl flags:
        if (idl_flag .gt. 0) then 
            call p2gtid('optwf:nopt_iter', nopt_iter, 6, 1)
            call p2gtfd('optwf:energy_tol', energy_tol, 1.d-3, 1)
            call p2gtfd('optwf:dparm_norm_min', dparm_norm_min, 1.0d0, 1)

            call p2gtfd('optwf:sr_tau', sr_tau, 0.02, 1)
            call p2gtfd('optwf:sr_adiag', sr_adiag, 0.01, 1)
            call p2gtfd('optwf:sr_eps', sr_eps, 0.001, 1)

            call p2gtfd('optwf:dl_mom', dl_mom, 0.0, 1)
            call p2gtad('optwf:dl_alg', dl_alg, 'nag', 1)
        end if

        if(method.eq.'linear'.and.MXREDUCED.ne.MXORBOP) 
     &    call fatal_error('READ_INPUT: MXREDUCED.ne.MXORBOP')
    !     if((method.eq.'sr_n'.or.method.eq.'lin_d').and.nstep*nblk_max.gt.MCONF)
    !  &    call fatal_error('READ_INPUT: nstep*nblk_max.gt.MCONF')
        endif

        if(ioptjas.eq.1.or.ioptorb.eq.1.or.ioptci.eq.1) 
     &  write(6,'(''Computing/writing quantities for optimization with method '',a10)') method

c Jastrow optimization flag (vmc/dmc only)
        if(ioptjas.gt.0) then
          write(6,'(''Jastrow derivatives are sampled'')') 
          write(6,'(''Number of Jastrow derivatives'',i5)') nparmj
          if(ijas.eq.4) then
            call cuspinit4(1)
           else
            call fatal_error('READ_INPUT: jasderiv only for ijas=4')
          endif
          call p2gtid('optwf:ngrad_jas_blocks',ngrad_jas_blocks,0,1)
         else
          nparmj=0
        endif

c ORB optimization flags (vmc/dmc only)
        call p2gtid('optwf:isample_cmat',isample_cmat,1,1)
CVARDOC flag: Correlation matrix will be sampled
        call p2gtid('optwf:save_blocks',isavebl,0,1)
CVARDOC Save block averages for error analysis
        call p2gtid('optwf:force_blocks',nefp_blocks,1,1)
CVARDOC Average orb-forces over n blocks.
        call p2gtid('optwf:iorbsample',iorbsample,1,1)
CVARDOC Sample frequency of orbital derivatives
        if(ioptorb.ne.0) then
          write(6,'(''Orbital derivatives are sampled'')')
          write(6,'(''ORB-PT blocks in force average = '',t30,i10)') nefp_blocks
          if(isavebl.ne.0)then
            write(6,'(''ORB-PT block averages will be saved'')')
            idump_blockav=43
            open(unit=idump_blockav,file='efpci_blockav.dat',status
     &         ='unknown',form='unformatted')
            write(idump_blockav) nefpterm
          endif
        endif

c CI flags
        call p2gtid('optwf:ioptci',ioptci,0,1)
CVARDOC flag: Linear CI will be sampled
        call p2gtid('ci:iciprt',iciprt,0,1)
CVARDOC flag: CI averages will be printed
        if(ioptci.ne.0) then
          write(6,'(''CI is sampled'')')
          write(6,'(''CI printout flag = '',t30,i10)') iciprt
          if(ioptjas.eq.0.and.ioptorb.eq.0.and.method.eq.'hessian') then
            method='linear'
            write(6,'(''Reset optimization method to linear'')')
          endif
c TMP due to changing kref -> also for ncsf=0, we need to have cxdet(i) carrying the phase
c         if(ncsf.eq.0) call fatal_error('ncsf.eq.0 - further changes needed due to kref')
          if(ncsf.gt.0) then
            nciterm=ncsf
           else
            nciterm=nciprim
          endif
         else
          nciprim=0
          nciterm=0
        endif
        write(6,'(''CI number of coefficients'',t30,i10)') nciterm
        if(ncsf.eq.0.and.nciprim.gt.MXCITERM) call fatal_error('INPUT: nciprim gt MXCITERM')
        if(nciterm.gt.MXCITERM) call fatal_error('INPUT: nciterm gt MXCITERM')


c Multiple states/efficiency/guiding flags
        call p2gtid('mstates:iguiding',iguiding,0,1)
CVARDOC flag: Use guiding wave function constructed from mstates
        if(iguiding.gt.0) then
          write(6,'(''Guiding function: square root of sum of squares'')')
          keyname='weights_guiding:'
          call get_weights(keyname,weights_g,iweight_g,nstates_g)
        endif
        call p2gtid('mstates:iefficiency',iefficiency,0,1)
CVARDOC flag: Efficiency for sampling states inputed in multiple_cistates
        if(iefficiency.gt.0) nstates_psig=nstates

       else
        ioptorb=0
        ioptci=0
        ioptjas=0
        iefficiency=0
        iguiding=0
        nstates=1
      endif

c QMMM classical potential
      call p2gtid('qmmm:iqmmm',iqmmm,0,1)
      if(iqmmm.gt.0) 
     & write(6,'(''QMMM external potential'')')
      if(iqmmm.gt.0) call qmmm_extpot_read

c Read in point charges
      call p2gtid('efield:iefield',iefield,0,1)
      if(iefield.gt.0) then
c External charges fetched in read_efield
        write(6,'(''Read in'',i4,'' external charges'')') ncharges
        call efield_compute_extint
      endif

c PCM polarization charges
c  ipcm=1 computes only the cavity (no qmc calculations)
c  ipcm=2 runs qmc and creates/updates polarization charges 
c  ipcm=3 runs qmc with fixed polarization charges
      call p2gtid('pcm:ipcm',ipcm,0,1)
      call p2gtid('pcm:ipcmprt',ipcmprt,0,1)

      isurf=0
      ncopcm=0
      nvopcm=0
      ichpol=0
      if(ipcm.ne.0) then
        if(ipcm.eq.2) ichpol=1
        if(ipcm.eq.1) isurf=1

        call p2gtad('pcm:file_cavity',pcmfile_cavity,'pcm000.dat',1)
        call p2gtad('pcm:file_chs',pcmfile_chs,'chsurf_old',1)
        call p2gtad('pcm:file_chv',pcmfile_chv,'chvol_old',1)

        nstep2=nstep/2
        call p2gtid('pcm:nblk_chv',nscv,nblk,1)
        call p2gtid('pcm:nstep_chv',iscov,nstep2,1)
        if(iscov.eq.0) call fatal_error('READ_INPUT: iscov eq 0')

        call p2gtf('pcm:eps_solv',eps_solv,1)
        call p2gtfd('pcm:fcol',fcol,1.d0,1)
        call p2gtfd('pcm:rcolv',rcolv,0.04d0,1)
        call p2gti('pcm:npmax',npmax,1)
        qfree=-nelec
        do i=1,ncent
          qfree=qfree+znuc(iwctype(i))
        enddo
        
        write(6,'(''PCM polarization charges '')')
        write(6,'(''pcm ipcm   =  '',t30,i3)') ipcm
        write(6,'(''pcm ichpol =  '',t30,i3)') ichpol
        write(6,'(''pcm isurf  =  '',t30,i3)') isurf
        write(6,'(''pcm file (cavity) ='',t30,a20)') pcmfile_cavity
        write(6,'(''pcm file (chs)    ='',t30,a20)') pcmfile_chs
        write(6,'(''pcm file (chv)    ='',t30,a20)') pcmfile_chv
        write(6,'(''pcm nconf sampled for chv ='',t30,i10)') nscv
        write(6,'(''pcm frequency for chv ='',t30,i10)') iscov
        write(6,'(''pcm epsilon_solvent ='',t30,f7.3)') eps_solv
        write(6,'(''pcm rcolv ='',t30,f7.3)') rcolv
        write(6,'(''pcm fcol  ='',t30,f7.3)') fcol
        write(6,'(''pcm npmax ='',t30,i10)') npmax

        call pcm_extpot_read(fcol,npmax)

        !  We use the UNDEFINED, IUNDEFINED from grid_mod (that are the same as pcm_3dgrid)
        call p2gtid('pcm:ipcm_3dgrid',ipcm_3dgrid,0,1)
        call p2gtid('pcm:nx_pcm',ipcm_nstep3d(1),IUNDEFINED,1)
        call p2gtid('pcm:ny_pcm',ipcm_nstep3d(2),IUNDEFINED,1)
        call p2gtid('pcm:nz_pcm',ipcm_nstep3d(3),IUNDEFINED,1)
 
        call p2gtfd('pcm:dx_pcm',pcm_step3d(1),UNDEFINED,1)
        call p2gtfd('pcm:dy_pcm',pcm_step3d(2),UNDEFINED,1)
        call p2gtfd('pcm:dz_pcm',pcm_step3d(3),UNDEFINED,1)
  
        call p2gtfd('pcm:x0_pcm',pcm_origin(1),UNDEFINED,1)
        call p2gtfd('pcm:y0_pcm',pcm_origin(2),UNDEFINED,1)
        call p2gtfd('pcm:z0_pcm',pcm_origin(3),UNDEFINED,1)
  
        call p2gtfd('pcm:xn_pcm',pcm_endpt(1),UNDEFINED,1)
        call p2gtfd('pcm:yn_pcm',pcm_endpt(2),UNDEFINED,1)
        call p2gtfd('pcm:zn_pcm',pcm_endpt(3),UNDEFINED,1)

        call p2gtfd('pcm:shift',PCM_SHIFT,4.d0,1)
        if(ipcm_3dgrid.gt.0) then
         if(ipcm.ne.3) call fatal('READ_INPUT:ipcm_3dgrid gt 0 & ipcm ne 3')
         call pcm_setup_grid
         call pcm_setup_3dspl
        endif
      endif

c QM-MMPOL  (fxed charges)
c  immpol=1 runs qmc (QM-MM)  and creates the first set of induced dipoles on MM sites
c  immpol=2 runs qmc (QM-MMPOL) and updates the set of induced dipoles on MM sites
c  immpol=3 runs qmc  (QM-MMPOL) with fixed induced dipoles
      call p2gtid('mmpol:immpol',immpol,0,1)
      call p2gtid('mmpol:immpolprt',immpolprt,0,1)

      isites_mmpol=0
      ich_mmpol=0
      if(immpol.ne.0) then
        if(immpol.eq.2) ich_mmpol=1

        call p2gtad('mmpol:file_sites',mmpolfile_sites,'mmpol000.dat',1)
        call p2gtad('mmpol:file_mmdipo',mmpolfile_chmm,'mmdipo_old',1)
        call p2gtfd('mmpol:a_cutoff',a_cutoff,2.5874d0,1)
        call p2gtfd('mmpol:rcolm',rcolm,0.04d0,1)

        write(6,'(''QM-MMPOL fixed charges and induced dipoles '')')
        write(6,'(''mmpol immpol   =  '',t30,i3)') immpol
        write(6,'(''mmpol ich_mmpol =  '',t30,i3)') ich_mmpol
        write(6,'(''mmpol isites  =  '',t30,i3)') isites_mmpol
        write(6,'(''mmpol file (sites)    ='',t30,a20)') mmpolfile_sites
        write(6,'(''mmpol file (chmm)    ='',t30,a20)') mmpolfile_chmm
        write(6,'(''mmpol a_cutoff ='',t30,f7.3)') a_cutoff
        write(6,'(''mmpol rcolm ='',t30,f7.3)') rcolm

        call mmpol_extpot_read
      endif

c Additional properties (<X> etc)
      call p2gtid('properties:sample',iprop,0,1)
CVARDOC flag: properties will be sampled
      call p2gtid('properties:print',ipropprt,0,1)
CVARDOC flag: properties will be printed

      if(iprop.ne.0) then
       nprop=MAXPROP
       write(6,'(''Properties will be sampled'')')
       write(6,'(''Properties printout flag = '',t30,i10)') ipropprt
       call prop_cc_nuc(znuc,cent,iwctype,nctype_tot,ncent_tot,ncent,cc_nuc)
      endif
      
c Pseudopotential section:
      call p2gtid('pseudo:nloc',nloc,0,1)


CVARDOC flag: type of pseudopotential (0: all electron)
      if(nloc.gt.0) then
        write(6,'(/,''pseudopotential calculation, nloc ='',t30,i10)') nloc
        if(nloc.eq.1) then
          call readps
         elseif(nloc.ge.2.and.nloc.le.5) then
          call p2gtid('pseudo:nquad',nquad,6,1)
CVARDOC number of quadrature points
          write(6,'(''nquad='',t30,i10)') nquad
          if(nquad.gt.MPS_QUAD) call fatal_error('INPUT: nquad > MPS_QUAD')
          if(nloc.eq.4)then
            call set_ps_gauss_filenames()
            call readps_gauss
           elseif(nloc.eq.5) then
            call set_ps_champ_filenames()
            call readps_champ
           else
            call set_ps_tm_filenames()
            call readps_tm
          endif
        endif
        call gesqua (nquad,xq,yq,zq,wq)
       else
        write(6,'(/,''all-electron calculation '')')
      endif

c Geometrical section
      call p2gti('atoms:nctype',nctype,1)
CVARDOC number of center types (atom species)
      call p2gti('atoms:natom',ncent,1)
CVARDOC number of centers (atoms)
      call p2gtid('atoms:addghostype',newghostype,0,1)
CVARDOC number of ghost atom types
      call p2gtid('atoms:nghostcent',nghostcent,0,1)
CVARDOC number of ghost centers
      if(max(3,nctype).gt.MCTYP3X) call fatal_error('INPUT: max(3,nctype) > MCTYP3X')
      ! if(nctype+newghostype.gt.MCTYPE) call fatal_error('INPUT: nctype+newghostype > MCTYPE')
      ! if(ncent+nghostcent.gt.MCENT) call fatal_error('INPUT: ncent+nghostcent > MCENT')

      write(6,'(/,''nctype,ncent ='',t30,2i5)') nctype,ncent
      if(newghostype+nghostcent.gt.0)
     & write(6,'(''newghostype,nghostcent ='',t30,2i5)') newghostype,nghostcent
      write(6,'(''iwctype ='',t30,20i3,(20i3))') (iwctype(i),i=1,ncent+nghostcent)

      if(iperiodic.ne.0) then
        call pw_setup_input

       else
c Center positions fetched in read_geometry
        write(6,'(/,''center positions'')')
        do 1 ic=1,ncent
    1     write(6,'(''center'',i4,1x,''('',3f9.5,'')'')') ic,(cent(k,ic),k=1,3)
        write(6,*)
        do 2 ic=1,nghostcent
    2     write(6,'(''ghost '',i4,1x,''('',3f9.5,'')'')') ic,(cent(k,ncent+ic),k=1,3)
        write(6,*)
      endif

c Wave function parameters
c Fetched in read_lcao, read_exponents, read_bas_num_info,read_determinants,read_jastrow_parameter

c Determinantal section
      if(ibasis.eq.1) then
        write(6,'(/,''Orbitals on localized basis'')')
        write(6,'(''total no. of basis ='',t30,i10)') nbasis
      endif
      if(ibasis.eq.2) write(6,'(/,''PW orbitals'')')
      if(numr.gt.0) write(6,'(''numerical basis used'')')

      if(ibasis.eq.1) then
        call write_orb_loc
        call set_bas_num_filenames()
        if(numr.gt.0) then
          do 10 iwft=1,nwftype
   10       call read_bas_num(iwft)
        ibas0(1)=1
        ibas1(1)=nbastyp(iwctype(1))
        do 15 ic=2,ncent
          ibas0(ic)=ibas1(ic-1)+1
   15     ibas1(ic)=ibas1(ic-1)+nbastyp(iwctype(ic))
        endif
      else if(ibasis.eq.2) then
        call read_orb_pw_tm
      endif

      write(6,'(/,''no. of determinants ='',t30,i10)') ndet

      write(6,'(/,''determinant coefficients'')')
      write(6,'(20f10.6)') (cdet(k,1,1),k=1,ndet)
      write(6,'(/,''orbitals in determinants'')')
      do 200 i=1,ndet
        do 200 j=1,nelec
  200     if(iworbd(j,i).gt.norb) call fatal_error('INPUT: iworbd > norb')
      do 230 i=1,ndet
        if(nup*ndn.gt.0) then
         write(fmt,'(''('',i3,''i3,3x,'',i3,''i3)'')') nup,ndn
        else
         write(fmt,'(''('',i3,''i3)'')') max(nup,ndn)
        endif
  230   write(6,fmt) (iworbd(j,i),j=1,nup),(iworbd(j+nup,i),j=1,ndn)

c Jastrow section
      call p2gtid('jastrow:ianalyt_lap',ianalyt_lap,1,1)
      write(6,'(''ianalyt_lap='',i3)') ianalyt_lap
      if(ianalyt_lap.eq.0) then
        if(nloc.gt.0)
     &  call fatal_error('No numerical jastrow derivatives with pseudopotentials')
        if(iperiodic.gt.0)
     &  call fatal_error('No numerical jastrow derivatives with periodic system: distances in jastrow_num not correct')
        if(ioptjas.gt.0)
     &  call fatal_error('No numerical jastrow derivatives and parms derivatives')
      endif

      if(ijas.ne.4 .and. iperiodic.gt.0)
     &call fatal_error('Only ijas=4 implemented for periodic systems')

      write(6,'(''ijas,isc,nspin1,nspin2='',5i4)')
     & ijas,isc,nspin1,nspin2

      if(ijas.eq.4) write(6,'(''new transferable standard form 4'')')
      if(ijas.eq.5) write(6,'(''new transferable standard form 5'')')
      if(ijas.eq.6) write(6,'(''new transferable standard form 6'')')

      if(isc.eq.2) write(6,'(
     &''dist scaled r=(1-exp(-scalek*r))/scalek'')')
      if(isc.eq.3) write(6,'(
     &''dist scaled r=(1-exp(-scalek*r-(scalek*r)**2/2))/scalek'')')
      if(isc.eq.4) write(6,'(
     &''dist scaled r=r/(1+scalek*r)'')')
      if(isc.eq.5) write(6,'(
     &''dist scaled r=r/(1+(scalek*r)**2)**.5'')')

      if(ijas.ge.4.and.ijas.le.6) then
        write(6,'(''norda,nordb,nordc='',3i5)') norda,nordb,nordc
        if(isc.ge.2) write(6,'(''scalek,a21='',2f10.5)') scalek(1),a21

        mparmja=2+max(0,norda-1)
        mparmjb=2+max(0,nordb-1)
        mparmjc=nterms4(nordc)
        write(6,'(''mparmja,mparmjb,mparmjc='',3i5)') mparmja,mparmjb,mparmjc
        do 301 it=1,nctype
  301     write(6,'(''a='',x,7f10.6,(8f10.6))')
     &                  (a4(iparm,it,1),iparm=1,mparmja)
        do 302 isp=nspin1,nspin2b
  302     write(6,'(''b='',x,7f10.6,(8f10.6))')
     &                 (b(iparm,isp,1),iparm=1,mparmjb)
        do 303 it=1,nctype
  303     write(6,'(''c='',x,7f10.6,(8f10.6))')
     &                 (c(iparm,it,1),iparm=1,mparmjc)
        
c Note: Fock terms yet to be put in ijas=4,5,6
      endif

c Call set_scale_dist to evaluate constants that need to be reset if
c scalek is being varied. If cutjas=0, then reset cutjas to infinity
c Warning: At present we are assuming that the same scalek is used for
c primary and secondary wavefns.  Otherwise c1_jas6i,c1_jas6,c2_jas6
c should be dimensioned to MWF
      if(isc.eq.6.or.isc.eq.7.or.isc.eq.16.or.isc.eq.17) then
        if(iperiodic.ne.0 .and. cutjas_tmp.gt.cutjas) then
          write(6,'(''**Warning: input cutjas > half shortest sim. cell lattice vector;
     &    cutjas reset from'',f9.5,'' to'',f9.5)') cutjas_tmp,cutjas
         else
          cutjas=cutjas_tmp
          write(6,'(''input cutjas='',d12.5)') cutjas_tmp
        endif
        if(cutjas.gt.0.d0) then
          cutjasi=1/cutjas
         else
          write(6,'(''cutjas reset to infinity'')')
          cutjas=1.d99
          cutjasi=0
        endif
        call set_scale_dist(1)
       else
        cutjas=1.d99
        cutjasi=0
        c1_jas6i=1
        c1_jas6=1
        c2_jas6=0
        asymp_r=0
        do 310 it=1,nctype
 310      asymp_jasa(it)=0
        do 320 i=1,2
 320      asymp_jasb(i)=0
      endif
      call set_scale_dist(1)

c Write out information about calculation of energy gradients 
c and Z matrix
      write(6,*)
      if(ngradnts.gt.0 .and. igrdtype.eq.1) call inpwrt_grdnts_cart()
      if(ngradnts.gt.0 .and. igrdtype.eq.2) call inpwrt_grdnts_zmat()
      if(izmatrix.eq.1) call inpwrt_zmatrix()
      write(6,*)

c get normalization for basis functions
      if(numr.eq.0) then
        do 330 iwft=1,nwftype
  330     call basis_norm(iwft,anorm,0)
      endif
c check if the orbitals coefficients are to be multiplied by a constant parameter
      call p2gtfd('general:scalecoef',scalecoef,1.0d0,1)
      ! call resize_tensor(coef, norb+nadorb, 2)
      if(scalecoef.ne.1.0d0) then
        do 340 iwft=1,nwftype
          do 340 iorb=1,norb+nadorb
	    do 340 j=1,nbasis
  340         coef(j,iorb,iwft)=coef(j,iorb,iwft)*scalecoef
        write(6,'(/,''Orbital coefficients scaled by a constant parameter = '',f7.4)') scalecoef
	write(6,*)
      endif

c get nuclear potential energy
c     call pot_nn(cent,znuc,iwctype,ncent,pecent)

c verify number of orbitals and setup optorb
      call verify_orbitals
      
c get parameters for the grid of the orbitals
      call p2gtid('general:i3dgrid',i3dgrid,0,1)
      call p2gtid('general:i3dsplorb',i3dsplorb,0,1)
      call p2gtid('general:i3dlagorb',i3dlagorb,0,1)
      call p2gtid('general:i3ddensity',i3ddensity,0,1)

      if((i3dsplorb.ge.1).or.(i3dlagorb.ge.1).or.(i3ddensity.ge.1))
     & i3dgrid=1

      if(i3dgrid.ge.1) then 

c Read the grid information: 
        call p2gtid('3dgrid:nstepx',nstep3d(1),IUNDEFINED,1)
        call p2gtid('3dgrid:nstepy',nstep3d(2),IUNDEFINED,1)
        call p2gtid('3dgrid:nstepz',nstep3d(3),IUNDEFINED,1)

        call p2gtfd('3dgrid:stepx',step3d(1),UNDEFINED,1)
        call p2gtfd('3dgrid:stepy',step3d(2),UNDEFINED,1)
        call p2gtfd('3dgrid:stepz',step3d(3),UNDEFINED,1)
 
        call p2gtfd('3dgrid:x0',origin(1),UNDEFINED,1)
        call p2gtfd('3dgrid:y0',origin(2),UNDEFINED,1)
        call p2gtfd('3dgrid:z0',origin(3),UNDEFINED,1)
 
        call p2gtfd('3dgrid:xn',endpt(1),UNDEFINED,1)
        call p2gtfd('3dgrid:yn',endpt(2),UNDEFINED,1)
        call p2gtfd('3dgrid:zn',endpt(3),UNDEFINED,1)

C Grid setup:
        call setup_grid

        if(i3dlagorb.ge.1) then
          write(6,'(/,''orbitals on a grid: splines interpolation'')')
          call setup_3dlagorb
          i3dsplorb=0
         elseif(i3dsplorb.ge.1) then
          write(6,'(/,''orbitals on a grid: Lagrange interpolation'')')
          call setup_3dsplorb
          i3dlagorb=0
        endif
      endif

      return
      end
 
c-----------------------------------------------------------------------
      subroutine flaginit
c Initialize flags used to identify presence/absence of blocks in input

      use inputflags, only: iznuc, igeometry, ibasis_num, ilcao, iexponents
      use inputflags, only: ideterminants, ijastrow_parameter, ioptorb_def, ilattice
      use inputflags, only: ici_def, iforces, icsfs, icharge_efield
      use inputflags, only: imultideterminants, imodify_zmat, izmatrix_check
      use inputflags, only: ihessian_zmat

      implicit real*8(a-h,o-z)

      iznuc=0
      igeometry=0
      ibasis_num=0
      ilcao=0
      iexponents=0
      ideterminants=0
      ijastrow_parameter=0
      ioptorb_def=0
      ilattice=0
      ici_def=0
      iforces=0
      icsfs=0
      icharge_efield=0
      imultideterminants=0
      izmatrix_check=0
      imodify_zmat=0
      ihessian_zmat=0

      return
      end

c-----------------------------------------------------------------------
      subroutine flagcheck
c Check that the required blocks are there in the input

      use force_mod, only: MFORCE, MWF
      use vmc_mod, only: MELEC, MORB
      use numbas, only: numr
      use optorb_mix, only: norbopt, norbvirt
      use efield, only: iefield
      use inputflags, only: iznuc, igeometry, ibasis_num, ilcao, iexponents
      use inputflags, only: ideterminants, ijastrow_parameter, ioptorb_def, ilattice
      use inputflags, only: ici_def, iforces, icsfs, igradients, icharge_efield
      use inputflags, only: imultideterminants, ioptorb_mixvirt, imodify_zmat, izmatrix_check
      use inputflags, only: ihessian_zmat
      use mstates_ctrl, only: iguiding
      ! might not be needed
      use mstates_mod, only: MSTATES
      use atom, only: znuc 
      use contrl_per, only: iperiodic, ibasis
      use force_analy, only: iforce_analy, iuse_zmat
      use forcepar, only: nforce
      use optwf_contrl, only: ioptci, ioptorb
      use optwf_contrl, only: no_active
      use wfsec, only: nwftype
      use orbval, only: ddorb, dorb, nadorb, ndetorb, orb
      use elec, only: ndn, nup
      use const, only: nelec 
      use coefs, only: norb, next_max 

      implicit real*8(a-h,o-z)

      call p2gti('electrons:nelec',nelec,1)
      call p2gti('electrons:nup',nup,1)
      call p2gtid('general:nwftype',nwftype,1,1)
      call p2gtid('general:nforce',nforce,1,1)
      ! if(nforce.gt.MFORCE) call fatal_error('INPUT: nforce > MFORCE')
      call p2gtid('general:nwftype',nwftype,1,1)
      !if(nwftype.gt.MWF) call fatal_error('INPUT: nwftype gt MWF')
      call p2gtid('general:iperiodic',iperiodic,0,1)
      call p2gtid('general:ibasis',ibasis,1,1)

      call p2gtid('optwf:ioptorb',ioptorb,0,1)
      call p2gtid('optwf:ioptci',ioptci,0,1)
      call p2gtid('optwf:no_active',no_active,0,1)

      call p2gtid('mstates:iguiding',iguiding,0,1)
      call p2gtid('efield:iefield',iefield,0,1)
      call p2gtid('optgeo:iforce_analy',iforce_analy,0,0)
      call p2gtid('optgeo:iuse_zmat',iuse_zmat,0,0)
      next_max=norb-ndetorb
      call p2gtid('optwf:nextorb',nadorb, next_max,1)
      if (nadorb.gt.norb) call fatal_error('nadorb > norb')
      
      

      if(iznuc.eq.0) call fatal_error('INPUT: block znuc missing')
      if(igeometry.eq.0) call fatal_error('INPUT: block geometry missing')
      if(ibasis.eq.1.and.numr.gt.0.and.ibasis_num.eq.0) call fatal_error('INPUT: block basis missing')
      if(iperiodic.eq.0.and.ilcao.eq.0) call fatal_error('INPUT: block lcao missing')
      if(iperiodic.gt.0.and.ilattice.eq.0) call fatal_error('INPUT: lattice vectors missing')
      if(ijastrow_parameter.eq.0) call fatal_error('INPUT: block jastrow_parameter missing')
      if(iefield.gt.0.and.icharge_efield.eq.0) call fatal_error('INPUT: block efield missing')
       
      write(6,'(''========================================'')')
      if(iexponents.eq.0) then 
        write(6,'(''INPUT: block exponents missing: all exponents set to 1'')')
        call inputzex
      endif
      if(icsfs.eq.0) then
        write(6,'(''INPUT: block csf missing: nstates set to 1'')')
        call inputcsf
      endif
      if(nforce.ge.1.and.iforces.eq.0.and.igradients.eq.0) then
        write(6,'(''INPUT: block forces_displace or gradients_* missing: geometries set equal to primary'')')
        call inputforces
      endif
      if(iforce_analy.gt.0) then
        if(iuse_zmat.gt.0.and.izmatrix_check.eq.0) call fatal_error('INPUT: block connectionzmatrix missing')
        if(imodify_zmat.eq.0) call modify_zmat_define
        if(ihessian_zmat.eq.0) call hessian_zmat_define
      endif
      if(imultideterminants.eq.0) then
        write(6,'(''INPUT: multideterminant bloc MISSING'')')
        call multideterminants_define(0,0)
      endif
      if(ioptorb.ne.0) then
        if(ioptorb_mixvirt.eq.0) then
          norbopt=0
          norbvirt=0
        endif
        if(ioptorb_def.eq.0) then
          write(6,'(''INPUT: definition of orbital variations missing'')')
          call optorb_define
        endif
      endif
      if(ioptci.ne.0.and.ici_def.eq.0) then
        write(6,'(''INPUT: definition of OPTCI operators missing'')')
        call optci_define
      endif

      if(nwftype.gt.1) then
        if(ijastrow_parameter.ne.nwftype) then
          write(6,'(''INPUT: block jastrow_parameter missing for one wave function'')')
          write(6,'(''INPUT: jastrow_parameter blocks equal for all wave functions'')')
          call inputjastrow(nwftype)
        endif
        if(iperiodic.eq.0.and.ilcao.ne.nwftype) then
          write(6,'(''Warning INPUT: block lcao missing for one wave function'')')
          write(6,'(''Warning INPUT: lcao blocks equal for all wave functions'')')
          call inputlcao(nwftype)
        endif
        if(ideterminants.ne.nwftype) then
          write(6,'(''Warning INPUT: block determinants missing for one wave function'')')
          write(6,'(''Warning INPUT: determinants blocks equal for all wave functions'')')
          call inputdet(nwftype)
        endif
        write(6,*)
      endif

      write(6,'(''========================================'')')
      return
      end

c-----------------------------------------------------------------------

      subroutine set_ps_gauss_filenames()
c ### Set name files of gaussian pseudopotentials.
      use general, only: pooldir, pp_id, atomtyp, filename, atomsymbol
      use general, only: filenames_ps_gauss
      use atom, only: nctype
      implicit real*8(a-h,o-z)
c Allocation of the array storing the filenames of gaussian basis: 
      allocate(filenames_ps_gauss(nctype))
      do ic=1,nctype
        if(ic.lt.10) then
          write(atomtyp,'(i1)') ic
        elseif(ic.lt.100) then
          write(atomtyp,'(i2)') ic
        endif
        if(pp_id.eq.'none') then
            call fatal_error('READ_INPUT: ECP name missing')
           else
            call p2gtad('atom_types:'//atomtyp(1:index(atomtyp,' ')-1)
     &         ,atomsymbol,'X',1)
             filename=pooldir(1:index(pooldir,' ')-1)//
     &         '/'//
     &         pp_id(1:index(pp_id,' ')-1)//
     &         '.gauss_ecp.dat.'//
     &         atomsymbol(1:index(atomsymbol,' ')-1)
         endif 
         filenames_ps_gauss(ic)=filename
        enddo
        end subroutine

c-----------------------------------------------------------------------

      subroutine set_ps_champ_filenames()
c ### Set name files of CHAMP-formatted pseudopotentials.
      use atom, only: nctype
      use general, only: pooldir, pp_id, atomtyp, filename, atomsymbol
      use general, only: filenames_ps_champ 
      implicit real*8(a-h,o-z)
c Allocation of the array storing the filenames of gaussian basis: 
      allocate(filenames_ps_champ(nctype))
      
      do ict=1, nctype
        if(ict.lt.10) then
          write(atomtyp,'(i1)') ict
         elseif(ict.lt.100) then
          write(atomtyp,'(i2)') ict
        endif
        if(pp_id.eq.'none') then
c old naming convention
          filename=pooldir(1:index(pooldir,' ')-1)//'/'//
     &               'pseudopot_champ'//atomtyp(1:index(atomtyp,' ')-1)
         else
c new naming convention
          call p2gtad('atom_types:'//atomtyp(1:index(atomtyp,' ')-1)
     &         ,atomsymbol,'X',1)
          filename=pooldir(1:index(pooldir,' ')-1)//
     &             '/'//
     &             pp_id(1:index(pp_id,' ')-1)//
     &             '.pseudopot_champ.'//
     &             atomsymbol(1:index(atomsymbol,' ')-1)
        endif      
        filenames_ps_champ(ict)=filename
      enddo
      end subroutine

c-----------------------------------------------------------------------

      subroutine set_ps_tm_filenames()
c ### Set name files of Troullier-Martins pseudopotentials.
      use atom, only: nctype
      use general, only: pooldir, pp_id, atomtyp, filename, atomsymbol
      use general, only: filenames_ps_tm
      use pseudo, only: nloc
      implicit real*8(a-h,o-z)

      allocate(filenames_ps_tm(nctype))
      do ic=1,nctype
        if(pp_id.eq.'none') then
c old naming convention
          if(nloc.eq.2) then
            filename=pooldir(1:index(pooldir,' ')-1)//'/'//
     &           'pseudo.dat.'//atomtyp(1:index(atomtyp,' ')-1)
           elseif(nloc.eq.3) then
            filename=pooldir(1:index(pooldir,' ')-1)//'/'//
     &           'pseudopot'//atomtyp(1:index(atomtyp,' ')-1)
          endif
        else
c new naming convention
          call p2gtad('atom_types:'//atomtyp(1:index(atomtyp,' ')-1)
     &     ,atomsymbol,'X',1)
          if(nloc.eq.2) then
            filename=pooldir(1:index(pooldir,' ')-1)//
     &           '/'//
     &           pp_id(1:index(pp_id,' ')-1)//
     &           '.pseudo.dat.'//
     &           atomsymbol(1:index(atomsymbol,' ')-1)
           elseif(nloc.eq.3) then
            filename=pooldir(1:index(pooldir,' ')-1)//
     &           '/'//
     &           pp_id(1:index(pp_id,' ')-1)//
     &           '.pseudopot.'//
     &           atomsymbol(1:index(atomsymbol,' ')-1)
          endif
        endif
        filenames_ps_tm(ic)=filename
      enddo
      end subroutine

c-----------------------------------------------------------------------

      subroutine set_bas_num_filenames()
c ### Set numerical num. orbital filenames.
      use atom, only: nctype
      use general, only: pooldir, pp_id, bas_id, atomtyp, filename, atomsymbol
      use general, only: filenames_bas_num, wforce 
      use ghostatom, only: newghostype
      implicit real*8(a-h,o-z)
c Allocation of the array storing the filenames of numerical basis: 
      allocate(filenames_bas_num(nctype+newghostype))
  
      do ic=1,nctype+newghostype
        if(ic.lt.10) then
          write(atomtyp,'(i1)') ic
         elseif(ic.lt.100) then
          write(atomtyp,'(i2)') ic
         elseif(iwf.lt.1000) then
           write(wforce,'(i3)') iwf
         endif
        if(bas_id.eq.'none') then
c old file name convention 
          filename=pooldir(1:index(pooldir,' ')-1)//'/'//
     &             'basis.'//atomtyp(1:index(atomtyp,' ')-1)
          if(iwf.ge.2) then
            filename=filename(1:index(filename,' ')-1)//'.'//wforce
          endif
        else
c new convention
          call p2gtad('atom_types:'//atomtyp(1:index(atomtyp,' ')-1)
     &              ,atomsymbol,'X',1)
          filename=pooldir(1:index(pooldir,' ')-1)//
     &           '/'//
     &           bas_id(1:index(bas_id,' ')-1)//
     &           '.basis.'//
     &           atomsymbol(1:index(atomsymbol,' ')-1)
        endif
        filenames_bas_num(ic)=filename
      enddo
      end subroutine