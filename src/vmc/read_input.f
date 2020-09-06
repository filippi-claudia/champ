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

c Transfer from lists to fortran variables, print out, check,
c and read in everything which is still in the old format
      call process_input

      if(index(mode,'mc').ne.0 ) call p2vin('input.log',1)

      return
      end

c-----------------------------------------------------------------------
      subroutine process_input
c Written by Cyrus Umrigar, Claudia Filippi, Friedemann Schautz,
c and Anthony Scemema
      use atom, only: znuc, cent, pecent, iwctype, nctype, ncent
      use jaspar, only: nspin1, nspin2, is
      use ghostatom, only: newghostype, nghostcent
      use const, only: pi, hb, etrial, delta, deltai, fbias, nelec, imetro, ipr
      use jaspar1, only: cjas1, cjas2
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
      use optwf_contrl, only: ioptci, ioptjas, ioptorb
      use optwf_parms, only: nparmj
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
      use dorb_m, only: iworbd
      use contrl_per, only: iperiodic, ibasis
      use force_analy, only: iforce_analy, iuse_zmat, alfgeo
      use pseudo, only: nloc
      use optorb_cblock, only: idump_blockav
      use gradjerrb, only: ngrad_jas_blocks
      use qua, only: nquad, wq, xq, yq, zq
      use mmpol_cntrl, only: ich_mmpol, immpol, immpolprt, isites_mmpol
      use mmpol_parms, only: chmm
      use mmpol_fdc, only: a_cutoff, rcolm
      use grid3dflag, only: i3ddensity, i3dgrid, i3dlagorb, i3dsplorb
      use efield, only: iefield, iscreen, ncharges
      use mstates_ctrl, only: iefficiency, iguiding, nstates_psig
      use mstates3, only: iweight_g, weights_g
      use ci000, only: iciprt, nciprim, nciterm
      use pcm_cntrl, only: icall, ichpol, ipcm, ipcmprt, isurf
      use pcm_unit, only: pcmfile_cavity, pcmfile_chs, pcmfile_chv
      use pcm_parms, only: ch, eps_solv, iscov, nch, nchs, nchs1, nchs2
      use pcm_parms, only: nchv, ncopcm, nesph, nscv, nvopcm, re, re2
      use pcm_parms, only: retk, surk, xe, xpol, ye, ze

      use prp000, only: iprop, ipropprt, nprop
      use pcm_fdc, only: feps, fs, qfree, qvol, rcol, rcolt, rcolv
      implicit real*8(a-h,o-z)



      parameter (zero=0.d0,one=1.d0,two=2.d0,four=4.d0)

      include 'vmc.h'
      include 'dmc.h'
      include 'force.h'
      include 'pseudo.h'
      include 'numbas.h'
      include 'optjas.h'
      include 'optci.h'
      include 'optorb.h'
      include 'sr.h'
      include 'ewald.h'
      include 'mmpol.h'
      include 'pcm.h'
      include 'pcm_3dgrid.h'
      include 'efield.h'
      include 'mstates.h'
      include 'properties.h'


      common /orbval/ orb(MELEC,MORB),dorb(3,MELEC,MORB),ddorb(MELEC,MORB),ndetorb,nadorb

      character*20 fmt
      character*32 keyname
      character*10 eunit
      character*16 cseed
      character*20 dl_alg
      dimension irn(4),cent_tmp(3),anorm(MBASIS)

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
      call stripquotes(title)

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
      if(nforce.gt.MFORCE) call fatal_error('INPUT: nforce > MFORCE')
      call p2gtid('general:nwftype',nwftype,1,1)
      write(6,'(/,''number of wave functions='',t30,i10)') nwftype
      if(nwftype.gt.nforce) call fatal_error('INPUT: nwftype gt nforce')
      if(nwftype.gt.MWF) call fatal_error('INPUT: nwftype gt MWF')

c Electron section
      call p2gti('electrons:nelec',nelec,1)
      write(6,'(''number of electrons ='',t30,i10)') nelec
      if(nelec.gt.MELEC) call fatal_error('INPUT: nelec exceeds MELEC')

      call p2gti('electrons:nup',nup,1)
      if(nup.gt.MELEC/2) call fatal_error('INPUT: nup exceeds MELEC/2')
      ndn=nelec-nup
      write(6,'(''number of up,dn electrons ='',t30,2i5)') nup,ndn

c VMC parameters
      if(index(mode,'vmc').ne.0) then

        call p2gtid('vmc:imetro',imetro,6,1)
        write(6,'(''version of Metropolis ='',t30,i10)') imetro

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
        if(nloc.eq.0) call fatal_error('INPUT: no all-electron calculations supported')

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
        call p2gtad('optwf:dl_alg',dl_alg,'nag',1)
        call p2gtid('optwf:nblk_max',nblk_max,nblk,1)

        if(method.eq.'linear'.and.MXREDUCED.ne.MXORBOP) 
     &    call fatal_error('READ_INPUT: MXREDUCED.ne.MXORBOP')
        if((method.eq.'sr_n'.or.method.eq.'lin_d').and.nstep*nblk_max.gt.MCONF)
     &    call fatal_error('READ_INPUT: nstep*nblk_max.gt.MCONF')
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

        call p2gtid('pcm:ipcm_3dgrid',ipcm_3dgrid,0,1)
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
       call prop_cc_nuc(znuc,cent,iwctype,MCTYPE,MCENT,ncent,cc_nuc)
      endif
      
c Pseudopotential section
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
            call readps_gauss
           elseif(nloc.eq.5) then
            call readps_champ
           else
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
      if(nctype+newghostype.gt.MCTYPE) call fatal_error('INPUT: nctype+newghostype > MCTYPE')
      if(ncent+nghostcent.gt.MCENT) call fatal_error('INPUT: ncent+nghostcent > MCENT')

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
      subroutine read_znuc(iu)
C$INPUT znuc inp
CKEYDOC nuclear charge for each atom type and ghost type

      use atom, only: znuc, nctype
      use inputflags, only: iznuc,igeometry,ibasis_num,ilcao,iexponents,
     &             ideterminants,ijastrow_parameter, ioptorb_def,ilattice,
     &             ici_def,iforces,icsfs,imstates,igradients,icharge_efield,
     &             imultideterminants,ioptorb_mixvirt,imodify_zmat,izmatrix_check,
     &             ihessian_zmat 

      implicit real*8(a-h,o-z)

      include 'vmc.h'

      call p2gti('atoms:nctype',nctype,1)
      call p2gtid('atoms:addghostype',newghostype,0,1)
      if(nctype+newghostype.gt.MCTYPE) call fatal_error('INPUT: nctype+newghostype > MCTYPE')

      call incpos(iu,itmp,1)
      read(iu,*) (znuc(i),i=1,nctype+newghostype)
      iznuc=1
      call p2chkend(iu, 'znuc')
      end

c-----------------------------------------------------------------------
      subroutine read_lcao(norb_tmp,nbasis_tmp,iwft,filename)
C$INPUT lcao i i i=1 a=<input> 
CKEYDOC Orbital coefficients wrt complete basis.
CKEYDOC Usage:  {\tt lcao  norb,nbasis,filename,norbv}
CKEYDOC norb: number of orbitals for trial wave function
CKEYDOC nbasis: number of basis functiobns
CKEYDOC iwft: wave function type (used when nforce>1 and wftype>1)
CKEYDOC filename: file containing orbitals coefficients

      use coefs, only: coef, nbasis, norb
      use inputflags, only: iznuc,igeometry,ibasis_num,ilcao,iexponents,
     &             ideterminants,ijastrow_parameter, ioptorb_def,ilattice,
     &             ici_def,iforces,icsfs,imstates,igradients,icharge_efield,
     &             imultideterminants,ioptorb_mixvirt,imodify_zmat,izmatrix_check,
     &             ihessian_zmat 

      use pcm_fdc, only: feps, fs, qfree, qvol, rcol, rcolt, rcolv
      implicit real*8(a-h,o-z)


      include 'vmc.h'
      include 'force.h'

c fs NOTE: additional variable norbv for efp orbitals removed 


      character filename*(*)

      call file(iu,filename,'old',1,0)
      nbasis=nbasis_tmp
      norb=norb_tmp
      nototal=norb
      if(nbasis.gt.MBASIS) call fatal_error('LCAO: nbasis > MBASIS')
      if(nototal.gt.MORB) call fatal_error('LCAO: number of orbitals > MORB')

      call p2gtid('general:nwftype',nwftype,1,1)
      if(iwft.gt.nwftype) call fatal_error('LCAO: wave function type > nwftype')

      do 20 i=1,nototal
        call incpos(iu,itmp,1)
   20   read(iu,*) (coef(j,i,iwft),j=1,nbasis)
      ilcao=ilcao+1
      if(filename.eq.'<input>') then
       call p2chkend(iu, 'lcao')
      endif
      end

c-----------------------------------------------------------------------
      subroutine read_geometry(iu)
C$INPUT geometry inp
CKEYDOC position and type for each atom and ghost atom

      use atom, only: cent, iwctype, ncent
      use inputflags, only: iznuc,igeometry,ibasis_num,ilcao,iexponents,
     &             ideterminants,ijastrow_parameter, ioptorb_def,ilattice,
     &             ici_def,iforces,icsfs,imstates,igradients,icharge_efield,
     &             imultideterminants,ioptorb_mixvirt,imodify_zmat,izmatrix_check,
     &             ihessian_zmat 

      implicit real*8(a-h,o-z)

      include 'vmc.h'

      call p2gti('atoms:natom',ncent,1)
      call p2gtid('atoms:nghostcent',nghostcent,0,1)
      if(ncent+nghostcent.gt.MCENT) call fatal_error('INPUT: ncent+nghostcent > MCENT')

      do 20 i=1,ncent+nghostcent
        call incpos(iu,itmp,1)
  20    read(iu,*) (cent(k,i),k=1,3),iwctype(i)

      igeometry=1
      call p2chkend(iu, 'geometry')
      end

c-----------------------------------------------------------------------

      subroutine read_exponents(iu,iwft)
C$INPUT exponents inp i=1
CKEYDOC Basis function exponents (only if no numerical basis)

      use coefs, only: nbasis
      use basis, only: zex
      use inputflags, only: iznuc,igeometry,ibasis_num,ilcao,iexponents,
     &             ideterminants,ijastrow_parameter, ioptorb_def,ilattice,
     &             ici_def,iforces,icsfs,imstates,igradients,icharge_efield,
     &             imultideterminants,ioptorb_mixvirt,imodify_zmat,izmatrix_check,
     &             ihessian_zmat 

      implicit real*8(a-h,o-z)

      include 'vmc.h'
      include 'force.h'

      call incpos(iu,itmp,1)
      read(iu,*) (zex(i,iwft),i=1,nbasis)
      iexponents=iexponents+1
      call p2chkend(iu, 'exponents')
      end

c-----------------------------------------------------------------------
      subroutine read_determinants(iu,nd,iwft)
C$INPUT determinants inp i i=1
CKEYDOC CI coefficients and occupation of determinants in wf

      use dets, only: cdet, ndet
      use dorb_m, only: iworbd
      use inputflags, only: iznuc,igeometry,ibasis_num,ilcao,iexponents,
     &             ideterminants,ijastrow_parameter, ioptorb_def,ilattice,
     &             ici_def,iforces,icsfs,imstates,igradients,icharge_efield,
     &             imultideterminants,ioptorb_mixvirt,imodify_zmat,izmatrix_check,
     &             ihessian_zmat 

      implicit real*8(a-h,o-z)

      include 'vmc.h'
      include 'force.h'
      include 'mstates.h'

      ndet=nd
      if(ndet.gt.MDET) then
        write (6,*)  "ndet=", ndet
        write (6,*)  "MDET=", MDET
        call fatal_error('DET: ndet > MDET')
      endif

      call p2gti('electrons:nelec',nelec,1)
      if(nelec.gt.MELEC) call fatal_error('INPUT: nelec exceeds MELEC')
      call incpos(iu,itmp,1)
      read(iu,*) (cdet(i,1,iwft),i=1,ndet)
c     if(iwft.eq.1) then
        do 20 i=1,ndet
          call incpos(iu,itmp,1)
   20     read(iu,*) (iworbd(j,i),j=1,nelec)
c     endif
      ideterminants=ideterminants+1
      call p2chkend(iu, 'determinants')
      end

c-----------------------------------------------------------------------
      subroutine read_multideterminants(iu,nd)
C$INPUT multideterminants inp i 
CKEYDOC CI coefficients and occupation of determinants in wf
      use dets, only: ndet
      use multidet, only: irepcol_det, ireporb_det, numrep_det
      use inputflags, only: iznuc,igeometry,ibasis_num,ilcao,iexponents,
     &             ideterminants,ijastrow_parameter, ioptorb_def,ilattice,
     &             ici_def,iforces,icsfs,imstates,igradients,icharge_efield,
     &             imultideterminants,ioptorb_mixvirt,imodify_zmat,izmatrix_check,
     &             ihessian_zmat 

      implicit real*8(a-h,o-z)

      include 'vmc.h'
      include 'force.h'
      include 'mstates.h'

      if(nd.ne.ndet-1) call fatal_error('INPUT: problem in multidet')

      call incpos(iu,itmp,1)
      do k=2,nd+1
        read(iu,*) (numrep_det(k,iab),iab=1,2)
        do iab=1,2
          do irep=1,numrep_det(k,iab)
            read(iu,*) irepcol_det(irep,k,iab),ireporb_det(irep,k,iab)
          enddo
        enddo
      enddo

      imultideterminants=imultideterminants+1
      call p2chkend(iu, 'multideterminants')
      end

c-----------------------------------------------------------------------
      subroutine read_jastrow_parameter(iu,iwft)
C$INPUT jastrow_parameter inp i=1
CKEYDOC Parameters of Jastrow factor (depends on value of ijas!)

      use jaspar, only: nspin1, nspin2
      use elec, only: ndn
      use jaspar3, only: a, b, c, scalek
      use jaspar4, only: a4, norda, nordb, nordc
      use jaspar6, only: cutjas
      use bparm, only: nocuspb, nspin2b
      use contr2, only: ifock, ijas
      use contr2, only: isc
      use inputflags, only: iznuc,igeometry,ibasis_num,ilcao,iexponents,
     &             ideterminants,ijastrow_parameter, ioptorb_def,ilattice,
     &             ici_def,iforces,icsfs,imstates,igradients,icharge_efield,
     &             imultideterminants,ioptorb_mixvirt,imodify_zmat,izmatrix_check,
     &             ihessian_zmat 

      implicit real*8(a-h,o-z)

      include 'vmc.h'
      include 'force.h'

      call p2gti('jastrow:ijas',ijas,1)
      call p2gti('jastrow:isc',isc,1)
      call p2gtid('jastrow:nspin1',nspin1,1,1)
      call p2gtid('jastrow:nspin2',nspin2,1,1)
      call p2gtid('jastrow:ifock',ifock,0,1)

      call p2gti('atoms:natom',ncent,1)
      call p2gti('atoms:nctype',nctype,1)

      if(ijas.lt.4.or.ijas.gt.6) call fatal_error('JASTROW: only ijas=4,5,6 implemented')
      if(ndn.eq.1.and.nspin2.eq.3) call fatal_error('JASTROW: 1 spin down and nspin2=3')
      if((ijas.eq.4.or.ijas.eq.5).and.
     &(isc.ne.2.and.isc.ne.4.and.isc.ne.6.and.isc.ne.7.and.
     &isc.ne.12.and.isc.ne.14.and.isc.ne.16.and.isc.ne.17))
     & call fatal_error('JASTROW: if ijas=4 or 5, isc must be one of 2,4,6,7,12,14,16,17')
      if((ijas.eq.6).and.(isc.ne.6.and.isc.ne.7))
     & call fatal_error('JASTROW: if ijas=6, isc must be 6 or 7')

      nspin2b=iabs(nspin2)
      nocuspb=0
      if(nspin2.lt.0) then
        if(nspin2.eq.-1) nocuspb=1
        nspin2=1
      endif

      if(ijas.ge.4.and.ijas.le.6) then
        if(ifock.gt.0) call fatal_error('JASTROW: fock not yet implemented for ijas=4,5,6')
        read(iu,*) norda,nordb,nordc
        if(isc.ge.2) read(iu,*) scalek(iwft),a21
        mparmja=2+max(0,norda-1)
        mparmjb=2+max(0,nordb-1)
        mparmjc=nterms4(nordc)
        do 70 it=1,nctype
          read(iu,*) (a4(iparm,it,iwft),iparm=1,mparmja)
   70     call incpos(iu,itmp,1)
        do 80 isp=nspin1,nspin2b
          read(iu,*) (b(iparm,isp,iwft),iparm=1,mparmjb)
   80     call incpos(iu,itmp,1)
        do 90 it=1,nctype
          read(iu,*) (c(iparm,it,iwft),iparm=1,mparmjc)
   90     call incpos(iu,itmp,1)
      endif
c Read cutoff for Jastrow4,5,6
      if(isc.eq.6.or.isc.eq.7) read(iu,*) cutjas

      ijastrow_parameter=ijastrow_parameter+1
      call p2chkend(iu, 'jastrow_parameter')
      end

c-----------------------------------------------------------------------
      subroutine read_bas_num_info(iu,numeric)
C$INPUT basis inp i
CKEYDOC Basis function types and pointers to radial parts tables
C$INPUT qmc_bf_info inp i
CKEYDOC alternative name for keyword basis because of GAMBLE input
      use numbas, only: iwrwf, numr
      use numbas1, only: iwlbas, nbastyp
      use basis, only: n1s, n2s, n2p, n3s, n3p, n3dzr, n3dx2, n3dxy, n3dxz, n3dyz
      use basis, only: n4s, n4p, n4fxxx, n4fyyy, n4fzzz, n4fxxy, n4fxxz, n4fyyx, n4fyyz
      use basis, only: n4fzzx, n4fzzy, n4fxyz, nsa, npa, ndzra, ndxya, ndxza, ndyza, ndx2a
      use inputflags, only: iznuc,igeometry,ibasis_num,ilcao,iexponents,
     &             ideterminants,ijastrow_parameter, ioptorb_def,ilattice,
     &             ici_def,iforces,icsfs,imstates,igradients,icharge_efield,
     &             imultideterminants,ioptorb_mixvirt,imodify_zmat,izmatrix_check,
     &             ihessian_zmat 

      implicit real*8(a-h,o-z)


      include 'vmc.h'
      include 'force.h'
      include 'numbas.h'

      call p2gti('atoms:nctype',nctype,1)
      call p2gtid('atoms:addghostype',newghostype,0,1)
      if(nctype+newghostype.gt.MCTYPE) call fatal_error('ATOMS: nctype+newghostype > MCTYPE')
      numr=numeric
      do 10 i=1,nctype+newghostype
        read(iu,*) n1s(i),n2s(i),(n2p(j,i),j=1,3)
     &  ,n3s(i),(n3p(j,i),j=1,3)
     &  ,n3dzr(i),n3dx2(i),n3dxy(i),n3dxz(i),n3dyz(i)
     &  ,n4s(i),(n4p(j,i),j=1,3)
     &  ,n4fxxx(i),n4fyyy(i),n4fzzz(i),n4fxxy(i),n4fxxz(i)
     &  ,n4fyyx(i),n4fyyz(i),n4fzzx(i),n4fzzy(i),n4fxyz(i)
     &  ,nsa(i),(npa(j,i),j=1,3)
     &  ,ndzra(i),ndx2a(i),ndxya(i),ndxza(i),ndyza(i)
      call incpos(iu,itmp,1)
        if(numr.gt.0) then
          if(n2s(i).ne.0.or.n3s(i).ne.0.or.n4s(i).ne.0.or.
     &      n3p(1,i).ne.0.or.n3p(2,i).ne.0.or.n3p(3,i).ne.0.or.
     &      n4p(1,i).ne.0.or.n4p(2,i).ne.0.or.n4p(3,i).ne.0.or.
     &      nsa(i).ne.0.or.npa(1,i).ne.0.or.npa(2,i).ne.0.or.
     &      npa(3,i).ne.0.or.ndzra(i).ne.0.or.ndx2a(i).ne.0.or.
     &      ndxya(i).ne.0.or.ndxza(i).ne.0.or.ndyza(i).ne.0)
     &      call fatal_error('BASIS: n1s,n2p,n3d only for numerical basis')

          nbastyp(i)=iabs(n1s(i))
     &           +iabs(n2p(1,i))+iabs(n2p(2,i))+iabs(n2p(3,i))
     &           +iabs(n3dzr(i))+iabs(n3dx2(i))
     &           +iabs(n3dxy(i))+iabs(n3dxz(i))+iabs(n3dyz(i))
     &           +iabs(n4fxxx(i))+iabs(n4fyyy(i))+iabs(n4fzzz(i))+iabs(n4fxxy(i))+iabs(n4fxxz(i))
     &           +iabs(n4fyyx(i))+iabs(n4fyyz(i))+iabs(n4fzzx(i))+iabs(n4fzzy(i))+iabs(n4fxyz(i))

          if(nbastyp(i).gt.MRWF) call fatal_error('BASIS: nbastyp > MRWF')

          read(iu,*) (iwrwf(ib,i),ib=1,nbastyp(i))
          call incpos(iu,itmp,1)
         else
          if(n4fxxx(i).ne.0.or.n4fyyy(i).ne.0.or.n4fzzz(i).ne.0.or.
     &       n4fxxy(i).ne.0.or.n4fxxz(i).ne.0.or.n4fyyx(i).ne.0.or.
     &       n4fyyz(i).ne.0.or.n4fzzx(i).ne.0.or.n4fzzy(i).ne.0.or.
     &       n4fxyz(i).ne.0) call fatal_error('BASIS: n4f only for numerical basis')
        endif
   10 continue

      if(numr.gt.0) then

        do 1000 i=1,nctype+newghostype
          jj=0
          do 20 j=1,iabs(n1s(i))
            jj=jj+1
   20       iwlbas(jj,i)=1
          do 30 j=1,iabs(n2p(1,i))
            jj=jj+1
   30       iwlbas(jj,i)=2
          do 40 j=1,iabs(n2p(2,i))
            jj=jj+1
   40       iwlbas(jj,i)=3
          do 50 j=1,iabs(n2p(3,i))
            jj=jj+1
   50       iwlbas(jj,i)=4
          do 60 j=1,iabs(n3dzr(i))
            jj=jj+1
   60       iwlbas(jj,i)=5
          do 70 j=1,iabs(n3dx2(i))
            jj=jj+1
   70       iwlbas(jj,i)=6
          do 80 j=1,iabs(n3dxy(i))
            jj=jj+1
   80       iwlbas(jj,i)=7
          do 90 j=1,iabs(n3dxz(i))
            jj=jj+1
   90       iwlbas(jj,i)=8
          do 100 j=1,iabs(n3dyz(i))
            jj=jj+1
  100       iwlbas(jj,i)=9
          do 110 j=1,iabs(n4fxxx(i))
            jj=jj+1
  110       iwlbas(jj,i)=10
          do 120 j=1,iabs(n4fyyy(i))
            jj=jj+1
  120       iwlbas(jj,i)=11
          do 130 j=1,iabs(n4fzzz(i))
            jj=jj+1
  130       iwlbas(jj,i)=12
          do 140 j=1,iabs(n4fxxy(i))
            jj=jj+1
  140       iwlbas(jj,i)=13
          do 150 j=1,iabs(n4fxxz(i))
            jj=jj+1
  150       iwlbas(jj,i)=14
          do 160 j=1,iabs(n4fyyx(i))
            jj=jj+1
  160       iwlbas(jj,i)=15
          do 170 j=1,iabs(n4fyyz(i))
            jj=jj+1
  170       iwlbas(jj,i)=16
          do 180 j=1,iabs(n4fzzx(i))
            jj=jj+1
  180       iwlbas(jj,i)=17
          do 190 j=1,iabs(n4fzzy(i))
            jj=jj+1
  190       iwlbas(jj,i)=18
          do 200 j=1,iabs(n4fxyz(i))
            jj=jj+1
  200       iwlbas(jj,i)=19
          
 1000   continue
      endif
c     write(6,*) 'HELLO_INPUT',numr,nctype
c     do i=1,nctype
c     write(6,*) 'HELLO_INPUT',(iwlbas(j,i),j=1,nbastyp(i))
c     enddo
      ibasis_num=1
      call p2chkend(iu, 'basis')
      end

c----------------------------------------------------------------------
      subroutine read_lattice(iu)
C$INPUT lattice inp
CKEYDOC Lattice vectors of primitive and simulation cell
      implicit real*8(a-h,o-z)
      call do_read_lattice(iu)
      end
c-----------------------------------------------------------------------
      subroutine read_forces(iu)
C$INPUT forces_displace inp
CKEYDOC Displacement parameters and wave function types

      use forcepar, only: nforce
      use forcestr, only: delc
      use wfsec, only: iwftype
      use inputflags, only: iznuc,igeometry,ibasis_num,ilcao,iexponents,
     &             ideterminants,ijastrow_parameter, ioptorb_def,ilattice,
     &             ici_def,iforces,icsfs,imstates,igradients,icharge_efield,
     &             imultideterminants,ioptorb_mixvirt,imodify_zmat,izmatrix_check,
     &             ihessian_zmat 

      implicit real*8(a-h,o-z)

      include 'vmc.h'
      include 'force.h'

      call p2gti('atoms:natom',ncent,1)
      if(ncent.gt.MCENT) call fatal_error('FORCES: ncent > MCENT')

      call p2gtid('general:nforce',nforce,1,1)
      if(nforce.gt.MFORCE) call fatal_error('FORCES: nforce > MFORCE')

      do 60 i=1,nforce
        do 60 ic=1,ncent
          call incpos(iu,itmp,1)
   60     read(iu,*)  (delc(k,ic,i),k=1,3)
      call incpos(iu,itmp,1)
      read(iu,*) (iwftype(i),i=1,nforce)
      if(iwftype(1).ne.1) call fatal_error('INPUT: iwftype(1) ne 1')

      iforces=1
      call p2chkend(iu, 'forces')

      return
      end
c-----------------------------------------------------------------------
      subroutine read_csf(ncsf_read,nstates_read,fn)
C$INPUT csf i i=1 a=<input>

      use csfs, only: ccsf, ncsf, nstates
      use inputflags, only: iznuc,igeometry,ibasis_num,ilcao,iexponents,
     &             ideterminants,ijastrow_parameter, ioptorb_def,ilattice,
     &             ici_def,iforces,icsfs,imstates,igradients,icharge_efield,
     &             imultideterminants,ioptorb_mixvirt,imodify_zmat,izmatrix_check,
     &             ihessian_zmat 

      implicit real*8(a-h,o-z)

      include 'vmc.h'
      include 'force.h'
      include 'mstates.h'
      include 'optci.h'

      character fn*(*)

      call ptfile(iu,fn,'old')

      ncsf=ncsf_read
      if(ncsf.gt.MDET) 
     $ call fatal_error('CSF: too many csf')

      nstates=nstates_read
      if(nstates.gt.MSTATES) 
     $ call fatal_error('CSF: too many states')

      do i=1,nstates
        read(iu,*) (ccsf(j,i,1),j=1,ncsf)
      enddo

      icsfs=1

      if(fn.eq.'<input>') then
       call p2chkend(iu, 'csf')
      endif

      end   
c-----------------------------------------------------------------------
      subroutine read_csfmap(fn)
C$INPUT csfmap a=<input>
CKEYDOC Read mapping between csf and determinants.
      use csfs, only: ccsf, cxdet, iadet, ibdet, icxdet, ncsf, nstates

      use dets, only: cdet, ndet
      implicit real*8(a-h,o-z)


      include 'vmc.h'
      include 'force.h'
      include 'mstates.h'
      character fn*(*)
c
      call ptfile(iu,fn,'old')
c
      read(iu,*) ncsf_check,ndet_check,nmap_check
      write(6,'(''csfmap'',3i4)') ncsf_check,ndet_check,nmap_check
      if(ndet_check.ne.ndet) 
     $ call fatal_error('CSFMAP: wrong number of determinants')
      if(ncsf_check.ne.ncsf) 
     $ call fatal_error('CSFMAP: wrong number of csf')
      if(nmap_check.gt.float(MDET)*MDET) 
     $ call fatal_error('CSFMAP: too many determinants in map list')
c
      nptr=1
      do 10 i=1,ncsf
       read(iu,*) nterm
       iadet(i)=nptr
       ibdet(i)=nptr+nterm-1
       do 12 j=1,nterm
        read(iu,*) id,c
        icxdet(nptr)=id
        cxdet(nptr)=c
        nptr=nptr+1
        if(nptr.gt.MDET*MDETCSFX)
     $    call fatal_error('CSFMAP: problem with nmap')
 12    enddo
 10   enddo
      if(nmap_check.ne.nptr-1)
     $ call fatal_error('CSFMAP: problem with nmap')
      nmap=nptr


      write(6,'(''Warning: det coef overwritten with csf'')') 
      do 30 k=1,nstates
        do 15 j=1,ndet
 15       cdet(j,k,1)=0
        do 30 icsf=1,ncsf
          do 20 j=iadet(icsf),ibdet(icsf)
            jx=icxdet(j)
            cdet(jx,k,1)=cdet(jx,k,1)+ccsf(icsf,k,1)*cxdet(j)
 20       continue
 30   continue

      if(fn.eq.'<input>') then
       call p2chkend(iu, 'csfmap')
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine flaginit
c Initialize flags used to identify presence/absence of blocks in input

      use inputflags, only: iznuc,igeometry,ibasis_num,ilcao,iexponents,
     &             ideterminants,ijastrow_parameter, ioptorb_def,ilattice,
     &             ici_def,iforces,icsfs,imstates,igradients,icharge_efield,
     &             imultideterminants,ioptorb_mixvirt,imodify_zmat,izmatrix_check,
     &             ihessian_zmat 

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

      use numbas, only: numr
      use optorb_mix, only: norbopt, norbvirt
      use efield, only: iefield, iscreen, ncharges
      use inputflags, only: iznuc,igeometry,ibasis_num,ilcao,iexponents,
     &             ideterminants,ijastrow_parameter, ioptorb_def,ilattice,
     &             ici_def,iforces,icsfs,imstates,igradients,icharge_efield,
     &             imultideterminants,ioptorb_mixvirt,imodify_zmat,izmatrix_check,
     &             ihessian_zmat 

      use mstates_ctrl, only: iefficiency, iguiding, nstates_psig
      implicit real*8(a-h,o-z)


      include 'vmc.h'
      include 'force.h'
      include 'numbas.h'
      include 'optci.h'
      include 'optorb.h'
      include 'efield.h'

      common /orbval/ orb(MELEC,MORB),dorb(3,MELEC,MORB),ddorb(MELEC,MORB),ndetorb,nadorb

      call p2gtid('general:nforce',nforce,1,1)
      if(nforce.gt.MFORCE) call fatal_error('INPUT: nforce > MFORCE')
      call p2gtid('general:nwftype',nwftype,1,1)
      if(nwftype.gt.MWF) call fatal_error('INPUT: nwftype gt MWF')
      call p2gtid('general:iperiodic',iperiodic,0,1)
      call p2gtid('general:ibasis',ibasis,1,1)
      call p2gtid('optwf:ioptorb',ioptorb,0,1)
      call p2gtid('optwf:ioptci',ioptci,0,1)
      call p2gtid('mstates:iguiding',iguiding,0,1)
      call p2gtid('efield:iefield',iefield,0,1)
      call p2gtid('optgeo:iforce_analy',iforce_analy,0,0)
      call p2gtid('optgeo:iuse_zmat',iuse_zmat,0,0)

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
      if(nforce.gt.1.and.iforces.eq.0.and.igradients.eq.0) then
        write(6,'(''INPUT: block forces_displace or gradients_* missing: geometries set equal to primary'')')
        call inputforces
      endif
      if(iforce_analy.gt.0) then
        if(iuse_zmat.gt.0.and.izmatrix_check.eq.0) call fatal_error('INPUT: block connectionzmatrix missing')
        if(imodify_zmat.eq.0) call modify_zmat_define
        if(ihessian_zmat.eq.0) call hessian_zmat_define
      endif
      if(imultideterminants.eq.0) then
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
      subroutine inputcsf
c Check that the required blocks are there in the input

      use csfs, only: ncsf, nstates
      use inputflags, only: iznuc,igeometry,ibasis_num,ilcao,iexponents,
     &             ideterminants,ijastrow_parameter, ioptorb_def,ilattice,
     &             ici_def,iforces,icsfs,imstates,igradients,icharge_efield,
     &             imultideterminants,ioptorb_mixvirt,imodify_zmat,izmatrix_check,
     &             ihessian_zmat 

      use ci000, only: iciprt, nciprim, nciterm

      implicit real*8(a-h,o-z)


      include 'vmc.h'
      include 'force.h'
      include 'optci.h'
      include 'mstates.h'

      nstates=1
      ncsf=0

      call p2gtid('optwf:ioptci',ioptci,0,1)
      if(ioptci.ne.0.and.ici_def.eq.1) nciterm=nciprim
      return
      end
c----------------------------------------------------------------------
      subroutine inputzex
c Set the exponents to one when using a numerical basis
      use numbas, only: numr
      use coefs, only: nbasis
      use basis, only: zex

      implicit real*8(a-h,o-z)

      include 'vmc.h'
      include 'force.h'
      include 'numbas.h'


      call p2gtid('general:nwftype',nwftype,1,1)
      call p2gtid('general:iperiodic',iperiodic,0,1)
      if(nwftype.gt.MWF) call fatal_error('WF: nwftype gt MWF')

      if(numr.eq.0.and.iperiodic.eq.0) 
     & call fatal_error('ZEX: numr=0 and iperiodic=0 but no zex are inputed')
      do 10 iwft=1,nwftype
        do 10 i=1,nbasis
   10     zex(i,iwft)=1

      end

c----------------------------------------------------------------------
      subroutine inputdet(nwftype)
c Set the cdet to be equal
      use dets, only: cdet, ndet
      implicit real*8(a-h,o-z)


      include 'vmc.h'
      include 'force.h'
      include 'mstates.h'

       do 10 iwft=2,nwftype
         do 10 k=1,ndet
   10      cdet(k,1,iwft)=cdet(k,1,1)

      end
c----------------------------------------------------------------------
      subroutine inputlcao(nwftype)
c Set the lcao to be equal
      use coefs, only: coef, nbasis, norb
      implicit real*8(a-h,o-z)


      include 'vmc.h'
      include 'force.h'


       do 10 iwft=2,nwftype
         do 10 i=1,norb
           do 10 j=1,nbasis
   10        coef(j,i,iwft)=coef(j,i,1)

      end
c----------------------------------------------------------------------
      subroutine inputjastrow(nwftype)
c Set the jastrow to be equal

      use jaspar, only: nspin1, nspin2
      use jaspar3, only: a, b, c, scalek
      use jaspar4, only: a4, norda, nordb, nordc
      use bparm, only: nspin2b
      use contr2, only: ifock, ijas
      use contr2, only: isc
      use inputflags, only: iznuc,igeometry,ibasis_num,ilcao,iexponents,
     &             ideterminants,ijastrow_parameter, ioptorb_def,ilattice,
     &             ici_def,iforces,icsfs,imstates,igradients,icharge_efield,
     &             imultideterminants,ioptorb_mixvirt,imodify_zmat,izmatrix_check,
     &             ihessian_zmat 

      implicit real*8(a-h,o-z)

      include 'vmc.h'
      include 'force.h'

      call p2gti('jastrow:ijas',ijas,1)
      call p2gti('jastrow:isc',isc,1)
      call p2gtid('jastrow:nspin1',nspin1,1,1)
      call p2gtid('jastrow:nspin2',nspin2,1,1)
      call p2gtid('jastrow:ifock',ifock,0,1)

      call p2gti('atoms:natom',ncent,1)
      call p2gti('atoms:nctype',nctype,1)

      if(ijas.ge.4.and.ijas.le.6) then
        mparmja=2+max(0,norda-1)
        mparmjb=2+max(0,nordb-1)
        mparmjc=nterms4(nordc)
        do 90 iwft=2,nwftype
          scalek(iwft)=scalek(1)
          do 70 it=1,nctype
            do 70 iparm=1,mparmja
   70         a4(iparm,it,iwft)=a4(iparm,it,1)
        do 80 isp=nspin1,nspin2b
          do 80 iparm=1,mparmjb
   80       b(iparm,isp,iwft)=b(iparm,isp,1)
        do 90 it=1,nctype
          do 90 iparm=1,mparmjc
   90       c(iparm,it,iwft)=c(iparm,it,1)
      endif

      end
c----------------------------------------------------------------------
      subroutine inputforces
c Set all force displacements to zero
      use forcepar, only: nforce
      use wfsec, only: iwftype, nwftype
      implicit real*8(a-h,o-z)



      include 'vmc.h'
      include 'force.h'


      call p2gti('atoms:natom',ncent,1)
      if(ncent.gt.MCENT) call fatal_error('FORCES: ncent > MCENT')

      call p2gtid('general:nforce',nforce,1,1)

      call set_displace_zero(nforce)

      call p2gtid('general:nwftype',nwftype,1,1)
      if(nwftype.gt.MWF) call fatal_error('FORCES: nwftype gt MWF')
      
      if(nwftype.eq.1) then
        do 70 i=1,nforce
   70     iwftype(i)=1
       elseif(nwftype.eq.nforce) then
        do 80 i=1,nforce
   80     iwftype(i)=i
       else
        call fatal_error('FORCES: need to specify iwftype')
      endif
   
      end
c-----------------------------------------------------------------------
      subroutine read_jasderiv(iu)
C$INPUT jasderiv inp
      use atom, only: nctype
      use jaspar, only: nspin1, is
      use jaspar4, only: norda, nordb, nordc
      use jaspointer, only: npoint, npointa
      use numbas, only: numr

      use optwf_nparmj, only: nparma, nparmb, nparmc, nparmf
      use optwf_parms, only: nparmj
      use optwf_wjas, only: iwjasa, iwjasb, iwjasc
      use bparm, only: nspin2b
      use contr2, only: ijas
      use contr2, only: isc
      implicit real*8(a-h,o-z)










      include 'vmc.h'
      include 'optjas.h'
      include 'numbas.h'
      include 'force.h'




      na1=1
      na2=nctype

      read(iu,*) (nparma(ia),ia=na1,na2),
     &(nparmb(isp),isp=nspin1,nspin2b),(nparmc(it),it=1,nctype),
     &(nparmf(it),it=1,nctype)

      if(ijas.ge.4.and.ijas.le.6) then
        do 5 it=1,nctype
          if(numr.eq.0) then
c All-electron with analytic slater basis
            if((nparma(it).gt.0.and.norda.eq.0).or.(nparma(it).gt.norda+1)) then
              write(6,'(''it,norda,nparma(it)'',3i5)') it,norda,nparma(it)
              stop 'nparma too large for norda'
            endif
           else
c Pseudopotential with numerical basis: cannot vary a(1) or a(2)
            if(norda.eq.1) stop 'makes no sense to have norda=1 for numr>0'
            if((norda.eq.0.and.nparma(it).gt.0).or.(norda.gt.0.and.nparma(it).gt.norda-1)) then
              write(6,'(''it,norda,nparma(it)'',3i5)') it,norda,nparma(it)
              stop 'nparma too large for norda'
            endif
          endif
          if(isc.le.7 .and.
     &       ((nordc.le.2.and.nparmc(it).gt.0)
     &    .or.(nordc.eq.3.and.nparmc(it).gt.2)
     &    .or.(nordc.eq.4.and.nparmc(it).gt.7)
     &    .or.(nordc.eq.5.and.nparmc(it).gt.15)
     &    .or.(nordc.eq.6.and.nparmc(it).gt.27)
     &    .or.(nordc.eq.7.and.nparmc(it).gt.43))) then
            write(6,'(''it,nordc,nparmc(it)'',3i5)') it,nordc,nparmc(it)
            stop 'nparmc too large for nordc in J_een with cusp conds'
          endif
          if(isc.gt.7 .and.
     &       ((nordc.le.1.and.nparmc(it).gt.0)
     &    .or.(nordc.eq.2.and.nparmc(it).gt.2)
     &    .or.(nordc.eq.3.and.nparmc(it).gt.6)
     &    .or.(nordc.eq.4.and.nparmc(it).gt.13)
     &    .or.(nordc.eq.5.and.nparmc(it).gt.23)
     &    .or.(nordc.eq.6.and.nparmc(it).gt.37)
     &    .or.(nordc.eq.7.and.nparmc(it).gt.55))) then
            write(6,'(''it,nordc,nparmc(it)'',3i5)') it,nordc,nparmc(it)
            stop 'nparmc too large for nordc without cusp conds'
          endif
    5   continue
c For the b coefs. we assume that b(1) is fixed by the cusp-cond.
        do 6 isp=1,nspin1,nspin2b
            if(nparmb(isp).gt.nordb) then
              write(6,'(''isp,nordb,nparmb(isp)'',3i5)') isp,nordb,nparmb(isp)
              stop 'nparmb too large for nordb'
            endif
    6   continue
      endif

c compute nparmj
      nparmj=0
      npointa(1)=0
      do 30 ia=na1,na2
        if(ia.gt.1) npointa(ia)=npointa(ia-1)+nparma(ia-1)
   30   nparmj=nparmj+nparma(ia)
      do 35 isp=nspin1,nspin2b
   35   nparmj=nparmj+nparmb(isp)
      npoint(1)=nparmj
      do 45 it=1,nctype
        if(it.gt.1) npoint(it)=npoint(it-1)+nparmc(it-1)
   45   nparmj=nparmj+nparmc(it)+nparmf(it)

      if(nparmj.gt.MPARMJ) call fatal_error('JASDERIV: MPARMJ too small') 

      do 60 it=1,nctype
   60   read(iu,*) (iwjasa(iparm,it),iparm=1,nparma(it))
      do 65 isp=nspin1,nspin2b
   65   read(iu,*) (iwjasb(iparm,isp),iparm=1,nparmb(isp))
      do 70 it=1,nctype
   70   read(iu,*) (iwjasc(iparm,it),iparm=1,nparmc(it))

c     ifitparms=1
      call p2chkend(iu, 'jasderiv')
      end
c-----------------------------------------------------------------------
      subroutine read_sym(nsym,mo,fn)
C$INPUT sym_labels i i a=<input>
CKEYDOC Read symmetry information
      use coefs, only: norb
      use optorb, only: irrep
      implicit real*8(a-h,o-z)


      include 'vmc.h'
      include 'force.h'


      character fn*(*)
      character atmp*80

      call ptfile(iu,fn,'old')
      nirrep=nsym
      if(norb.ne.0.and.norb.ne.mo) then
       write(6,'(2i5)') norb,mo
       call fatal_error('READSYM: wrong number of orbitals') 
      else
       norb=mo
      endif
c Ignore irrep text labels
      read(iu,'(a80)') atmp
      read(iu,*) (irrep(io),io=1,norb)

      if(fn.eq.'<input>') then
       call p2chkend(iu, 'sym_labels')
      endif
      end
c-----------------------------------------------------------------------
      subroutine read_optorb_mixvirt(moopt,movirt,fn)
C$INPUT optorb_mixvirt i i a=<input>
CKEYDOC Read which virtual orbitals are mixed with the occupied ones

      use optorb_mix, only: iwmix_virt, norbopt, norbvirt
      use coefs, only: norb
      use inputflags, only: iznuc,igeometry,ibasis_num,ilcao,iexponents,
     &             ideterminants,ijastrow_parameter, ioptorb_def,ilattice,
     &             ici_def,iforces,icsfs,imstates,igradients,icharge_efield,
     &             imultideterminants,ioptorb_mixvirt,imodify_zmat,izmatrix_check,
     &             ihessian_zmat 

      implicit real*8(a-h,o-z)

      include 'vmc.h'
      include 'force.h'

      character fn*(*)
      character atmp*80

      norbopt=moopt
      norbvirt=movirt
      call ptfile(iu,fn,'old')
      if(norb.ne.0.and.norbopt.gt.norb) then
       write(6,'(3i5)') norb,moopt,movirt
       call fatal_error('READMIXVIRT: wrong number of orbitals')
      endif

      do 50 io=1,norbopt
  50    read(iu,*) (iwmix_virt(io,jo),jo=1,norbvirt)

      ioptorb_mixvirt=1

      if(fn.eq.'<input>') then
       call p2chkend(iu, 'optorb_mixvirt')
      endif
      end
c-----------------------------------------------------------------------
      subroutine read_energies(mo,fn)
C$INPUT energies i a=<input>
C$INPUT eigenvalues i a=<input>
CKEYDOC Read orbital energies 
      use coefs, only: norb
      use optorb, only: orb_energy
      implicit real*8(a-h,o-z)


      include 'vmc.h'
      include 'force.h'


      character fn*(*)

      call ptfile(iu,fn,'old')
      if(norb.ne.0.and.norb.ne.mo) then
        write(6,'(2i5)') norb,mo
        call fatal_error('READEIG: wrong number of orbitals') 
      endif
      read(iu,*) (orb_energy(io),io=1,norb)

      if(fn.eq.'<input>') then
       call p2chkend(iu, 'energies')
      endif
      end
c----------------------------------------------------------------------
      subroutine read_dmatrix(no,ns,fn)
C$INPUT dmatrix i i a=<input> 
CKEYDOC Read diagonal density matrix information.
      use sa_weights, only: iweight, nweight, weights
      use coefs, only: norb
      use optorb, only: dmat_diag
      implicit real*8(a-h,o-z)




      include 'vmc.h'
      include 'force.h'
      include 'mstates.h'

      character fn*(*)

      dimension dmat(MORB),iwdmat(MSTATES)

      call p2gtid('general:ipr',ipr,-1,1)
      call ptfile(iu,fn,'old')

      ndetorb=no
      if(ndetorb.gt.norb) call fatal('READ_DMATRIX: wrong number of orbitals')

      call get_weights('weights:',weights,iweight,nweight)
      if(ns.ne.nweight) call fatal('READ_DMATRIX: wrong number of dmatrices')

      read(iu,*) (iwdmat(i),i=1,nweight)
      do 10 iw=1,nweight
  10    if(iwdmat(iw).ne.iweight(iw)) call fatal('READ_DMATRIX: iwdmat')

      do 20 i=1,norb
  20    dmat_diag(i)=0.d0

      do 25 iw=1,nweight
        read(iu,*) (dmat(j),j=1,ndetorb)
        do 25 j=1,ndetorb
  25      dmat_diag(j)=dmat_diag(j)+weights(iw)*dmat(j)
       do 30 i=1,ndetorb 
  30     if(dabs(dmat_diag(i)-1.d0).lt.1.d-6) dmat_diag(i)=1.d0

      if(ipr.gt.2) then 
       write(6,'(''diagonal elements of the density matrix'')')
       write(6,'(100f10.6)') (dmat_diag(i),i=1,ndetorb)
      endif

      if(fn.eq.'<input>') then
       call p2chkend(iu, 'dmatrix')
      endif

      return
      end
c----------------------------------------------------------------------
      subroutine get_weights(field,weights,iweight,nweight)
      use csfs, only: nstates

      implicit real*8(a-h,o-z)

      include 'vmc.h'
      include 'force.h'
      include 'mstates.h'
c weights for state averaging


      dimension weights(MSTATES),iweight(MSTATES)

      character vname*(32)
      character field*(32)
      wsum=0.d0
      nweight=0

      write(6,*) field,field(1:index(field,' '))
      do 10 i=1,nstates
        wdef=0.d0
        call append_number(field(1:index(field,' ')-1),i,vname,nv,0)
        call p2gtfd(vname(1:nv),w,wdef,0)
CVARDOC Input of weights for individual states.
        w=dabs(w)
        if(w.gt.1d-6) then
          nweight=nweight+1
          iweight(nweight)=i
          weights(nweight)=w
          wsum=wsum+w
        endif
   10 continue

      do 20 i=1,nweight
   20  weights(i)=weights(i)/wsum

      if(nweight.eq.0) then
        nweight=1
        iweight(1)=1
        weights(1)=1.d0
      endif

c TEMPORARY
      if(nweight.ne.nstates) 
     & call fatal_error('GET_WEIGHTS: problems with nweight')

      end
c----------------------------------------------------------------------
      subroutine read_cavity_spheres(iu,nspheres)
C$INPUT cavity_spheres inp i 
CKEYDOC Read centers of cavity spheres and radii
      use pcm_parms, only: ch, eps_solv, iscov, nch, nchs, nchs1, nchs2
      use pcm_parms, only: nchv, ncopcm, nesph, nscv, nvopcm, re, re2
      use pcm_parms, only: retk, surk, xe, xpol, ye, ze
      implicit real*8(a-h,o-z)


      include 'vmc.h'
      include 'pcm.h'

      nesph=nspheres 
      do i=1,nesph
        call incpos(iu,itmp,1)
        read(iu,*) xe(i),ye(i),ze(i),re(i)
        re2(i)=re(i)**2.0d0
      enddo
      call p2chkend(iu, 'cavity_spheres')

      return
      end
c----------------------------------------------------------------------
      subroutine read_gradnts_cart(iu)
C$INPUT gradients_cartesian inp
CKEYDOC Read for which x,y,z cartesian coordiantes of 
CKEYDOC atoms energy gradients are to be calculated for.

c     Written by Omar Valsson

      use forcepar, only: nforce
      use forcestr, only: delc
      use grdntsmv, only: igrdaidx, igrdcidx, igrdmv
      use grdntspar, only: delgrdxyz, igrdtype, ngradnts
      use wfsec, only: iwftype
      use inputflags, only: iznuc,igeometry,ibasis_num,ilcao,iexponents,
     &             ideterminants,ijastrow_parameter, ioptorb_def,ilattice,
     &             ici_def,iforces,icsfs,imstates,igradients,icharge_efield,
     &             imultideterminants,ioptorb_mixvirt,imodify_zmat,izmatrix_check,
     &             ihessian_zmat 

      implicit real*8(a-h,o-z)

      include 'vmc.h'
      include 'force.h'


      call p2gti('atoms:natom',ncent,1)
      if(ncent.gt.MCENT) call fatal_error('GRADIENTS_CARTESIAN: ncent > MCENT')


      call p2gtfd('gradients:delgrdxyz',delgrdxyz,0.001d0,1)

      call p2gtid('gradients:igrdtype',igrdtype,1,1)
      if(igrdtype.ne.1) call fatal_error('GRADIENTS_CARTESIAN: igrdtype /= 1')

      call p2gtid('general:nforce',nforce,1,1)      
      call p2gtid('gradients:ngradnts',ngradnts,0,1)      
      if( (2*ngradnts+1).ne.nforce) call 
     &  fatal_error('GRADIENTS_CARTESIAN: (2*ngradnts+1)  /=  nforce')


      do 60 i=1,nforce
        iwftype(i)=1
        do 60 ic=1,ncent
          do 60 k=1,3
            igrdmv(k,ic)=0
   60       delc(k,ic,i)=0.0d0

      ia=2
      do 70 ic=1,ncent
        call incpos(iu,itmp,1)
        read(iu,*)  (igrdmv(k,ic),k=1,3)
          do 70 k=1,3            
            if(igrdmv(k,ic).lt.0 .or. igrdmv(k,ic).gt.1) 
     &        call fatal_error('GRADIENTS_CARTESIAN: igrdmv \= 0,1')
            if(igrdmv(k,ic).eq.1) then
              igrdaidx(ia/2)=ic
              igrdcidx(ia/2)=k
              delc(k,ic,ia)=delgrdxyz
              delc(k,ic,ia+1)=-delgrdxyz
              ia=ia+2
            endif      
   70     continue

      igradients=1

      call p2chkend(iu, 'gradients_cartesian')

      return
      end
c-----------------------------------------------------------------------
      subroutine read_gradnts_zmat(iu)
C$INPUT gradients_zmatrix inp
CKEYDOC Read for which Z matrix (internal) coordiantes of 
CKEYDOC atoms energy gradients are to be calculated for.

c      Written by Omar Valsson.

      use forcepar, only: nforce
      use forcestr, only: delc
      use grdntsmv, only: igrdaidx, igrdcidx, igrdmv
      use grdntspar, only: delgrdba, delgrdbl, delgrdda, igrdtype, ngradnts
      use zmatrix, only: izmatrix
      use wfsec, only: iwftype
      use inputflags, only: iznuc,igeometry,ibasis_num,ilcao,iexponents,
     &             ideterminants,ijastrow_parameter, ioptorb_def,ilattice,
     &             ici_def,iforces,icsfs,imstates,igradients,icharge_efield,
     &             imultideterminants,ioptorb_mixvirt,imodify_zmat,izmatrix_check,
     &             ihessian_zmat 

      implicit real*8(a-h,o-z)

      include 'vmc.h'
      include 'force.h'

      call p2gti('atoms:natom',ncent,1)
      if(ncent.gt.MCENT) call fatal_error('GRADIENTS_ZMATRIX: ncent > MCENT')

      call p2gtfd('gradients:delgrdbl',delgrdbl,0.001d0,1)
      call p2gtfd('gradients:delgrdba',delgrdba,0.01d0,1)
      call p2gtfd('gradients:delgrdda',delgrdda,0.01d0,1)

      call p2gtid('gradients:igrdtype',igrdtype,2,1)
      if(igrdtype.ne.2) call fatal_error('GRADIENTS_ZMATRIX: igrdtype /= 2')

      if(izmatrix.ne.1) call fatal_error('GRADIENTS_ZMATRIX: No Z matrix connection matrix')

      call p2gtid('general:nforce',nforce,1,1)      
      call p2gtid('gradients:ngradnts',ngradnts,0,1)      
      if( (2*ngradnts+1).ne.nforce) call 
     &  fatal_error('GRADIENTS_ZMATRIX: (2*ngradnts+1)  /=  nforce')

      do 60 i=1,nforce
        iwftype(i)=1
        do 60 ic=1,ncent
          do 60 k=1,3
            igrdmv(k,ic)=0
   60       delc(k,ic,i)=0.0d0

      ia=2

      do 70 ic=1,ncent
        call incpos(iu,itmp,1)
        read(iu,*)  (igrdmv(k,ic),k=1,3)
          do 70 k=1,3            
            if(igrdmv(k,ic).lt.0 .or. igrdmv(k,ic).gt.1) 
     &        call fatal_error('GRADIENTS_ZMATRIX: igrdmv \= 0,1')
            if(igrdmv(k,ic).eq.1) then
              igrdaidx(ia/2)=ic
              igrdcidx(ia/2)=k
              call grdzmat_displ(k,ic,ia  ,+1.0d0)
              call grdzmat_displ(k,ic,ia+1,-1.0d0)
              ia=ia+2
            endif      
   70     continue

 
      igradients=1

      call p2chkend(iu, 'gradients_zmatrix')

      return
      end
c-----------------------------------------------------------------------

      subroutine read_modify_zmat(iu)
C$INPUT modify_zmatrix inp
CKEYDOC Read for which Z matrix (internal) coordiantes of 
CKEYDOC atoms energy gradients are to be calculated for.

      use grdntsmv, only: igrdmv
      use inputflags, only: iznuc,igeometry,ibasis_num,ilcao,iexponents,
     &             ideterminants,ijastrow_parameter, ioptorb_def,ilattice,
     &             ici_def,iforces,icsfs,imstates,igradients,icharge_efield,
     &             imultideterminants,ioptorb_mixvirt,imodify_zmat,izmatrix_check,
     &             ihessian_zmat 

      implicit real*8(a-h,o-z)

      include 'vmc.h'
      include 'force.h'


      call p2gti('atoms:natom',ncent,1)
      if(ncent.gt.MCENT) call fatal_error('MODIFY_ZMATRIX: ncent > MCENT')

      do 70 ic=1,ncent
        call incpos(iu,itmp,1)
        read(iu,*)  (igrdmv(k,ic),k=1,3)
          do 70 k=1,3            
            if(igrdmv(k,ic).lt.0 .or. igrdmv(k,ic).gt.1) 
     &        call fatal_error('MODIFY_ZMATRIX: igrdmv \= 0,1')
   70     continue

      imodify_zmat=1
      call p2chkend(iu, 'modify_zmatrix')

      return
      end
c-----------------------------------------------------------------------

      subroutine read_hessian_zmat(iu)
C$INPUT hessian_zmatrix inp
CKEYDOC Read for which Z matrix (internal) coordiantes of 
CKEYDOC atoms energy gradients are to be calculated for.

      use grdnthes, only: hessian_zmat
      use inputflags, only: iznuc,igeometry,ibasis_num,ilcao,iexponents,
     &             ideterminants,ijastrow_parameter, ioptorb_def,ilattice,
     &             ici_def,iforces,icsfs,imstates,igradients,icharge_efield,
     &             imultideterminants,ioptorb_mixvirt,imodify_zmat,izmatrix_check,
     &             ihessian_zmat 

      implicit real*8(a-h,o-z)

      include 'vmc.h'
      include 'force.h'



      call p2gti('atoms:natom',ncent,1)
      if(ncent.gt.MCENT) call fatal_error('HESSIAN_ZMATRIX: ncent > MCENT')

      do 70 ic=1,ncent
        call incpos(iu,itmp,1)
        read(iu,*)  (hessian_zmat(k,ic),k=1,3)
          do 70 k=1,3            
            if(hessian_zmat(k,ic).le.0 ) 
     &        call fatal_error('HESSIAN_ZMATRIX: hess <=  0')
   70     continue

      ihessian_zmat=1
      call p2chkend(iu, 'hessian_zmatrix')

      return
      end
c-----------------------------------------------------------------------

      subroutine read_zmat_conn(iu)
C$INPUT zmatrix_connectionmatrix inp
CKEYDOC Read the atom connection matrix for the Z matrix.
CKEYDOC It is need when calculating forces in Z matrix
CKEYDOC coordinates.

c      Written by Omar Valsson

      use atom, only: cent, ncent
      use zmatrix, only: czcart, czint, czcart_ref, izcmat, izmatrix
      use inputflags, only: iznuc,igeometry,ibasis_num,ilcao,iexponents,
     &             ideterminants,ijastrow_parameter, ioptorb_def,ilattice,
     &             ici_def,iforces,icsfs,imstates,igradients,icharge_efield,
     &             imultideterminants,ioptorb_mixvirt,imodify_zmat,izmatrix_check,
     &             ihessian_zmat 

      implicit real*8(a-h,o-z)

      include 'vmc.h'
      include 'force.h'

      
      do 10 ic=1,3
        do 10 k=1,3
   10     czcart_ref(k,ic)=cent(k,ic)

      do 20 ic=1,ncent
        do 20 k=1,3
          izcmat(k,ic)=0
          czint(k,ic)=0.0d0
   20     czcart(k,ic)=cent(k,ic)


      do 30 ic=1,ncent
        call incpos(iu,itmp,1)
        read(iu,*)  (izcmat(k,ic),k=1,3)
        do 30 k=1,3
   30     if(izcmat(k,ic).ge.ic) call fatal_error('ZMATRIX: Error in connection matrix')
 
      call cart2zmat(MCENT,czcart,izcmat,czint)
      call zmat2cart_rc(MCENT,izcmat,czint,czcart,czcart_ref)

      izmatrix=1
      izmatrix_check=1

      call p2chkend(iu, 'zmatrix_connectionmatrix')

      return
      end
c-----------------------------------------------------------------------
      subroutine read_efield(ncharges_tmp,iscreen_tmp,filename)
C$INPUT efield i i a=<input>

      use efield_blk, only: ascreen, bscreen, qcharge, xcharge, ycharge, zcharge
      use efield, only: iefield, iscreen, ncharges
      use inputflags, only: iznuc,igeometry,ibasis_num,ilcao,iexponents,
     &             ideterminants,ijastrow_parameter, ioptorb_def,ilattice,
     &             ici_def,iforces,icsfs,imstates,igradients,icharge_efield,
     &             imultideterminants,ioptorb_mixvirt,imodify_zmat,izmatrix_check,
     &             ihessian_zmat 

      implicit real*8(a-h,o-z)

      include 'vmc.h'
      include 'force.h'
      include 'efield.h'

      character filename*(*)

      call file(iu,filename,'old',1,0)
      ncharges=ncharges_tmp
      iscreen=iscreen_tmp
      write(6,*) 'reading in',ncharges,' charges!'

      if(ncharges.gt.MCHARGES) call fatal_error('EFIELD: ncharges > MCHARGES')

      do 20 i=1,ncharges
        call incpos(iu,itmp,1)
   20   read(iu,*) xcharge(i),ycharge(i),zcharge(i),qcharge(i),ascreen(i),bscreen(i)
      icharge_efield=icharge_efield+1
      write(6,*) 'icharge_efield=',icharge_efield

      if(filename.eq.'<input>') then
       call p2chkend(iu, 'efield')
      endif
      end
c-----------------------------------------------------------------------
      subroutine set_displace_zero(nforce_tmp)
      use forcestr, only: delc
      use pcm_force, only: sch_s
      use pcm_cntrl, only: icall, ichpol, ipcm, ipcmprt, isurf
      use pcm_parms, only: ch, eps_solv, iscov, nch, nchs, nchs1, nchs2
      use pcm_parms, only: nchv, ncopcm, nesph, nscv, nvopcm, re, re2
      use pcm_parms, only: retk, surk, xe, xpol, ye, ze
      implicit real*8(a-h,o-z)





      include 'vmc.h'
      include 'force.h'
      include 'pcm.h'


      call p2gti('atoms:natom',ncent,1)
      if(ncent.gt.MCENT) call fatal_error('FORCES: ncent > MCENT')

      do 60 i=1,nforce_tmp
        do 60 ic=1,ncent
          do 60 k=1,3
   60       delc(k,ic,i)=0.d0

      if(ipcm.eq.3) then
        do 65 i=1,nforce_tmp
          do 65 j=1,nchs
   65       sch_s(j,i)=ch(j)
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine modify_zmat_define

      use grdntsmv, only: igrdmv
      use inputflags, only: iznuc,igeometry,ibasis_num,ilcao,iexponents,
     &             ideterminants,ijastrow_parameter, ioptorb_def,ilattice,
     &             ici_def,iforces,icsfs,imstates,igradients,icharge_efield,
     &             imultideterminants,ioptorb_mixvirt,imodify_zmat,izmatrix_check,
     &             ihessian_zmat 

      implicit real*8(a-h,o-z)

      include 'vmc.h'
      include 'force.h'

      call p2gti('atoms:natom',ncent,1)
      if(ncent.gt.MCENT) call fatal_error('MODIFY_ZMATRIX: ncent > MCENT')

      do 10 ic=1,ncent
        do 10 k=1,3
 10       igrdmv(k,ic)=1

      return
      end
c-----------------------------------------------------------------------
      subroutine hessian_zmat_define

      use grdnthes, only: hessian_zmat
      use inputflags, only: iznuc,igeometry,ibasis_num,ilcao,iexponents,
     &             ideterminants,ijastrow_parameter, ioptorb_def,ilattice,
     &             ici_def,iforces,icsfs,imstates,igradients,icharge_efield,
     &             imultideterminants,ioptorb_mixvirt,imodify_zmat,izmatrix_check,
     &             ihessian_zmat 

      implicit real*8(a-h,o-z)

      include 'vmc.h'
      include 'force.h'

      call p2gti('atoms:natom',ncent,1)
      if(ncent.gt.MCENT) call fatal_error('HESSIAN_ZMATRIX: ncent > MCENT')

      do 10 ic=1,ncent
        do 10 k=1,3
 10       hessian_zmat(k,ic)=1.d0

      return
      end
