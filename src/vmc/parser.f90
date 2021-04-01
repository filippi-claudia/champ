!
!     Shows FDF capabilities..
!
subroutine parser
  use fdf 
  use prec

! CHAMP modules
  use contr3,         only: mode
  use allocation_mod, only: allocate_vmc, allocate_dmc
  use periodic_table

! in the replacement of preprocess input
  use elec,           only: ndn, nup
  use const,          only: nelec
  use ghostatom,      only: newghostype, nghostcent
  use atom,           only: nctype, ncent
  use contrl,         only: nstep, nblk, nblk_max
  use contr3,         only: mode
  use wfsec,          only: nwftype
  use forcepar,       only: nforce
  use force_mod,      only: MFORCE
  use method_opt,     only: method

! variables from process input
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

! Note the additions: Ravindra
  use control_vmc, only: vmc_idump,  vmc_irstar, vmc_isite, vmc_nconf, vmc_nblk, vmc_nblk_max
  use control_vmc, only: vmc_nblkeq, vmc_nconf_new, vmc_nstep, vmc_icharged_atom, vmc_nblk_ci  
! Note the additions: Ravindra
  use control_dmc, only: dmc_idump, dmc_irstar, dmc_isite, dmc_nconf, dmc_nblk, dmc_nblk_max
  use control_dmc, only: dmc_nblkeq, dmc_nconf_new, dmc_nstep, dmc_icharged_atom, dmc_nblk_ci  

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
  use dmc_mod, only: MWALK  


! Note the following modules are new additions


!  
  implicit none
!--------------------------------------------------------------- Local Variables
  integer, parameter         :: maxa = 100
  logical                    :: doit, debug

  character(len=72)          :: fname, filename, fmt, title

  type(block_fdf)            :: bfdf
  type(parsed_line), pointer :: pline

  character(len=20)          :: real_format    = '(A, T20, F14.8)'
  character(len=20)          :: int_format     = '(A, T20, I8)'
  character(len=80)          :: string_format  = '(A, T40, A)'  
  character(len=132)         :: file_basis, file_pseudo, path_pool
! for determinants sections
  ! integer                    :: nelectrons, nexcitation, iostat
  ! integer, allocatable       :: det_alpha(:), det_beta(:)
  ! real(selected_real_kind(6,15)), allocatable :: det_coeff(:)
  ! character(len=20)          :: temp1, temp2, temp3, temp4, temp5
!------------------------------------------------------------------------- BEGIN
! debug purpose only
  character(len=72)          :: optwf, blocking_vmc, blocking_dmc


! from process input subroutine

  character(len=20)          :: fmt
  character(len=32)          :: keyname
  character(len=10)          :: eunit
  character(len=16)          :: cseed
  integer                    :: irn(4), cent_tmp(3)
  integer, allocatable       :: anorm(:) ! dimensions = nbasis



! Initialize # get the filenames from the commandline arguments
  call fdf_init('test-champ.inp', 'test-champ.out')

! Headers in general
  title       = fdf_get('title', 'Untitled')
  pool_dir    = fdf_get('pool', '.')
  pp_id       = fdf_get('pseudopot', 'none')  
  bas_id      = fdf_get('basis', 'none')  

! some units  
  eunit       = fdf_get('unit', 'Hartrees')
  hb          = fdf_get('mass', 0.5d0)

! %module electrons
  nelec       = fdf_get('nelec', 1)
  nup         = fdf_get('nup', 1)  
  ndn         = nelec-nup

! %module atoms
  nctype      = fdf_get('nctype', 1)    
  ncent       = fdf_get('ncent', 1)     
  newghostype = fdf_get('newghostype', 0)       
  nghostcent  = fdf_get('nghostcent', 0)       

! %module general
  nforce      = fdf_get('nforce', 1)    
  MFORCE      = nforce 
  nwftype     = fdf_get('nwftype', 1)      

! %module jastrow
  ijas        = fdf_get('ijas', 1)      
  isc         = fdf_get('isc', 1)
  nspin1      = fdf_get('nspin1', 1)
  nspin2      = fdf_get('nspin2', 1)  
  ifock       = fdf_get('ifock', 0)

  ! general:
  iperiodic   = fdf_get('iperiodic', 0)  
  ibasis      = fdf_get('ibasis', 1)    
  

  ! mstates:
  iguiding    = fdf_get('iguiding', 0)      
  ! efield:  
  iefield     = fdf_get('iefield', 0)       

  ! optgeo:  
  iforce_analy= fdf_get('iforce_analy', 0)  
  iuse_zmat   = fdf_get('iuse_zmat', 0)     
  ! optwf:    
  nadorb      = fdf_get('nextorb', next_max)


  ! module vmc
  imetro      = fdf_get('imetro', 6)     
  node_cutoff = fdf_get('node_cutoff', 0)       
  eps_node_cutoff = fdf_get('eps_node_cutoff', 1.0d-7)         
  
  delta       = fdf_get('delta', 1)       
  deltar      = fdf_get('deltar', 1)         
  deltat      = fdf_get('deltat', 1)         
  fbias       = fdf_get('fbias', 1.0d0)           


  !module dmc
  idmc        = fdf_get('idmc', 2)       
  ipq         = fdf_get('ipq', 1)       
  itau_eff    = fdf_get('itau_eff', 1)       
  iacc_rej    = fdf_get('iacc_rej', 1)       
  icross      = fdf_get('icross', 1)       
  icuspg      = fdf_get('icuspg', 0)
  idiv_v      = fdf_get('idiv_v', 0)
  icut_br     = fdf_get('icut_br', 0)  
  icut_e      = fdf_get('icut_e', 0)    

! attention conflict dmc and vmc variable names match
  dmc_node_cutoff = fdf_get('dmc_node_cutoff', 0)       
  dmc_eps_node_cutoff = fdf_get('dmc_eps_node_cutoff', 1.0d-7)           
!  call p2gtid('dmc:node_cutoff',node_cutoff,0,1)
!  call p2gtfd('dmc:enode_cutoff',eps_node_cutoff,1.d-7,1)

  nfprod      = fdf_get('nfprod', 1)    
  tau         = fdf_get('tau', 1)      
  rttau=dsqrt(tau)

  etrial      = fdf_get('etrial', 1)    
  nfprod      = fdf_get('nfprod', 200)      
  itausec     = fdf_get('itausec', 1)        

  icasula     = fdf_get('icasula', 0)      
  nloc        = fdf_get('nloc', 0)        

! vmc
  vmc_nstep     = fdf_get('vmc_nstep', 1)    
  vmc_nblk      = fdf_get('vmc_nblk', 1)      
  vmc_nblkeq    = fdf_get('vmc_nblkeq', 2)        
  vmc_nblk_max  = fdf_get('vmc_nblk_max', 1)      
  vmc_nconf     = fdf_get('vmc_nconf', 1)        
  vmc_nconf_new = fdf_get('vmc_nconf_new', 1)      
  vmc_idump     = fdf_get('vmc_idump', 1)        
  vmc_irstar    = fdf_get('vmc_irstar', 0)        
  vmc_isite     = fdf_get('vmc_isite', 1)          
  vmc_icharged_atom     = fdf_get('vmc_icharged_atom', 0)            


!dmc
  dmc_nstep     = fdf_get('dmc_nstep', 1)    
  dmc_nblk      = fdf_get('dmc_nblk', 1)      
  dmc_nblkeq    = fdf_get('dmc_nblkeq', 2)          
  dmc_nblk_max  = fdf_get('dmc_nblk_max', 1)      
  dmc_nconf     = fdf_get('dmc_nconf', 1)        
  dmc_nconf_new = fdf_get('dmc_nconf_new', 1)      
  dmc_idump     = fdf_get('dmc_idump', 1)        
  dmc_irstar    = fdf_get('dmc_irstar', 0)        
  dmc_isite     = fdf_get('dmc_isite', 1)          
  dmc_icharged_atom     = fdf_get('dmc_icharged_atom', 0)            

!forces
  izvzb         = fdf_get('izvzb', 0)          
  alfgeo        = fdf_get('alfgeo', 1.0d0)            

!optimization flags vmc/dmc
  ioptwf        = fdf_get('ioptwf', 0)    
  method        = fdf_get('method', 'linear')            
  idl_flag      = fdf_get('idl_flag', 0)      
  ilbfgs_flag   = fdf_get('ilbfgs_flag', 0)        
  ilbfgs_m      = fdf_get('ilbfgs_m', 5)        
  i_sr_rescale  = fdf_get('i_sr_rescale', 0)        

  ibeta         = fdf_get('ibeta', -1)    
  ratio         = fdf_get('ratio', ratio_j)    
  iapprox       = fdf_get('iapprox', 0)    
  ncore         = fdf_get('ncore', 0)      
  iuse_orbeigv  = fdf_get('iuse_orbeigv', 0)    
  
  ioptjas       = fdf_get('ioptjas', 0)    
  ioptorb       = fdf_get('ioptorb', 0)    
  ioptci        = fdf_get('ioptci', 0)    
  no_active     = fdf_get('no_active', 0)        

! attention check the keyword nblk_max. it appears in opt/vmc/dmc
  nblk_max      = fdf_get('nblk_max', nblk) !or vmc_nblk      

  energy_tol    = fdf_get('energy_tol', 1.d-3)
  dparm_norm_min = fdf_get('dparm_norm_min', 1.0d0)
  add_diag(1)   = fdf_get('add_diag',1.d-6)
  nopt_iter     = fdf_get('nopt_iter',6)
  micro_iter_sr = fdf_get('micro_iter_sr', 1)

  ifunc_omega   = fdf_get('func_omega', 0)
  omega0        = fdf_get('omega', 0.d0)
  n_omegaf      = fdf_get('n_omegaf', nopt_iter)
  n_omegat      = fdf_get('n_omegat', 0)

  nvec          = fdf_get('lin_nvec', 5)
  nvecx         = fdf_get('lin_nvecx', MVEC)
  alin_adiag    = fdf_get('lin_adiag', 0.01)
  alin_eps      = fdf_get('lin_eps', 0.001)
  lin_jdav      = fdf_get('lin_jdav',0)
  multiple_adiag = fdf_get('multiple_adiag',0)

  sr_tau        = fdf_get('sr_tau', 0.02)
  sr_adiag      = fdf_get('sr_adiag', 0.01)
  sr_eps        = fdf_get('sr_eps', 0.001)

  iroot_geo     = fdf_get('iroot_geo', 0)
  nblk_ci       = fdf_get('nblk_ci', nblk)  !attention
  ilastvmc      = fdf_get('ilastvmc',1)
  dl_mom        = fdf_get('dl_mom', 0.0)
  dl_alg        = fdf_get('dl_alg','nag')


  ngrad_jas_blocks = fdf_get('ngrad_jas_blocks',0)
  isample_cmat  = fdf_get('isample_cmat', 1)
  isavebl       = fdf_get('save_blocks', 0)
  nefp_blocks   = fdf_get('force_blocks',1)
  iorbsample    = fdf_get('iorbsample',1)
  !ci
  iciprt        = fdf_get('ci:iciprt',0)  
  !mstates
  iguiding      = fdf_get('iguiding',0)
  iefficiency   = fdf_get('iefficiency',0)
  !pcm
  iefield       = fdf_get('efield:iefield',0)
  ipcm          = fdf_get('ipcm',0)
  ipcmprt       = fdf_get('ipcmprt',0)

  pcmfile_cavity = fdf_get('file_cavity','pcm000.dat')
  pcmfile_chs   = fdf_get('file_chs','chsurf_old')
  pcmfile_chv   = fdf_get('file_chv','chvol_old')
  nscv          = fdf_get('nblk_chv',nblk)
  iscov         = fdf_get('nstep_chv',nstep2)
  eps_solv      = fdf_get('eps_solv',1)
  fcol          = fdf_get('fcol',1.d0)
  rcolv         = fdf_get('rcolv',0.04d0)
  npmax         = fdf_get('npmax',1)
  ipcm_3dgrid   = fdf_get('ipcm_3dgrid',ipcm_3dgrid,0)
  ipcm_nstep3d(1) = fdf_get('nx_pcm',IUNDEFINED)
  ipcm_nstep3d(2) = fdf_get('ny_pcm',IUNDEFINED)
  ipcm_nstep3d(3) = fdf_get('nz_pcm',IUNDEFINED)
  pcm_step3d(1) = fdf_get('dx_pcm',UNDEFINED)
  pcm_step3d(2) = fdf_get('dy_pcm',UNDEFINED)
  pcm_step3d(3) = fdf_get('dz_pcm',UNDEFINED)
  pcm_origin(1) = fdf_get('x0_pcm',UNDEFINED)
  pcm_origin(2) = fdf_get('y0_pcm',UNDEFINED)
  pcm_origin(3) = fdf_get('z0_pcm',UNDEFINED)
  pcm_endpt(1)  = fdf_get('xn_pcm',UNDEFINED)
  pcm_endpt(2)  = fdf_get('yn_pcm',UNDEFINED)
  pcm_endpt(3)  = fdf_get('zn_pcm',UNDEFINED)
  PCM_SHIFT     = fdf_get('shift',4.d0)


  immpol        = fdf_get('immpol',0)
  immpolprt     = fdf_get('immpolprt',0)
  mmpolfile_sites = fdf_get('file_sites','mmpol000.dat')
  mmpolfile_chmm = fdf_get('file_mmdipo','mmdipo_old')
  a_cutoff      = fdf_get('a_cutoff',2.5874d0)
  rcolm         = fdf_get('rcolm',0.04d0)
  iprop         = fdf_get('sample',0)
  ipropprt      = fdf_get('print',0)
  nloc          = fdf_get('nloc',0)
  nquad         = fdf_get('nquad',6)    

  ! file_basis = fdf_load_filename('basis', 'default.bas')
  ! write(6,fmt=string_format) 'filename basis :: ', trim(file_basis)

! module dependent processing . These will be replaced by inliners

! %module blocking_vmc  
  if (fdf_defined("blocking_vmc")) then
    mode      = fdf_get('mode_dmc', 'dmc_one_mpi1')        
    vmc_nstep     = fdf_get('vmc_nstep', 1)    
    vmc_nblk      = fdf_get('vmc_nblk', 1)      
    ! set variable default from %optwf module
    vmc_nblk_max  = fdf_get('vmc_nblk_max', 1)      
  endif 

  if (fdf_defined("blocking_dmc")) then
    dmc_nstep     = fdf_get('dmc_nstep', 1)    
    dmc_nblk      = fdf_get('dmc_nblk', 1)      
    ! set variable default from %optwf module
    dmc_nblk_max  = fdf_get('dmc_nblk_max', 1)      
  endif 
  
! %module optwf
  if (fdf_defined("optwf")) then
    method        = fdf_get('method', 'linear')        
    nwftype = 3; MFORCE = 3
  endif

  call compute_mat_size_new()


  call allocate_vmc()
  call allocate_dmc()






! Some sanity check  !! Make sure that all the variables are parsed before this line

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
    if(ijastrow_parameter .ne. nwftype) then
      write(6,'(''INPUT: block jastrow_parameter missing for one wave function'')')
      write(6,'(''INPUT: jastrow_parameter blocks equal for all wave functions'')')
      call inputjastrow(nwftype)
    endif

    if(iperiodic .eq. 0 .and. ilcao .ne. nwftype) then
      write(6,'(''Warning INPUT: block lcao missing for one wave function'')')
      write(6,'(''Warning INPUT: lcao blocks equal for all wave functions'')')
      call inputlcao(nwftype)
    endif

    if(ideterminants .ne. nwftype) then
      write(6,'(''Warning INPUT: block determinants missing for one wave function'')')
      write(6,'(''Warning INPUT: determinants blocks equal for all wave functions'')')
      call inputdet(nwftype)
    endif

    write(6,*)
endif









! ! &optwf ioptwf 1 ioptci 1 ioptjas 1 ioptorb 1
!   optimize_wavefunction = fdf_boolean("optimize_wavefunction", .false.)
!   write(6,*) ' optimize_wavefunction = ', optimize_wavefunction

!   optimize_ci = fdf_boolean('optimize_ci', .false.)
!   write(6,*) ' optimize_ci = ', optimize_ci

!   optimize_jastrow = fdf_boolean('optimize_jastrow', .false.)
!   write(6,*) ' optimize_jastrow = ', optimize_jastrow

!   optimize_orbitals = fdf_boolean('optimize_orbitals', .false.)
!   write(6,*) ' optimize_orbitals = ', optimize_orbitals

!   write(6,'(A)')  
!   write(6,*) '------------------------------------------------------'


! !Integer numbers (keyword, default_value). The variable is assigned default_value when keyword is not present
!   ! &optwf ncore 0 nextorb 280 no_active 0
!   ! &optwf nblk_max 200 nopt_iter 2
!   ncore = fdf_integer('ncore', 0)
!   write(6,fmt=int_format) 'NCore =', ncore

!   nextorb = fdf_integer('nextorb', 0)
!   write(6,fmt=int_format) 'Next Orb =', nextorb

!   no_active = fdf_integer('no_active', 0)
!   write(6,fmt=int_format) 'no_active =', no_active

!   nblk_max = fdf_integer('nblk_max', 0)
!   write(6,fmt=int_format) 'nblk max =', nblk_max

!   nopt_iter = fdf_integer('nopt_iter', 0)
!   write(6,fmt=int_format) 'nopt_iter =', nopt_iter


! ! floats (keyword, default_value) variable is assigned default_value when keyword is not present

!   ! &optwf sr_tau 0.025 sr_eps 0.001 sr_adiag 0.01
!   ! &optwf isample_cmat 0 energy_tol 0.0
  
!   sr_tau = fdf_get('sr_tau', 0.025d0)
!   write(6,fmt=real_format) 'sr_tau:', sr_tau

!   sr_eps = fdf_get('sr_eps', 0.001d0)
!   write(6,fmt=real_format) 'sr_eps:', sr_eps

!   sr_adiag = fdf_get('sr_adiag', 0.01d0)
!   write(6,fmt=real_format) 'sr_adiag:', sr_adiag

!   energy_tol = fdf_get('energy_tol', 0.00001d0)
!   write(6,fmt=real_format) 'energy_tol:', energy_tol

!   ! &optwf method sr_n multiple_adiag 0

!   opt_method = fdf_get('opt_method', "sr_n")
!   write(6,*) 'Optimization method ', opt_method

!   multiple_adiag = fdf_get('multiple_adiag', .false.)
!   write(6,*) 'multiple_adiag:', multiple_adiag


!   ! logical :: true, .true., yes, T, and TRUE are equivalent
!   debug = fdf_boolean('Debug', .TRUE.)
!   write(6,'(A, L2)') 'Debug:', debug


! ! ianalyt_lap 1 isc 2 nspin1 1 nspin2 1 ifock 0
!   analytic_laplacian = fdf_get('ianalyt_lap', 1)
!   write(6,*) 'analytic laplacian from global.fdf pointer explained ', ianalyt_lap

!   nspin1 = fdf_get('nspin1', 1)
!   write(6,*) 'nspin1 from global.fdf ', nspin1

!   nspin2 = fdf_get('nspin2', 1)
!   write(6,*) 'nspin2 from global.fdf ', nspin2

!   ifock = fdf_get('ifock', 1)
!   write(6,*) 'ifock from global.fdf ', ifock


!   ! mixed types in one line (for example, reading a number with units)
  ! tau = fdf_get('tau', 0.05)
  ! write(6,fmt=real_format) 'DMC tau = ', tau

  ! etrial = fdf_physical('etrial', -20.d0, 'eV')
  ! write(6,fmt=real_format) 'Energy CutOff in eV :: ', etrial

!   write(6,'(A)')  

!   write(6,*) '------------------------------------------------------'



! !  Additional keywords. check if they clash with existing

!   excess_charge = fdf_integer('excess_charge', 0)
!   write(6,fmt=int_format) 'Excess charges =', excess_charge

!   multiplicity = fdf_integer('multiplicity', 1)   ! default multiplicity singlet. An assertion is needed
!   write(6,fmt=int_format) 'multiplicity =', multiplicity

  


! !  write(6,'(A,4X)') 'optimize_wavefunction using bline', (subblock(i), i = 1, 4)

!   if (fdf_block('general', bfdf)) then
!     write(*,*) "inside general block"
!     i = 1
!     do while(fdf_bline(bfdf, pline))    
!       doit = fdf_bsearch(pline, "pool")    
!       write(*,*) "pool found", doit      
!       i = i + 1
!     enddo
!   endif

!   write(6,'(A)')  

!   write(6,*) '------------------------------------------------------'
  
  
  

!   if (.not. fdf_block('molecule', bfdf)) then
!       !   External file reading
!           write(6,*) 'Reading coordinates of the molecule from an external file'
!           ia = 1

!           open (unit=12,file=file_molecule, iostat=iostat, action='read' )
!           if (iostat .ne. 0) stop "Problem in opening the molecule file"
!           read(12,*) natoms
!           print*, "natoms ", natoms
!           if (.not. allocated(cent)) allocate(cent(3,natoms))
          
!           read(12,'(A)')  key
!           print*, "Comment :: ", trim(key)
!           do i = 1, natoms
!             read(12,*) symbol(i), cent(1,i), cent(2,i), cent(3,i)
!           enddo
!           close(12)

!           write(6,*) 'Coordinates from Molecule load construct: '
!           do ia= 1, natoms
!             write(6,'(A4,3F10.6)') symbol(ia), (cent(i,ia),i=1,3)
!           enddo
!   endif
 
!   write(6,'(A)')  
!   write(6,*) '------------------------------------------------------'




!   if (fdf_block('molecule', bfdf)) then
!     !   External file reading
!         write(6,*) 'beginning of external file coordinates block  '
!         ia = 1
! !        write(*,*) "linecount", fdf_block_linecount("molecule")
    
!         do while((fdf_bline(bfdf, pline)))
! !         get the integer from the first line 
!           if ((pline%id(1) .eq. "i") .and. (pline%ntokens .eq. 1)) then        ! check if it is the only integer present in a line
!             natoms = fdf_bintegers(pline, 1)
!             write(*,*) "Number of atoms = ", natoms
!           endif

!           if (.not. allocated(cent)) allocate(cent(3,natoms))
        
!           if (pline%ntokens == 4) then
!             symbol(ia) = fdf_bnames(pline, 1)
!             do i= 1, 3
!               cent(i,ia) = fdf_bvalues(pline, i)
!             enddo
!             ia = ia + 1
!           endif
!         enddo

!         write(6,*) 'Coordinates from single line Molecule block: '
!         do ia= 1, natoms
!           write(6,'(A4,3F10.6)') symbol(ia), (cent(i,ia),i=1,3)
!         enddo
!       endif

!   write(6,'(A)')  

!   write(6,*) '------------------------------------------------------'


! !  Molecule coordinate block begins here  for demonstration

!   if (fdf_block('Coordinates', bfdf)) then
!     ia = 1
!     do while(fdf_bline(bfdf, pline))
!       symbol(ia) = fdf_bnames(pline, 1)
!       do i= 1, 3
!         xa(i,ia) = fdf_bvalues(pline, i)
!       enddo
!       ia = ia + 1
!     enddo
!     write(6,*) 'Coordinates from explicit data block:'
!     do j = 1, ia
!       write(6,'(A, 4x, 3F10.6)') symbol(j), (xa(i,j),i=1,3) 
!     enddo
!   endif

!   write(6,*) '------------------------------------------------------'


!   if (fdf_block('inline_xyz', bfdf)) then
! !   Forward reading 
!     write(6,*) 'Reading an inline_xyz block  '
!     ia = 1

!     do while((fdf_bline(bfdf, pline)))

!       if (pline%ntokens == 1) then
!         number_of_atoms = fdf_bintegers(pline, 1)
!         write(*,*) "Number of atoms", number_of_atoms
!       endif
!       na = number_of_atoms

!       if (pline%ntokens == 4) then
!         symbol(ia) = fdf_bnames(pline, 1)
!         do i= 1, 3
!           xa(i,ia) = fdf_bvalues(pline, i)
!         enddo
!         ia = ia + 1
!       endif
!     enddo

!     write(6,*) 'Inline XYZ Coordinates block:'
!     do ia= 1, na
!       write(6,'(A4,3F10.6)') symbol(ia), (xa(i,ia),i=1,3)
!     enddo
!   endif

!   write(6,'(A)')  

!   write(6,*) '------------------------------------------------------'


!   !  Determinants as a block. read directly from the input file
! !    under construction
!     if (fdf_block('determinants', bfdf)) then
!       ia = 1
!       do while(fdf_bline(bfdf, pline))
!         symbol(ia) = fdf_bnames(pline, 1)
!         do i= 1, 3
!           xa(i,ia) = fdf_bvalues(pline, i)
!         enddo
!         ia = ia + 1
!       enddo
!       na = ia - 1 

!     endif

!     ! if (fdf_block('Coordinates', bfdf)) then
!     !   write(6,*) 'Coordinates:'
!     !   do ia = 1, na
!     !     write(6,'(A, 4x, 3F10.6)') symbol(ia), (xa(i,ia),i=1,3) 
!     !   enddo
!     ! endif


!     write(6,*) '------------------------------------------------------'



!   if (.not. fdf_block('determinants', bfdf)) then
!     if ( fdf_load_defined('determinants') ) then
!       !   External file reading
!         write(6,'(A)')  " Determinants Block"

!         write(6,*) '------------------------------------------------------'      

!         write(6,*) 'Reading the determinants block from an external file '

!         open (unit=11,file=file_determinants, iostat=iostat, action='read' )
!         if (iostat .ne. 0) stop "Problem in opening the determinant file"
!         read(11,*) temp1, temp2, nelectrons, temp3, nalpha

!         read(11,*)  temp1, ndeterminants, nexcitation
!         if (.not. allocated(det_coeff)) allocate(det_coeff(ndeterminants))           

!         read(11,*) (det_coeff(i), i=1,ndeterminants)
!         write(fmt,*)  '(', ndeterminants, '(f11.8,1x))'
!         write(*,fmt) (det_coeff(i), i=1,ndeterminants)
! !        write(*,'(<ndeterminants>(f11.8, 1x))') (det_coeff(i), i=1,ndeterminants)    ! for Intel Fortran  

!         nbeta       = nelectrons - nalpha        
! !       allocate the orbital mapping array        
!         if (.not. allocated(iworbd)) allocate(iworbd(nelectrons, ndeterminants))
        
!         write(*,*) "total number of       electrons ", nelectrons
!         write(*,*) "      number of alpha electrons ", nalpha        
!         write(*,*) "      number of beta  electrons ", nbeta

!         do i = 1, ndeterminants
!           read(11,*) (iworbd(j,i), j=1,nelectrons)
!         enddo

!         write(fmt,*)  '(i4,1x)'        
!         do i = 1, ndeterminants
!           write(*,'(<nelectrons>(i4, 1x))') (iworbd(j,i), j=1,nelectrons)
!         enddo
     
!         read(11,*) temp1
!         if (temp1 == "end" ) write(*,*) "Determinant File read successfully "
!         close(11)

!     endif ! condition if load determinant is present

!   endif ! condition determinant block not present

!   write(6,'(A)')  

!   write(6,*) '------------------------------------------------------'


  call fdf_shutdown()

!----------------------------------------------------------------------------END
end subroutine parser


subroutine flaginit_new
! Flags used to identify the presence or absence of certain blocks/modules
  
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
end subroutine flaginit_new

subroutine compute_mat_size_new()
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

end subroutine compute_mat_size_new