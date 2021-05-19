
subroutine parser
  use fdf     ! modified libfdf
  use prec    ! modified libfdf

  use, intrinsic :: iso_fortran_env, only : iostat_end

! CHAMP modules
  use contr3,         	only: mode
  use contrl_file,    	only: file_input, file_output, file_error
  use contrl_file,    	only: iunit, ounit, errunit
  use allocation_mod, 	only: allocate_vmc, allocate_dmc
  use periodic_table, 	only: atom_t, element

! in the replacement of preprocess input
  use elec,           	only: ndn, nup
  use const,          	only: nelec
  use atom,           	only: nctype, ncent
  use contrl,         	only: nstep, nblk, nblk_max
  use wfsec,          	only: nwftype
  use forcepar,       	only: nforce
  use force_mod,      	only: MFORCE
  use method_opt,     	only: method

! variables from process input
  use sr_mod,         	only: MCONF, MVEC
  use pseudo_mod,     	only: MPS_QUAD
  use properties,     	only: MAXPROP
  use optorb_mod,     	only: MXORBOP, MXREDUCED
  use optci,          	only: MXCITERM
  use mmpol_mod,      	only: mmpolfile_sites, mmpolfile_chmm
  use force_mod,      	only: MFORCE, MWF
  use vmc_mod, 			only: MELEC, MORB, MBASIS, MCENT, MCTYPE, MCTYP3X
  use atom, 			only: znuc, cent, pecent, iwctype, nctype, ncent, ncent_tot, nctype_tot, symbol, atomtyp
  use jaspar, 			only: nspin1, nspin2, is
  use ghostatom, 		only: newghostype, nghostcent
  use const, 			only: pi, hb, etrial, delta, deltai, fbias, nelec, imetro, ipr
  use jaspar1, 			only: cjas1, cjas2
  use general, 			only: pooldir, pp_id, bas_id
  use general, 			only: filenames_bas_num, wforce
  use csfs, 			only: cxdet, ncsf, nstates
  use dets, 			only: cdet, ndet
  use elec, 			only: ndn, nup
  use forcepar, 		only: nforce
  use grdntspar, 		only: igrdtype, ngradnts
  use header, 			only: title
  use jaspar2, 			only: a1, a2
  use jaspar3, 			only: a, b, c, nord, scalek
  use jaspar4, 			only: a4, norda, nordb, nordc
  use jaspar6, 			only: asymp_jasa, asymp_jasb, asymp_r, c1_jas6, c1_jas6i, c2_jas6
  use jaspar6, 			only: cutjas, cutjasi
  use numbas, 			only: numr
  use numbas1, 			only: nbastyp
  use numbas2, 			only: ibas0, ibas1
  use optwf_contrl, 	only: ioptci, ioptjas, ioptorb, ioptwf
  use optwf_contrl, 	only: idl_flag, ilbfgs_flag, ilbfgs_m, dl_mom, dl_alg
  use optwf_contrl, 	only: ibeta, ratio_j, iapprox, ncore
  use optwf_contrl, 	only: iuse_orbeigv
  use optwf_contrl, 	only: no_active
  use optwf_parms, 		only: nparmj
  use optwf_sr_mod, 	only: i_sr_rescale, izvzb
  use pars, 			only: Z, a20, a21
  use rlobxy, 			only: rlobx
  use sa_weights, 		only: iweight, nweight, weights
  use wfsec, 			only: nwftype
  use zmatrix, 			only: izmatrix
  use bparm, 			only: nocuspb, nspin2b
  use casula, 			only: i_vpsp, icasula
  use coefs, 			only: coef, nbasis, norb
  use const2, 			only: deltar, deltat
  use contr2, 			only: ianalyt_lap, ijas
  use contr2, 			only: isc
  use contrldmc, 		only: iacc_rej, icross, icuspg, icut_br, icut_e, idiv_v, idmc, ipq
  use contrldmc, 		only: itau_eff, nfprod, rttau, tau

! Note the additions: Ravindra
  use control_vmc, 		only: vmc_idump,  vmc_irstar, vmc_isite, vmc_nconf, vmc_nblk, vmc_nblk_max
  use control_vmc, 		only: vmc_nblkeq, vmc_nconf_new, vmc_nstep, vmc_icharged_atom, vmc_nblk_ci
! Note the additions: Ravindra
  use control_dmc, 		only: dmc_idump, dmc_irstar, dmc_isite, dmc_nconf, dmc_nblk, dmc_nblk_max
  use control_dmc, 		only: dmc_nblkeq, dmc_nconf_new, dmc_nstep, dmc_icharged_atom, dmc_nblk_ci

  use dorb_m, 			only: iworbd
  use contrl_per, 		only: iperiodic, ibasis
  use force_analy, 		only: iforce_analy, iuse_zmat, alfgeo
  use force_dmc, 		only: itausec, nwprod
  use forcestr,     only: delc
  use wfsec, only: iwftype
  use pseudo, 			only: nloc
  use optorb_cblock, 	only: idump_blockav
  use gradjerrb, 		only: ngrad_jas_blocks
  use qua, 				only: nquad, wq, xq, yq, zq
  use mmpol_cntrl, 		only: ich_mmpol, immpol, immpolprt, isites_mmpol
  use mmpol_parms, 		only: chmm
  use mmpol_fdc, 		only: a_cutoff, rcolm
  use grid3dflag, 		only: i3ddensity, i3dgrid, i3dlagorb, i3dsplorb
  use grid_mod, 		only: UNDEFINED, IUNDEFINED
  use efield, 			only: iefield, ncharges
  use mstates_ctrl, 	only: iefficiency, iguiding, nstates_psig
  use mstates3, 		only: iweight_g, weights_g
  use ci000, 			only: iciprt, nciprim, nciterm
  use pcm_cntrl, 		only: ichpol, ipcm, ipcmprt, isurf
  use pcm_unit, 		only: pcmfile_cavity, pcmfile_chs, pcmfile_chv
  use pcm_parms, 		only: eps_solv, iscov
  use pcm_parms, 		only: ncopcm, nscv, nvopcm
  use prp000, 			only: iprop, ipropprt, nprop
  use pcm_fdc, 			only: qfree, rcolv
  use pcm_grid3d_contrl,only: ipcm_3dgrid
  use pcm_grid3d_param, only: ipcm_nstep3d, pcm_step3d, pcm_origin, pcm_endpt, allocate_pcm_grid3d_param
  use pcm_3dgrid, 		only: PCM_SHIFT, PCM_UNDEFINED, PCM_IUNDEFINED
  use prp003, 			only: cc_nuc
  use method_opt, 		only: method
  use optorb_cblock, 	only: nefp_blocks, isample_cmat, iorbsample
  use orbval, 			only: ddorb, dorb, nadorb, ndetorb, orb
  use array_resize_utils, only: resize_tensor
  use grid3d_param, 	only: endpt, nstep3d, origin, step3d
  use inputflags, 		only: node_cutoff, eps_node_cutoff, dmc_node_cutoff, dmc_eps_node_cutoff, iqmmm, scalecoef
  use optwf_contrl, 	only: energy_tol, dparm_norm_min, nopt_iter, micro_iter_sr
  use optwf_contrl, 	only: nvec, nvecx, alin_adiag, alin_eps, lin_jdav, multiple_adiag
  use optwf_contrl, 	only: ilastvmc, iroot_geo
  use optwf_contrl, 	only: sr_tau , sr_adiag, sr_eps
  use optwf_func, 		only: ifunc_omega, omega0, n_omegaf, n_omegat
  use optwf_corsam, 	only: add_diag
  use dmc_mod, 			only: MWALK

  use optorb_mix,   only: norbopt, norbvirt

  use grdntspar, 		only: delgrdxyz, igrdtype, ngradnts
  use grdntspar, 		only: delgrdba, delgrdbl, delgrdda, ngradnts

  use inputflags, 		only: iznuc, igeometry, ibasis_num, ilcao, iexponents
  use inputflags, 		only: ideterminants, ijastrow_parameter, ioptorb_def, ilattice
  use inputflags, 		only: ici_def, iforces, icsfs, icharge_efield
  use inputflags, 		only: imultideterminants, imodify_zmat, izmatrix_check
  use inputflags, 		only: ihessian_zmat
  use basis,          only: zex

! Note the following modules are new additions


!
  implicit none
!--------------------------------------------------------------- Local Variables
  integer, parameter         :: maxa = 100
  logical                    :: doit, debug

  character(len=72)          :: fname, key
  character(len=20)          :: temp1, temp2, temp3, temp4, temp5
  integer                    :: ifock , ratio, isavebl

  type(block_fdf)            :: bfdf
  type(parsed_line), pointer :: pline

  character(len=100)         :: real_format    = '(A, T40, F14.8)'
  character(len=100)         :: int_format     = '(A, T40, I8)'
  character(len=100)         :: string_format  = '(A, T40, A)'
  character(len=100)         :: fmt32          = '(A, T32, A)'

!------------------------------------------------------------------------- BEGIN
! debug purpose only
  character(len=72)          :: optwf, blocking_vmc, blocking_dmc
  character(len=72)          :: file_basis, 			file_molecule
  character(len=72)          ::	file_determinants, 		file_symmetry
  character(len=72)          :: file_jastrow, 			file_jastrow_der
  character(len=72)          ::	file_orbitals, 			file_pseudo
  character(len=72)          :: file_exponents, 		file_optorb_mixvirt
  character(len=72)          ::	file_multideterminants, file_forces
  character(len=72)          :: file_eigenvalues, 		file_basis_num_info
  character(len=72)          :: file_dmatrix,			file_cavity_spheres
  character(len=72)          :: file_gradients_zmatrix, file_gradients_cartesian
  character(len=72)          :: file_modify_zmatrix,	file_hessian_zmatrix
  character(len=72)          :: file_efield,			file_zmatrix_connection


! from process input subroutine

  character(len=20)          :: fmt
  character(len=32)          :: keyname
  character(len=10)          :: eunit
  character(len=16)          :: cseed
  integer                    :: irn(4), cent_tmp(3)
  integer, allocatable       :: anorm(:) ! dimensions = nbasis

! local counter variables
  integer                    :: i,j,k, iostat
  type(atom_t)               :: atoms
  character(len=2), allocatable   :: unique(:)


! Initialize # get the filenames from the commandline arguments
  call fdf_init(file_input, 'parser.log')


  call flaginit_new()
  !! Number of input variables found so far :: 171


! %module general (complete)
  mode        = fdf_get('mode', 'vmc_one_mpi')
  title       = fdf_get('title', 'Untitled')
  pooldir     = fdf_get('pool', '.')
  pp_id       = fdf_get('pseudopot', 'none')
  bas_id      = fdf_get('basis', 'none')
  nforce      = fdf_get('nforce', 1)
  MFORCE      = nforce
  nwftype     = fdf_get('nwftype', 1)
  iperiodic   = fdf_get('iperiodic', 0)
  ibasis      = fdf_get('ibasis', 1)
!  cseed       = fdf_get('seed', 1)
  ipr         = fdf_get('ipr', -1)
  eunit       = fdf_get('unit', 'Hartrees')
  hb          = fdf_get('mass', 0.5d0)
  scalecoef   = fdf_get('scalecoef',1.0d0)
  i3dgrid     = fdf_get('i3dgrid',0)
  i3dsplorb   = fdf_get('i3dsplorb',0)
  i3dlagorb   = fdf_get('i3dlagorb',0)
  i3ddensity  = fdf_get('i3ddensity',0)

! %module electrons (complete)
  nelec       = fdf_get('nelec', 1)
  nup         = fdf_get('nup', 1)
  ndn         = nelec-nup

! %module atoms (complete)
!  nctype      = fdf_get('nctype', 1)
!  ncent       = fdf_get('natom', 1)
  newghostype = fdf_get('newghostype', 0)
  nghostcent  = fdf_get('nghostcent', 0)

! %module jastrow (complete)
  ijas        = fdf_get('ijas', 1)
  isc         = fdf_get('isc', 1)
  nspin1      = fdf_get('nspin1', 1)
  nspin2      = fdf_get('nspin2', 1)
  ifock       = fdf_get('ifock', 0)
  ianalyt_lap = fdf_get('ianalyt_lap',1)

! %module optgeo (complete)
  iforce_analy= fdf_get('iforce_analy', 0)
  iuse_zmat   = fdf_get('iuse_zmat', 0)
  izvzb       = fdf_get('izvzb', 0)
  alfgeo      = fdf_get('alfgeo', 1.0d0)
  iroot_geo   = fdf_get('iroot_geo', 0)

! %module gradients
  delgrdxyz   = fdf_get('delgrdxyz', 0.001d0)
  igrdtype    = fdf_get('igrdtype', 1)
  ngradnts    = fdf_get('ngradnts', 0)
  delgrdbl    = fdf_get('delgrdbl', 0.001d0)
  delgrdba    = fdf_get('delgrdba', 0.01d0)
  delgrdda    = fdf_get('delgrdda', 0.01d0)
  igrdtype    = fdf_get('igrdtype', 2)
  ngradnts    = fdf_get('ngradnts', 0)

! %module iguiding (complete)
  iguiding      = fdf_get('iguiding',0)
  iefficiency   = fdf_get('iefficiency',0)

! %module efield (complete)
  iefield     = fdf_get('iefield', 0)

! module vmc (complete)
  imetro      = fdf_get('imetro', 6)
  node_cutoff = fdf_get('node_cutoff', 0)
  eps_node_cutoff = fdf_get('enode_cutoff', 1.0d-7)
  delta       = fdf_get('delta', 1.)
  deltar      = fdf_get('deltar', 1.)
  deltat      = fdf_get('deltat', 1.)
  fbias       = fdf_get('fbias', 1.0d0)

! %module vmc / blocking_vmc (complete)
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
  vmc_nblk_ci       = fdf_get('vmc_nblk_ci', vmc_nblk)

!module dmc (complete)
  idmc        = fdf_get('idmc', 2)
  ipq         = fdf_get('ipq', 1)
  itau_eff    = fdf_get('itau_eff', 1)
  iacc_rej    = fdf_get('iacc_rej', 1)
  icross      = fdf_get('icross', 1)
  icuspg      = fdf_get('icuspg', 0)
  idiv_v      = fdf_get('idiv_v', 0)
  icut_br     = fdf_get('icut_br', 0)
  icut_e      = fdf_get('icut_e', 0)
  dmc_node_cutoff = fdf_get('dmc_node_cutoff', 0)
  dmc_eps_node_cutoff = fdf_get('dmc_enode_cutoff', 1.0d-7)
  nfprod      = fdf_get('nfprod', 1)
  tau         = fdf_get('tau', 1.)
  rttau=dsqrt(tau)
  etrial      = fdf_get('etrial', 1.)
  nfprod      = fdf_get('nfprod', 200)
  itausec     = fdf_get('itausec', 1)
  icasula     = fdf_get('icasula', 0)

! %module dmc / blocking_dmc (complete)
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
  dmc_nblk_ci       = fdf_get('dmc_nblk_ci', dmc_nblk)




!optimization flags vmc/dmc
  ioptwf        = fdf_get('ioptwf', 0)
  method        = fdf_get('method', 'linear')
  idl_flag      = fdf_get('idl_flag', 0)
  ilbfgs_flag   = fdf_get('ilbfgs_flag', 0)
  ilbfgs_m      = fdf_get('ilbfgs_m', 5)
  i_sr_rescale  = fdf_get('sr_rescale', 0)
  ibeta         = fdf_get('ibeta', -1)
  ratio         = fdf_get('ratio', ratio_j)
  iapprox       = fdf_get('iapprox', 0)
  ncore         = fdf_get('ncore', 0)
  iuse_orbeigv  = fdf_get('iuse_orbeigv', 0)
  ioptjas       = fdf_get('ioptjas', 0)
  ioptorb       = fdf_get('ioptorb', 0)
  ioptci        = fdf_get('ioptci', 0)
  no_active     = fdf_get('no_active', 0)
  energy_tol    = fdf_get('energy_tol', 1.d-3)
  dparm_norm_min = fdf_get('dparm_norm_min', 1.0d0)
!  add_diag(1)   = fdf_get('add_diag',1.d-6)
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
  ilastvmc      = fdf_get('ilastvmc',1)
  dl_mom        = fdf_get('dl_mom', 0.0)
  dl_alg        = fdf_get('dl_alg','nag')
  ngrad_jas_blocks = fdf_get('ngrad_jas_blocks',0)
  isample_cmat  = fdf_get('isample_cmat', 1)
  isavebl       = fdf_get('save_blocks', 0)
  nefp_blocks   = fdf_get('force_blocks',1)
  iorbsample    = fdf_get('iorbsample',1)
! attention please
!  nadorb        = fdf_get('nextorb', next_max)



! %module ci (complete)
  iciprt        = fdf_get('ci:iciprt',0)


!%module pcm (complete)

  ipcm          = fdf_get('ipcm',0)
  ipcmprt       = fdf_get('ipcmprt',0)
  pcmfile_cavity = fdf_get('file_cavity','pcm000.dat')
  pcmfile_chs   = fdf_get('file_chs','chsurf_old')
  pcmfile_chv   = fdf_get('file_chv','chvol_old')
!  nscv          = fdf_get('nblk_chv',nblk)   ! attention
!  iscov         = fdf_get('nstep_chv',nstep2) ! attention
  eps_solv      = fdf_get('eps_solv',1)
!  fcol          = fdf_get('fcol',1.d0)
  rcolv         = fdf_get('rcolv',0.04d0)
!  npmax         = fdf_get('npmax',1)


  call allocate_pcm_grid3d_param()
  ipcm_3dgrid   = fdf_get('ipcm_3dgrid',0)
  ipcm_nstep3d(1) 	= fdf_get('nx_pcm',PCM_IUNDEFINED)
  ipcm_nstep3d(2) 	= fdf_get('ny_pcm',PCM_IUNDEFINED)
  ipcm_nstep3d(3) 	= fdf_get('nz_pcm',PCM_IUNDEFINED)
  pcm_step3d(1) 	= fdf_get('dx_pcm',PCM_UNDEFINED)
  pcm_step3d(2) 	= fdf_get('dy_pcm',PCM_UNDEFINED)
  pcm_step3d(3) 	= fdf_get('dz_pcm',PCM_UNDEFINED)
  pcm_origin(1) 	= fdf_get('x0_pcm',PCM_UNDEFINED)
  pcm_origin(2) 	= fdf_get('y0_pcm',PCM_UNDEFINED)
  pcm_origin(3) 	= fdf_get('z0_pcm',PCM_UNDEFINED)
  pcm_endpt(1)  	= fdf_get('xn_pcm',PCM_UNDEFINED)
  pcm_endpt(2)  	= fdf_get('yn_pcm',PCM_UNDEFINED)
  pcm_endpt(3)  	= fdf_get('zn_pcm',PCM_UNDEFINED)
  PCM_SHIFT     	= fdf_get('shift',4.d0)

! %module mmpol (complete)
  immpol        = fdf_get('immpol',0)
  immpolprt     = fdf_get('immpolprt',0)
  mmpolfile_sites = fdf_get('file_sites','mmpol000.dat')
  mmpolfile_chmm = fdf_get('file_mmdipo','mmdipo_old')
  a_cutoff      = fdf_get('a_cutoff',2.5874d0)
  rcolm         = fdf_get('rcolm',0.04d0)

! %module properties (complete)
  iprop         = fdf_get('sample',0)
  ipropprt      = fdf_get('print',0)
  nquad         = fdf_get('nquad',6)

! %module pseudo (complete)
  nloc          = fdf_get('nloc',0)

! %module qmmm (complete)
!  iqmm          = fdf_get('iqmm',0)


  ! Filenames parsing
  file_basis        		    = fdf_load_filename('basis', 			'default.bas')
  file_molecule     		    = fdf_load_filename('molecule', 		'default.xyz')
  file_determinants 		    = fdf_load_filename('determinants', 	'default.det')
  file_symmetry     		    = fdf_load_filename('symmetry', 		'default.sym')
  file_jastrow      		    = fdf_load_filename('jastrow', 			'default.jas')
  file_jastrow_der  		    = fdf_load_filename('jastrow_der', 		'default.jasder')
  file_orbitals     		    = fdf_load_filename('orbitals', 		'default.orb')
  file_exponents    		    = fdf_load_filename('exponents', 		'exponents.exp')
  file_pseudo 			        = fdf_load_filename('pseudo', 			'default.psp')
  file_optorb_mixvirt       = fdf_load_filename('optorb_mixvirt', 	'default.mix')
  file_multideterminants    = fdf_load_filename('multideterminants', 'default.mdet')
  file_forces       		    = fdf_load_filename('forces', 			'default.for')
  file_eigenvalues	       	= fdf_load_filename('eigenvalues', 		'default.eig')
  file_basis_num_info       = fdf_load_filename('basis_num_info', 	'default.bni')
  file_dmatrix		       	  = fdf_load_filename('dmatrix', 			'default.dmat')
  file_cavity_spheres       = fdf_load_filename('cavity_spheres', 	'default.cav')
  file_gradients_zmatrix    = fdf_load_filename('gradients_zmatrix','default.gzmat')
  file_gradients_cartesian  = fdf_load_filename('gradients_cartesian', 'default.gcart')
  file_modify_zmatrix       = fdf_load_filename('modify_zmatrix', 	'default.mzmat')
  file_hessian_zmatrix      = fdf_load_filename('hessian_zmatrix', 	'default.hzmat')
  file_zmatrix_connection   = fdf_load_filename('zmatrix_connection', 'default.zmcon')
  file_efield	       		    = fdf_load_filename('efield', 			'default.efield')

  call header_printing()

! Processing of data read from the parsed files or setting them with defaults

! (1) Molecular geometry file exclusively in .xyz format [#####]

  if ( fdf_load_defined('molecule') ) then
    call read_molecule_file(file_molecule)
  elseif (fdf_block('molecule', bfdf)) then
    call fdf_read_molecule_block(bfdf)
  else
    write(errunit,'(a)') "Error:: No information about molecular coordiates provided."
    write(errunit,'(3a,i6)') "Stats for nerds :: in file ",__FILE__, " at line ", __LINE__
    error stop
  endif

! (2) Determinants (excluding csf and csfmap) [#####]

  if ( fdf_load_defined('determinants') ) then
    call read_determinants_file(file_determinants)
  elseif ( fdf_block('determinants', bfdf)) then
  ! call fdf_read_determinants_block(bfdf)
    write(errunit,'(a)') "Error:: No information about determinants provided."
    write(errunit,'(3a,i6)') "Stats for nerds :: in file ",__FILE__, " at line ", __LINE__
    error stop
  else
    call inputdet()
  endif

! (3) CSF [#####]

  if ( fdf_load_defined('determinants') .and. ndet .gt. 1 ) then
    call read_csf_file(file_determinants)
  elseif (fdf_block('csf', bfdf)) then
    call fdf_read_csf_block(bfdf)
  else
    ! No csf present; set default values (in replacement of inputcsf)
    nstates = 1 ; ncsf = 0
    if (ioptci .ne. 0) nciterm = nciprim
  endif

! (4) CSFMAP [#####]

  if ( fdf_load_defined('determinants') .and. ndet .gt. 1 ) then
    call read_csfmap_file(file_determinants)
  elseif (fdf_block('determinants', bfdf)) then
  ! call fdf_read_csfmap_block(bfdf)
    write(errunit,'(a)') "Error:: No information about csfmaps provided."
    write(errunit,'(3a,i6)') "Stats for nerds :: in file ",__FILE__, " at line ", __LINE__
    error stop
  endif

! (5) Jastrow Parameters (either block or from a file)

  if ( fdf_load_defined('jastrow') ) then
      call read_jastrow_file(file_jastrow)
  elseif (fdf_block('jastrow', bfdf)) then
  !call fdf_read_jastrow_block(bfdf)
    write(errunit,'(a)') "Error:: No information about jastrow provided in the block."
    write(errunit,'(3a,i6)') "Stats for nerds :: in file ",__FILE__, " at line ", __LINE__
    error stop
  else
    ! no information about jastrow present. Set some values
    call inputjastrow()
  endif

! (6) LCAO orbitals

  if ( fdf_load_defined('orbitals') ) then
    call read_orbitals_file(file_orbitals)
  elseif ( fdf_block('orbitals', bfdf)) then
  ! call fdf_read_orbitals_block(bfdf)
    write(errunit,'(a)') "Error:: No information about orbitals provided."
    write(errunit,'(3a,i6)') "Stats for nerds :: in file ",__FILE__, " at line ", __LINE__
    error stop
  else
    call inputlcao()
  endif

! (7) exponents

  if ( fdf_load_defined('exponents') ) then
    call read_exponents_file(file_exponents)
  elseif ( fdf_block('exponents', bfdf)) then
  ! call fdf_read_exponents_block(bfdf)
    write(errunit,'(a)') "Error:: No information about exponents provided."
    write(errunit,'(3a,i6)') "Stats for nerds :: in file ",__FILE__, " at line ", __LINE__
    error stop
  else
    zex = 1   ! debug check condition about numr == 0
  endif

! (8) Jastrow derivative Parameters (either block or from a file)

  if ( fdf_load_defined('jastrow_der') ) then
    call read_jasderiv_file(file_jastrow_der)
  elseif ( fdf_block('jastrow_der', bfdf)) then
  !call fdf_read_jastrow_derivative_block(bfdf)
    write(errunit,'(a)') "Error:: No information about jastrow derivatives provided in the block."
    write(errunit,'(3a,i6)') "Stats for nerds :: in file ",__FILE__, " at line ", __LINE__
    error stop
  else
    write(errunit,'(a)') "Error:: No information about jastrow derivatives provided in the block."
    write(errunit,'(3a,i6)') "Stats for nerds :: in file ",__FILE__, " at line ", __LINE__
    error stop
  endif

! (9) Symmetry information of orbitals (either block or from a file)

  if ( fdf_load_defined('symmetry') ) then
    call read_symmetry_file(file_symmetry)
  elseif ( fdf_block('symmetry', bfdf)) then
  ! call fdf_read_symmetry_block(bfdf)
    write(errunit,'(a)') "Error:: No information about orbital symmetries provided in the block."
    write(errunit,'(3a,i6)') "Stats for nerds :: in file ",__FILE__, " at line ", __LINE__
    error stop
  else
    write(errunit,'(a)') "Error:: No information about orbital symmetries provided in the block."
    write(errunit,'(3a,i6)') "Stats for nerds :: in file ",__FILE__, " at line ", __LINE__
    error stop
  endif

! (10) optorb_mixvirt information of orbitals (either block or from a file)

  if(ioptorb .ne. 0) then   ! read next file only if orb optimization is requested
    if ( fdf_load_defined('optorb_mixvirt') ) then
      call read_optorb_mixvirt_file(file_optorb_mixvirt)
    elseif ( fdf_block('optorb_mixvirt', bfdf)) then
    ! call fdf_read_optorb_mixvirt_block(bfdf)
      write(errunit,'(a)') "Error:: No information about optorb_mixvirt provided in the block."
      write(errunit,'(3a,i6)') "Stats for nerds :: in file ",__FILE__, " at line ", __LINE__
      error stop
    else
      !use optorb_mix, only: norbopt, norbvirt
      norbopt   = 0
      norbvirt  = 0
    endif

    if(ioptorb_def.eq.0) then
      write(6,'(''INPUT: definition of orbital variations missing'')')
      call optorb_define
    endif
  endif

! (11) Eigenvalues information of orbitals (either block or from a file)

  if ( fdf_load_defined('eigenvalues') ) then
    call read_eigenvalues_file(file_eigenvalues)
  elseif ( fdf_block('eigenvalues', bfdf)) then
  ! call fdf_read_eigenvalues_block(bfdf)
    write(errunit,'(a)') "Error:: No information about eigenvalues provided in the block."
    write(errunit,'(3a,i6)') "Stats for nerds :: in file ",__FILE__, " at line ", __LINE__
    error stop
  else
    write(errunit,'(a)') "Error:: No information about eigenvalues provided in the block."
    write(errunit,'(3a,i6)') "Stats for nerds :: in file ",__FILE__, " at line ", __LINE__
!    error stop
  endif

! (12) Basis num information (either block or from a file)

  if ( fdf_load_defined('basis_num_info') ) then
    call read_basis_num_info_file(file_basis_num_info)
  elseif (.not. fdf_block('basis_num_info', bfdf)) then
  ! call fdf_read_eigenvalues_block(bfdf)
    write(errunit,'(a)') "Error:: No information about eigenvalues provided in the block."
    write(errunit,'(3a,i6)') "Stats for nerds :: in file ",__FILE__, " at line ", __LINE__
    error stop
  else
    write(errunit,'(a)') "Error:: No information about eigenvalues provided in the block."
    write(errunit,'(3a,i6)') "Stats for nerds :: in file ",__FILE__, " at line ", __LINE__
    error stop
  endif

! (13) Forces information (either block or from a file) [#####]

  if (fdf_load_defined('forces') ) then
    call read_forces_file(file_forces)
  elseif (fdf_block('forces', bfdf)) then
    call fdf_read_forces_block(bfdf)
  else
    call inputforces()
  endif

! (14) Dmatrix information (either block or from a file)

  if ( fdf_load_defined('dmatrix') ) then
    call read_dmatrix_file(file_dmatrix)
  elseif (fdf_block('dmatrix', bfdf)) then
  ! call fdf_read_dmatrix_block(bfdf)
    write(errunit,'(a)') "Error:: No information about dmatrix provided in the block."
    write(errunit,'(3a,i6)') "Stats for nerds :: in file ",__FILE__, " at line ", __LINE__
    error stop
  else
    write(errunit,'(a)') "Error:: No information about dmatrix provided in the input."
    write(errunit,'(3a,i6)') "Stats for nerds :: in file ",__FILE__, " at line ", __LINE__
!    error stop
  endif

! (15) pseudo information (either block or from a file)

  if ( fdf_defined('pseudopot') ) then
    call readps_gauss()
  elseif (nloc .eq. 0) then
    write(ounit,'(a)') "Warning:: Is this an all electron calculation?"
  else
    write(errunit,'(a)') "Error:: No information about pseudo provided in the input."
    write(errunit,'(3a,i6)') "Stats for nerds :: in file ",__FILE__, " at line ", __LINE__
  endif

! (16) basis information (either block or from a file)

  if ( fdf_defined('basis') ) then
    call read_bas_num(1)   ! i == iwf debug
  else
    write(errunit,'(a)') "Error:: No information about basis provided in the block."
    write(errunit,'(3a,i6)') "Stats for nerds :: in file ",__FILE__, " at line ", __LINE__
    error stop
  endif

! (17) multideterminants information (either block or from a file)

  if ( fdf_load_defined('multideterminants') ) then
    call read_multideterminants_file(file_multideterminants)
  elseif (fdf_block('multideterminants', bfdf)) then
  ! call fdf_read_multideterminants_block(bfdf)
    write(errunit,'(a)') "Error:: No information about multideterminants provided in the block."
    write(errunit,'(3a,i6)') "Stats for nerds :: in file ",__FILE__, " at line ", __LINE__
    error stop
  else
    call multideterminants_define(0,0)
  endif

! (18) cavity_spheres information (either block or from a file)

  if ( fdf_load_defined('cavity_spheres') ) then
    call read_cavity_spheres_file(file_cavity_spheres)
  elseif ( fdf_block('cavity_spheres', bfdf)) then
  ! call fdf_read_cavity_spheres_block(bfdf)
    write(errunit,'(a)') "Error:: No information about cavity_spheres provided in the block."
    write(errunit,'(3a,i6)') "Stats for nerds :: in file ",__FILE__, " at line ", __LINE__
    error stop
  else
    write(errunit,'(a)') "Error:: No information about cavity_spheres provided in the block."
    write(errunit,'(3a,i6)') "Stats for nerds :: in file ",__FILE__, " at line ", __LINE__
!    error stop
  endif

! (19) gradients_zmatrix information (either block or from a file)

  if ( fdf_load_defined('gradients_zmatrix') ) then
    call read_gradients_zmatrix_file(file_gradients_zmatrix)
  elseif ( fdf_block('gradients_zmatrix', bfdf)) then
  ! call fdf_read_gradients_zmatrix_block(bfdf)
    write(errunit,'(a)') "Error:: No information about gradients_zmatrix provided in the block."
    write(errunit,'(3a,i6)') "Stats for nerds :: in file ",__FILE__, " at line ", __LINE__
    error stop
  else
    write(errunit,'(a)') "Error:: No information about gradients_zmatrix provided in the block."
    write(errunit,'(3a,i6)') "Stats for nerds :: in file ",__FILE__, " at line ", __LINE__
!    error stop
  endif

! (20) gradients_cartesian information (either block or from a file)

  if ( fdf_load_defined('gradients_cartesian') ) then
    call read_gradients_cartesian_file(file_gradients_cartesian)
  elseif ( fdf_block('gradients_cartesian', bfdf)) then
  ! call fdf_read_gradients_cartesian_block(bfdf)
    write(errunit,'(a)') "Error:: No information about gradients_cartesian provided in the block."
    write(errunit,'(3a,i6)') "Stats for nerds :: in file ",__FILE__, " at line ", __LINE__
    error stop
  else
    write(errunit,'(a)') "Error:: No information about gradients_cartesian provided in the block."
    write(errunit,'(3a,i6)') "Stats for nerds :: in file ",__FILE__, " at line ", __LINE__
!    error stop
  endif

! (21) modify_zmatrix information (either block or from a file)

  if(iforce_analy.gt.0) then
    if ( fdf_load_defined('modify_zmatrix') ) then
      call read_modify_zmatrix_file(file_modify_zmatrix)
    elseif ( fdf_block('modify_zmatrix', bfdf)) then
    ! call fdf_read_modify_zmatrix_block(bfdf)
      write(errunit,'(a)') "Error:: No information about modify_zmatrix provided in the block."
      write(errunit,'(3a,i6)') "Stats for nerds :: in file ",__FILE__, " at line ", __LINE__
      error stop
    else
      call modify_zmat_define
    endif
  endif

! (22) hessian_zmatrix information (either block or from a file)
  if(iforce_analy.gt.0) then
    if ( fdf_load_defined('hessian_zmatrix') ) then
      call read_hessian_zmatrix_file(file_hessian_zmatrix)
    elseif ( fdf_block('hessian_zmatrix', bfdf)) then
    ! call fdf_read_hessian_zmatrix_block(bfdf)
      write(errunit,'(a)') "Error:: No information about hessian_zmatrix provided in the block."
      write(errunit,'(3a,i6)') "Stats for nerds :: in file ",__FILE__, " at line ", __LINE__
      error stop
    else
      call hessian_zmat_define
    endif
  endif

! (23) zmatrix_connection information (either block or from a file)

  if ( fdf_load_defined('zmatrix_connection') ) then
    call read_zmatrix_connection_file(file_zmatrix_connection)
  elseif ( fdf_block('zmatrix_connection', bfdf)) then
  ! call fdf_read_zmatrix_connection_block(bfdf)
    write(errunit,'(a)') "Error:: No information about zmatrix_connection provided in the block."
    write(errunit,'(3a,i6)') "Stats for nerds :: in file ",__FILE__, " at line ", __LINE__
    error stop
  else
    write(errunit,'(a)') "Error:: No information about zmatrix_connection provided in the block."
    write(errunit,'(3a,i6)') "Stats for nerds :: in file ",__FILE__, " at line ", __LINE__
!    error stop
  endif

! (24) efield information (either block or from a file)

  if ( fdf_load_defined('efield') ) then
    call read_efield_file(file_efield)
  elseif ( fdf_block('efield', bfdf)) then
  ! call fdf_read_efield_block(bfdf)
    write(errunit,'(a)') "Error:: No information about efield provided in the block."
    write(errunit,'(3a,i6)') "Stats for nerds :: in file ",__FILE__, " at line ", __LINE__
    error stop
  else
    write(errunit,'(a)') "Error:: No information about efield provided in the block."
    write(errunit,'(3a,i6)') "Stats for nerds :: in file ",__FILE__, " at line ", __LINE__
!    error stop
  endif

! Done reading all the files


! %module optwf
  if (fdf_defined("optwf")) then
    nwftype = 3; MFORCE = 3
  endif

  call compute_mat_size_new()
  call allocate_vmc()
  call allocate_dmc()

  call fdf_shutdown()

  ! The following portion can be shifted to another subroutine.
  ! It does the processing of the input read so far and initializes some
  ! arrays if something is missing.

  if(ioptci.ne.0 .and. ici_def.eq.0) then
    write(6,'(''INPUT: definition of OPTCI operators missing'')')
    call optci_define
  endif

!----------------------------------------------------------------------------END
  contains

  !! Here all the subroutines that handle the block data are written

  subroutine fdf_read_molecule_block(bfdf)
    implicit none

    type(block_fdf)            :: bfdf
    type(parsed_line), pointer :: pline

    write(ounit,*) ' Molecular Coordinates from molecule block '

	  j = 1 !local counter
    do while((fdf_bline(bfdf, pline)))
!     get the integer from the first line
      if ((pline%id(1) .eq. "i") .and. (pline%ntokens .eq. 1)) then  ! check if it is the only integer present in a line
        ncent = fdf_bintegers(pline, 1)
        write(ounit,fmt=int_format) " Number of atoms ::  ", ncent
      endif

      if (.not. allocated(cent)) allocate(cent(3,ncent))
      if (.not. allocated(symbol)) allocate(symbol(ncent))
      if (.not. allocated(iwctype)) allocate(iwctype(ncent))
      if (.not. allocated(unique)) allocate(unique(ncent))

      if (pline%ntokens .ne. 4) then  ! check if it is the only integer present in a line
        write(ounit,*) " Comment from the file ::  ", trim(pline%line)
      endif


      if (pline%ntokens == 4) then
        symbol(j) = fdf_bnames(pline, 1)
        do i= 1, 3
          cent(i,j) = fdf_bvalues(pline, i)
        enddo
        j = j + 1
      endif
    enddo


    ncent_tot = ncent + nghostcent

    ! Count unique type of elements
    nctype = 1
    unique(1) = symbol(1)

    do j= 2, ncent
        if (any(unique == symbol(j) ))  cycle
        nctype = nctype + 1
        unique(nctype) = symbol(j)
    enddo

    write(ounit,*) " Number of distinct types of elements (nctype) :: ", nctype
    write(ounit,*)

    if (.not. allocated(atomtyp)) allocate(atomtyp(nctype))
    if (.not. allocated(znuc)) allocate(znuc(nctype))

    ! get the correspondence for each atom according to the rule defined for atomtypes
    do j = 1, ncent
        do k = 1, nctype
            if (symbol(j) == unique(k))   iwctype(j) = k
        enddo
    enddo

    ! Get the correspondence rule
    do k = 1, nctype
        atomtyp(k) = unique(k)
    enddo
    if (allocated(unique)) deallocate(unique)

    ! Get the znuc for each unique atom
    do j = 1, nctype
        atoms = element(atomtyp(j))
        znuc(j) = atoms%nvalence
    enddo

    nctype_tot = nctype + newghostype

    write(ounit,*) 'Atomic symbol, coordinates, and iwctype from the molecule coordinates file '
    write(ounit,*)
    do j= 1, ncent
        write(ounit,'(A4, 2x, 3F12.6, 2x, i3)') symbol(j), (cent(i,j),i=1,3), iwctype(j)
    enddo

    write(ounit,*)
    write(ounit,*) " Values of znuc (number of valence electrons) "
    write(ounit,'(10F12.6)') (znuc(j), j = 1, nctype)
    write(ounit,*)

    write(ounit,*) '------------------------------------------------------'
  end subroutine fdf_read_molecule_block

  subroutine fdf_read_forces_block(bfdf)
    implicit none
    type(block_fdf)            :: bfdf
    type(parsed_line), pointer :: pline
    integer                    :: i,j,k

    if (.not. allocated(delc)) allocate (delc(3, ncent, nforce))
    if (.not. allocated(iwftype)) allocate (iwftype(nforce))

    i = 1; j = 1
    do while((fdf_bline(bfdf, pline)))
      if (pline%ntokens == 1) i = fdf_bintegers(pline, 1)
      if (pline%ntokens == 3) then
        do k = 1, 3
          delc(k, j, i) = fdf_bvalues(pline, k)
        enddo ! xyz
        j = j + 1
      endif ! expect only three values in a line
    enddo ! parse entire file

    write(ounit,*) 'Force displacements from the %block forces  '
    write(ounit,*)
    do i = 1, nforce
      write(ounit,'(a,i4)') 'Number ::',i
      write(ounit,*) '-----------------------------------------------------------------------'
      write(ounit,'(a, t15, a, t27, a, t39, a, t45)') 'Symbol', 'x', 'y', 'z'
      write(ounit,'(t14, a, t26, a, t38, a )') '(A)', '(A)', '(A)'
      write(ounit,*) '-----------------------------------------------------------------------'
      do j= 1, ncent
        write(ounit,'(A4, 2x, 3F12.6)') symbol(j), (delc(k, j, i),k=1,3)
      enddo
    enddo
    write(ounit,*)
    write(ounit,*) '-----------------------------------------------------------------------'
  end subroutine fdf_read_forces_block

  subroutine fdf_read_csf_block(bfdf)
    use csfs, only: ccsf, ncsf, nstates
    implicit none

    type(block_fdf)            :: bfdf
    type(parsed_line), pointer :: pline
    integer                    :: i,j,k

!   Format of the CSF block
!   %block csf
!    ncsf 200     nstates 2
!   # State 1           State 2           State N     ...
!   csf(1,1,1)          csf(1,2,1)        csf(1,N,1)
!   csf(2,1,1)          csf(2,2,1)        csf(2,N,1)
!   .                   .                 .
!   csf(ncsf,1,1)       csf(ncsf,2,1)     csf(ncsf,N,1)
!   %endblock

!     First get the two numbers required for allocations
    ncsf    = fdf_bintegers(bfdf%mark%pline, 1) ! 1st integer in the line
    nstates = fdf_bintegers(bfdf%mark%pline, 2) ! 2nd integer in the line
    if (.not. allocated(ccsf)) allocate(ccsf(ncsf, nstates, nwftype))

    j = 1
    do while((fdf_bline(bfdf, pline)))

      if (pline%ntokens == nstates) then
        do k = 1, nstates
          ccsf(j,k,1) = fdf_bvalues(pline, k)
        enddo
        j = j + 1
      endif
    enddo

    write(ounit,*)
    write(ounit,*) " CSF coefficients from %block csf"

    write(ounit,'(10(1x, a9, i3, 1x))') ((" State: ", i), i =1, nstates)
    do j = 1, ncsf
      write(ounit,'(10(1x, f12.8, 1x))') (ccsf(j,i,1), i=1,nstates)
    enddo

  end subroutine fdf_read_csf_block

  subroutine fdf_read_jastrow_block(bfdf)

    implicit none

    type(block_fdf)            :: bfdf
    type(parsed_line), pointer :: pline
    integer                    :: i,j,k, iwft
    integer                    :: mparmja, mparmjb, mparmjc, nterms4
!   Format of the jastrow block
!   %block csf
!   jastrow_parameter 1
!    5  5  0           norda,nordb,nordc
!     0.60000000   0.00000000     scalek,a21
!     0.00000000   0.00000000   0.05946443  -0.68575835   0.42250502  -0.10845009 (a(iparmj),iparmj=1,nparma)
!     0.00000000   0.00000000  -0.13082284  -0.06620300   0.18687803  -0.08503472 (a(iparmj),iparmj=1,nparma)
!     0.50000000   0.38065787   0.16654238  -0.05430118   0.00399345   0.00429553 (b(iparmj),iparmj=1,nparmb)
!   %endblock


    allocate (scalek(nwftype))


    j = 1
    do while((fdf_bline(bfdf, pline)))
!     First get the first number required for allocations
      if ((pline%id(2) .eq. "i") .and. (pline%ntokens .eq. 2)) then  ! check if it is the only integer present in a line
        iwft  = fdf_bintegers(pline, 1) ! 1st integer in the line
      endif

      if ((pline%id(2) .eq. "i") .and. (pline%id(4) .eq. "i") .and. (pline%id(6) .eq. "i")) then
        norda = fdf_bintegers(pline, 1) ! 1st integer in the line
        nordb = fdf_bintegers(pline, 2) ! 2nd integer in the line
        nordc = fdf_bintegers(pline, 3) ! 3rd integer in the line

        mparmja = 2 + max(0, norda - 1)
        mparmjb = 2 + max(0, nordb - 1)
        mparmjc = nterms4(nordc)

        print*, "mparmja ", mparmja
        print*, "mparmjb ", mparmjb
        print*, "mparmjc ", mparmjc
      endif

      if ((pline%id(2) .eq. "r") .and. (pline%id(4) .eq. "r") .and. (pline%id(1) .eq. "n") .and. (pline%id(3) .eq. "n")) then  ! check if it is the only integer present in a line
        scalek(iwft) = fdf_bvalues(pline, 1) ! 1st integer in the line
        print*, "scalek ", scalek(iwft)
        a21          = fdf_bvalues(pline, 2) ! 2nd integer in the line
        print*, "a21 ", a21
      endif

      !print*, "check all the ids ", pline%id(1:6)
      ! fix this part

      ! if ((pline%id(1) .eq. "n") .and. (pline%ntokens .eq. mparmja+1) ) then  ! check if it is the only integer present in a line
      !   do j = 1, nctype
      !     do i = 1, mparmja
      !       a4(i,j,iwft) = fdf_bvalues(pline, i)
      !       print *, a4(i,j,iwft)
      !     enddo
      !   enddo
      ! endif



    enddo
  end subroutine fdf_read_jastrow_block
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
