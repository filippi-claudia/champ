module parser_mod
      use error,   only: fatal_error
implicit none
contains
subroutine parser
  !> This subroutine parses the input file using the modified libfdf parser
  !> and assigns the values to variables and arrays of different modules.
  !> @author Ravindra Shinde
  !> @email  r.l.shinde@utwente.nl
  !> @date   11-08-2021
  !> @version 1.2

      use fdf               ! modified libfdf
      use allocation_mod, only: allocate_dmc,allocate_vmc
      use array_resize_utils, only: resize_tensor
      use basis,   only: zex
      use bparm,   only: nocuspb,nspin2b
      use casula,  only: i_vpsp,icasula
      use ci000,   only: iciprt,nciprim,nciterm
      use coefs,   only: nbasis,next_max
      use const,   only: etrial, esigmatrial
      use constants, only: hb,pi
      use contrl_file, only: errunit,file_error,file_input
      use contrl_file, only: file_output,iunit,ounit
      use contrl_per, only: ibasis,iperiodic
      use contrldmc, only: iacc_rej,icross,icuspg,icut_br,icut_e,idiv_v, ibranching_c
      use contrldmc, only: idmc,ipq,itau_eff,nfprod,rttau,tau,limit_wt_dmc
      use control, only: ipr,mode
      use control_dmc, only: dmc_idump,dmc_irstar,dmc_isite,dmc_nblk
      use control_dmc, only: dmc_nblkeq,dmc_nconf,dmc_nconf_new
      use control_dmc, only: dmc_nstep
      use control_vmc, only: vmc_icharged_atom,vmc_idump,vmc_irstar
      use control_vmc, only: vmc_isite,vmc_nblk,vmc_nblk_ci,vmc_nblk_max
      use control_vmc, only: vmc_nblkeq,vmc_nconf,vmc_nconf_new
      use control_vmc, only: vmc_nstep
      use csfs,    only: anormo,ccsf,cxdet,iadet,ibdet,icxdet,maxcsf,ncsf,nstates
      use cuspinit4_mod, only: cuspinit4
      use custom_broadcast, only: bcast
      use dets,    only: nmap
      use dmc_mod, only: mwalk,set_mwalk
      use dorb_m,  only: iworbd
      use efield,  only: iefield,ncharges
      use efield_f_mod, only: efield_compute_extint
      use fragments, only: nfrag, ifragcent, ibranching_cfrag, etrialfrag
      use general, only: bas_id,filenames_bas_num,pooldir,pp_id,write_walkalize
      use get_norbterm_mod, only: get_norbterm
      use gradjerrb, only: ngrad_jas_blocks
      use grdntspar, only: delgrdba,delgrdbl,delgrdda,delgrdxyz,igrdtype
      use grdntspar, only: ngradnts
      use grid3d,  only: setup_grid
      use grid3d_orbitals, only: setup_3dlagorb,setup_3dsplorb
      use grid3d_param, only: endpt,nstep3d,origin,step3d
      use grid3dflag, only: i3ddensity,i3dgrid,i3dlagorb,i3dsplorb
      use grid_mod, only: IUNDEFINED,UNDEFINED
      use header,  only: title
      use inputflags, only: dmc_eps_node_cutoff,dmc_node_cutoff
      use inputflags, only: eps_node_cutoff,ibasis_num,icharge_efield
      use inputflags, only: ici_def,icsfs,ideterminants,iexponents
      use inputflags, only: iforces,igeometry,igradients,ihessian_zmat
      use inputflags, only: ijastrow_parameter,ilattice,ilcao
      use inputflags, only: imodify_zmat,imultideterminants,ioptorb_def
      use inputflags, only: ioptorb_mixvirt,iqmmm,izmatrix_check,iznuc
      use inputflags, only: node_cutoff,scalecoef
      use jastrow, only: norda,nordb,nordc
      use jastrow, only: asymp_r
      use jastrow, only: a4,allocate_jasasymp,asymp_jasa,asymp_jasb,b,c
      use jastrow, only: ijas,ijas_lr,is,isc,neqsx,nordj,nordj1
      use jastrow, only: nspin1,nspin2,scalek
      use jastrow4_mod, only: nterms4
      use m_force_analytic, only: alfgeo,iforce_analy,iuse_zmat
      use metropolis, only: imetro, vmc_tau
      use metropolis, only: delta,deltai,deltar,deltat,fbias
      use misc_grdnts, only: inpwrt_grdnts_cart,inpwrt_grdnts_zmat
      use misc_grdnts, only: inpwrt_zmatrix
      use mmpol_cntrl, only: ich_mmpol,immpol,immpolprt,isites_mmpol
      use mmpol_fdc, only: a_cutoff,rcolm
      use mmpol_mod, only: mmpolfile_chmm,mmpolfile_sites
      use mmpol_parms, only: chmm
      use mpiconf, only: wid, idtask
      use mpitimer, only: elapsed_time
      use mstates3, only: iweight_g,weights_g
      use mstates_ctrl, only: iefficiency,iguiding,nstates_psig
      use mstates_mod, only: MSTATES
      use multidet, only: kref_fixed
      use multiple_geo, only: alfstr,istrech
      use multiple_geo, only: MFORCE,MWF,delc,itausec,iwftype,nforce
      use multiple_geo, only: nwftype,nwprod,pecent
      use numbas,  only: numr
      use numbas1, only: nbastyp
      use numbas2, only: ibas0,ibas1
      use optci,   only: mxciterm
      use optci_mod, only: optci_define
      use optorb,  only: irrep
      use optorb_cblock, only: idump_blockav,iorbsample,isample_cmat
      use optorb_cblock, only: nefp_blocks,norbterm
      use optorb_f_mod, only: optorb_define
      use optorb_mix, only: norbopt,norbvirt
      use optorb_mod, only: mxreduced
      use optwf_control, only: alin_adiag,alin_eps,dl_alg,dl_mom
      use optwf_control, only: dparm_norm_min,energy_tol,iapprox,ibeta
      use optwf_control, only: idl_flag,ilastvmc,ilbfgs_flag,ilbfgs_m
      use optwf_control, only: ioptci,ioptjas,ioptorb,ioptwf,iroot_geo
      use optwf_control, only: iuse_orbeigv,lin_jdav,method
      use optwf_control, only: micro_iter_sr,multiple_adiag,ncore
      use optwf_control, only: no_active,nopt_iter,nparm,nvec,nvecx
      use optwf_control, only: orbitals_ortho,ratio_j,sr_adiag,sr_eps,sr_tau
      use optwf_corsam, only: add_diag
      use optwf_func, only: ifunc_omega,n_omegaf,n_omegat,omega0
      use optwf_handle_wf, only: set_nparms_tot
      use optwf_parms, only: nparmj
      use orbval,  only: ddorb,dorb,nadorb,ndetorb,orb
      use mpi
      use pathak_mod, only: ipathak, eps_max, deps
      use pathak_mod, only: init_pathak
      use parser_read_data, only: header_printing
      use parser_read_data, only: read_basis_num_info_file,read_csf_file
      use parser_read_data, only: read_csfmap_file
      use parser_read_data, only: read_determinants_file
      use parser_read_data, only: read_dmatrix_file,read_efield_file
      use parser_read_data, only: read_eigenvalues_file
      use parser_read_data, only: read_exponents_file,read_forces_file
      use parser_read_data, only: read_gradients_cartesian_file
      use parser_read_data, only: read_gradients_zmatrix_file
      use parser_read_data, only: read_hessian_zmatrix_file
      use parser_read_data, only: read_jasderiv_file,read_jastrow_file
      use parser_read_data, only: read_modify_zmatrix_file
      use parser_read_data, only: read_molecule_file
      use parser_read_data, only: read_multideterminants_file
      use parser_read_data, only: read_optorb_mixvirt_file
      use parser_read_data, only: read_orbitals_file,read_symmetry_file
      use parser_read_data, only: read_zmatrix_connection_file
      use pcm,     only: MCHS
      use pcm_3dgrid, only: PCM_IUNDEFINED,PCM_SHIFT,PCM_UNDEFINED
      use pcm_cntrl, only: ichpol,ipcm,ipcmprt,isurf
      use pcm_fdc, only: qfree,rcolv
      use pcm_grid3d_contrl, only: ipcm_3dgrid
      use pcm_grid3d_param, only: allocate_pcm_grid3d_param,ipcm_nstep3d
      use pcm_grid3d_param, only: pcm_endpt,pcm_origin,pcm_step3d
      use pcm_parms, only: eps_solv,iscov,ncopcm,nscv,nvopcm
      use pcm_unit, only: pcmfile_cavity,pcmfile_chs,pcmfile_chv
      use periodic, only: ngnorm, ngvec
      use periodic_table, only: atom_t,element
      use pot,     only: pot_nn
      use precision_kinds, only: dp
      use properties_mod, only: prop_cc_nuc
      use prp000,  only: iprop,ipropprt,nprop
      use prp003,  only: cc_nuc
      use pseudo,  only: nloc
      use pseudo_mod, only: MPS_QUAD
      use qua,     only: nquad,wq,xq,yq,zq
      use random_mod, only: setrn, jumprn
      use read_bas_num_mod, only: read_bas_num,readps_gauss
      use sa_weights, only: iweight,nweight,weights
      use scale_dist_mod, only: set_scale_dist
      use set_input_data, only: hessian_zmat_define,inputdet,inputforces
      use set_input_data, only: inputjastrow,inputlcao
      use set_input_data, only: modify_zmat_define
      use set_input_data, only: multideterminants_define
      use slater,  only: cdet,coef,ndet,norb
      use sr_mat_n, only: isr_lambda,ortho,sr_lambda
      use sr_mod,  only: i_sr_rescale,izvzb,mconf,mparm
      use system,  only: atomtyp,cent,iwctype,ncent,ncent_tot,nctype
      use system,  only: nctype_tot,ndn,nelec,newghostype,nghostcent,nup
      use system,  only: symbol,znuc
      use vd_mod,  only: dmc_ivd
      use verify_orbitals_mod, only: verify_orbitals
      use vmc_mod, only: mterms,norb_tot
      use vmc_mod, only: nwftypejas,stoj,jtos,nstoj_tot,nstojmax,extraj
      use vmc_mod, only: nwftypeorb,stoo,otos,nstoo_tot,nstoomax,extrao
      use vmc_mod, only: nbjx,stobjx,bjxtoj,bjxtoo,nstoj_tot,nstoo_tot
      use write_orb_loc_mod, only: write_orb_loc
      use zmatrix, only: izmatrix
#if defined(TREXIO_FOUND)
      use trexio_read_data, only: read_trexio_basis_file
      use trexio_read_data, only: read_trexio_determinant_file
      use trexio_read_data, only: read_trexio_ecp_file
      use trexio_read_data, only: read_trexio_molecule_file
      use trexio_read_data, only: read_trexio_orbitals_file
      use trexio_read_data, only: read_trexio_symmetry_file
      use trexio_read_data, only: write_trexio_basis_num_info_file
      use trexio_read_data, only: file_trexio_path, file_trexio_new
      use verify_orbitals_mod, only: verify_orbitals
      use write_orb_loc_mod, only: write_orb_loc
      use zmatrix, only: izmatrix
      use contrl_file,       only: backend
      use trexio            ! trexio library for reading and writing hdf5 files
#endif

#if defined(TREXIO_FOUND) && defined(QMCKL_FOUND)
      use qmckl_data
      use jastrow_qmckl_mod, only: jastrow_init_qmckl
      use orbitals_qmckl_mod, only: init_orbitals_qmckl
#endif

  use, intrinsic :: iso_fortran_env, only : iostat_end

  !! Allocate_periodic
  use periodic,         only: npoly,np_coul, np_jas,cutg, cutg_big, alattice
  use periodic,         only: rlatt, rlatt_inv, n_images, ell
  use ewald_breakup,    only: set_ewald
  use periodic,         only: allocate_periodic
  use ewald_test,       only: allocate_ewald_test, deallocate_ewald_test
  use m_pseudo,         only: allocate_m_pseudo

! CHAMP modules

! in the replacement of preprocess input

! variables from process input

! Note the additions: Ravindra
! Note the additions: Ravindra

! Note the following modules are new additions

!
  implicit none

!--------------------------------------------------------------- Local Variables
  integer, parameter         :: maxa = 100
  logical                    :: doit, debug

  character(len=72)          :: fname, key
  character(len=20)          :: temp
  integer                    :: ierr, ratio, isavebl

  real(dp)                   :: wsum

  type(block_fdf)            :: bfdf
  type(parsed_line), pointer :: pline

  character(len=100)         :: real_format    = '(A, T40, ":: ", T42, F25.16)'
  character(len=100)         :: int_format     = '(A, T40, ":: ", T50, I0)'
  character(len=100)         :: string_format  = '(A, T40, ":: ", T50, A)'
  character(len=100)         :: array_format   = '(A, "(",I0,")", T40, ":: ", T42, F25.16)'

!------------------------------------------------------------------------- BEGIN

  character(:), allocatable  :: optwf, blocking_vmc, blocking_dmc
  character(:), allocatable  :: file_basis
  character(:), allocatable  :: file_molecule
  character(:), allocatable  :: file_determinants
  character(:), allocatable  :: file_symmetry
  character(:), allocatable  :: file_jastrow
  character(:), allocatable  :: file_jastrow_der
  character(:), allocatable  :: file_orbitals
  character(:), allocatable  :: file_pseudo
  character(:), allocatable  :: file_exponents
  character(:), allocatable  :: file_optorb_mixvirt
  character(:), allocatable  :: file_efield
  character(:), allocatable  :: file_zmatrix_connection
  character(:), allocatable  :: file_eigenvalues
  character(:), allocatable  :: file_basis_num_info
  character(:), allocatable  :: file_dmatrix
  character(:), allocatable  :: file_cavity_spheres
  character(:), allocatable  :: file_modify_zmatrix
  character(:), allocatable  :: file_hessian_zmatrix
  character(:), allocatable  :: file_gradients_zmatrix
  character(:), allocatable  :: file_gradients_cartesian
  character(:), allocatable  :: file_multideterminants
  character(:), allocatable  :: file_forces
  character(:), allocatable  :: file_trexio
  character(:), allocatable  :: trex_backend
  character(:), allocatable  :: file_lattice

! from process input subroutine

  character(len=20)          :: fmt
  character(len=32)          :: keyname
  character(len=10)          :: eunit
  character(len=32)          :: cseed
  integer                    :: irn(8), cent_tmp(3), nefpterm, nstates_g
  real(dp), allocatable       :: anorm(:) ! dimensions = nbasis

! local counter variables
  integer                    :: i,j,k, n, iostat
  integer                    :: ic, iwft, istate, imax
  type(atom_t)               :: atoms
  real(dp)                   :: acsfmax,acsfnow
  character(len=2), allocatable   :: unique(:)

  real(dp), parameter        :: zero = 0.d0
  real(dp), parameter        :: one  = 1.d0
  real(dp), parameter        :: two  = 2.d0

#if defined(TREXIO_FOUND) && defined(QMCKL_FOUND)
  integer(qmckl_exit_code)   :: rc
  integer*8                  :: n8
  integer*8                  :: ncheck, ictx
  integer*8                  :: norb_qmckl(qmckl_no_ctx_max)
  integer, allocatable       :: keep(:)
  character*(1024)           :: err_message = ''

  logical                    :: do_nucl_fitcusp
  real(dp), allocatable      :: nucl_fitcusp_radius(:)
  real(dp), parameter        :: a_cusp = 1.74891d0
  real(dp), parameter        :: b_cusp = 0.126057d0
#endif

! Initialize # get the filenames from the commandline arguments
  call fdf_init(file_input, 'parser.log')

  call flaginit_new()
  !! Number of input variables found so far :: 171

! %module general (complete)
  mode        = fdf_get('mode', 'vmc_one_mpi')
  title       = adjustl(fdf_get('title', 'Untitled'))
  pooldir     = fdf_get('pool', './')
  pp_id       = fdf_get('pseudopot', 'none')
  bas_id      = fdf_get('basis', 'none')
  nforce      = fdf_get('nforce', 1)
  MFORCE      = nforce
  nwftype     = fdf_get('nwftype', 1)
  nwftypejas  = fdf_get('nwftypejas', 1)
  nwftypeorb  = fdf_get('nwftypeorb', 1)
  iperiodic   = fdf_get('iperiodic', 0)
  nstates     = fdf_get('nstates', 1)
  ibasis      = fdf_get('ibasis', 1)
  cseed       = fdf_string('seed', "18123437465549275475255231234865")
  ipr         = fdf_get('ipr', -1)
  eunit       = fdf_get('unit', 'Hartrees')
  !hb          = fdf_get('mass', 0.5d0) ! Always 0.5
  scalecoef   = fdf_get('scalecoef',1.0d0)
  i3dgrid     = fdf_get('i3dgrid',0)
  i3dsplorb   = fdf_get('i3dsplorb',0)
  i3dlagorb   = fdf_get('i3dlagorb',0)
  i3ddensity  = fdf_get('i3ddensity',0)
  istrech     = fdf_get('istrech',0)
  alfstr      = fdf_get('alfstr',4.0d0)
  write_walkalize  = fdf_get('write_walkalize', .false.)

  ! trexio
  trex_backend = fdf_get('backend', 'hdf5')
#if defined(TREXIO_FOUND)
  if (trex_backend == "hdf5") backend = TREXIO_HDF5
  if (trex_backend == "text") backend = TREXIO_TEXT
#endif

  ! Ewald module for periodic
  ! npoly, order polynimial split ewald-breakup
  npoly  = fdf_get('npoly', 8)
  ! polynomial order of the cuttoff better even value
  np_coul= fdf_get('np_coul', 2)
  np_jas = fdf_get('np_jas', 3)
  ! cutoffs in reciprocal space
  cutg   = fdf_get('cutg', 1.0d0)
  cutg_big   = fdf_get('cutg_big', 1.d0)
  ! number of images for ao's evaluation in PBC
  n_images  = fdf_get('n_images', 1)
  !alattice = fdf_get('alattice', 1.0d0)

! %module electrons (complete)
  nelec       = fdf_get('nelec', 1)
  nup         = fdf_get('nup', 1)
  ndn         = nelec-nup

! %module atoms (complete)
  newghostype = fdf_get('newghostype', 0)
  nghostcent  = fdf_get('nghostcent', 0)

! %module jastrow (complete)
  ijas        = fdf_get('ijas', 4)
  isc         = fdf_get('isc', 2)
  nspin1      = fdf_get('nspin1', 1)
  nspin2      = fdf_get('nspin2', 1)

  ijas_lr     = fdf_get('ijas_lr', 0)

! %module optgeo (complete)
  iforce_analy= fdf_get('iforce_analy', 0)
  iuse_zmat   = fdf_get('iuse_zmat', 0)
  izvzb       = fdf_get('izvzb', 0)
  if (iforce_analy .gt. 0) then
    alfgeo      = fdf_get('alfgeo', 1.0d0)
  endif
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
  iguiding    = fdf_get('iguiding',0)
  iefficiency = fdf_get('iefficiency',0)

! %module efield (complete)
  iefield     = fdf_get('iefield', 0)

! module vmc (complete)
  imetro      = fdf_get('imetro', 6)
  node_cutoff = fdf_get('node_cutoff', 0)
  eps_node_cutoff = fdf_get('enode_cutoff', 1.0d-7)
  delta       = fdf_get('delta', 1.0d0)
  deltar      = fdf_get('deltar', 5.0d0)
  deltat      = fdf_get('deltat', 1.0d0)
  fbias       = fdf_get('fbias', 1.0d0)
  vmc_tau      = fdf_get('vmc_tau', 0.5d0)

! %module vmc / blocking_vmc (complete)
  vmc_nstep     = fdf_get('vmc_nstep', 1)
  vmc_nblk      = fdf_get('vmc_nblk', 1)
  vmc_nblkeq    = fdf_get('vmc_nblkeq', 2)
  vmc_nblk_max  = fdf_get('nblk_max', vmc_nblk)
  vmc_nconf     = fdf_get('vmc_nconf', 1)
  vmc_nconf_new = fdf_get('vmc_nconf_new', 0)
  vmc_idump     = fdf_get('vmc_idump', 1)
  vmc_irstar    = fdf_get('vmc_irstar', 0)
  vmc_isite     = fdf_get('vmc_isite', 1)
  vmc_icharged_atom     = fdf_get('vmc_icharged_atom', 0)
  vmc_nblk_ci   = fdf_get('vmc_nblk_ci', vmc_nblk)

  kref_fixed    = fdf_get('kref_fixed', 1)

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
  ibranching_c   = fdf_get('ibranching_c', 0.0d0)
  nfrag       = fdf_get('nfrag', 0)
  limit_wt_dmc= fdf_get('limit_wt_dmc', 0)
  dmc_node_cutoff = fdf_get('dmc_node_cutoff', 0)
  dmc_eps_node_cutoff = fdf_get('dmc_enode_cutoff', 1.0d-7)
  tau         = fdf_get('tau', 1.0d0)
  etrial      = fdf_get('etrial', 1.0d0)
  esigmatrial = fdf_get('esigmatrial', 1.0d0)
  nfprod      = fdf_get('nfprod', 100)
  nwprod      = fdf_get('nwprod', 1)
  itausec     = fdf_get('itausec', 1)
  icasula     = fdf_get('icasula', 0)
  dmc_ivd     = fdf_get('dmc_ivd', 0)
  ipathak     = fdf_get('ipathak', 0)
  call init_pathak()
  eps_max     = fdf_get('eps_max', 0.d0)
  deps        = fdf_get('deps', 0.d0)

! %module dmc / blocking_dmc (complete)
  dmc_nstep     = fdf_get('dmc_nstep', 1)
  dmc_nblk      = fdf_get('dmc_nblk', 1)
  dmc_nblkeq    = fdf_get('dmc_nblkeq', 2)
  dmc_nconf     = fdf_get('dmc_nconf', 1)
  dmc_nconf_new = fdf_get('dmc_nconf_new', 0)
  dmc_idump     = fdf_get('dmc_idump', 1)
  dmc_irstar    = fdf_get('dmc_irstar', 0)
  dmc_isite     = fdf_get('dmc_isite', 1)


!optimization flags vmc/dmc
! %module optwf

  ioptwf        = fdf_get('ioptwf', 0)
  method        = fdf_get('method', 'sr_n')
  ioptjas       = fdf_get('ioptjas', 0)
  ioptorb       = fdf_get('ioptorb', 0)
  ioptci        = fdf_get('ioptci', 0)
  nopt_iter     = fdf_get('nopt_iter',6)
  micro_iter_sr = fdf_get('micro_iter_sr', 1)
  isample_cmat  = fdf_get('isample_cmat', 1)
  energy_tol    = fdf_get('energy_tol', 1.d-3)

  if (fdf_defined("optwf")) then
    if ( method .eq. 'linear' ) then
      MFORCE = 3  ! Only set MFORCE here. nwftype=3 is set just before the allocation
    endif
    if ( method .eq. 'sr_n' .or. method .eq. 'lin_dav' .or. method .eq. 'mix_n') then
      isample_cmat=0
      energy_tol=0.d0
    endif
  endif

  dparm_norm_min = fdf_get('dparm_norm_min', 1.0d0)
  ilastvmc      = fdf_get('ilastvmc',1)

  idl_flag      = fdf_get('idl_flag', 0)
  dl_mom        = fdf_get('dl_mom', 0.0d0)
  dl_alg        = fdf_get('dl_alg','nag')

  ilbfgs_flag   = fdf_get('ilbfgs_flag', 0)
  ilbfgs_m      = fdf_get('ilbfgs_m', 5)

  ibeta         = fdf_get('ibeta', -1)
  ratio         = fdf_get('ratio', ratio_j)
  iapprox       = fdf_get('iapprox', 0)
  iuse_orbeigv  = fdf_get('iuse_orbeigv', 0)

  no_active     = fdf_get('no_active', 0)
  ncore         = fdf_get('ncore', 0)
  orbitals_ortho = fdf_get('orbitals_ortho', .false.)

  multiple_adiag = fdf_get('multiple_adiag',0)
! attention needed here.
  if (.not. allocated(add_diag)) allocate (add_diag(MFORCE))
  add_diag(1)   = fdf_get('add_diag',1.d-6)

  ifunc_omega   = fdf_get('func_omega', 0)
  if (ifunc_omega .gt. 0) then
    omega0        = fdf_get('omega', 0.d0)
    n_omegaf      = fdf_get('n_omegaf', nopt_iter)
    n_omegat      = fdf_get('n_omegat', 0)
  endif
  nvec          = fdf_get('lin_nvec', 5)
  nvecx         = fdf_get('lin_nvecx', 160)
  alin_adiag    = fdf_get('lin_adiag', 1.0d-2)
  alin_eps      = fdf_get('lin_eps', 1.0d-3)
  lin_jdav      = fdf_get('lin_jdav',0)

  sr_tau        = fdf_get('sr_tau', 2.0d-2)
  sr_adiag      = fdf_get('sr_adiag', 1.0d-2)
  sr_eps        = fdf_get('sr_eps', 1.0d-3)
  i_sr_rescale  = fdf_get('sr_rescale', 0)

  ngrad_jas_blocks = fdf_get('ngrad_jas_blocks',0)
  isavebl       = fdf_get('save_blocks', 0)
  nefp_blocks   = fdf_get('force_blocks',1)
  iorbsample    = fdf_get('iorbsample',1)

! %module ci (complete)
  iciprt        = fdf_get('iciprt',0)

!%module pcm (complete)

!   ipcm          = fdf_get('ipcm',0)
!   ipcmprt       = fdf_get('ipcmprt',0)
!   pcmfile_cavity = fdf_get('file_cavity','pcm000.dat')
!   pcmfile_chs   = fdf_get('file_chs','chsurf_old')
!   pcmfile_chv   = fdf_get('file_chv','chvol_old')
! !  nscv          = fdf_get('nblk_chv',nblk)   ! attention
! !  iscov         = fdf_get('nstep_chv',nstep2) ! attention
!   eps_solv      = fdf_get('eps_solv',1)
! !  fcol          = fdf_get('fcol',1.d0)
!   rcolv         = fdf_get('rcolv',0.04d0)
!  npmax         = fdf_get('npmax',1)


  ! call allocate_pcm_grid3d_param()
  ! ipcm_3dgrid   = fdf_get('ipcm_3dgrid',0)
  ! ipcm_nstep3d(1)  = fdf_get('nx_pcm',PCM_IUNDEFINED)
  ! ipcm_nstep3d(2)  = fdf_get('ny_pcm',PCM_IUNDEFINED)
  ! ipcm_nstep3d(3)  = fdf_get('nz_pcm',PCM_IUNDEFINED)
  ! pcm_step3d(1)  = fdf_get('dx_pcm',PCM_UNDEFINED)
  ! pcm_step3d(2)  = fdf_get('dy_pcm',PCM_UNDEFINED)
  ! pcm_step3d(3)  = fdf_get('dz_pcm',PCM_UNDEFINED)
  ! pcm_origin(1)  = fdf_get('x0_pcm',PCM_UNDEFINED)
  ! pcm_origin(2)  = fdf_get('y0_pcm',PCM_UNDEFINED)
  ! pcm_origin(3)  = fdf_get('z0_pcm',PCM_UNDEFINED)
  ! pcm_endpt(1)   = fdf_get('xn_pcm',PCM_UNDEFINED)
  ! pcm_endpt(2)   = fdf_get('yn_pcm',PCM_UNDEFINED)
  ! pcm_endpt(3)   = fdf_get('zn_pcm',PCM_UNDEFINED)
  ! PCM_SHIFT      = fdf_get('shift',4.d0)

! %module mmpol (complete)
  ! immpol        = fdf_get('immpol',0)
  ! immpolprt     = fdf_get('immpolprt',0)
  ! mmpolfile_sites = fdf_get('file_sites','mmpol000.dat')
  ! mmpolfile_chmm = fdf_get('file_mmdipo','mmdipo_old')
  ! a_cutoff      = fdf_get('a_cutoff',2.5874d0)
  ! rcolm         = fdf_get('rcolm',0.04d0)

! %module properties (complete)
  iprop         = fdf_get('sample',0)
  ipropprt      = fdf_get('print',0)

! %module pseudo (complete)
  nloc          = fdf_get('nloc',4)  ! for pseudo in Gauss format
  nquad         = fdf_get('nquad',6)
! %module qmmm (complete)
!  iqmm          = fdf_get('iqmm',0)

! attention please. The following line moved here because next_max was not defined yet.
  nadorb        = fdf_get('nextorb', -1)  ! the default should be next_max


  ! Filenames parsing
  file_trexio            = fdf_load_filename('trexio',   'default.hdf5')
  file_basis             = fdf_load_filename('basis',   'default.bas')
  file_molecule          = fdf_load_filename('molecule',   'default.xyz')
  file_determinants      = fdf_load_filename('determinants',  'default.det')
  file_symmetry          = fdf_load_filename('symmetry',   'default.sym')
  file_jastrow           = fdf_load_filename('jastrow',   'default.jas')
  file_jastrow_der       = fdf_load_filename('jastrow_der',   'default.jasder')
  file_orbitals          = fdf_load_filename('orbitals',   'default.orb')
  file_exponents         = fdf_load_filename('exponents',   'exponents.exp')
  file_pseudo       = fdf_load_filename('pseudo',   'default.psp')
  file_optorb_mixvirt       = fdf_load_filename('optorb_mixvirt',  'default.mix')
  file_multideterminants    = fdf_load_filename('multideterminants',    'default.mdet')
  file_forces            = fdf_load_filename('forces',   'default.for')
  file_eigenvalues     = fdf_load_filename('eigenvalues',   'default.eig')
  file_basis_num_info       = fdf_load_filename('basis_num_info',  'default.bni')
  file_dmatrix      = fdf_load_filename('dmatrix',   'default.dmat')
  file_cavity_spheres       = fdf_load_filename('cavity_spheres',  'default.cav')
  file_gradients_zmatrix    = fdf_load_filename('gradients_zmatrix',    'default.gzmat')
  file_gradients_cartesian  = fdf_load_filename('gradients_cartesian',  'default.gcart')
  file_modify_zmatrix       = fdf_load_filename('modify_zmatrix',  'default.mzmat')
  file_hessian_zmatrix      = fdf_load_filename('hessian_zmatrix',  'default.hzmat')
  file_zmatrix_connection   = fdf_load_filename('zmatrix_connection',   'default.zmcon')
  file_efield             = fdf_load_filename('efield',   'default.efield')
  file_lattice              = fdf_load_filename('lattice',              'lattice.txt')

  call header_printing()

! to be moved in a separate subroutine
! printing some information about calculation setup.

  write(ounit,*) '____________________________________________________________________'
  write(ounit,*)
  write(ounit,'(a)') " General input parameters :: "
  write(ounit,*) '____________________________________________________________________'
  write(ounit,*)

  write(ounit,*) " Modules present in the input file :: "
  do i = 1, number_of_modules
    write(ounit,'(a,i0,a,a)') '  (', i, ')', modulenames(i)
  enddo
  write(ounit,*)

  select case (mode)
  case ('vmc_all_mpi')
    write(ounit,'(a,a)') " Calculation mode :: Variational MC all-electron move for ", title
  case ('vmc_one_mpi')
    write(ounit,'(a,a)') " Calculation mode :: Variational MC one-electron move mpi for ",  title
  case ('dmc_all_mpi1')
    write(ounit,'(a,a)') " Calculation mode :: Diffusion MC all-electron move, mpi no global pop for ", title
    call fatal_error('INPUT: This calculations mode not implemented yet')
  case ('dmc_one_mpi1')
    write(ounit,'(a,a)') " Calculation mode :: Diffusion MC one-electron move, mpi no global pop for ", title
  case ('dmc_one_mpi2')
    write(ounit,'(a,a)') " Calculation mode :: Diffusion MC one-electron move, mpi global pop comm for ",  title
  case default
    write(ounit,'(a,a)') " Calculation mode :: ",  title
    call fatal_error('INPUT: This calculations mode not implemented yet')
  end select

  write(ounit,*)
  pooldir = trim(pooldir)
  write(ounit,'(a,a)') " Pool directory for common input files :: ",  pooldir

  write(ounit,*)
  if (len_trim(cseed) < 32) then
    call fatal_error('INPUT: Random number seed is shorter than 32 characters.')
  endif
  if (wid) read(cseed,'(8i4)') irn
  call bcast(irn)

  write(ounit,'(a25,t40,8i5)') " Random number seed root", irn
  call setrn(irn)
  do n=1,idtask ! unique 2^128 non-overlapping sequences for each process
    call jumprn()
  end do

  write(ounit,*)
  write(ounit, string_format) " All energies are in units of ", eunit
  write(ounit, real_format) " hbar**2/(2.*m) ",  hb
  write(ounit, int_format) " Number of geometries ", nforce
  write(ounit, int_format) " Number of wave functions ", nwftype
  if(nwftype.gt.nforce) call fatal_error('INPUT: nwftype gt nforce')
  if(nwftype.gt.MWF) call fatal_error('INPUT: MWF needs to be still modified in module - to fix')
  write(ounit,*)

  call elapsed_time ( "Parsing input file and printing headers : " )
  ! Printing header information and common calculation parameters ends here

! Molecular geometry file in .xyz format [#####]
  write(ounit,*)
  write(ounit,'(a)') " System Information :: Geometry : "
  write(ounit,*) '____________________________________________________________________'
  write(ounit,*)

  write(ounit,*)
  if (fdf_block('molecule', bfdf)) then
    call fdf_read_molecule_block(bfdf)
  elseif ( fdf_load_defined('molecule') ) then
    call read_molecule_file(file_molecule)
  elseif ( fdf_load_defined('trexio') ) then
#if defined(TREXIO_FOUND)
    call read_trexio_molecule_file(file_trexio)
#else
  write(errunit,'(a)') "Error:: Not compiled with TREXIO support but trexio file is present in input file"
  error stop
#endif
  else
    write(errunit,'(a)') "Error:: No information about molecular coordiates provided."
    !write(errunit,'(3a,i6)') "Stats for nerds :: in file ",__FILE__, " at line ", __LINE__
    error stop
  endif

  call elapsed_time ( "Reading molecular coordinates file : " )


! By this point all information about geometry and znuc should be present.
  iznuc     = 1
  igeometry = 1
! Molecular geometry section ends here

! Electrons
  write(ounit,*)
  write(ounit,'(a)') " System Information :: Electrons : "
  write(ounit,*) '____________________________________________________________________'
  write(ounit,*)

  !checks
  if(nup.lt.nelec/2) call fatal_error('INPUT: nelec/2 exceeds nup')

  write(ounit,*)
  write(ounit,int_format) " Number of total electrons = ", nelec
  write(ounit,int_format) " Number of alpha electrons = ", nup
  write(ounit,int_format) " Number of beta  electrons = ", ndn
  write(ounit,*)
! Electrons section ends here


! Pseudopotential section:
! Pseudopotential information (either block or from a file)

  write(ounit,*)
  write(ounit,'(a)') " Calculation Parameters :: Pseudopotential / All Electron : "
  write(ounit,*) '____________________________________________________________________'
  write(ounit,*)

  if ( fdf_defined('pseudopot') ) then
    call readps_gauss()
    ! nquad :: number of quadrature points
    write(ounit,*)
    write(ounit,int_format ) " number of quadrature points (nquad) = ", nquad
#if defined(TREXIO_FOUND)
  elseif ( fdf_load_defined('trexio') .and. nloc .ne. 0) then
    call read_trexio_ecp_file(file_trexio)
    write(ounit,*)
    write(ounit,int_format ) " number of quadrature points (nquad) = ", nquad
#endif
  elseif (nloc .eq. 0) then
    write(ounit,'(a)') "Warning:: Is this an all electron calculation?"
  else
    write(errunit,'(a)') "Error:: No information about pseudo provided in the input."
    !write(errunit,'(3a,i6)') "Stats for nerds :: in file ",__FILE__, " at line ", __LINE__
  endif
  write(ounit,*)

  if (nloc .ne. 0) call elapsed_time ( "Reading ECP files : " )

  ! Pseudopotential section ends here



! More about calculation parameters :: VMC / DMC settings

  if( mode(1:3) == 'vmc') then
    write(ounit,*)
    write(ounit,'(a)') " Calculation Parameters :: VMC : "
    write(ounit,*) '____________________________________________________________________'
    write(ounit,*)

    write(ounit,int_format)  " Version of Metropolis = ", imetro
    if(node_cutoff.gt.0) then
      write(ounit,real_format) " Sampling finite guiding wave function at nodes "
      write(ounit,real_format) " VMC eps node cutoff   = ", eps_node_cutoff
    endif

    if (imetro.eq.1) then
      write(ounit,real_format)  " Drift-diffusion tau = ", vmc_tau
    else
      if(deltar .lt. one) then
        write(ounit,'(a)') '**Warning value of deltar reset to 2.'
        deltar = two
      endif
      if(deltat .lt. zero .or. deltat .gt. two) then
        write(ounit,'(a)') '**Warning value of deltat reset to 2.'
        deltat = two
      endif
      write(ounit,real_format)  " Radial step multiplier  ", deltar
      write(ounit,real_format)  " cos(theta) step size  ",   deltat
    endif

    ! Truncate fbias so that fbias and the sampled quantity are never negative
    fbias=dmin1(two,dmax1(zero,fbias))
    write(ounit,real_format)  " Force bias  ",   fbias

    write(ounit,int_format) " Number of VMC steps/block  ", vmc_nstep
    write(ounit,int_format) " Number of VMC blocks after eq.  ", vmc_nblk
    write(ounit,int_format) " Number of VMC blocks before eq.  ", vmc_nblkeq
    write(ounit,int_format) " Number of VMC configurations saved  ", vmc_nconf_new
  endif

  ! DMC
  if( mode(1:3) == 'dmc' ) then
    write(ounit,*)
    write(ounit,'(a)') " Calculation Parameters :: DMC : "
    write(ounit,*) '____________________________________________________________________'
    write(ounit,*)

    rttau=dsqrt(tau)

    write(ounit,int_format) " Version of DMC ",  idmc
    if( dmc_node_cutoff.gt.0 ) then
      write(ounit,real_format) " Sampling finite guiding wave function at nodes "
      write(ounit,real_format) " DMC eps node cutoff   = ", dmc_eps_node_cutoff
    endif
    write(ounit,int_format) " nfprod ",  nfprod
    write(ounit,real_format) " tau ", tau

    write(ounit,int_format) " ipq ", ipq
    write(ounit,int_format) " itau_eff ", itau_eff
    write(ounit,int_format) " iacc_rej ", iacc_rej
    write(ounit,int_format) " icross ", icross
    write(ounit,int_format) " icuspg ", icuspg
    write(ounit,int_format) " idiv_v ", idiv_v
    write(ounit,int_format) " icut_br ", icut_br
    write(ounit,int_format) " icut_e ", icut_e
    write(ounit,int_format) " ipq ", ipq
    write(ounit,real_format) " etrial ", etrial
    write(ounit,int_format)  " casula ",  icasula
    write(ounit,int_format) " node_cutoff ", dmc_node_cutoff

    if (dmc_node_cutoff.gt.0) write(ounit,real_format) " enode cutoff = ", dmc_eps_node_cutoff

    if (icasula.eq.-1.and.dmc_irstar.eq.1) write(ounit,*) 'Restart is not consistent with icasula = -1'

    if (iabs(idmc).ne.2) call fatal_error('INPUT: only idmc=2 supported')

    if (nloc.eq.0) call fatal_error('INPUT: no all-electron DMC calculations supported')

    if ((iforce_analy.gt.0.and.dmc_ivd.gt.0).or.nforce.gt.1) write(ounit,int_format) " nwprod", nwprod
    if (.not. fdf_defined('etrial')) call fatal_error("etrial required for DMC calculations")

  else
    icasula=0
  endif

  ! Inizialized to zero for call to hpsi in vmc or dmc with no casula or/and in acuest
  i_vpsp=0

  ! Reduce printing in case of a large calculation

  if(vmc_nstep*(vmc_nblk+2*vmc_nblkeq) .gt. 104000) ipr=-1
  if(vmc_irstar.eq.1) vmc_nblkeq=0

  if(dmc_nstep*(dmc_nblk+2*dmc_nblkeq) .gt. 104000) ipr=-1
  if(dmc_irstar.eq.1) dmc_nblkeq=0

  if( mode(1:3) == 'dmc' ) then
    write(ounit,int_format) " Number of DMC steps/block = ", dmc_nstep
    write(ounit,int_format) " Number of DMC blocks after eq. = ", dmc_nblk
    write(ounit,int_format) " Number of DMC blocks before eq. = ", dmc_nblkeq

    write(ounit,int_format) " Target walker population = ", dmc_nconf
    call set_mwalk()  ! to set up maximum allowed number of walkers
    if(dmc_nconf .le. 0) call fatal_error('INPUT: target population <= 0')
    if(dmc_nconf .gt. mwalk) call fatal_error('INPUT: target population > mwalk')
    write(ounit,int_format) " Number of configurations saved = ", dmc_nconf_new
  endif

  call elapsed_time ("Setting VMC/DMC parameters : ")

  ! VMC/DMC calculation parameters settings ends here.


! LCAO orbitals (must be loaded before reading basis )

  write(ounit,*)
  write(ounit,'(a)') " Calculation Parameters :: Orbital Information : "
  write(ounit,*) '____________________________________________________________________'
  write(ounit,*)

  if ( fdf_load_defined('orbitals') ) then
    call read_orbitals_file(file_orbitals)
  elseif ( fdf_load_defined('trexio') ) then
#if defined(TREXIO_FOUND)
    call read_trexio_orbitals_file(file_trexio, .false.)
#endif

  elseif ( fdf_block('orbitals', bfdf)) then
  ! call fdf_read_orbitals_block(bfdf)
    write(errunit,'(a)') "Error:: No information about orbitals provided."
    !write(errunit,'(3a,i6)') "Stats for nerds :: in file ",__FILE__, " at line ", __LINE__
    error stop
  else
    if(nwftype.gt.1) then
      if(iperiodic.eq.0.and.ilcao.ne.nwftype) then
        write(ounit,*) "Warning INPUT: block lcao missing for one wave function"
        write(ounit,*) "Warning INPUT: lcao blocks equal for all wave functions"
        call inputlcao
      endif
    endif
  endif

  call elapsed_time ("Reading molecular coefficients file : ")

! (9) Symmetry information of orbitals (either block or from a file)

  if ( fdf_load_defined('symmetry') ) then
    call read_symmetry_file(file_symmetry)
  elseif ( fdf_load_defined('trexio') ) then
#if defined(TREXIO_FOUND)
    call read_trexio_symmetry_file(file_trexio)
#endif
  elseif ( fdf_block('symmetry', bfdf)) then
  ! call fdf_read_symmetry_block(bfdf)
    write(errunit,'(a)') "Warning:: No information about orbital symmetries provided in the block."
    write(errunit,'(3a,i6)') "Stats for nerds :: in file ",__FILE__, " at line ", __LINE__
!   if( mode(1:3) == 'vmc' ) error stop
  else
    write(errunit,'(a)') "Warning:: No information about orbital symmetries provided in the block."
    write(errunit,'(3a,i6)') "Stats for nerds :: in file ",__FILE__, " at line ", __LINE__
!    if( mode(1:3) == 'vmc' ) error stop
    ! if no symmetry file present, assume same symmetry for all orbitals
    if (.not. allocated(irrep)) allocate (irrep(norb_tot))
    irrep(1:norb_tot) = 1
    write(ounit,*)
    write(ounit,*) '____________________________________________________________________'
    write(ounit, *) " Orbital symmetries are set to default "
    write(ounit, '(10(1x, i3))') (irrep(i), i=1, norb_tot)
    write(ounit,*) '____________________________________________________________________'
    write(ounit,*)
  endif

  call elapsed_time ("Reading/setting symmetry file : ")

! Basis num information (either block or from a file)

  write(ounit,*)
  write(ounit,'(a)') " Calculation Parameters :: Numerical Basis Information : "
  write(ounit,*) '____________________________________________________________________'
  write(ounit,*)

  if ( fdf_load_defined('basis_num_info') ) then
    call read_basis_num_info_file(file_basis_num_info)
  elseif ( fdf_load_defined('trexio') ) then
#if defined(TREXIO_FOUND)
    call write_trexio_basis_num_info_file(file_trexio)
#endif
  elseif (.not. fdf_block('basis_num_info', bfdf)) then
  ! call fdf_read_eigenvalues_block(bfdf)
    write(errunit,'(a)') "Error :: No information about basis num info provided in the block."
    !write(errunit,'(3a,i6)') "Stats for nerds :: in file ",__FILE__, " at line ", __LINE__
    error stop
  else
    write(errunit,'(a)') "Error :: No information about basis num info provided in the block."
    !write(errunit,'(3a,i6)') "Stats for nerds :: in file ",__FILE__, " at line ", __LINE__
    error stop
  endif

  call elapsed_time ("Reading basis information file : ")

  write(ounit,*)
  write(ounit,'(a)') " Calculation Parameters :: Jastrow and Jastrow Derivatives : "
  write(ounit,*) '____________________________________________________________________'
  write(ounit,*)

  if(ijas.ne.1.and.ijas.ne.4) call fatal_error('Only ijas 1 and 4 implemented')
  if(iperiodic.eq.0.and.ijas_lr.gt.0) call fatal_error('No long-range Jastrow for non-periodic system')

! Jastrow Parameters (either block or from a file)

  if ( fdf_load_defined('jastrow') ) then
    call read_jastrow_file(file_jastrow)
  elseif (fdf_block('jastrow', bfdf)) then
  !call fdf_read_jastrow_block(bfdf)
    write(errunit,'(a)') "Error:: No information about jastrow provided in the block."
    !write(errunit,'(3a,i6)') "Stats for nerds :: in file ",__FILE__, " at line ", __LINE__
    error stop
  else
    ! no information about jastrow present.
    if(nwftype.gt.1) then
      if(ijastrow_parameter.ne.nwftype) then
        write(ounit,*) "INPUT: block jastrow_parameter missing for one wave function"
        write(ounit,*) "INPUT: jastrow_parameter blocks equal for all wave functions"
        call inputjastrow
      endif
    endif
  endif

  call elapsed_time ("Reading Jastrow file : ")

! Jastrow derivative Parameters (either block or from a file)

  if ( fdf_load_defined('jastrow_der') ) then
    call read_jasderiv_file(file_jastrow_der)
  elseif ( fdf_block('jastrow_der', bfdf)) then
  !call fdf_read_jastrow_derivative_block(bfdf)
    write(errunit,'(a)') "Error:: No information about jastrow derivatives provided in the block."
    !write(errunit,'(3a,i6)') "Stats for nerds :: in file ",__FILE__, " at line ", __LINE__
    if( mode(1:3) == 'vmc' ) error stop
  elseif (ioptjas .ne. 0) then
    write(errunit,'(a)') "Error:: No information about jastrow derivatives provided in the block."
    !write(errunit,'(3a,i6)') "Stats for nerds :: in file ",__FILE__, " at line ", __LINE__
    if( mode(1:3) == 'vmc' ) error stop
  endif

  call elapsed_time ("Reading Jastrow derivative file : ")

!printing some information and warnings and checks about Jastrow
  write(ounit, * )

  write(ounit, int_format ) " ijas = ", ijas
  write(ounit, int_format ) " isc = ",  isc
  write(ounit, int_format ) " nspin1 = ", nspin1
  write(ounit, int_format ) " nspin2 = ", nspin2

  if(ijas.ne.1.and.iperiodic.gt.0) &
    call fatal_error('Only ijas=1 implemented for periodic systems')

  if(ijas.eq.4) then
    asymp_r=0

    write(ounit,'(a)') " transferable standard form 4"

    if(isc.eq.2) write(ounit,'(a)') " dist scaled r=(1-exp(-scalek*r))/scalek"
    if(isc.eq.3) write(ounit,'(a)') " dist scaled r=(1-exp(-scalek*r-(scalek*r)**2/2))/scalek"
    if(isc.eq.4) write(ounit,'(a)') " dist scaled r=r/(1+scalek*r)"
    if(isc.eq.5) write(ounit,'(a)') " dist scaled r=r/(1+(scalek*r)**2)**.5"

    call allocate_jasasymp()  ! Needed for the following two arrays
    do j=1,nwftypejas
      call set_scale_dist(j,ipr)
    enddo
  endif
  call elapsed_time ("Setting Jastrow parameters : ")

! Determinants (only)

  write(ounit,*)
  write(ounit,'(a)') " Calculation Parameters :: Determinants : "
  write(ounit,*) '____________________________________________________________________'
  write(ounit,*)

  if ( fdf_load_defined('determinants') ) then
    call read_determinants_file(file_determinants)
    if (ioptci .ne. 0) mxciterm = ndet
  elseif ( fdf_block('determinants', bfdf)) then
    if (ioptci .ne. 0) mxciterm = ndet
  ! call fdf_read_determinants_block(bfdf)
#if defined(TREXIO_FOUND)
  elseif ( fdf_load_defined('trexio') ) then
    call read_trexio_determinant_file(file_trexio)
    if (ioptci .ne. 0) mxciterm = ndet
#endif
  elseif(nwftype.gt.1) then
      if(ideterminants.ne.nwftype) then
        write(ounit,*) "Warning INPUT: block determinants missing for one wave function"
        write(ounit,*) "Warning INPUT: determinants blocks equal for all wave functions"
        call inputdet
      endif
  else
    write(errunit,'(a)') "Error:: No information about determinants provided."
    write(errunit,'(3a,i6)') "Stats for nerds :: in file ",__FILE__, " at line ", __LINE__
    error stop
  endif

  ! allocation after determinants and basis
  if (fdf_defined("optwf")) then
    if ( method .eq. 'linear' ) then
      nwftype = 3
      nforce = 3
    endif
  endif

  call elapsed_time ("Reading determinants only from a file : ")

! (3) CSF only

  if ( fdf_load_defined('determinants') .and. ndet .gt. 1 ) then
    call read_csf_file(file_determinants)
    if (ioptci .ne. 0 .and. ncsf .gt. 0 ) mxciterm = ncsf
  elseif (fdf_block('csf', bfdf)) then
    call fdf_read_csf_block(bfdf)
    if (ioptci .ne. 0) mxciterm = ncsf
  else
    ! No csf present; set default values; This replaces inputcsf
    nstates = 1
    ncsf = 0
    if (ioptci .ne. 0 .and. ici_def .eq. 1) nciterm = nciprim

    if( (method(1:3) == 'lin')) then
      if (.not. allocated(ccsf)) allocate(ccsf(ndet, nstates, 3))
    else
      if (.not. allocated(ccsf)) allocate(ccsf(ndet, nstates, nwftype))
    endif
    do j = 1, ndet
      ccsf(j,1,1) = cdet(j,1,1)
    enddo
  endif

  if(.not. allocated(maxcsf)) allocate(maxcsf(nstates))

  if(ncsf.gt.0 .and. method(1:3) .ne. 'lin') then
      do istate=1,nstates
        acsfmax=dabs(ccsf(1,istate,1))
        maxcsf(istate)=1
        do j=2,ncsf
          acsfnow=dabs(ccsf(j,istate,1))
          if(acsfnow.gt.acsfmax) then
            acsfmax=acsfnow
            maxcsf(istate)=j
          endif
        enddo
        write(ounit,'(''Saving max scf:state,csf,value'',2i5,f20.15)') &
                          istate,maxcsf(istate),acsfmax
      enddo
    else
      maxcsf(1:nstates)=1
   endif

! (4) CSFMAP [#####]

  if ( fdf_load_defined('determinants') .and. ndet .gt. 1 ) then
    call read_csfmap_file(file_determinants)
  elseif (fdf_block('determinants', bfdf)) then
  ! call fdf_read_csfmap_block(bfdf)
    write(errunit,'(a)') "Error:: No information about csfmaps provided."
    !write(errunit,'(3a,i6)') "Stats for nerds :: in file ",__FILE__, " at line ", __LINE__
    error stop
  else
   nmap = ndet
   if (.not. allocated(cxdet)) allocate (cxdet(nmap))
   if (.not. allocated(iadet)) allocate (iadet(ndet))
   if (.not. allocated(ibdet)) allocate (ibdet(ndet))
   if (.not. allocated(icxdet)) allocate (icxdet(nmap))

   do i = 1, ndet
     iadet(i) = i
     ibdet(i) = i
     icxdet(i) = i
     cxdet(i) = 1.0d0
   enddo

   write(ounit,*) " Determinant - CSF has one-to-one mapping  "
  endif

  call elapsed_time ("Reading CSF and CSFMAP file : ")


  if (mode(1:3) == 'vmc') write(ounit, *) "nstoj_tot, nstoo_tot, nstates", nstoj_tot, nstoo_tot, nstates
  if (extraj.eq.1.and.nstoj_tot.ne.nstates) &
          call fatal_error('Some states have not been assigned a jastrow type')
  if (extrao.eq.1.and.nstoo_tot.ne.nstates) &
          call fatal_error('Some states have not been assigned an orbital set')
  if (method.eq.'sr_n'.and.extraj.eq.0.and.ioptwf.ge.1.and.nstates.gt.1) &
          call fatal_error('For multistate sr_n optimization must use jastrows_to_states in jastrow input')
  if (method.eq.'sr_n'.and.extrao.eq.0.and.ioptwf.ge.1.and.nstates.gt.1) &
          call fatal_error('For multistate sr_n optimization must use orbitals_to_states in lcao input')

  write(ounit,'(A)') 'Mapping of jastrow, orbital, and mixed quantities to each state.'
  allocate(stoj(nstates))
  allocate(stoo(nstates))
  do istate=1,nstates
    stoj(istate) = 0
    stoo(istate) = 0
  enddo

! Jastrow mapping
  if (extraj.eq.0) then
    do istate=1,nstates
      stoj(istate)=1
      if (mode(1:3) == 'vmc') write(ounit,'(A)') "State  -->  Jastrow #"
      if (mode(1:3) == 'vmc') write(ounit,'(i4,A,i4)') istate, '   -->', stoj(istate)
    enddo
  else
    do i=1,nwftypejas
      do j=1,nstojmax
        !write(ounit,*) "nstojmax", nstojmax
        !write(ounit,*) "i,j,jtos(i,j)", i, j, jtos(i,j)
        istate=jtos(i,j)
        if (istate.ne.0) stoj(istate)=i
        if (istate.ne.0) then
          if (mode(1:3) == 'vmc') write(ounit,'(A)') "State  -->  Jastrow #"
          if (mode(1:3) == 'vmc') write(ounit,'(i4,A,i4)') istate, '   -->', stoj(istate)
        endif
      enddo
    enddo
  endif

  do istate=1,nstates
    if (stoj(istate) .eq. 0) then
      if (mode(1:3) == 'vmc') write(ounit,'(A,i4,A)') " State ", istate, " has not been assigned a jastrow type. "
      call fatal_error('JASTROW INPUT: a state has not been assigned a jastrow type.')
    endif
  enddo



! Orbital mapping
  if (extrao.eq.0 .and. .not. fdf_defined("trexio") ) then
    do istate=1,nstates
      stoo(istate)=1
      if (mode(1:3) == 'vmc') write(ounit,'(A)') "State  -->  Orbital set"
      if (mode(1:3) == 'vmc') write(ounit,'(i4,A,i4)') istate, '   -->', stoo(istate)
    enddo
  else
    do i=1,nwftypeorb
      do j=1,nstoomax
        !write(ounit,*) "nstoomax", nstoomax
        !write(ounit,*) "i,j,otos(i,j)", i, j, otos(i,j)
        istate=otos(i,j)
        if (istate.ne.0) stoo(istate)=i
        if (istate.ne.0) then
          if (mode(1:3) == 'vmc') write(ounit,'(A)') "State  -->  Orbital set "
          if (mode(1:3) == 'vmc') write(ounit,'(i4,A,i4)') istate, '   -->', stoo(istate)
       endif
      enddo
    enddo
  endif
  call bcast(stoo)
  call bcast(stoj)

  do istate=1,nstates
    if (stoo(istate) .eq. 0) then
      if (mode(1:3) == 'vmc') write(ounit,'(A,i4,A)') " State ", istate, " has not been assigned an orbital set . "
      call fatal_error('LCAO INPUT: a state has not been assigned an orbital set.')
    endif
  enddo

  allocate(stobjx(nstates))
  if (nwftypejas.eq.1.and.nwftypeorb.eq.1) then
    nbjx = 1
    allocate(bjxtoo(1))
    allocate(bjxtoj(1))
    bjxtoo(1)=1
    bjxtoj(1)=1
    do istate=1,nstates
      stobjx(istate) = 1
    enddo
  else
    nbjx = 1
    stobjx(1)=1
    do istate=2,nstates
      j=0
      do k=1,istate-1
        if (stoo(istate).eq.stoo(k).and.stoj(istate).eq.stoj(k)) then
          stobjx(istate)=k
          j=1
          exit
        endif
      enddo
      if (j.eq.0) then
        stobjx(istate)=istate
        nbjx = nbjx + 1
      endif
    enddo
    allocate(bjxtoo(nbjx))
    allocate(bjxtoj(nbjx))
    do istate=1,nstates
      if (istate.eq.1) then
        imax=stobjx(istate)
        bjxtoo(1)=stoo(1)
        bjxtoj(1)=stoj(1)
      else
        if (stobjx(istate).gt.imax) then
          bjxtoo(stobjx(istate))=stoo(istate)
          bjxtoj(stobjx(istate))=stoj(istate)
        endif
        imax=max(imax,stobjx(istate))
      endif
    enddo
  endif

  do istate=1,nstates
    if (mode(1:3) == 'vmc') write(ounit,'(A)') "State  -->  Mixed Quantity #  <--  Jastrow #, Orbital set"
    if (mode(1:3) == 'vmc') write(ounit,'(i4,A,i4,A,2i4)') istate, '   -->', stobjx(istate), '   <--', bjxtoj(stobjx(istate)), bjxtoo(stobjx(istate))
  enddo

  ! Know the number of orbitals for optimization.
  if (ioptorb .ne. 0) call get_norbterm()

  ! Jastrow mterms needed for allocations
  mterms = nterms4(nordc)
  ! Add up all the parameters. It will be used to allocate arrays.
  nciterm = mxciterm
  call set_nparms_tot()
  ! Set maximum number of parameters. For multistate orbital optimization
  ! the following additional terms will be present. The last +1 is failsafe mechanism.
  if (method.eq.'sr_n') then
    mparm = nparm*nstates + 2 !necessary because storing 2 quantities in sr_o, and atimesn increments ddot by mparm.
  else
    mparm = nparm + (nstates-1)*(norbterm) + 1
  endif

  if(vmc_nblk.gt.vmc_nblk_max) then
    write(ounit,*) "Warning, vmc_nblk gt vmc_nblk_max. Setting vmc_nblk_max = vmc_nblk"
    vmc_nblk_max=vmc_nblk
    ! or we force a fatal error.
  endif

  ! allocate ewald module and initialize the module
  if (iperiodic.gt.0) then
     call allocate_periodic()
     call allocate_ewald_test()
  ! for periodic calculations
     if ( fdf_load_defined('lattice') ) then
        call read_lattice_file(file_lattice)
     endif
     call set_ewald
     call deallocate_ewald_test()
  endif

! Additional Properties
! properties will be sampled iprop
! properties will be printed ipropprt
  nprop=1
  if(iprop.ne.0) then
     if (iperiodic.gt.0) then
        !        nprop=5+ngnorm
        nprop=6+(ngvec-1)+2*(ngvec-1)
     else
        nprop=5
     endif

     write(ounit,'(a)' ) " Properties will be sampled "
     write(ounit,*) " NPROP ", nprop
     write(ounit,int_format ) " Properties printout flag = ", ipropprt
!    call prop_cc_nuc(znuc,cent,iwctype,nctype_tot,ncent_tot,ncent,cc_nuc)
  endif
  
  
  call compute_mat_size_new()
  call allocate_vmc()
  call allocate_dmc()

  ! read fragment indeces
  if ( (nfrag.gt.1).and.(ndet.gt.1) ) call fatal_error('READ_INPUT: Fragments not implemented for multideterminant wavefunctions')
  if (nfrag.eq.1) call fatal_error("READ_INPUT: nfrag=1, do not use nfrag if you dont have any fragments")
  if(nfrag.gt.0) then
    write(ounit, *) ""

    if ( fdf_islist('ifragcent') .and. fdf_islinteger('ifragcent') ) then
      i = -1
      call fdf_list('ifragcent',i,ifragcent)
      write(ounit,'(a)' )
      write(ounit,'(tr1,a,i0,a)') ' ifragcent has ',i,' entries'
      if(i.ne.ncent) call fatal_error('READ_INPUT: ifragcent array must contain ncent entries.')
      ! The ifragcent array is already allocated in allocate_dmc so this does not have to be done again.
      call fdf_list('ifragcent',i,ifragcent)
      write(temp, '(a,i0,a)') '(a,', ncent, '(i4))'
      !write(ounit, '(a, <ncent>i)') 'Fragmentation indices : ', (ifragcent(i), i=1,ncent)
      write(ounit, temp) 'Fragmentation indices : ', (ifragcent(i), i=1,ncent) ! GNU version

      if (maxval(ifragcent) .ne. nfrag) call fatal_error("READ_INPUT: The number of fragments (nfrag) must match the maximal fragment index in ifragcent")
    else
      call fatal_error("READ_INPUT: ifragcent must be defined if nfrag is used")
    endif
  endif
  ! read ibranching_cfrag
  if(nfrag.gt.1.and.(icut_e.eq.-5.or.icut_e.eq.-6)) then
    write(ounit, *) ""

    if ( fdf_islist('ibranching_cfrag') .and. fdf_islreal('ibranching_cfrag') ) then
      i = -1
      call fdf_list('ibranching_cfrag',i,ibranching_cfrag)
      write(ounit,'(a)' )
      write(ounit,'(tr1,a,i0,a)') ' ibranching_cfrag has ',i,' entries'
      if(i.ne.nfrag) call fatal_error('READ_INPUT: ibranching_cfrag array must contain nfrag entries.')
      ! The ifragcent array is already allocated in allocate_dmc so this does not have to be done again.
      call fdf_list('ibranching_cfrag',i,ibranching_cfrag)
      write(temp, '(a,i0,a)') '(a,', nfrag, '(f12.6))'
      !write(ounit, '(a, <nfrag>f12.6)') 'branching c: ', (ibranching_cfrag(i), i=1,nfrag)
      write(ounit, temp) 'branching c: ', (ibranching_cfrag(i), i=1,nfrag) ! GNU version
    else
      call fatal_error("READ_INPUT: chosen icut_e requires ibranching_cfrag")
    endif
  endif
  ! read etrialfrag
  if(nfrag.gt.1.and.(icut_e.eq.-5.or.icut_e.eq.-6)) then
    write(ounit, *) ""

    if ( fdf_islist('etrialfrag') .and. fdf_islreal('etrialfrag') ) then
      i = -1
      call fdf_list('etrialfrag',i,etrialfrag)
      write(ounit,'(a)' )
      write(ounit,'(tr1,a,i0,a)') ' etrialfrag has ',i,' entries'
      if(i.ne.nfrag) call fatal_error('READ_INPUT: etrialfrag array must contain nfrag entries.')
      ! The etrialfrag array is already allocated in allocate_dmc so this does not have to be done again.
      call fdf_list('etrialfrag',i,etrialfrag)
      !write(ounit, '(a, <nfrag>f12.6)') 'etrialfrag: ', (etrialfrag(i), i=1,nfrag)
      write(ounit, temp) 'etrialfrag: ', (etrialfrag(i), i=1,nfrag) ! GNU version
    else
      call fatal_error("READ_INPUT: chosen icut_e requires etrialfrag")
    endif
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
    if(imultideterminants.eq.0) then
      write(errunit,*) "INPUT: multideterminant bloc MISSING"
      call multideterminants_define(0)
    endif
  endif
  imultideterminants = 1


! (7) exponents

  if ( fdf_load_defined('exponents') ) then
    call read_exponents_file(file_exponents)
  elseif ( fdf_block('exponents', bfdf)) then
  ! call fdf_read_exponents_block(bfdf)
    write(errunit,'(a)') "Error:: No information about exponents provided."
    write(errunit,'(3a,i6)') "Stats for nerds :: in file ",__FILE__, " at line ", __LINE__
    error stop
  else
    if( (method(1:3) == 'lin')) then
      if (.not. allocated(zex)) allocate (zex(nbasis, 3))
    else
      if (.not. allocated(zex)) allocate (zex(nbasis, nwftype))
    endif
    zex = 1   ! debug check condition about numr == 0
  endif
  iexponents = iexponents + 1



! (11) Eigenvalues information of orbitals (either block or from a file)

  if ( fdf_load_defined('eigenvalues') ) then
    call read_eigenvalues_file(file_eigenvalues)
  elseif ( fdf_block('eigenvalues', bfdf)) then
  ! call fdf_read_eigenvalues_block(bfdf)
    write(errunit,'(a)') "Error:: No information about eigenvalues provided in the block."
    write(errunit,'(3a,i6)') "Stats for nerds :: in file ",__FILE__, " at line ", __LINE__
    error stop
  else
    write(errunit,'(a)') "Warning:: No information about eigenvalues provided in the block."
    write(errunit,'(3a,i6)') "Stats for nerds :: in file ",__FILE__, " at line ", __LINE__
!    error stop
  endif

  call elapsed_time ("Major allocations for VMC/DMC : ")

! basis information (either block or from a file)

  write(ounit,*)
  write(ounit,'(a)') " Calculation Parameters :: Basis : "
  write(ounit,*) '____________________________________________________________________'
  write(ounit,*)


  if(numr.gt.0)   write(ounit,'(a)') " Numerical basis used"

  if(ibasis.eq.1) then
    write(ounit,'(a)') " Orbitals on localized basis "
    write(ounit, int_format) " Total no. of basis = ", nbasis
    ! call write_orb_loc

    if ( fdf_defined('basis') ) then
      if(numr.gt.0) then
        do iwft=1,nwftype
          call read_bas_num(iwft)
        enddo
        ! See if this is really allocated at this point
        if (.not. allocated(ibas0)) allocate (ibas0(ncent_tot))
        if (.not. allocated(ibas1)) allocate (ibas1(ncent_tot))
        ibas0(1)=1
        ibas1(1)=nbastyp(iwctype(1))
        do ic=2,ncent
          ibas0(ic)=ibas1(ic-1)+1
          ibas1(ic)=ibas1(ic-1)+nbastyp(iwctype(ic))
        enddo
      endif
    ! elseif (fdf_block('basis', bfdf)) then
    !   call fdf_read_basis_block(bfdf)
#if defined(TREXIO_FOUND)
    elseif ( fdf_load_defined('trexio') ) then
      call read_trexio_basis_file(file_trexio)
      ! See if this is really allocated at this point
     if (.not. allocated(ibas0)) allocate (ibas0(ncent_tot))
     if (.not. allocated(ibas1)) allocate (ibas1(ncent_tot))
       ibas0(1)=1
       ibas1(1)=nbastyp(iwctype(1))
       do ic=2,ncent
         ibas0(ic)=ibas1(ic-1)+1
         ibas1(ic)=ibas1(ic-1)+nbastyp(iwctype(ic))
       enddo
#endif
    else
      write(errunit,'(a)') "Error:: No information about basis provided in the block."
      write(errunit,'(3a,i6)') "Stats for nerds :: in file ",__FILE__, " at line ", __LINE__
      error stop
    endif
 endif

  call elapsed_time ("Reading basis file : ")

  ! Basis information section ends here

! check if the orbitals coefficients are to be multiplied by a constant parameter
  if(scalecoef.ne.1.0d0) then
    if((method.eq.'sr_n'.and.nwftypeorb.gt.1)) then
      k = nwftype
      nwftype = nstates
    endif
    do  iwft=1,nwftype
      do  i=1,norb+nadorb
        do  j=1,nbasis
            coef(j,i,iwft)=coef(j,i,iwft)*scalecoef
        enddo
      enddo
    enddo
    write(ounit, real_format) " Orbital coefficients scaled by a constant parameter = ",  scalecoef
    write(ounit,*)
    if((method.eq.'sr_n'.and.nwftypeorb.gt.1)) then
      nwftype = k
    endif
  endif

! verify number of orbitals and setup optorb
! verification already handeled in read_data file.



!! Grid information

  if((i3dsplorb.ge.1).or.(i3dlagorb.ge.1).or.(i3ddensity.ge.1))  i3dgrid=1

  if(i3dgrid.ge.1) then

    ! Grid setup:
            call setup_grid

            if(i3dlagorb.ge.1) then
              write(ounit, '(a)') " Orbitals on a grid: splines interpolation"
              call setup_3dlagorb
              i3dsplorb=0
             elseif(i3dsplorb.ge.1) then
              write(ounit, '(a)') " Orbitals on a grid: Lagrange interpolation"
              call setup_3dsplorb
              i3dlagorb=0
            endif
  endif

  ! Analytical forces flags (vmc only)
  if( mode(1:3) == 'vmc' ) then
    if(iforce_analy.gt.0) then
      write(ounit,'(a)' ) " Geometry optimization with analytic gradients"
      if(iuse_zmat.gt.0) write(ounit,'(a)' ) " Using internal coordinates "
      write(ounit,'(a,t36,f12.6)') " starting alfgeo = ", alfgeo
    endif
  endif

! Contents from flagcheck. Moved here because ndet and norb should be defined by now
! CF: norb is set to max # occupied orb in WF in optorb_define/verify_orbitals
  if(ioptorb.eq.0) then
    call verify_orbitals()
   else
    if(ioptorb_mixvirt.eq.0) then
      norbopt=0
      norbvirt=0
    endif
    if(ioptorb_def.eq.0) then
      write(ounit,*) "INPUT: definition of orbital variations missing"
      call optorb_define(0)
    endif
  endif
  if(ioptci.ne.0.and.ici_def.eq.0) then
    write(ounit,*) "INPUT: definition of OPTCI operators missing"
    call optci_define
  endif


! Optimization flags WF (vmc/dmc only)
  if( (mode(1:3) == 'vmc') .or. (mode(1:3) == 'dmc') ) then

    write(ounit,*)
    write(ounit,'(a)') " Calculation Parameters :: optimizations : "
    write(ounit,*) '____________________________________________________________________'
    write(ounit,*)

    if(ioptwf.gt.0) then
      write(ounit,'(a)' ) " Perform wave function optimization in vmc/dmc"
      write(ounit,'(a,a)' ) " Computing/writing quantities for optimization with method = ", method
      if(nstates.gt.1 .and. ioptwf.gt.0 .and. method.eq.'sr_n') then
      !    call fatal_error('READ_INPUT: nstates>1 and sr_n')
           write(ounit,'(a)' ) " Performing multi-state sr_n, setting ortho=1"
           ortho=1
      endif
    elseif((ioptjas.eq.1) .or. (ioptorb.eq.1) .or. (ioptci.eq.1) ) then
      write(ounit,'(a)' ) " Only sample derivatives of wave function for external use"
      write(ounit,'(a,a)' ) " Computing/writing quantities for optimization with method = ", method
    endif

    if(ioptwf.gt.0.or.ioptjas+ioptorb+ioptci.ne.0) then
      if(method.eq.'lin_d' .or. method.eq.'mix_n') then
        if(lin_jdav.eq.0) then
                write(ounit,'(a)' ) " Use old Regterg"
        elseif(lin_jdav.eq.1) then
                write(ounit,'(a)' ) " Use new Davidson"
        else
                write(ounit,'(a)' ) " Use new Jacobi-Davidson"
        endif
      endif
      if(method.eq.'linear' .and. mxreduced.ne.norbterm )  &
          call fatal_error('READ_INPUT: mxreduced .ne. norbterm')
    endif

! Optimization flag Jastrow  (vmc/dmc only)
    write(ounit,*)
    if(ioptjas.gt.0) then
      write(ounit,'(a)' ) " Jastrow derivatives are sampled "
      write(ounit,int_format) " Number of Jastrow derivatives ",  nparmj
      call cuspinit4(1)
    else
      nparmj=0
    endif

! ORB optimization flags (vmc/dmc only)
    if(ioptorb.ne.0) then
      write(ounit,'(a)' ) " Orbital derivatives are sampled"
      write(ounit,int_format)  " ORB-PT blocks in force average = ", nefp_blocks
      if(isavebl.ne.0)then
        write(ounit,'(a)' ) " ORB-PT block averages will be saved "
        if (wid) then
          idump_blockav=43
          open(unit=idump_blockav,file='efpci_blockav.dat',status='unknown',form='unformatted')
          write(idump_blockav) nefpterm
        endif
      endif
    endif

! CI optimization flags
    if(ioptci.ne.0) then
      write(ounit,'(a)' ) " CI is sampled "
      write(ounit,int_format)  " CI printout flag = ", iciprt
      if( (ioptjas.eq.0) .and. (ioptorb.eq.0) .and. (method.eq.'hessian') ) then
        method='linear'
        write(ounit,'(a)' ) " Reset optimization method to linear"
      endif
      if(ncsf.gt.0) then
        nciterm=ncsf
      else
  ! TMP due to changing kref -> also for ncsf=0, we need to have cxdet(i) carrying the phase
  !     if(kref_fix.eq.0) call fatal_error('ncsf.eq.0 - further changes needed due to kref')
        nciterm=nciprim
      endif
    else
      nciprim=0
      nciterm=0
    endif

    ! write(ounit,int_format)  " CI number of coefficients ", nciterm
    ! write(ounit,int_format)  " nciprim ", nciprim
    mxciterm = nciprim  ! validate this change debug ravindra

    if((ncsf.eq.0) .and. (nciprim.gt.mxciterm) ) call fatal_error('INPUT: nciprim gt mxciterm')
    if(nciterm.gt.mxciterm) call fatal_error('INPUT: nciterm gt mxciterm')


! Multiple states/efficiency/guiding flags
    ! Use guiding wave function constructed from mstates
    if(iguiding.gt.0) then
      if(node_cutoff.gt.0)  call fatal_error('INPUT: guiding wave function AND node_cutoff > 0')

      write(ounit, *) "Guiding function: square root of sum of squares"

      ! Part which handles the guiding weights
      if (.not. allocated(weights_g)) allocate (weights_g(MSTATES))
      if (.not. allocated(iweight_g)) allocate (iweight_g(MSTATES))

      if ( fdf_islreal('weights_guiding') .and. fdf_islist('weights_guiding') &
          .and. (.not. fdf_islinteger('weights_guiding')) ) then
        i = -1
        call fdf_list('weights_guiding',i,weights_g)
        write(ounit,'(a)' )
        write(ounit,'(tr1,a,i0,a)') ' Guiding weights has ',i,' entries'
        call fdf_list('weights_guiding',i,weights_g)
        write(temp, '(a,i0,a)') '(a,', MSTATES, '(f12.6))'
        !write(ounit, '(a,<MSTATES>(f12.6))') ' Weights_guiding : ', weights_g(1:i) ! Intel version
        write(ounit, temp) ' Weights_guiding : ', weights_g(1:i)                   ! GNU version
      else
        weights_g = 1.0d0
        write(ounit,'(a,t40, 10f12.6)') 'Default weights_guiding ', weights_g(1:nstates)
      end if

      wsum = 0.d0
      nweight = 0
      do i = 1, nstates
        if (weights_g(i) .gt. 1d-6) then
            nweight = nweight + 1
            iweight_g(nweight) = i
            wsum = wsum + weights_g(i)
        endif
      enddo

      do i = 1, nweight
        weights_g(i) = weights_g(i)/wsum
      enddo

      if (nweight .eq. 0) then
          nweight = 1
          iweight_g(1) = 1
          weights_g(1) = 1.d0
      endif
    ! The above part should be moved to get_weights subroutine
    endif

    ! Efficiency for sampling states inputed in multiple_cistates
    if(iefficiency.gt.0) nstates_psig=nstates
  else
    ioptorb=0
    ioptci=0
    ioptjas=0
    iefficiency=0
    iguiding=0
    nstates=1
  endif ! if loop of condition of either vmc/dmc ends here
! Read in anormo if multi-state sr_n
  if(iguiding.gt.0) then
    write(ounit, *) "ANORMO: Determining normalization constants for guiding wave function."

    if (.not. allocated(anormo)) allocate (anormo(MSTATES))

    if ( fdf_islreal('anorm') .and. fdf_islist('anorm') &
        .and. (.not. fdf_islinteger('anorm')) ) then
      i = -1
      call fdf_list('anorm',i,anormo)
      write(ounit,'(a)' )
      write(ounit,'(tr1,a,i0,a)') ' anorm has ',i,' entries'
      if(i.ne.nstates) call fatal_error('READ_INPUT: anorm array must contain nstate entries.')
      call fdf_list('anorm',i,anormo)
      write(temp, '(a,i0,a)') '(a,', MSTATES, '(f12.6))'
      !write(ounit, '(a,<MSTATES>(f12.6))') ' anormo : ', anormo(1:i) ! Intel version
      write(ounit, temp) ' anorm : ', anormo(1:i)                   ! GNU version
    else
      anormo = 1.0d0
      write(ounit,'(a,t40, 10f12.6)') 'Using default anorm correction ', anormo(1:nstates)
    endif
  endif

! Read in sr_lambda if multi-state sr_n
  if(nstates.gt.1.and.method.eq.'sr_n'.and.ioptwf.eq.1) then
    write(ounit, *) "SR lambda: imposing overlap penalties"

    ! Part which handles the overlap penalty factors
    if (.not. allocated(isr_lambda)) allocate (isr_lambda(MSTATES*(MSTATES-1)/2))
    if (.not. allocated(sr_lambda)) allocate (sr_lambda(MSTATES,MSTATES))

    if ( fdf_islreal('sr_lambda') .and. fdf_islist('sr_lambda') &
        .and. (.not. fdf_islinteger('sr_lambda')) ) then
      i = -1
      call fdf_list('sr_lambda',i,isr_lambda)
      write(ounit,'(a)' )
      write(ounit,'(tr1,a,i0,a)') ' SR lambda has ',i,' entries'
      if(i.ne.nstates*(nstates-1)/2) call fatal_error('READ_INPUT: sr_lambda array, &
          &must contain nstates*(nstates-1)/2 entries, [ 1-2, 1-3, ..., 1-nstates, &
          &2-3, 2-4, ..., 2-nstates, ..., (nstates-1)-nstates ]')
      call fdf_list('sr_lambda',i,isr_lambda)
      write(temp, '(a,i0,a)') '(a,', MSTATES*(MSTATES-1)/2, '(f12.6))'
      !write(ounit, '(a,<MSTATES*(MSTATES-1)/2>(f12.6))') ' SR lambda : ', sr_lambda(1:i) ! Intel version
      write(ounit, temp) ' SR lambda : ', isr_lambda(1:i)                   ! GNU version
    else
      isr_lambda = 0.0d0
      write(ounit,'(a,t40, 10f12.6)') 'Default SR lambda ', isr_lambda(1:nstates)
    endif
    ! reshaping into symmetric matrix, its less messy to reference
    do i = 1, nstates
      do j = 1, nstates
        if(j.gt.i) then
          sr_lambda(i,j) = isr_lambda((i-1)*nstates-i*(i-1)/2+j-i)
!         sr_lambda(j,i) = sr_lambda(i,j)
! TMP -> to change optwf_sr/compute_grad
          sr_lambda(j,i) = 0.d0
        endif
      enddo
    enddo
  endif

! QMMM classical potential
  ! if(iqmmm.gt.0) then
  !   write(ounit,'(a)' ) "QMMM external potential "
  !   call qmmm_extpot_read
  ! endif

! Read in point charges
  if(iefield.gt.0) then
! External charges fetched in read_efield
    write(ounit,'(a,i4,a)')  " Read in ", ncharges, " external charges"
    call efield_compute_extint
  endif

! Additional Properties
! properties will be sampled iprop
! properties will be printed ipropprt
! if(iprop.ne.0) then
!   nprop=MAXPROP
!   write(ounit,'(a)' ) " Properties will be sampled "
!   write(ounit,int_format ) " Properties printout flag = ", ipropprt
!   call prop_cc_nuc(znuc,cent,iwctype,nctype_tot,ncent_tot,ncent,cc_nuc)
! endif

! (13) Forces information (either block or from a file) [#####]

  if (fdf_block('forces', bfdf)) then
    call fdf_read_forces_block(bfdf)
  elseif (fdf_load_defined('forces') ) then
    call read_forces_file(file_forces)
  else
    if(nforce.ge.1.and.iforces.eq.0.and.igradients.eq.0) then
      write(errunit,*) "INPUT: block forces_displace or gradients_* missing: geometries set equal to primary"
      call inputforces()
    endif
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
    write(errunit,'(a)') "Warning:: No information about dmatrix provided in the input."
    write(errunit,'(3a,i6)') "Stats for nerds :: in file ",__FILE__, " at line ", __LINE__
!    error stop
  endif


! Part which handles the weights. needs modifications for guiding

  if (.not. allocated(weights)) allocate (weights(MSTATES))
  if (.not. allocated(iweight)) allocate (iweight(MSTATES))


  if ( fdf_islreal('weights') .and. fdf_islist('weights') &
      .and. (.not. fdf_islinteger('weights')) ) then
    i = -1
    call fdf_list('weights',i,weights)
    write(ounit,'(a)' )
    write(ounit,'(tr1,a,i0,a)') ' Weights has ',i,' entries'
    call fdf_list('weights',i,weights)
    write(temp, '(a,i0,a)') '(a,', MSTATES, '(f12.6))'
    !write(ounit, '(a,<MSTATES>(f12.6))') 'weights : ', weights(1:i)  ! Intel version
    write(ounit, temp) 'weights : ', weights(1:i)                    ! GNU version
  else
    weights = 1.0d0
    write(ounit,'(a,t40, 10f12.6)') 'Default weights ', weights(1:nstates)
  end if

  wsum = 0.d0
  nweight = 0
  do i = 1, nstates
    if (weights(i) .gt. 1d-6) then
        nweight = nweight + 1
        iweight(nweight) = i
        wsum = wsum + weights(i)
    endif
  enddo

  do i = 1, nweight
    weights(i) = weights(i)/wsum
  enddo

  if (nweight .eq. 0) then
      nweight = 1
      iweight(1) = 1
      weights(1) = 1.d0
  endif

! The above part should be moved to get_weights subroutine


! Processing of data read from the parsed files or setting them with defaults


! (10) optorb_mixvirt information of orbitals (either block or from a file)

  if ( fdf_load_defined('optorb_mixvirt') ) then
    call read_optorb_mixvirt_file(file_optorb_mixvirt)
  elseif ( fdf_block('optorb_mixvirt', bfdf)) then
  ! call fdf_read_optorb_mixvirt_block(bfdf)
    write(errunit,'(a)') "Error:: No information about optorb_mixvirt provided in the block."
    write(errunit,'(3a,i6)') "Stats for nerds :: in file ",__FILE__, " at line ", __LINE__
    error stop
  else
    write(errunit,'(a)') "Warning:: No information about optorb_mixvirt provided in the block."
    write(errunit,'(3a,i6)') "Stats for nerds :: in file ",__FILE__, " at line ", __LINE__
  endif


! (18) cavity_spheres information (either block or from a file)

!   if ( fdf_load_defined('cavity_spheres') ) then
!     call read_cavity_spheres_file(file_cavity_spheres)
!   elseif ( fdf_block('cavity_spheres', bfdf)) then
!   ! call fdf_read_cavity_spheres_block(bfdf)
!     write(errunit,'(a)') "Error:: No information about cavity_spheres provided in the block."
!     !write(errunit,'(3a,i6)') "Stats for nerds :: in file ",__FILE__, " at line ", __LINE__
!     error stop
!   else
!     write(errunit,'(a)') "Error:: No information about cavity_spheres provided in the block."
!     !write(errunit,'(3a,i6)') "Stats for nerds :: in file ",__FILE__, " at line ", __LINE__
! !    error stop
!   endif


! ZMATRIX begins here
! gradients_zmatrix information (either block or from a file)

  if ( fdf_load_defined('gradients_zmatrix') ) then
    call read_gradients_zmatrix_file(file_gradients_zmatrix)
  elseif ( fdf_block('gradients_zmatrix', bfdf)) then
  ! call fdf_read_gradients_zmatrix_block(bfdf)
    write(errunit,'(a)') "Error:: No information about gradients_zmatrix provided in the block."
    write(errunit,'(3a,i6)') "Stats for nerds :: in file ",__FILE__, " at line ", __LINE__
    error stop
  else
    write(errunit,'(a)') "Warning:: No information about gradients_zmatrix provided in the block."
    write(errunit,'(3a,i6)') "Stats for nerds :: in file ",__FILE__, " at line ", __LINE__
!    error stop
  endif

! gradients_cartesian information (either block or from a file)

  if ( fdf_load_defined('gradients_cartesian') ) then
    call read_gradients_cartesian_file(file_gradients_cartesian)
  elseif ( fdf_block('gradients_cartesian', bfdf)) then
  ! call fdf_read_gradients_cartesian_block(bfdf)
    write(errunit,'(a)') "Error:: No information about gradients_cartesian provided in the block."
    write(errunit,'(3a,i6)') "Stats for nerds :: in file ",__FILE__, " at line ", __LINE__
    error stop
  else
    write(errunit,'(a)') "Warning:: No information about gradients_cartesian provided in the block."
    write(errunit,'(3a,i6)') "Stats for nerds :: in file ",__FILE__, " at line ", __LINE__
!    error stop
  endif

! modify_zmatrix information (either block or from a file)

  if(iforce_analy.gt.0) then
    if ( fdf_load_defined('modify_zmatrix') ) then
      call read_modify_zmatrix_file(file_modify_zmatrix)
    elseif ( fdf_block('modify_zmatrix', bfdf)) then
    ! call fdf_read_modify_zmatrix_block(bfdf)
      write(errunit,'(a)') "Error:: No information about modify_zmatrix provided in the block."
      write(errunit,'(3a,i6)') "Stats for nerds :: in file ",__FILE__, " at line ", __LINE__
      error stop
    else
      if(imodify_zmat.eq.0) call modify_zmat_define
    endif
  endif

! hessian_zmatrix information (either block or from a file)
  if(iforce_analy.gt.0) then
    if ( fdf_load_defined('hessian_zmatrix') ) then
      call read_hessian_zmatrix_file(file_hessian_zmatrix)
    elseif ( fdf_block('hessian_zmatrix', bfdf)) then
    ! call fdf_read_hessian_zmatrix_block(bfdf)
      write(errunit,'(a)') "Error:: No information about hessian_zmatrix provided in the block."
      write(errunit,'(3a,i6)') "Stats for nerds :: in file ",__FILE__, " at line ", __LINE__
      error stop
    else
      if(ihessian_zmat.eq.0) call hessian_zmat_define
    endif
  endif

! zmatrix_connection information (either block or from a file)

  if ( fdf_load_defined('zmatrix_connection') ) then
    call read_zmatrix_connection_file(file_zmatrix_connection)
  elseif ( fdf_block('zmatrix_connection', bfdf)) then
  ! call fdf_read_zmatrix_connection_block(bfdf)
    write(errunit,'(a)') "Error:: No information about zmatrix_connection provided in the block."
    write(errunit,'(3a,i6)') "Stats for nerds :: in file ",__FILE__, " at line ", __LINE__
    error stop
  else
    write(errunit,'(a)') "Warning:: No information about zmatrix_connection provided in the block."
    write(errunit,'(3a,i6)') "Stats for nerds :: in file ",__FILE__, " at line ", __LINE__
    if(iuse_zmat.gt.0.and.izmatrix_check.eq.0) call fatal_error('INPUT: block connectionzmatrix missing')
  endif

! Some checks on Z Matrixs.
! Write out information about calculation of energy gradients and Z matrix
  write(ounit,*)
  if(ngradnts.gt.0 .and. igrdtype.eq.1) call inpwrt_grdnts_cart()
  if(ngradnts.gt.0 .and. igrdtype.eq.2) call inpwrt_grdnts_zmat()
  if(izmatrix.eq.1) call inpwrt_zmatrix()
  write(ounit,*)

! ZMATRIX section ends here

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
    write(errunit,'(3a,i6)') "Stats for nerds :: in file ",trim(__FILE__), " at line ", __LINE__
!    error stop
  endif

! Done reading all the files

! ! Not sure if this line should be here or not.
! call pot_nn(cent,znuc,iwctype,ncent,pecent)

! Make sure that all the blocks are read. Use inputflags here to check

  if(iznuc.eq.0) call fatal_error('INPUT: block znuc missing')
  if(igeometry.eq.0) call fatal_error('INPUT: block geometry missing')
  if(ijastrow_parameter.eq.0) call fatal_error('INPUT: block jastrow_parameter missing')
  if(iefield.gt.0.and.icharge_efield.eq.0) call fatal_error('INPUT: block efield missing')

  call fdf_shutdown()

  ! The following portion can be shifted to another subroutine.
  ! It does the processing of the input read so far and initializes some
  ! arrays if something is missing.

  !qmckl initialization

#if defined(TREXIO_FOUND) && defined(QMCKL_FOUND)
  if (use_qmckl) then

     if (nwftypeorb.gt.1) call fatal_error('Error: QMCKL does not yet support multi-orbital calculations. ')

     qmckl_no_ctx = 2
     if(ioptorb.gt.0) qmckl_no_ctx = 3

     ! Create a new QMCkl context
     do ictx=1,qmckl_no_ctx
       qmckl_ctx(ictx) = qmckl_context_create()
       write(ounit, *) " QMCkl initial context created  " , qmckl_ctx(ictx) , " successfully "
     enddo

     if(ioptorb.gt.0) then

       file_trexio_new = file_trexio(1:index(file_trexio,'.hdf5')-1)//'_orbchanged.hdf5'
       if((file_trexio_new(1:6) == '$pool/') .or. (file_trexio_new(1:6) == '$POOL/')) then
           file_trexio_path = pooldir // file_trexio_new(7:)
       else
           file_trexio_path = file_trexio_new
       endif

       if(wid) then
         if (trexio_inquire(file_trexio_path) .eq. TREXIO_SUCCESS) then
           write(ounit,'(a)') "Removing existing " // file_trexio_path // " file"
           call system('rm -v ' // file_trexio_path)
         endif

         rc = trexio_cp(file_trexio, file_trexio_path)
         if (rc .ne. TREXIO_SUCCESS) call fatal_error('INPUT: QMCkl error: Unable to copy trexio file')
       endif
       call MPI_Barrier( MPI_COMM_WORLD, ierr )

       do ictx=1,qmckl_no_ctx
         rc = qmckl_trexio_read(qmckl_ctx(ictx), file_trexio_path, 1_8*len(trim(file_trexio_path)))
         write(ounit, *) "Status QMCKl trexio read file_trexio_path", rc
         if (rc /= QMCKL_SUCCESS) call fatal_error('PARSER: QMCkl error: Unable to read TREXIO file')
       enddo

      else ! ioptorb.eq.0
        do ictx=1,qmckl_no_ctx
          rc = qmckl_trexio_read(qmckl_ctx(ictx), file_trexio, 1_8*len(trim(file_trexio)))
          write(ounit, *) "Status QMCKl trexio read file_trexio_path", rc
          if (rc /= QMCKL_SUCCESS) call fatal_error('PARSER: QMCkl error: Unable to read TREXIO file')
        enddo
     endif

     ! get mo's number should correspond to norb_tot

     call jastrow_init_qmckl(qmckl_no_ctx)

     call init_orbitals_qmckl(.True.)



  endif
#endif

  !----------------------------------------------------------------------------END
  contains

  !! Here all the subroutines that handle the block data are written
  !! Quick Tutorial about FDF syntax
  !!   'v' ('value') is matched by both an 'integer' and a 'real'.
  !!   'j' is matched by both an 'integer' and a 'name'.
  !!   's' is matched by an 'integer', a 'real', and a 'name'.
  !!   'x' is matched by any kind of token.
  !!   'a' is matched by a list with integers
  !!   'c' is matched by a list with reals
  !!   'e' is matched by a list with integers or reals
  !!   'd' is reserved for future dictionaries...

!! for periodic assuming square cell still
    subroutine read_lattice_file(file_lattice)

      use contrl_file, only: errunit,ounit
      use m_string_operations, only: wordcount
      use custom_broadcast, only: bcast
      use mpiconf, only: wid
      use periodic, only: alattice
      use periodic, only: rlatt, rlatt_inv
      use precision_kinds, only: dp

      implicit none

      character(len=72), intent(in)   :: file_lattice
      character(len=90)               :: file_lattice_path, line
      real(dp) :: latt
      integer                         :: iunit, iostat, count
      logical                         :: exist
      integer :: i,k

      !   External file reading

      if((file_lattice(1:6) == '$pool/') .or. (file_lattice(1:6) == '$POOL/')) then
         file_lattice_path = pooldir // file_lattice(7:)
      else
         file_lattice_path = file_lattice
      endif

      write(ounit,*) '-----------------------------------------------------------------------'
      write(ounit,string_format)  " Reading Lattice Parameters from the file :: ",  trim(file_lattice_path)
      write(ounit,*)

      if (wid) then
         inquire(file=file_lattice_path, exist=exist)
         if (exist) then
            open (newunit=iunit,file=file_lattice_path, iostat=iostat, action='read' )
            if (iostat .ne. 0) call fatal_error("Problem opening the Super-cell file")
         else
            call fatal_error (" Super-cell file "// trim(file_lattice) // " does not exist.")
         endif
      endif

      !Initialization to avoid garbage

      rlatt = 0.d0
      rlatt_inv = 0.d0
      count=0

      if (wid) then
         do i = 1, 3
            !! Reading each row as a vector of the box
            read(iunit,*, iostat=iostat) rlatt(i,1), rlatt(i,2), rlatt(i,3)
            if(iostat.ne.0.and.count.eq.0) then
               if(rlatt(1,1).gt.0.d0) then
                  alattice=rlatt(1,1)
                  count=count+1
                  exit
               else
                  write(ounit, *) rlatt(i,1), rlatt(i,2), rlatt(i,3)
                  call fatal_error("Error in reading lattice parameters file")
               endif
            endif
            count=count+1
            if(iostat.ne.0.and.count.ne.1) then
               write(ounit, *) rlatt(i,1), rlatt(i,2), rlatt(i,3)
               call fatal_error("Error in reading lattice parameters file")
            endif
         enddo

         if(count.eq.1) then
            write(ounit,*) 'This is a cubic cell'
            write(ounit,*) 'The lattice constant is', alattice

            rlatt(1,1) = alattice
            rlatt(2,2) = alattice
            rlatt(3,3) = alattice
         else if(count.eq.3) then
            write(ounit,*) "The simulation cell is:"

            write(ounit,*) "a", rlatt(1,1), rlatt(1,2), rlatt(1,3)
            write(ounit,*) "b", rlatt(2,1), rlatt(2,2), rlatt(2,3)
            write(ounit,*) "c", rlatt(3,1), rlatt(3,2), rlatt(3,3)

            !! assuming still rectangular box
            alattice=rlatt(1,1)
            do i = 2, 3
               if(rlatt(i,i).lt.alattice) alattice=rlatt(i,i)
            enddo

            if(alattice.le.0.d0) call fatal_error("Wrong lattice parameter")
         else
            call fatal_error("Error reading lattice file")
         endif

         !   regarding Ewald Breakup assume column not row vectors for the lattice
         rlatt=TRANSPOSE(rlatt)
      endif

      call bcast(rlatt)
      if (wid) close(iunit)

end subroutine read_lattice_file

  subroutine fdf_read_molecule_block(bfdf)
    implicit none

    type(block_fdf)                 :: bfdf
    type(parsed_line), pointer      :: pline
    double precision, allocatable   :: nval(:)
    integer                         :: count = 0
    ! %block molecule
    ! 4
    ! some comment (symbol, x,y,z)
    ! C   -3.466419  0.298187  0
    ! C    3.466419 -0.298187  0
    ! H   -3.706633  2.326423  0
    ! H    3.706633 -2.326423  0
    ! %endblock

    ! Example of znuc assignment in the block (last column)
    ! %block molecule
    ! 4
    ! some comment (symbol, x,y,z, znuc)
    ! C1   -3.466419  0.298187  0   4.0
    ! C2    3.466419 -0.298187  0   4.0
    ! H1   -3.706633  2.326423  0   1.0
    ! H2    3.706633 -2.326423  0   1.0
    ! %endblock

    ! Keep compiler happy
    if (.not.allocated(nval)) allocate(nval(1))
    nval = 0._8

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
      if (.not. allocated(nval)) allocate(nval(ncent))

      count = pline%ntokens
      ! get the coordinates: 4 tokens per line; first char (n) and three (r)reals or (i)ints.
      if ((pline%ntokens==4).and.((pline%id(1).eq."n").and.((any(pline%id(2:4).eq."r")) .or. (any(pline%id(2:4).eq.("i"))) ))) then
        symbol(j) = fdf_bnames(pline, 1)
        do i= 1, 3
          cent(i,j) = fdf_bvalues(pline, i)
        enddo
        j = j + 1
      ! get the coordinates: 5 tokens per line; first char (n) and three (r)reals or (i)ints for coords and 4th for nvalence/znuc
      elseif ((pline%ntokens==5).and.((pline%id(1).eq."n").and.((any(pline%id(2:4).eq."r")) .or. (any(pline%id(2:4).eq.("i"))) ))) then
        symbol(j) = fdf_bnames(pline, 1)
        do i= 1, 3
          cent(i,j) = fdf_bvalues(pline, i)
        enddo
        nval(j) = fdf_bvalues(pline, 4)
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
            if (symbol(j) == unique(k))  then
              iwctype(j) = k
              if (count .gt. 4) znuc(k) = nval(j)
          endif
        enddo
    enddo

    ! Get the correspondence rule
    do k = 1, nctype
        atomtyp(k) = unique(k)
    enddo
    if (allocated(unique)) deallocate(unique)

    if (count == 4) then
      ! Get the znuc for each unique atom
      do j = 1, nctype
          atoms = element(atomtyp(j))
          znuc(j) = atoms%nvalence
      enddo
    endif

    ncent_tot = ncent + nghostcent
    nctype_tot = nctype + newghostype

    write(ounit,*) '-----------------------------------------------------------------------'
    write(ounit,'(a, t15, a, t27, a, t39, a, t45, a)') 'Symbol', 'x', 'y', 'z', 'Type'
    write(ounit,'(t14, a, t26, a, t38, a )') '(bohr)', '(bohr)', '(bohr)'
    write(ounit,*) '-----------------------------------------------------------------------'

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
    integer                    :: i,j,k,l

!  Format of the forces block
!
!    %block forces
!    1 1
!   -0.733652000000     0.1   -1.157935000000
!    0.733652000000     0.1    1.157935000000
!    0.298187000000     0.1   -3.466419000000
!   -0.298187000000     0.1    3.466419000000
!    2.326423000000     0.1   -3.706633000000
!   -2.326423000000     0.1    3.706633000000
!   -2.772181000000     0.1   -0.963193000000
!    2.772181000000     0.1    0.963193000000
!   -0.855551000000     0.1   -5.146950000000
!    0.855551000000     0.1    5.146950000000
!    #
!   -3.733652000000     2.1   -1.157935000000
!    3.733652000000     2.1    1.157935000000
!    3.298187000000     2.1   -3.466419000000
!   -3.298187000000     2.1    3.466419000000
!    3.326423000000     2.1   -3.706633000000
!   -3.326423000000     2.1    3.706633000000
!   -3.772181000000     2.1   -0.963193000000
!    3.772181000000     2.1    0.963193000000
!   -3.855551000000     2.1   -5.146950000000
!    3.855551000000     2.1    5.146950000000
!    %endblock

    if (.not. allocated(delc)) allocate (delc(3, ncent, nforce))
    if (.not. allocated(iwftype)) allocate (iwftype(nforce))

    i = 1; j = 1
    do while((fdf_bline(bfdf, pline)))
      if (pline%ntokens == nforce .and. (pline%id(1) .eq. "i") ) then
        do l = 1, nforce
          iwftype(l) = fdf_bintegers(pline, l)
        enddo
      endif
      if (pline%ntokens == 3 .and. (any(pline%id(1:3).eq."r")) ) then
        do k = 1, 3
          delc(k, j, i) = fdf_bvalues(pline, k)
        enddo ! xyz

        if (mod(j, ncent) == 0) then
          i = i + 1
          j = 1
        else
          j = j + 1
        endif

      endif ! expect only three values in a line
    enddo ! parse entire file


    write(ounit,*) 'Force displacements from the %block forces  '
    write(ounit,*)
    do i = 1, nforce
      write(ounit,*) '-----------------------------------------------------------------------'
      write(ounit,'(a,i4)') 'Number (iwftype) :: ',i
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

    if( (method(1:3) == 'lin')) then
      if (.not. allocated(ccsf)) allocate(ccsf(ncsf, nstates, 3))
    else
      if (.not. allocated(ccsf)) allocate(ccsf(ncsf, nstates, nwftype))
    endif

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

    write(ounit,'(10(1x, a9, i3, 1x))') (" State: ", i, i =1, nstates)
    do j = 1, ncsf
      write(ounit,'(10(1x, f12.8, 1x))') (ccsf(j,i,1), i=1,nstates)
    enddo

  end subroutine fdf_read_csf_block

  subroutine fdf_read_jastrow_block(bfdf)
    use jastrow4_mod,   only: nterms4

    implicit none

    type(block_fdf)            :: bfdf
    type(parsed_line), pointer :: pline
    integer                    :: i,j,k, iwft
    integer                    :: mparmja, mparmjb, mparmjc
!   Format of the jastrow block
!   %block csf
!   jastrow_parameter 1
!    5  5  0           norda,nordb,nordc
!     0.60000000   0.00000000     scalek
!     0.00000000   0.00000000   0.05946443  -0.68575835   0.42250502  -0.10845009 (a(iparmj),iparmj=1,nparma)
!     0.00000000   0.00000000  -0.13082284  -0.06620300   0.18687803  -0.08503472 (a(iparmj),iparmj=1,nparma)
!     0.50000000   0.38065787   0.16654238  -0.05430118   0.00399345   0.00429553 (b(iparmj),iparmj=1,nparmb)
!   %endblock


    if( (method(1:3) == 'lin')) then
      allocate (scalek(3))
    else
      allocate (scalek(nwftype))
    endif

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

        nordj = max(norda, nordb, nordc)
        nordj1 = nordj + 1
        neqsx = 6*nordj

        mparmja = 2 + max(0, norda - 1)
        mparmjb = 2 + max(0, nordb - 1)
        mparmjc = nterms4(nordc)

        print*, "mparmja ", mparmja
        print*, "mparmjb ", mparmjb
        print*, "mparmjc ", mparmjc
      endif

      ! check if it is the only integer present in a line
      if ((pline%id(2) .eq. "r") .and. (pline%id(4) .eq. "r") .and. (pline%id(1) .eq. "n") .and. (pline%id(3) .eq. "n")) then
        scalek(iwft) = fdf_bvalues(pline, 1) ! 1st integer in the line
        print*, "scalek ", scalek(iwft)
      endif

      !print*, "check all the ids ", pline%id(1:6)
      ! fix this part

      ! if ((pline%id(1) .eq. "n") .and. (pline%ntokens .eq. mparmja+1) ) then  !check if it is the only integer present in a line
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

        implicit none

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
  ! use vmc_mod, only: nmat_dim, nmat_dim2
  ! use const, only: nelec
  ! use system, only: nctype_tot, ncent_tot

  use sr_mod, only: mparm, mobs, mconf
  use control, only: mode
  use control_vmc, only: vmc_nstep, vmc_nblk_max
  use control_dmc, only: dmc_nstep
  use vmc_mod, only: set_vmc_size
  use optci, only: set_optci_size
  use optorb_mod, only: set_optorb_size
  use gradhess_all, only: set_gradhess_all_size
  ! use sr_mod, only: set_sr_size

  implicit none

  ! leads to circular dependecy of put in sr_mod ..
  mobs = 10 + 6*mparm

  if( mode(1:3) == 'vmc' ) then
     mconf = vmc_nstep * vmc_nblk_max
  else
     mconf = dmc_nstep * vmc_nblk_max
  endif

  call set_vmc_size
  call set_optci_size
  call set_optorb_size
  call set_gradhess_all_size
  ! call set_sr_size

end subroutine compute_mat_size_new
end module
