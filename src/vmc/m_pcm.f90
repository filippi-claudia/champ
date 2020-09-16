module pcm 
  integer, parameter :: MCHS=1
  integer, parameter :: MCHV=1
  integer, parameter :: MSPHERE=30
  private 
  public :: MCHS, MCHV, MSPHERE
  save 
 end module pcm

 module pcm_3dgrid
  !     flags and dimensions for the 3d grid objects
  use precision_kinds, only: dp
  integer, parameter :: MGRID_PCM=1
  integer, parameter :: IUNDEFINED = -1234567890
  integer, parameter :: MGRID_PCM2=MGRID_PCM*MGRID_PCM
  integer, parameter :: MGRID_PCM3=MGRID_PCM2*MGRID_PCM
  real(dp), parameter :: UNDEFINED = -1234567890.d0

  private 
  public :: MGRID_PCM, MGRID_PCM2, MGRID_PCM3
  public :: UNDEFINED, IUNDEFINED 
  save 
 end module pcm_3dgrid
 
 module pcm_ah
   !> Arguments: ahca, bh
   use pcm, only: MCHS
   use precision_kinds, only: dp
 
   real(dp) :: ahca(MCHS,MCHS)
   real(dp) :: bh(MCHS)

   private
   public :: ahca, bh
   save
 end module pcm_ah

 module pcm_ameta
   !> Arguments: amdlg, eta
   use pcm, only: MCHS
   use precision_kinds, only: dp
 
   real(dp) :: amdlg(MCHS)
   real(dp) :: eta(3,MCHS)

   private
   public :: amdlg, eta
   save
 end module pcm_ameta

 module pcm_averages
   !> Arguments: spcmsum, spcmcum, spcmcm2, vpcmsum, vpcmcum, vpcmcm2
   ! qopcm_sum, qopcm_cum, qopcm_cm2, 
   ! enfpcm_sum(MCHS), enfpcm_cum(MCHS), enfpcm_cm2(MCHS)
   use pcm, only: MCHS
   use precision_kinds, only: dp
 
    real(dp) :: spcmsum
    real(dp) :: spcmcum
    real(dp) :: spcmcm2
    real(dp) :: vpcmsum
    real(dp) :: vpcmcum
    real(dp) :: vpcmcm2
    real(dp) :: qopcm_sum
    real(dp) :: qopcm_cum
    real(dp) :: qopcm_cm2
    real(dp) :: enfpcm_sum(MCHS)
    real(dp) :: enfpcm_cum(MCHS)
    real(dp) :: enfpcm_cm2(MCHS)
 
    private
    public :: spcmsum, spcmcum, spcmcm2, vpcmsum, vpcmcum, vpcmcm2
    public :: qopcm_sum, qopcm_cum, qopcm_cm2
    public :: enfpcm_sum, enfpcm_cum, enfpcm_cm2
    save
 end module pcm_averages

 module pcm_cntrl
    !> Arguments: ichpol, ipcm, ipcmprt, icall, isurf
 
    integer  :: icall
    integer  :: ichpol
    integer  :: ipcm
    integer  :: ipcmprt
    integer  :: isurf

    private
    public :: ichpol, ipcm, ipcmprt, icall, isurf
    save
 end module pcm_cntrl
 
 module pcm_fdc
   !> Arguments: fs, rcol, feps, rcolt, rcolv, qfree, qvol
   use precision_kinds, only: dp

   real(dp) :: feps
   real(dp) :: fs
   real(dp) :: qfree
   real(dp) :: qvol
   real(dp) :: rcol
   real(dp) :: rcolt
   real(dp) :: rcolv

   private
   public :: fs, rcol, feps, rcolt, rcolv, qfree, qvol
   save
 end module pcm_fdc

 module pcm_force
   !> Arguments: sch_s
   use pcm, only: MCHS
   use force_mod, only: MFORCE
   use precision_kinds, only: dp

   real(dp) :: sch_s(MCHS,MFORCE)

   private
   public :: sch_s
   save
 end module pcm_force

 module pcm_grid3d_contrl
   !> Arguments: ipcm_3dgrid
 
   integer  :: ipcm_3dgrid
 
   private
   public :: ipcm_3dgrid
   save
 end module pcm_grid3d_contrl

 module pcm_grid3d_array
   !> Arguments: pcm_cart_from_int
   use pcm_3dgrid, only: MGRID_PCM
   use precision_kinds, only: dp

   real(dp) :: pcm_cart_from_int(MGRID_PCM,3)

   private
   public :: pcm_cart_from_int
   save
 end module pcm_grid3d_array
 
 module pcm_grid3d_param
   !> Arguments: pcm_endpt, pcm_origin, ipcm_nstep3d, pcm_step3d
   use precision_kinds, only: dp
 
   integer  :: ipcm_nstep3d(3)
   real(dp) :: pcm_endpt(3)
   real(dp) :: pcm_origin(3)
   real(dp) :: pcm_step3d(3)
 
   private
   public :: pcm_endpt, pcm_origin, ipcm_nstep3d, pcm_step3d
   save
 end module pcm_grid3d_param

 module pcm_hpsi
   !> Arguments: enfpcm, pepcms, pepcmv, qopcm
   use pcm, only: MCHS
   use precision_kinds, only: dp

   real(dp) :: enfpcm(MCHS)
   real(dp) :: pepcms
   real(dp) :: pepcmv
   real(dp) :: qopcm

   private
   public :: enfpcm, pepcms, pepcmv, qopcm
   save
 end module pcm_hpsi

 module pcm_inda
   use pcm, only: MCHS
   !> Arguments: inda
 
   integer  :: inda(MCHS)

   private
   public :: inda
   save
 end module pcm_inda

 module m_pcm_num_spl
   !> Arguments: pcm_num_spl
   use pcm_3dgrid, only: MGRID_PCM
   use precision_kinds, only: dp
 
   real(dp) :: pcm_num_spl(8,MGRID_PCM,MGRID_PCM,MGRID_PCM)

   private
   public :: pcm_num_spl
   save
 end module m_pcm_num_spl

 module pcm_num_spl2
   !> Arguments: bc, wk
   use precision_kinds, only: dp

   real(dp) :: bc
   real(dp) :: wk

   private
   public :: bc, wk
   save
 end module pcm_num_spl2

 module pcm_parms
   !> Arguments: re, nchv, nesph, ze, iscov, eps_solv, xpol, 
   !             retk, ch, xe, nvopcm, nch, re2, ncopcm, surk, nscv, nchs, ye, nchs2, nchs1
   use pcm, only: MCHV, MSPHERE
   use precision_kinds, only: dp

   real(dp) :: ch(MCHV)
   real(dp) :: eps_solv
   integer  :: iscov
   integer  :: nch
   integer  :: nchs
   integer  :: nchs1
   integer  :: nchs2
   integer  :: nchv
   integer  :: ncopcm
   integer  :: nesph
   integer  :: nscv
   integer  :: nvopcm
   real(dp) :: re(MSPHERE)
   real(dp) :: re2(MSPHERE)
   real(dp) :: retk
   real(dp) :: surk
   real(dp) :: xe(MSPHERE)
   real(dp) :: xpol(3,MCHV)
   real(dp) :: ye(MSPHERE)
   real(dp) :: ze(MSPHERE)

   private
   public :: re, nchv, nesph, ze, iscov, eps_solv, xpol
   public :: ch, xe, nvopcm, nch, re2, ncopcm, surk
   public :: retk, nscv, nchs, ye, nchs2, nchs1
   save
 end module pcm_parms

 module pcm_pot
   !> Arguments: penupol, penups, penupv
   use precision_kinds, only: dp
 
   real(dp) :: penupol
   real(dp) :: penups
   real(dp) :: penupv
   private
 
   public :: penupol, penups, penupv
   save
 end module pcm_pot

 module pcm_xv_new
   !> Arguments: xv_new
   use pcm, only: MCHV
   use precision_kinds, only: dp

   real(dp) :: xv_new(3,MCHV)

   private
   public :: xv_new
   save
 end module pcm_xv_new

 module pcm_unit
   !> Arguments: pcmfile_cavity, pcmfile_chv, pcmfile_chs
 
    character*80 :: pcmfile_cavity
    character*80 :: pcmfile_chs
    character*80 :: pcmfile_chv

   private
   public :: pcmfile_cavity, pcmfile_chv, pcmfile_chs
   save
 end module pcm_unit

 module pcmo
   !> Arguments: enfpcmo, qopcmo, spcmo, vpcmo
   use pcm, only: MCHS
   use precision_kinds, only: dp

   real(dp) :: enfpcmo(MCHS)
   real(dp) :: qopcmo
   real(dp) :: spcmo
   real(dp) :: vpcmo

   private
   public :: enfpcmo, qopcmo, spcmo, vpcmo
   save
 end module pcmo