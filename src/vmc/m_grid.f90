module grid_mod
  !> Arguments
  ! flags and dimensions for the 3d grid objects
  use precision_kinds, only: dp, sp
  
  integer, parameter :: MXNSTEP = 1
  ! integer, parameter :: MXNSTEP = 50
  integer, parameter :: MXNSTEP2 = MXNSTEP*MXNSTEP
  integer, parameter :: MXNSTEP3 = MXNSTEP2*MXNSTEP
  
  integer, parameter :: IUNDEFINED = -1234567890
  real(dp), parameter :: UNDEFINED = -1234567890.d0, SHIFT = 2.d0
  
  real(sp) grid3d(MXNSTEP,MXNSTEP,MXNSTEP), cart_from_int(MXNSTEP,3)

  private 
  public :: MXNSTEP, MXNSTEP2, MXNSTEP3
  public :: IUNDEFINED, UNDEFINED, SHIFT 
  public :: grid3d, cart_from_int
  save 

 end module grid_mod

 module grid_spline_mod
  !> Arguments
  use precision_kinds, only: sp
  use vmc, only: MELEC
  use grid_mod, only: MXNSTEP

  integer, parameter :: MORB_OCC = MELEC/2+3
  real(sp) orb_num_spl(8,MXNSTEP,MXNSTEP,MXNSTEP,MORB_OCC)

  private 
  public :: MORB_OCC
  public :: orb_num_spl
  save 

 end module grid_spline_mod

 module grid_lagrange_mod
  !> argument
  use precision_kinds, only: sp
  use grid_mod, only: MXNSTEP
  use vmc, only: MELEC
  ! Number of Lagrange interpolation points/axis
  integer, parameter :: LAGMAX=4
  integer, parameter :: LAGSTART = -LAGMAX/2, LAGEND=LAGSTART+LAGMAX-1
  integer, parameter ::  MORB_OCC = MELEC/2
  
  !  Spline fits of the orbitals
  ! and boundary conditions (for the creation of the fit)
  real(sp) orb_num_lag(5,MXNSTEP,MXNSTEP,MXNSTEP,MORB_OCC)

  private 
  public :: LAGMAX, LAGSTART, LAGEND, MORB_OCC
  public :: orb_num_lag
  save
 end module grid_lagrange_mod

 module grid3d_param
   !> Arguments: nstep3d, endpt, origin, step3d
   use precision_kinds, only: dp

   real(dp) :: endpt(3)
   integer  :: nstep3d(3)
   real(dp) :: origin(3)
   real(dp) :: step3d(3)

   private
   public :: nstep3d, endpt, origin, step3d
   save
 end module grid3d_param

 module grid3dflag
   !> Arguments: i3dsplorb, i3dlagorb, i3dgrid, i3ddensity

   integer  :: i3ddensity
   integer  :: i3dgrid
   integer  :: i3dlagorb
   integer  :: i3dsplorb

   private
   public :: i3dsplorb, i3dlagorb, i3dgrid, i3ddensity
   save
 end module grid3dflag

  module orbital_num_lag
   !> Arguments: denom
   use precision_kinds, only: dp
   use grid_lagrange_mod, only: LAGSTART, LAGEND

   real(dp) :: denom(LAGSTART:LAGEND,3)
   real(dp) :: step_inv(3,3)

   private
   public :: denom, step_inv
   save
end module orbital_num_lag