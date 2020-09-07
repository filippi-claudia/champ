!> \brief File collecting all modules related to the grid 
!>
!> \author  P. Lopez-Tarifa & F. Zapata NLeSC(2019)

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

  PRIVATE
  PUBLIC :: MORB_OCC
  PUBLIC :: orb_num_spl
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