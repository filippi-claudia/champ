!> @brief Module for 3D grid infrastructure and spatial discretization.
!> @author CHAMP developers
!> @date 2022
!>
!> @details This module provides the basic infrastructure for 3D spatial grids
!> used to evaluate orbitals, densities, and other quantities on a regular mesh.
!> The grid is defined by MXNSTEP points per dimension, with utilities for
!> coordinate transformation between integer indices and Cartesian positions.
!>
!> Key components:
!> - Grid dimensions: MXNSTEP (points/axis), MXNSTEP2 (2D), MXNSTEP3 (3D total)
!> - Grid storage: grid3d (orbital/density values on mesh)
!> - Coordinate mapping: cart_from_int (index to Cartesian transformation)
!> - Undefined values: IUNDEFINED, UNDEFINED (sentinel values for uninitialized data)
!>
!> @note MXNSTEP is currently set to 1 for minimal memory usage; increase for production.
!> @note Grid is stored in single precision (sp) for memory efficiency.
module grid_mod
    
    use precision_kinds, only: dp,sp

    implicit none

    !> Maximum number of grid points per axis (currently set to 1)
    integer, parameter :: MXNSTEP = 1
    ! integer, parameter :: MXNSTEP = 50
    
    !> Total grid points in 2D slice (MXNSTEP^2)
    integer, parameter :: MXNSTEP2 = MXNSTEP*MXNSTEP
    
    !> Total grid points in 3D volume (MXNSTEP^3)
    integer, parameter :: MXNSTEP3 = MXNSTEP2*MXNSTEP

    !> Integer sentinel value for undefined/uninitialized grid indices
    integer, parameter :: IUNDEFINED = -1234567890
    
    !> Real sentinel value for undefined/uninitialized grid data
    real(dp), parameter :: UNDEFINED = -1234567890.d0
    
    !> Shift parameter for grid coordinate transformations
    real(dp), parameter :: SHIFT = 2.d0

    !> 3D grid for storing orbital or density values (MXNSTEP, MXNSTEP, MXNSTEP)
    real(sp), dimension(:, :, :), allocatable :: grid3d
    
    !> Mapping from integer grid indices to Cartesian coordinates (MXNSTEP, 3)
    real(sp), dimension(:, :), allocatable :: cart_from_int

    private
    public :: MXNSTEP, MXNSTEP2, MXNSTEP3
    public :: IUNDEFINED, UNDEFINED, SHIFT
    public :: grid3d, cart_from_int
    public :: allocate_grid_mod, deallocate_grid_mod
    save

contains
    !> @brief Allocate 3D grid and coordinate mapping arrays.
    !> @details Allocates grid3d(MXNSTEP, MXNSTEP, MXNSTEP) for storing values on the mesh
    !> and cart_from_int(MXNSTEP, 3) for index-to-Cartesian transformation.
    !> Both arrays are initialized to zero.
    subroutine allocate_grid_mod()
        if (.not. allocated(grid3d)) allocate (grid3d(MXNSTEP, MXNSTEP, MXNSTEP), source=0.0_sp)
        if (.not. allocated(cart_from_int)) allocate (cart_from_int(MXNSTEP, 3), source=0.0_sp)
    end subroutine allocate_grid_mod

    !> @brief Deallocate 3D grid and coordinate mapping arrays.
    !> @details Frees memory for grid3d and cart_from_int arrays.
    subroutine deallocate_grid_mod()
        if (allocated(cart_from_int)) deallocate (cart_from_int)
        if (allocated(grid3d)) deallocate (grid3d)
    end subroutine deallocate_grid_mod

end module grid_mod

!> @brief Module for tricubic spline representation of orbitals on 3D grids.
!> @author CHAMP developers
!> @date 2022
!>
!> @details This module stores tricubic spline coefficients for numerical orbital
!> representations on 3D grids. Tricubic splines provide smooth, continuous orbital
!> interpolation with continuous first and second derivatives, essential for accurate
!> force and energy calculations.
!>
!> The spline representation uses 8 coefficients per grid cell (for the 8 vertices)
!> and covers MORB_OCC occupied orbitals.
!>
!> Key components:
!> - MORB_OCC: Number of occupied orbitals to store
!> - orb_num_spl: 5D array of spline coefficients (8, MXNSTEP, MXNSTEP, MXNSTEP, MORB_OCC)
!>
!> @note The 8 coefficients per cell correspond to tricubic interpolation polynomial terms.
module grid_spline_mod
    use precision_kinds, only: sp
    use system, only: nelec
    use grid_mod, only: MXNSTEP

    implicit none

    !> Number of occupied orbitals to store on grid
    integer :: MORB_OCC
    
    !> Tricubic spline coefficients for numerical orbitals (8 coefs, MXNSTEP^3 grid, MORB_OCC orbitals)
    real(sp), dimension(:, :, :, :, :), allocatable :: orb_num_spl

    private
    public :: MORB_OCC
    public :: orb_num_spl
    public :: allocate_grid_spline_mod, deallocate_grid_spline_mod
    save

contains
    !> @brief Allocate tricubic spline coefficient array for numerical orbitals.
    !> @details Allocates orb_num_spl(8, MXNSTEP, MXNSTEP, MXNSTEP, MORB_OCC) with
    !> 8 spline coefficients per grid cell for MORB_OCC = nelec/2 + 3 orbitals.
    !> Initialized to zero.
    subroutine allocate_grid_spline_mod()
        use system, only: nelec
        MORB_OCC = nelec/2 + 3
        if (.not. allocated(orb_num_spl)) allocate (orb_num_spl(8, MXNSTEP, MXNSTEP, MXNSTEP, MORB_OCC), source=0.0_sp)
    end subroutine allocate_grid_spline_mod

    !> @brief Deallocate tricubic spline coefficient array.
    !> @details Frees memory for orb_num_spl array.
    subroutine deallocate_grid_spline_mod()
        if (allocated(orb_num_spl)) deallocate (orb_num_spl)
    end subroutine deallocate_grid_spline_mod

end module grid_spline_mod

!> @brief Module for Lagrange polynomial interpolation of orbitals on 3D grids.
!> @author CHAMP developers
!> @date 2022
!>
!> @details This module stores data for Lagrange polynomial interpolation of
!> numerical orbitals on 3D grids. Lagrange interpolation uses LAGMAX points
!> per axis (default 4) to construct piecewise polynomial approximations.
!>
!> The interpolation stencil is centered at each grid point with indices from
!> LAGSTART to LAGEND, providing local high-order interpolation without global
!> coupling. This is computationally cheaper than splines for online evaluation.
!>
!> Key components:
!> - LAGMAX: Number of Lagrange points per axis (4-point stencil)
!> - LAGSTART, LAGEND: Index range for stencil (-2 to 1 for LAGMAX=4)
!> - MORB_OCC: Number of occupied orbitals (nelec/2)
!> - orb_num_lag: Orbital data and derivatives for Lagrange interpolation
!>
!> @note 4-point Lagrange provides cubic interpolation (degree 3).
module grid_lagrange_mod
    use precision_kinds, only: sp
    use grid_mod, only: MXNSTEP
    use system, only: nelec

    implicit none

    !> Number of Lagrange interpolation points per axis (4-point stencil)
    integer, parameter :: LAGMAX = 4
    
    !> Starting index for Lagrange stencil (-LAGMAX/2 = -2 for LAGMAX=4)
    integer, parameter :: LAGSTART = -LAGMAX/2
    
    !> Ending index for Lagrange stencil (LAGSTART + LAGMAX - 1 = 1 for LAGMAX=4)
    integer, parameter :: LAGEND = LAGSTART + LAGMAX - 1
    
    !> Number of occupied orbitals to store on grid (nelec/2 for closed-shell)
    integer ::  MORB_OCC

    !> Orbital values and derivatives for Lagrange interpolation (5, MXNSTEP, MXNSTEP, MXNSTEP, MORB_OCC)
    real(sp), dimension(:, :, :, :, :), allocatable :: orb_num_lag

    private
    public :: LAGMAX, LAGSTART, LAGEND, MORB_OCC
    public :: orb_num_lag
    public :: allocate_grid_lagrange_mod, deallocate_grid_lagrange_mod
    save
contains
    !> @brief Allocate Lagrange interpolation data array for numerical orbitals.
    !> @details Allocates orb_num_lag(5, MXNSTEP, MXNSTEP, MXNSTEP, MORB_OCC) for
    !> storing orbital value, gradients, and Laplacian at grid points for
    !> MORB_OCC = nelec/2 occupied orbitals. Initialized to zero.
    subroutine allocate_grid_lagrange_mod()
        use system, only: nelec
        use grid_mod, only: MXNSTEP
        MORB_OCC = nelec/2
        if (.not. allocated(orb_num_lag)) allocate (orb_num_lag(5, MXNSTEP, MXNSTEP, MXNSTEP, MORB_OCC), source=0.0_sp)
    end subroutine allocate_grid_lagrange_mod

    !> @brief Deallocate Lagrange interpolation data array.
    !> @details Frees memory for orb_num_lag array.
    subroutine deallocate_grid_lagrange_mod()
        if (allocated(orb_num_lag)) deallocate (orb_num_lag)
    end subroutine deallocate_grid_lagrange_mod

end module grid_lagrange_mod

!> @brief Module for 3D grid geometric parameters and spatial configuration.
!> @author CHAMP developers
!> @date 2022
!>
!> @details This module defines the geometric properties of 3D grids including
!> origin, endpoints, grid spacing, and number of steps along each axis. These
!> parameters fully specify the spatial extent and resolution of the computational
!> grid used for orbital and density evaluation.
!>
!> Grid specification:
!> - origin(3): Starting point (x0, y0, z0) in Cartesian space
!> - endpt(3): Ending point (x1, y1, z1) in Cartesian space
!> - nstep3d(3): Number of grid points along each axis
!> - step3d(3): Grid spacing Δx, Δy, Δz computed as (endpt-origin)/nstep3d
!>
module grid3d_param
      use precision_kinds, only: dp

    implicit none

    !> Grid endpoint coordinates in Cartesian space (x1, y1, z1)
    real(dp), dimension(:), allocatable :: endpt
    
    !> Number of grid steps along each axis (nx, ny, nz)
    integer, dimension(:), allocatable :: nstep3d
    
    !> Grid origin coordinates in Cartesian space (x0, y0, z0)
    real(dp), dimension(:), allocatable :: origin
    
    !> Grid spacing along each axis (Δx, Δy, Δz)
    real(dp), dimension(:), allocatable :: step3d

    private
    public :: nstep3d, endpt, origin, step3d
    public :: allocate_grid3d_param, deallocate_grid3d_param
    save
contains
    !> @brief Allocate 3D grid parameter arrays.
    !> @details Allocates origin(3), endpt(3), step3d(3), and nstep3d(3) arrays.
    !> The nstep3d array is initialized to zero; others are uninitialized.
    subroutine allocate_grid3d_param()
        if (.not. allocated(endpt)) allocate (endpt(3))
        if (.not. allocated(nstep3d)) allocate (nstep3d(3), source=0)
        if (.not. allocated(origin)) allocate (origin(3))
        if (.not. allocated(step3d)) allocate (step3d(3))
    end subroutine allocate_grid3d_param

    !> @brief Deallocate 3D grid parameter arrays.
    !> @details Frees memory for step3d, origin, nstep3d, and endpt arrays.
    subroutine deallocate_grid3d_param()
        if (allocated(step3d)) deallocate (step3d)
        if (allocated(origin)) deallocate (origin)
        if (allocated(nstep3d)) deallocate (nstep3d)
        if (allocated(endpt)) deallocate (endpt)
    end subroutine deallocate_grid3d_param

end module grid3d_param

!> @brief Module for 3D grid flags.
!> @author CHAMP developers
!> @date 2020
!>
!> @details This module contains variables for various 3D grid arrays. 
!>
!> Available variables:
!> - i3dgrid: 
!> - i3dsplorb: 
!> - i3dlagorb: 
!> - i3ddensity: 
!>
module grid3dflag

    implicit none

    !> i3ddensity
    integer :: i3ddensity
    
    !> i3dgrid
    integer :: i3dgrid
    
    !> i3dlagorb
    integer :: i3dlagorb
    
    !> i3dsplorb
    integer :: i3dsplorb

    private
    public :: i3dsplorb, i3dlagorb, i3dgrid, i3ddensity
    save
end module grid3dflag

!> @brief Module for Lagrange interpolation denominators and inverse grid spacings.
!> @author CHAMP developers
!> @date 2025
!>
!> @details This module stores precomputed quantities for efficient Lagrange
!> polynomial interpolation of numerical orbitals.
!>
!> For Lagrange interpolation L_i(x) = Π_{j≠i} (x-x_j)/(x_i-x_j), the denominators
!> are the products Π_{j≠i} (x_i-x_j) for each basis point i in the stencil.
!>
!> Key components:
!> - denom(LAGSTART:LAGEND, 3): Lagrange basis denominators for each axis
!> - step_inv(3, 3): Inverse grid spacings for coordinate transformation
!>
!> @note Precomputation improves efficiency for repeated orbital evaluations.
module orbital_num_lag
      use grid_lagrange_mod, only: LAGEND,LAGSTART
      use precision_kinds, only: dp

    implicit none

    !> Lagrange polynomial denominators for each stencil point and axis (LAGSTART:LAGEND, 3)
    real(dp), dimension(:, :), allocatable :: denom
    
    !> Inverse grid spacings for coordinate transformation (3, 3)
    real(dp), dimension(:, :), allocatable :: step_inv

    private
    public :: denom, step_inv
    public :: allocate_orbital_num_lag, deallocate_orbital_num_lag
    save
contains
    !> @brief Allocate Lagrange interpolation denominator and inverse spacing arrays.
    !> @details Allocates denom(LAGSTART:LAGEND, 3) for Lagrange basis denominators
    !> and step_inv(3, 3) for inverse grid spacings. Arrays are uninitialized.
    subroutine allocate_orbital_num_lag()
      use grid_lagrange_mod, only: LAGEND,LAGSTART
        if (.not. allocated(denom)) allocate (denom(LAGSTART:LAGEND, 3))
        if (.not. allocated(step_inv)) allocate (step_inv(3, 3))
    end subroutine allocate_orbital_num_lag

    !> @brief Deallocate Lagrange interpolation arrays.
    !> @details Frees memory for step_inv and denom arrays.
    subroutine deallocate_orbital_num_lag()
        if (allocated(step_inv)) deallocate (step_inv)
        if (allocated(denom)) deallocate (denom)
    end subroutine deallocate_orbital_num_lag

end module orbital_num_lag

!> @brief Master module for coordinating 3D grid memory allocation and deallocation.
!> @author CHAMP developers
!> @date 2022
!>
!> @details This module provides top-level routines for managing memory allocation
!> and deallocation for all grid-related data structures. It coordinates allocation
!> across multiple sub-modules to ensure consistent initialization and cleanup.
!>
!> Sub-modules managed:
!> - grid_mod: Basic 3D grid infrastructure (grid3d, cart_from_int)
!> - grid_spline_mod: Tricubic spline coefficients (orb_num_spl)
!> - grid_lagrange_mod: Lagrange interpolation data (orb_num_lag)
!> - grid3d_param: Grid geometric parameters (origin, endpt, step3d, nstep3d)
!> - orbital_num_lag: Lagrange denominators and inverse spacings (denom, step_inv)
!>
!> @note Call allocate_m_grid() before any grid operations.
!> @note Call deallocate_m_grid() to free all grid memory at program termination.
module m_grid
contains
!> @brief Allocate all grid-related data structures.
!> @details Calls allocation routines for all grid sub-modules in proper order:
!> 1. allocate_grid_mod(): Basic grid infrastructure
!> 2. allocate_grid_spline_mod(): Tricubic spline coefficients
!> 3. allocate_grid_lagrange_mod(): Lagrange interpolation data
!> 4. allocate_grid3d_param(): Grid geometric parameters
!> 5. allocate_orbital_num_lag(): Lagrange denominators
subroutine allocate_m_grid()
      use grid3d_param, only: allocate_grid3d_param
      use grid_lagrange_mod, only: allocate_grid_lagrange_mod
      use grid_mod, only: allocate_grid_mod
      use grid_spline_mod, only: allocate_grid_spline_mod
      use orbital_num_lag, only: allocate_orbital_num_lag

    implicit none

    call allocate_grid_mod()
    call allocate_grid_spline_mod()
    call allocate_grid_lagrange_mod()
    call allocate_grid3d_param()
    call allocate_orbital_num_lag()
end subroutine allocate_m_grid

!> @brief Deallocate all grid-related data structures.
!> @details Calls deallocation routines for all grid sub-modules:
!> 1. deallocate_grid_mod(): Basic grid arrays
!> 2. deallocate_grid_spline_mod(): Spline coefficients
!> 3. deallocate_grid_lagrange_mod(): Lagrange data
!> 4. deallocate_grid3d_param(): Grid parameters
!> 5. deallocate_orbital_num_lag(): Lagrange denominators
subroutine deallocate_m_grid()
      use grid3d_param, only: deallocate_grid3d_param
      use grid_lagrange_mod, only: deallocate_grid_lagrange_mod
      use grid_mod, only: deallocate_grid_mod
      use grid_spline_mod, only: deallocate_grid_spline_mod
      use orbital_num_lag, only: deallocate_orbital_num_lag

    implicit none

    call deallocate_grid_mod()
    call deallocate_grid_spline_mod()
    call deallocate_grid_lagrange_mod()
    call deallocate_grid3d_param()
    call deallocate_orbital_num_lag()
end subroutine deallocate_m_grid
end module 
