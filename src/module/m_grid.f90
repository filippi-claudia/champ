module grid_mod
    !> Arguments
    ! flags and dimensions for the 3d grid objects
      use precision_kinds, only: dp,sp

    implicit none

    integer, parameter :: MXNSTEP = 1
    ! integer, parameter :: MXNSTEP = 50
    integer, parameter :: MXNSTEP2 = MXNSTEP*MXNSTEP
    integer, parameter :: MXNSTEP3 = MXNSTEP2*MXNSTEP

    integer, parameter :: IUNDEFINED = -1234567890
    real(dp), parameter :: UNDEFINED = -1234567890.d0, SHIFT = 2.d0

    real(sp), dimension(:, :, :), allocatable :: grid3d !(MXNSTEP, MXNSTEP, MXNSTEP)
    real(sp), dimension(:, :), allocatable :: cart_from_int !(MXNSTEP, 3)

    private
    public :: MXNSTEP, MXNSTEP2, MXNSTEP3
    public :: IUNDEFINED, UNDEFINED, SHIFT
    public :: grid3d, cart_from_int
    public :: allocate_grid_mod, deallocate_grid_mod
    save

contains
    subroutine allocate_grid_mod()
        if (.not. allocated(grid3d)) allocate (grid3d(MXNSTEP, MXNSTEP, MXNSTEP), source=0.0_sp)
        if (.not. allocated(cart_from_int)) allocate (cart_from_int(MXNSTEP, 3), source=0.0_sp)
    end subroutine allocate_grid_mod

    subroutine deallocate_grid_mod()
        if (allocated(cart_from_int)) deallocate (cart_from_int)
        if (allocated(grid3d)) deallocate (grid3d)
    end subroutine deallocate_grid_mod

end module grid_mod

module grid_spline_mod
    !> Arguments
    use precision_kinds, only: sp
    use system, only: nelec
    use grid_mod, only: MXNSTEP

    implicit none

    integer :: MORB_OCC
    real(sp), dimension(:, :, :, :, :), allocatable :: orb_num_spl !(8, MXNSTEP, MXNSTEP, MXNSTEP, MORB_OCC)

    private
    public :: MORB_OCC
    public :: orb_num_spl
    public :: allocate_grid_spline_mod, deallocate_grid_spline_mod
    save

contains
    subroutine allocate_grid_spline_mod()
        use system, only: nelec
        MORB_OCC = nelec/2 + 3
        if (.not. allocated(orb_num_spl)) allocate (orb_num_spl(8, MXNSTEP, MXNSTEP, MXNSTEP, MORB_OCC), source=0.0_sp)
    end subroutine allocate_grid_spline_mod

    subroutine deallocate_grid_spline_mod()
        if (allocated(orb_num_spl)) deallocate (orb_num_spl)
    end subroutine deallocate_grid_spline_mod

end module grid_spline_mod

module grid_lagrange_mod
    !> argument
    use precision_kinds, only: sp
    use grid_mod, only: MXNSTEP
    use system, only: nelec

    implicit none

    ! Number of Lagrange interpolation points/axis
    integer, parameter :: LAGMAX = 4
    integer, parameter :: LAGSTART = -LAGMAX/2, LAGEND = LAGSTART + LAGMAX - 1
    integer ::  MORB_OCC

    !  Spline fits of the orbitals
    ! and boundary conditions (for the creation of the fit)
    real(sp), dimension(:, :, :, :, :), allocatable :: orb_num_lag !(5, MXNSTEP, MXNSTEP, MXNSTEP, MORB_OCC)

    private
    public :: LAGMAX, LAGSTART, LAGEND, MORB_OCC
    public :: orb_num_lag
    public :: allocate_grid_lagrange_mod, deallocate_grid_lagrange_mod
    save
contains
    subroutine allocate_grid_lagrange_mod()
        use system, only: nelec
        use grid_mod, only: MXNSTEP
        MORB_OCC = nelec/2
        if (.not. allocated(orb_num_lag)) allocate (orb_num_lag(5, MXNSTEP, MXNSTEP, MXNSTEP, MORB_OCC), source=0.0_sp)
    end subroutine allocate_grid_lagrange_mod

    subroutine deallocate_grid_lagrange_mod()
        if (allocated(orb_num_lag)) deallocate (orb_num_lag)
    end subroutine deallocate_grid_lagrange_mod

end module grid_lagrange_mod

module grid3d_param
    !> Arguments: nstep3d, endpt, origin, step3d
      use precision_kinds, only: dp

    implicit none

    real(dp), dimension(:), allocatable :: endpt !(3)
    integer, dimension(:), allocatable :: nstep3d !(3)
    real(dp), dimension(:), allocatable :: origin !(3)
    real(dp), dimension(:), allocatable :: step3d !(3)

    private
    public :: nstep3d, endpt, origin, step3d
    public :: allocate_grid3d_param, deallocate_grid3d_param
    save
contains
    subroutine allocate_grid3d_param()
        if (.not. allocated(endpt)) allocate (endpt(3))
        if (.not. allocated(nstep3d)) allocate (nstep3d(3), source=0)
        if (.not. allocated(origin)) allocate (origin(3))
        if (.not. allocated(step3d)) allocate (step3d(3))
    end subroutine allocate_grid3d_param

    subroutine deallocate_grid3d_param()
        if (allocated(step3d)) deallocate (step3d)
        if (allocated(origin)) deallocate (origin)
        if (allocated(nstep3d)) deallocate (nstep3d)
        if (allocated(endpt)) deallocate (endpt)
    end subroutine deallocate_grid3d_param

end module grid3d_param

module grid3dflag
    !> Arguments: i3dsplorb, i3dlagorb, i3dgrid, i3ddensity

    implicit none

    integer :: i3ddensity
    integer :: i3dgrid
    integer :: i3dlagorb
    integer :: i3dsplorb

    private
    public :: i3dsplorb, i3dlagorb, i3dgrid, i3ddensity
    save
end module grid3dflag

module orbital_num_lag
    !> Arguments: denom
      use grid_lagrange_mod, only: LAGEND,LAGSTART
      use precision_kinds, only: dp

    implicit none

    real(dp), dimension(:, :), allocatable :: denom !(LAGSTART:LAGEND, 3)
    real(dp), dimension(:, :), allocatable :: step_inv !(3, 3)

    private
    public :: denom, step_inv
    public :: allocate_orbital_num_lag, deallocate_orbital_num_lag
    save
contains
    subroutine allocate_orbital_num_lag()
      use grid_lagrange_mod, only: LAGEND,LAGSTART
        if (.not. allocated(denom)) allocate (denom(LAGSTART:LAGEND, 3))
        if (.not. allocated(step_inv)) allocate (step_inv(3, 3))
    end subroutine allocate_orbital_num_lag

    subroutine deallocate_orbital_num_lag()
        if (allocated(step_inv)) deallocate (step_inv)
        if (allocated(denom)) deallocate (denom)
    end subroutine deallocate_orbital_num_lag

end module orbital_num_lag

module m_grid
contains
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
