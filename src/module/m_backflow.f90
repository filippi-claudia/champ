!> @brief Module for backflow inplementation
!> @date November 2025
!> @author Emiel Slootman
!> @details This module manages arrays and data structures for backflow calculations
module m_backflow
    use precision_kinds, only: dp
    use system, only: nelec, ncent_tot
    use vmc_mod, only: nmat_dim, nwftypeorb, norb_tot
    
    implicit none

    !> Input parameter flag to enable backflow (default: 0 = disabled)
    integer :: ibackflow

    !> Quasicoordinate arrays (3, nelec)
    real(dp), dimension(:, :), allocatable :: quasi_x
    !> Quasicoordinate derivatives arrays (3, nelec, 3, nelec)
    real(dp), dimension(:, :, :, :), allocatable :: dquasi_dx
    !> Quasicoordinate second derivatives arrays (3, nelec, nelec)
    real(dp), dimension(:, :, :), allocatable :: d2quasi_dx2

    !> Distance vector (needed if run without qmckl)
    real(dp), dimension(:, :, :), allocatable :: rvec_en_bf
    !> Distances (needed if run without qmckl)
    real(dp), dimension(:, :), allocatable :: r_en_bf


    !> Backflow derivatives of slater matrix (3, nmat_dim,2,nwftypeorb)
    real(dp), dimension(:, :, :, :), allocatable :: dslm
    !> Backflow second derivatives of slater matrix (3,3,nmat_dim,2,nwftypeorb)
    real(dp), dimension(:,:,:, :, :), allocatable :: d2slm
 
    !> d2orb
    real(dp), dimension(:, :, :, :, :), allocatable :: d2orb

    
    private
    public :: ibackflow
    public :: quasi_x, dquasi_dx, d2quasi_dx2
    public :: rvec_en_bf, r_en_bf
    public :: allocate_m_backflow, deallocate_m_backflow
    public :: dslm, d2slm, d2orb

contains
    !> Allocates memory for backflow arrays.
    subroutine allocate_m_backflow()
      if (ibackflow.gt.0) then
        if (.not. allocated(quasi_x)) allocate (quasi_x(3, nelec))
        if (.not. allocated(dquasi_dx)) allocate (dquasi_dx(3, nelec, 3, nelec))
        if (.not. allocated(d2quasi_dx2)) allocate (d2quasi_dx2(3, nelec, nelec))
        if (.not. allocated(rvec_en_bf)) allocate (rvec_en_bf(3, nelec, ncent_tot))
        if (.not. allocated(r_en_bf)) allocate (r_en_bf(nelec, ncent_tot))
        if (.not. allocated(dslm)) allocate (dslm(3, nmat_dim, 2, nwftypeorb))
        if (.not. allocated(d2slm)) allocate (d2slm(3, 3, nmat_dim, 2, nwftypeorb))
        if (.not. allocated(d2orb)) allocate (d2orb(3,3,norb_tot, nelec, nwftypeorb))
      endif
    end subroutine allocate_m_backflow
  
    !> Deallocates memory for backflow arrays.
    subroutine deallocate_m_backflow()
      if (ibackflow.gt.0) then
        if (allocated(quasi_x)) deallocate(quasi_x)
        if (allocated(dquasi_dx)) deallocate(dquasi_dx)
        if (allocated(d2quasi_dx2)) deallocate(d2quasi_dx2)
        if (allocated(rvec_en_bf)) deallocate(rvec_en_bf)
        if (allocated(r_en_bf)) deallocate(r_en_bf)
        if (allocated(dslm)) deallocate(dslm)
        if (allocated(d2slm)) deallocate(d2slm)
        if (allocated(d2orb)) deallocate(d2orb)
      endif
    end subroutine deallocate_m_backflow
  
end module m_backflow