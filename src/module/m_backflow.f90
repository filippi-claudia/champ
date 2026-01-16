!> @brief Module for backflow inplementation
!> @date November 2025
!> @author Emiel Slootman
!> @details This module manages arrays and data structures for backflow calculations
module m_backflow
    use precision_kinds, only: dp
    use system, only: nelec, ncent_tot, nup, nctype
    use vmc_mod, only: nmat_dim, nwftypeorb, norb_tot
    use qua, only: nquad
    
    implicit none

    !> Input parameter flag to enable backflow (default: 0 = disabled)
    integer :: ibackflow

    !> Quasicoordinate arrays (3, nelec)
    real(dp), dimension(:, :), allocatable :: quasi_x
    !> Quasicoordinate derivatives arrays (3, nelec, 3, nelec)
    real(dp), dimension(:, :, :, :), allocatable :: dquasi_dx
    !> Quasicoordinate second derivatives arrays (3, nelec, nelec)
    real(dp), dimension(:, :, :), allocatable :: d2quasi_dx2
    !> Parameter derivatives of quasicoordinates (3, nelec, nparmbf)
    real(dp), dimension(:, :, :), allocatable :: dquasi_dp

      !> Quasicoordinate arrays (3, nelec)
    real(dp), dimension(:, :), allocatable :: quasi_x_new
    !> Quasicoordinate derivatives arrays (3, nelec, 3, nelec)
    real(dp), dimension(:, :, :, :), allocatable :: dquasi_dx_new
    !> Quasicoordinate second derivatives arrays (3, nelec, nelec)
    real(dp), dimension(:, :, :), allocatable :: d2quasi_dx2_new

    !> Distance vector (needed if run without qmckl)
    real(dp), dimension(:, :, :), allocatable :: rvec_en_bf
    !> Distances (needed if run without qmckl)
    real(dp), dimension(:, :), allocatable :: r_en_bf


    !> Backflow derivatives of slater matrix (3, nmat_dim,2,nwftypeorb)
    real(dp), dimension(:, :, :, :), allocatable :: dslm
      !> Backflow derivatives of slater matrix (3, nmat_dim,2,nwftypeorb)
    real(dp), dimension(:, :, :, :), allocatable :: dslm_bf
    !> Backflow second derivatives of slater matrix (3,3,nmat_dim,2,nwftypeorb)
    real(dp), dimension(:,:,:, :, :), allocatable :: d2slm
 
    !> d2orb
    real(dp), dimension(:, :, :, :, :), allocatable :: d2orb

    !> orbn
    real(dp), dimension(:, :, :), allocatable :: orbn_bf
    !> dorbn
    real(dp), dimension(:, :, :, :), allocatable :: dorbn_bf
    !> slmin
    real(dp), dimension(:, :, :), allocatable :: slmin_bf
    !> detn_bf
    real(dp), dimension(:, :), allocatable :: detn_bf

    real(dp), dimension(:, :, :), allocatable :: r_ee
    real(dp), dimension(:, :, :, :), allocatable :: r_ee_gl
    real(dp), dimension(:, :, :), allocatable :: rvec_ee
    real(dp), dimension(:, :, :, :), allocatable :: r_en
    real(dp), dimension(:, :, :, :, :), allocatable :: r_en_gl
    real(dp), dimension(:, :, :), allocatable :: rvec_en
    real(dp), dimension(:, :), allocatable :: cutoff_deriv
    real(dp), dimension(:, :, :, :), allocatable :: p
    real(dp), dimension(:, :, :, :, :), allocatable :: d_p

    !> nl_slm (nmat_dim,2)
    real(dp), dimension(:, :), allocatable :: nl_slm

    ! Array for cusp dependency (phi, theta)
    real(dp), dimension(:, :, :), allocatable :: B
    real(dp), dimension(:, :), allocatable :: cusp_parameters, inv_cusp_parameters
    integer, dimension(:, :), allocatable :: cusp_indices, inv_cusp_indices
    ! Store (k, l, m) indices for each basis function column
    integer, dimension(:, :), allocatable :: basis_klm

    real(dp), dimension(:), allocatable :: parm_bf
    real(dp), dimension(:), allocatable :: deriv_parm_bf
    integer :: norda_bf, nordb_bf, nordc_bf, nparm_bf, cutoff_scale, maxord, ncparm_bf, c_cuspconst

    
    private
    public :: ibackflow
    public :: quasi_x, dquasi_dx, d2quasi_dx2, dquasi_dp
    public :: rvec_en_bf, r_en_bf
    public :: allocate_m_backflow, deallocate_m_backflow
    public :: dslm, d2slm, d2orb, nl_slm, nparm_bf, parm_bf, deriv_parm_bf, dslm_bf
    public :: orbn_bf, dorbn_bf, slmin_bf, detn_bf, norda_bf, nordb_bf, nordc_bf, cutoff_scale, ncparm_bf, c_cuspconst
    public :: quasi_x_new, dquasi_dx_new, d2quasi_dx2_new, maxord
    public :: r_ee, rvec_ee, r_en, rvec_en, r_ee_gl, r_en_gl, p, d_p, cutoff_deriv
    public :: B, cusp_parameters, cusp_indices, inv_cusp_parameters, inv_cusp_indices, basis_klm

contains
    !> Allocates memory for backflow arrays.
    subroutine allocate_m_backflow()
      if (ibackflow.gt.0) then
        if (.not. allocated(quasi_x)) allocate (quasi_x(3, nelec))
        if (.not. allocated(dquasi_dx)) allocate (dquasi_dx(3, nelec, 3, nelec))
        if (.not. allocated(d2quasi_dx2)) allocate (d2quasi_dx2(3, nelec, nelec))
        if (.not. allocated(rvec_en_bf)) allocate (rvec_en_bf(3, nelec, ncent_tot))
        if (.not. allocated(r_en_bf)) allocate (r_en_bf(nelec, ncent_tot))
        if (.not. allocated(dslm)) allocate (dslm(3, nup*nup, 2, nwftypeorb))
        if (.not. allocated(d2slm)) allocate (d2slm(3, 3, nup*nup, 2, nwftypeorb))
        if (.not. allocated(d2orb)) allocate (d2orb(3,3,norb_tot, nelec, nwftypeorb))
        if (.not. allocated(nl_slm)) allocate (nl_slm(nup*nup, 2))
        if (.not. allocated(orbn_bf)) allocate (orbn_bf(nelec, norb_tot, nwftypeorb))
        if (.not. allocated(dorbn_bf)) allocate (dorbn_bf(norb_tot, nelec, 3, nwftypeorb))
        if (.not. allocated(slmin_bf)) allocate (slmin_bf(nup*nup,2,nwftypeorb))
        if (.not. allocated(detn_bf)) allocate (detn_bf(2, nwftypeorb))
        if (.not. allocated(parm_bf)) allocate (parm_bf(nparm_bf))
        if (.not. allocated(dquasi_dp)) allocate (dquasi_dp(3, nelec, nparm_bf))
        if (.not. allocated(deriv_parm_bf)) allocate (deriv_parm_bf(nparm_bf))
        if (.not. allocated(quasi_x_new)) allocate (quasi_x_new(3, nelec))
        if (.not. allocated(dquasi_dx_new)) allocate (dquasi_dx_new(3, nelec, 3, nelec))
        if (.not. allocated(d2quasi_dx2_new)) allocate (d2quasi_dx2_new(3, nelec, nelec))
        if (.not. allocated(dslm_bf)) allocate (dslm_bf(3, nup*nup, 2, nwftypeorb))
        if (.not. allocated(r_ee)) allocate (r_ee(nelec, nelec, 0:maxord))
        if (.not. allocated(rvec_ee)) allocate (rvec_ee(3, nelec, nelec))
        if (.not. allocated(r_en)) allocate (r_en(nelec, ncent_tot, 0:maxord, 2))
        if (.not. allocated(rvec_en)) allocate (rvec_en(3, nelec, ncent_tot))
        if (.not. allocated(r_ee_gl)) allocate (r_ee_gl(nelec, 4, nelec, 0:maxord))
        if (.not. allocated(r_en_gl)) allocate (r_en_gl(nelec, 4, ncent_tot, 0:maxord, 2))
        if (.not. allocated(p)) allocate (p(nelec, ncent_tot, nordc_bf, nordc_bf))
        if (.not. allocated(d_p)) allocate (d_p(nelec, 4, ncent_tot, nordc_bf, nordc_bf))
        if (.not. allocated(cutoff_deriv)) allocate (cutoff_deriv(nelec, ncent_tot))
        if (.not. allocated(B)) allocate (B(c_cuspconst, ncparm_bf,2*nctype))
        if (.not. allocated(cusp_parameters)) allocate (cusp_parameters(c_cuspconst*nctype, ncparm_bf))
        if (.not. allocated(cusp_indices)) allocate (cusp_indices(c_cuspconst*nctype, ncparm_bf))
        if (.not. allocated(inv_cusp_parameters)) allocate (inv_cusp_parameters(nparm_bf, ncparm_bf))
        if (.not. allocated(inv_cusp_indices)) allocate (inv_cusp_indices(nparm_bf, ncparm_bf))
        if (.not. allocated(basis_klm)) allocate (basis_klm(ncparm_bf, 3))
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
        if (allocated(nl_slm)) deallocate(nl_slm)
        if (allocated(orbn_bf)) deallocate(orbn_bf)
        if (allocated(dorbn_bf)) deallocate(dorbn_bf)
        if (allocated(slmin_bf)) deallocate(slmin_bf)
        if (allocated(detn_bf)) deallocate(detn_bf)
        if (allocated(parm_bf)) deallocate(parm_bf)
        if (allocated(dquasi_dp)) deallocate(dquasi_dp)
        if (allocated(deriv_parm_bf)) deallocate(deriv_parm_bf)
        if (allocated(quasi_x_new)) deallocate(quasi_x_new)
        if (allocated(dquasi_dx_new)) deallocate(dquasi_dx_new)
        if (allocated(d2quasi_dx2_new)) deallocate(d2quasi_dx2_new)
        if (allocated(dslm_bf)) deallocate(dslm_bf)
        if (allocated(r_ee)) deallocate(r_ee)
        if (allocated(rvec_ee)) deallocate(rvec_ee)
        if (allocated(r_en)) deallocate(r_en)
        if (allocated(rvec_en)) deallocate(rvec_en)
        if (allocated(r_ee_gl)) deallocate(r_ee_gl)
        if (allocated(r_en_gl)) deallocate(r_en_gl)
        if (allocated(p)) deallocate(p)
        if (allocated(d_p)) deallocate(d_p)
        if (allocated(cutoff_deriv)) deallocate(cutoff_deriv)
        if (allocated(B)) deallocate(B)
        if (allocated(cusp_parameters)) deallocate(cusp_parameters)
        if (allocated(cusp_indices)) deallocate(cusp_indices)
        if (allocated(inv_cusp_parameters)) deallocate(inv_cusp_parameters)
        if (allocated(inv_cusp_indices)) deallocate(inv_cusp_indices)
        if (allocated(basis_klm)) deallocate(basis_klm)
      endif
    end subroutine deallocate_m_backflow
  
end module m_backflow