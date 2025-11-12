!> @brief Module for combined gradient, Hessian, and overlap matrices in wavefunction optimization.
!> @author CHAMP developers
!> @date 2020
!>
!> @details This module stores the complete gradient vector, Hessian matrix (h), and
!> overlap matrix (s) for all optimizable wavefunction parameters combined. This includes
!> Jastrow parameters, CI coefficients, and orbital parameters.
!>
!> The total parameter count nparmall = nparmj + mxcireduced + mxreduced, where:
!> - nparmj: Number of Jastrow parameters
!> - mxcireduced: Number of CI coefficients
!> - mxreduced: Number of orbital parameters
!> - For linear method optimization, nparmall and nparmjp1 are incremented by 1
!>
module gradhess_all
      use optwf_parms, only: nparmj
      use precision_kinds, only: dp

    implicit none

    !> Total number of all optimizable parameters (Jastrow + CI + orbital)
    integer :: nparmall
    
    !> Number of Jastrow parameters plus 1 for linear method
    integer :: nparmjp1
    
    !> Gradient vector for all parameters (nparmall)
    real(dp), dimension(:), allocatable :: grad
    
    !> Hessian matrix for all parameters (nparmall, nparmall)
    real(dp), dimension(:, :), allocatable :: h
    
    !> Overlap matrix for all parameters (nparmall, nparmall)
    real(dp), dimension(:, :), allocatable :: s

    private
    public :: nparmall, nparmjp1, grad, h, s
    public :: allocate_gradhess_all, deallocate_gradhess_all, set_gradhess_all_size
    save
contains
    !> @brief Set the size of gradient and Hessian arrays based on optimization type.
    !> @details Computes nparmall = nparmj + mxcireduced + mxreduced, the total number
    !> of parameters to optimize. For linear method optimization with orbital or Jastrow
    !> parameters, both nparmall and nparmjp1 are incremented by 1.
    !>
    subroutine set_gradhess_all_size()
      use optwf_control, only: method, ioptjas, ioptorb
      use optci, only: mxcireduced
      use optorb_mod, only: mxreduced
      use optwf_parms, only: nparmj

      nparmall = nparmj + mxcireduced + mxreduced
      nparmjp1 = nparmj
      if(method.eq.'linear'.and.(ioptorb.gt.0.or.ioptjas.gt.0)) then
        nparmall=nparmall+1
        nparmjp1=nparmjp1+1
      endif

    end subroutine set_gradhess_all_size

    !> @brief Allocate gradient vector, Hessian, and overlap matrices for all parameters.
    !> @details Allocates grad(nparmall), h(nparmall, nparmall), and s(nparmall, nparmall).
    !> Arrays are uninitialized. Must call set_gradhess_all_size() first.
    subroutine allocate_gradhess_all()
        if (.not. allocated(grad)) allocate (grad(nparmall))
        if (.not. allocated(h)) allocate (h(nparmall, nparmall))
        if (.not. allocated(s)) allocate (s(nparmall, nparmall))
    end subroutine allocate_gradhess_all

    !> @brief Deallocate gradient, Hessian, and overlap matrices.
    !> @details Frees memory for s, h, and grad arrays.
    subroutine deallocate_gradhess_all()
        if (allocated(s)) deallocate (s)
        if (allocated(h)) deallocate (h)
        if (allocated(grad)) deallocate (grad)
    end subroutine deallocate_gradhess_all

end module gradhess_all

!> @brief Module for Jastrow parameter gradient and Hessian components.
!> @author CHAMP developers
!> @date 2022
!>
!> @details This module stores first and second derivatives of the Jastrow factor (J) 
!> with respect to Jastrow parameters. These quantities are building
!> blocks for constructing the gradient vector and Hessian matrix used in optimization.
!>
!> Notation:
!> - dj: First derivative of Jastrow factor
!> - d2j: Second derivative of Jastrow factor
!> - de: 
!> - de_de: 
!> - dj_e: 
!> - dj_de: 
!> - e2: 
!>
!> All quantities are computed for MSTATES states (usually 1 for ground state).
!>
module gradhessj
    use mstates_mod, only: MSTATES
    use optwf_parms, only: nparmj
    use precision_kinds, only: dp

    implicit none

    !> Second derivative of Jastrow factor (nparmj, nparmj, MSTATES)
    real(dp), dimension(:, :, :), allocatable :: d2j
    
    !> d2j_e  (nparmj, nparmj, MSTATES).
    real(dp), dimension(:, :, :), allocatable :: d2j_e
    
    !> de  (nparmj, MSTATES).
    real(dp), dimension(:, :), allocatable :: de
    
    !> de_de (nparmj, nparmj, MSTATES).
    real(dp), dimension(:, :, :), allocatable :: de_de
    
    !> de_e (nparmj, MSTATES).
    real(dp), dimension(:, :), allocatable :: de_e
    
    !> First derivative of Jastrow factor (nparmj, MSTATES)
    real(dp), dimension(:, :), allocatable :: dj
    
    !> dj_de (nparmj, nparmj, MSTATES)
    real(dp), dimension(:, :, :), allocatable :: dj_de
    
    !> dj_dj (nparmj, nparmj, MSTATES)
    real(dp), dimension(:, :, :), allocatable :: dj_dj
    
    !> dj_dj_e (nparmj, nparmj, MSTATES)
    real(dp), dimension(:, :, :), allocatable :: dj_dj_e
    
    !> dj_e (nparmj, MSTATES)
    real(dp), dimension(:, :), allocatable :: dj_e
    
    !> dj_e2 (nparmj, MSTATES)
    real(dp), dimension(:, :), allocatable :: dj_e2
    
    !> e2 (nparmj, MSTATES)
    real(dp), dimension(:, :), allocatable :: e2

    private
    public :: d2j, d2j_e, de, de_de, de_e, dj, dj_de, dj_dj, dj_dj_e, dj_e, dj_e2, e2
    public :: allocate_gradhessj, deallocate_gradhessj
    save
contains
    !> @brief Allocate Jastrow gradient and Hessian component arrays.
    !> @details Allocates 12 arrays for storing first/second derivatives and products
    !> of Jastrow factor and energy with respect to Jastrow parameters. All arrays
    !> dimensioned for nparmj parameters and MSTATES states. Arrays are uninitialized.
    subroutine allocate_gradhessj()
      use mstates_mod, only: MSTATES
      use optwf_parms, only: nparmj
        if (.not. allocated(d2j)) allocate (d2j(nparmj, nparmj, MSTATES))
        if (.not. allocated(d2j_e)) allocate (d2j_e(nparmj, nparmj, MSTATES))
        if (.not. allocated(de)) allocate (de(nparmj, MSTATES))
        if (.not. allocated(de_de)) allocate (de_de(nparmj, nparmj, MSTATES))
        if (.not. allocated(de_e)) allocate (de_e(nparmj, MSTATES))
        if (.not. allocated(dj)) allocate (dj(nparmj, MSTATES))
        if (.not. allocated(dj_de)) allocate (dj_de(nparmj, nparmj, MSTATES))
        if (.not. allocated(dj_dj)) allocate (dj_dj(nparmj, nparmj, MSTATES))
        if (.not. allocated(dj_dj_e)) allocate (dj_dj_e(nparmj, nparmj, MSTATES))
        if (.not. allocated(dj_e)) allocate (dj_e(nparmj, MSTATES))
        if (.not. allocated(dj_e2)) allocate (dj_e2(nparmj, MSTATES))
        if (.not. allocated(e2)) allocate (e2(nparmj, MSTATES))
    end subroutine allocate_gradhessj

    !> @brief Deallocate Jastrow gradient and Hessian component arrays.
    !> @details Frees memory for all 12 derivative and product arrays.
    subroutine deallocate_gradhessj()
        if (allocated(e2)) deallocate (e2)
        if (allocated(dj_e2)) deallocate (dj_e2)
        if (allocated(dj_e)) deallocate (dj_e)
        if (allocated(dj_dj_e)) deallocate (dj_dj_e)
        if (allocated(dj_dj)) deallocate (dj_dj)
        if (allocated(dj_de)) deallocate (dj_de)
        if (allocated(dj)) deallocate (dj)
        if (allocated(de_e)) deallocate (de_e)
        if (allocated(de_de)) deallocate (de_de)
        if (allocated(de)) deallocate (de)
        if (allocated(d2j_e)) deallocate (d2j_e)
        if (allocated(d2j)) deallocate (d2j)
    end subroutine deallocate_gradhessj

end module gradhessj

!> @brief Module for storing old Jastrow gradient values for finite-difference calculations.
!> @author CHAMP developers
!> @date 2022
!>
!> @details This module stores previous values of Jastrow derivatives and energies
!> needed for numerical finite-difference gradient calculations. By comparing current
!> and old values after small parameter perturbations, numerical derivatives can be
!> computed.
!>
!> Key stored quantities:
!> - denergy_old: Previous energy derivatives for comparison
!> - gvalue_old: Previous Jastrow gradient values
!> - d1d2a_old, d1d2b_old: Previous first-derivative intermediates (type a/b)
!> - d2d2a_old, d2d2b_old: Previous second-derivative intermediates (type a/b)
!>
!> The type 'a' arrays are dimensioned by nctype_tot (number of atom types),
!> while type 'b' arrays have dimension 2 for electron-electron terms.
!>
!> @note Used in finite-difference optimization methods.
!> @note "_old" suffix indicates these store previous iteration values.
module gradhessjo
    use mstates_mod, only: MSTATES
    use optwf_parms, only: nparmj
    use precision_kinds, only: dp

    implicit none

    !> Previous first-derivative intermediate values for atom types (nctype_tot)
    real(dp), dimension(:), allocatable :: d1d2a_old
    
    !> Previous first-derivative intermediate values for electron pairs (2)
    real(dp), dimension(:), allocatable :: d1d2b_old
    
    !> Previous second-derivative intermediate values for atom types (nctype_tot)
    real(dp), dimension(:), allocatable :: d2d2a_old
    
    !> Previous second-derivative intermediate values for electron pairs (2)
    real(dp), dimension(:), allocatable :: d2d2b_old
    
    !> Previous energy derivatives (nparmj, MSTATES)
    real(dp), dimension(:, :), allocatable :: denergy_old
    
    !> Previous Jastrow gradient values (nparmj, nwftypejas)
    real(dp), dimension(:, :), allocatable :: gvalue_old

    private
    public   ::  d1d2a_old, d1d2b_old, d2d2a_old, d2d2b_old, denergy_old, gvalue_old
    public :: allocate_gradhessjo, deallocate_gradhessjo
    save
contains
    !> @brief Allocate arrays for storing old Jastrow gradient values.
    !> @details Allocates 6 arrays for storing previous values of Jastrow derivatives
    !> and energies: d1d2a_old(nctype_tot), d1d2b_old(2), d2d2a_old(nctype_tot),
    !> d2d2b_old(2), denergy_old(nparmj, MSTATES), gvalue_old(nparmj, nwftypejas).
    !> Arrays are uninitialized.
    subroutine allocate_gradhessjo()
      use mstates_mod, only: MSTATES
      use optwf_parms, only: nparmj
      use system, only: nctype_tot
      use vmc_mod, only: nwftypejas
        if (.not. allocated(d1d2a_old)) allocate (d1d2a_old(nctype_tot))
        if (.not. allocated(d1d2b_old)) allocate (d1d2b_old(2))
        if (.not. allocated(d2d2a_old)) allocate (d2d2a_old(nctype_tot))
        if (.not. allocated(d2d2b_old)) allocate (d2d2b_old(2))
        if (.not. allocated(denergy_old)) allocate (denergy_old(nparmj, MSTATES))
        if (.not. allocated(gvalue_old)) allocate (gvalue_old(nparmj,nwftypejas))
    end subroutine allocate_gradhessjo

    !> @brief Deallocate old Jastrow gradient value arrays.
    !> @details Frees memory for all 6 old-value arrays.
    subroutine deallocate_gradhessjo()
        if (allocated(gvalue_old)) deallocate (gvalue_old)
        if (allocated(denergy_old)) deallocate (denergy_old)
        if (allocated(d2d2b_old)) deallocate (d2d2b_old)
        if (allocated(d2d2a_old)) deallocate (d2d2a_old)
        if (allocated(d1d2b_old)) deallocate (d1d2b_old)
        if (allocated(d1d2a_old)) deallocate (d1d2a_old)
    end subroutine deallocate_gradhessjo

end module gradhessjo

!> @brief Module for CI coefficient gradient, Hessian, and overlap matrices.
!> @author CHAMP developers
!> @date 2022
!>
!> @details This module stores the gradient vector, Hessian matrix, and overlap
!> matrix specifically for Configuration Interaction (CI) coefficient optimization.
!> These matrices are used to optimize the linear combination of determinants
!> in multi-determinant wavefunctions.
!>
!> Dimensions:
!> - mxciterm: Total number of CI terms in the expansion
!> - mxcireduced: Number of independent (non-redundant) CI parameters to optimize
!>
!> The CI coefficients c_i multiply Slater determinants |D_i⟩ in the wavefunction.
!>
module gradhess_ci
      use optci,   only: mxcireduced,mxciterm
      use precision_kinds, only: dp

    implicit none

    !> Gradient vector for CI coefficients (mxciterm)
    real(dp), dimension(:), allocatable :: grad_ci
    
    !> Hessian matrix for CI coefficients (mxciterm, mxcireduced)
    real(dp), dimension(:, :), allocatable :: h_ci
    
    !> Overlap matrix for CI coefficients (mxciterm, mxcireduced)
    real(dp), dimension(:, :), allocatable :: s_ci

    private
    public   ::  grad_ci, h_ci, s_ci
    public :: allocate_gradhess_ci, deallocate_gradhess_ci
    save
contains
    !> @brief Allocate CI coefficient gradient, Hessian, and overlap matrices.
    !> @details Allocates grad_ci(mxciterm), h_ci(mxciterm, mxcireduced), and
    !> s_ci(mxciterm, mxcireduced) for CI optimization. Arrays are uninitialized.
    subroutine allocate_gradhess_ci()
      use optci, only: mxciterm, mxcireduced
        if (.not. allocated(grad_ci)) allocate (grad_ci(mxciterm))
        if (.not. allocated(h_ci)) allocate (h_ci(mxciterm, mxcireduced))
        if (.not. allocated(s_ci)) allocate (s_ci(mxciterm, mxcireduced))
    end subroutine allocate_gradhess_ci

    !> @brief Deallocate CI coefficient matrices.
    !> @details Frees memory for s_ci, h_ci, and grad_ci arrays.
    subroutine deallocate_gradhess_ci()
        if (allocated(s_ci)) deallocate (s_ci)
        if (allocated(h_ci)) deallocate (h_ci)
        if (allocated(grad_ci)) deallocate (grad_ci)
    end subroutine deallocate_gradhess_ci

end module gradhess_ci

!> @brief Module for Jastrow parameter gradient, Hessian, and overlap matrices.
!> @author CHAMP developers
!> @date 2022
!>
!> @details This module stores the gradient vector, Hessian matrix, and overlap
!> matrix specifically for Jastrow factor optimization. The Jastrow factor J(r)
!> introduces explicit electron-electron, electron-nucleus, and electron-electron-nucleus
!> correlation into the wavefunction.
!>
!> The optimizable Jastrow parameters control polynomial/spline coefficients in:
!> J = J_en(r_iα) + J_ee(r_ij) + J_een(r_iα, r_jα, r_ij)
!>
!>
!> @note Dimension nparmjp1 = nparmj or nparmj+1 depending on linear method.
module gradhess_jas
      use optwf_parms, only: nparmj
      use precision_kinds, only: dp

    implicit none

    !> Gradient vector for Jastrow parameters (nparmjp1)
    real(dp), dimension(:), allocatable :: grad_jas
    
    !> Hessian matrix for Jastrow parameters (nparmjp1, nparmjp1)
    real(dp), dimension(:, :), allocatable :: h_jas
    
    !> Overlap matrix for Jastrow parameters (nparmjp1, nparmjp1)
    real(dp), dimension(:, :), allocatable :: s_jas

    private
    public   ::  grad_jas, h_jas, s_jas
    public :: allocate_gradhess_jas, deallocate_gradhess_jas
    save
contains
    !> @brief Allocate Jastrow gradient, Hessian, and overlap matrices.
    !> @details Allocates grad_jas(nparmjp1), h_jas(nparmjp1, nparmjp1), and
    !> s_jas(nparmjp1, nparmjp1) using nparmjp1 from gradhess_all module.
    !> Arrays are uninitialized.
    subroutine allocate_gradhess_jas()
    use gradhess_all, only: nparmjp1

        if (.not. allocated(grad_jas)) allocate (grad_jas(nparmjp1))
        if (.not. allocated(h_jas)) allocate (h_jas(nparmjp1, nparmjp1))
        if (.not. allocated(s_jas)) allocate (s_jas(nparmjp1, nparmjp1))
    end subroutine allocate_gradhess_jas

    !> @brief Deallocate Jastrow matrices.
    !> @details Frees memory for s_jas, h_jas, and grad_jas arrays.
    subroutine deallocate_gradhess_jas()
        if (allocated(s_jas)) deallocate (s_jas)
        if (allocated(h_jas)) deallocate (h_jas)
        if (allocated(grad_jas)) deallocate (grad_jas)
    end subroutine deallocate_gradhess_jas

end module gradhess_jas

!> @brief Module for mixed Jastrow-CI Hessian and overlap matrix blocks.
!> @author CHAMP developers
!> @date 2022
!>
!> @details This module stores the off-diagonal blocks of Hessian and overlap
!> matrices coupling Jastrow parameters with CI coefficients. These are needed
!> for simultaneous optimization of both parameter types.
!>
module gradhess_mix_jas_ci
    use optci, only: mxciterm
    use optwf_parms, only: nparmj
    use precision_kinds, only: dp

    implicit none

    !> Mixed Hessian block coupling Jastrow and CI (2*nparmj, mxciterm)
    real(dp), dimension(:, :), allocatable :: h_mix_jas_ci
    
    !> Mixed overlap block coupling Jastrow and CI (nparmj, mxciterm)
    real(dp), dimension(:, :), allocatable :: s_mix_jas_ci

    private
    public   ::  h_mix_jas_ci, s_mix_jas_ci
    public :: allocate_gradhess_mix_jas_ci, deallocate_gradhess_mix_jas_ci
    save
contains
    !> @brief Allocate mixed Jastrow-CI Hessian and overlap matrices.
    !> @details Allocates h_mix_jas_ci(2*nparmj, mxciterm) and
    !> s_mix_jas_ci(nparmj, mxciterm) for simultaneous Jastrow+CI optimization.
    !> Arrays are uninitialized.
    subroutine allocate_gradhess_mix_jas_ci()
      use optci, only: mxciterm
      use optwf_parms, only: nparmj
        if (.not. allocated(h_mix_jas_ci)) allocate (h_mix_jas_ci(2*nparmj, mxciterm))
        if (.not. allocated(s_mix_jas_ci)) allocate (s_mix_jas_ci(nparmj, mxciterm))
    end subroutine allocate_gradhess_mix_jas_ci

    !> @brief Deallocate mixed Jastrow-CI matrices.
    !> @details Frees memory for s_mix_jas_ci and h_mix_jas_ci arrays.
    subroutine deallocate_gradhess_mix_jas_ci()
        if (allocated(s_mix_jas_ci)) deallocate (s_mix_jas_ci)
        if (allocated(h_mix_jas_ci)) deallocate (h_mix_jas_ci)
    end subroutine deallocate_gradhess_mix_jas_ci

end module gradhess_mix_jas_ci

!> @brief Module for mixed Jastrow-orbital Hessian and overlap matrix blocks.
!> @author CHAMP developers
!> @date 2022
!>
!> @details This module stores the off-diagonal blocks of Hessian and overlap
!> matrices coupling Jastrow parameters with orbital parameters. These are needed
!> for simultaneous optimization of Jastrow factor and orbital coefficients.
!>
module gradhess_mix_jas_orb
      use optorb_mod, only: mxreduced
      use optwf_parms, only: nparmj
      use precision_kinds, only: dp

    implicit none

    !> Mixed Hessian block coupling Jastrow and orbitals (2*nparmj, mxreduced)
    real(dp), dimension(:, :), allocatable :: h_mix_jas_orb
    
    !> Mixed overlap block coupling Jastrow and orbitals (nparmj, mxreduced)
    real(dp), dimension(:, :), allocatable :: s_mix_jas_orb

    private
    public   ::  h_mix_jas_orb, s_mix_jas_orb
    public :: allocate_gradhess_mix_jas_orb, deallocate_gradhess_mix_jas_orb
    save
contains
    !> @brief Allocate mixed Jastrow-orbital Hessian and overlap matrices.
    !> @details Allocates h_mix_jas_orb(2*nparmj, mxreduced) and
    !> s_mix_jas_orb(nparmj, mxreduced) for simultaneous Jastrow+orbital optimization.
    !> Arrays are uninitialized.
    subroutine allocate_gradhess_mix_jas_orb()
      use optorb_mod, only: mxreduced
      use optwf_parms, only: nparmj
        if (.not. allocated(h_mix_jas_orb)) allocate (h_mix_jas_orb(2*nparmj, mxreduced))
        if (.not. allocated(s_mix_jas_orb)) allocate (s_mix_jas_orb(nparmj, mxreduced))
    end subroutine allocate_gradhess_mix_jas_orb

    !> @brief Deallocate mixed Jastrow-orbital matrices.
    !> @details Frees memory for s_mix_jas_orb and h_mix_jas_orb arrays.
    subroutine deallocate_gradhess_mix_jas_orb()
        if (allocated(s_mix_jas_orb)) deallocate (s_mix_jas_orb)
        if (allocated(h_mix_jas_orb)) deallocate (h_mix_jas_orb)
    end subroutine deallocate_gradhess_mix_jas_orb

end module gradhess_mix_jas_orb

!> @brief Module for mixed CI-orbital Hessian and overlap matrix blocks.
!> @author CHAMP developers
!> @date 2025
!>
!> @details This module stores the off-diagonal blocks of Hessian and overlap
!> matrices coupling CI coefficients with orbital parameters. These are needed
!> for simultaneous optimization of multi-determinant expansion coefficients
!> and the molecular orbitals themselves.
!>
module gradhess_mix_orb_ci
    use optci, only: mxciterm
    use optorb_mod, only: mxreduced
    use precision_kinds, only: dp

    implicit none

    !> Mixed Hessian block coupling CI and orbitals (2*mxciterm, mxreduced)
    real(dp), dimension(:, :), allocatable :: h_mix_ci_orb
    
    !> Mixed overlap block coupling CI and orbitals (mxciterm, mxreduced)
    real(dp), dimension(:, :), allocatable :: s_mix_ci_orb

    private
    public   ::  h_mix_ci_orb, s_mix_ci_orb
    public :: allocate_gradhess_mix_orb_ci, deallocate_gradhess_mix_orb_ci
    save
contains
    !> @brief Allocate mixed CI-orbital Hessian and overlap matrices.
    !> @details Allocates h_mix_ci_orb(2*mxciterm, mxreduced) and
    !> s_mix_ci_orb(mxciterm, mxreduced) for simultaneous CI+orbital optimization.
    !> Arrays are uninitialized.
    subroutine allocate_gradhess_mix_orb_ci()
      use optci, only: mxciterm
      use optorb_mod, only: mxreduced
        if (.not. allocated(h_mix_ci_orb)) allocate (h_mix_ci_orb(2*mxciterm, mxreduced))
        if (.not. allocated(s_mix_ci_orb)) allocate (s_mix_ci_orb(mxciterm, mxreduced))
    end subroutine allocate_gradhess_mix_orb_ci

    !> @brief Deallocate mixed CI-orbital matrices.
    !> @details Frees memory for s_mix_ci_orb and h_mix_ci_orb arrays.
    subroutine deallocate_gradhess_mix_orb_ci()
        if (allocated(s_mix_ci_orb)) deallocate (s_mix_ci_orb)
        if (allocated(h_mix_ci_orb)) deallocate (h_mix_ci_orb)
    end subroutine deallocate_gradhess_mix_orb_ci

end module gradhess_mix_orb_ci

!> @brief Module for Jastrow gradient error estimation and block averaging.
!> @author CHAMP developers
!> @date 2022
!>
!> @details This module stores quantities for computing statistical errors on
!> Jastrow gradients using block averaging. Block averaging divides the Monte Carlo
!> run into blocks and computes statistics within and between blocks to estimate
!> the standard error.
!>
!> Key arrays for block statistics:
!> - dj_bsum, dj_e_bsum: 
!> - dj_save, dj_e_save: 
!> - grad_jas_bcum: Cumulative gradient sum across blocks
!> - grad_jas_bcm2: Cumulative sum of squared gradients for variance
!> - e_bsum: 
!>
!>
module gradjerr
    use mstates_mod, only: MSTATES
    use optwf_parms, only: nparmj
    use precision_kinds, only: dp

    implicit none

    !> dj_bsum (nparmj, MSTATES)
    real(dp), dimension(:, :), allocatable :: dj_bsum
    
    !> dj_e_bsum (nparmj, MSTATES)
    real(dp), dimension(:, :), allocatable :: dj_e_bsum
    
    !> dj_e_save (nparmj, MSTATES)
    real(dp), dimension(:, :), allocatable :: dj_e_save
    
    !> dj_save (nparmj, MSTATES)
    real(dp), dimension(:, :), allocatable :: dj_save
    
    !> e_bsum (MSTATES)
    real(dp), dimension(:),    allocatable :: e_bsum
    
    !> grad_jas_bcm2 (nparmj, MSTATES)
    real(dp), dimension(:, :), allocatable :: grad_jas_bcm2
    
    !> grad_jas_bcum (nparmj, MSTATES)
    real(dp), dimension(:, :), allocatable :: grad_jas_bcum

    private
    public   ::  dj_bsum, dj_e_bsum, dj_e_save, dj_save, e_bsum, grad_jas_bcm2, grad_jas_bcum
    public :: allocate_gradjerr, deallocate_gradjerr
    save
contains
    !> @brief Allocate arrays for Jastrow gradient error estimation.
    !> @details Allocates 7 arrays for block averaging and error estimation:
    !> dj_bsum, dj_e_bsum, dj_e_save, dj_save (nparmj, MSTATES),
    !> e_bsum(MSTATES), grad_jas_bcm2, grad_jas_bcum (nparmj, MSTATES).
    !> Arrays are uninitialized.
    subroutine allocate_gradjerr()
      use mstates_mod, only: MSTATES
      use optwf_parms, only: nparmj
        if (.not. allocated(dj_bsum)) allocate (dj_bsum(nparmj, MSTATES))
        if (.not. allocated(dj_e_bsum)) allocate (dj_e_bsum(nparmj, MSTATES))
        if (.not. allocated(dj_e_save)) allocate (dj_e_save(nparmj, MSTATES))
        if (.not. allocated(dj_save)) allocate (dj_save(nparmj, MSTATES))
        if (.not. allocated(e_bsum)) allocate (e_bsum(MSTATES))
        if (.not. allocated(grad_jas_bcm2)) allocate (grad_jas_bcm2(nparmj, MSTATES))
        if (.not. allocated(grad_jas_bcum)) allocate (grad_jas_bcum(nparmj, MSTATES))
    end subroutine allocate_gradjerr

    !> @brief Deallocate gradient error estimation arrays.
    !> @details Frees memory for all 7 block averaging and error arrays.
    subroutine deallocate_gradjerr()
        if (allocated(grad_jas_bcum)) deallocate (grad_jas_bcum)
        if (allocated(grad_jas_bcm2)) deallocate (grad_jas_bcm2)
        if (allocated(e_bsum)) deallocate (e_bsum)
        if (allocated(dj_save)) deallocate (dj_save)
        if (allocated(dj_e_save)) deallocate (dj_e_save)
        if (allocated(dj_e_bsum)) deallocate (dj_e_bsum)
        if (allocated(dj_bsum)) deallocate (dj_bsum)
    end subroutine deallocate_gradjerr

end module gradjerr

!> @brief Master module for coordinating gradient and Hessian memory management.
!> @author CHAMP developers
!> @date 2022
!>
!> @details This module provides top-level routines for managing memory allocation
!> and deallocation for all gradient and Hessian-related data structures used in
!> wavefunction optimization. It coordinates allocation across 9 sub-modules covering
!> Jastrow, CI, orbital, and mixed parameter types.
!>
!> Sub-modules managed:
!> - gradhess_all: Combined gradient/Hessian/overlap for all parameters
!> - gradhessj: Jastrow parameter derivative components
!> - gradhessjo: Old Jastrow values for finite differences
!> - gradhess_ci: CI coefficient gradient/Hessian/overlap
!> - gradhess_jas: Jastrow parameter gradient/Hessian/overlap
!> - gradhess_mix_jas_ci: Jastrow-CI coupling blocks
!> - gradhess_mix_jas_orb: Jastrow-orbital coupling blocks
!> - gradhess_mix_orb_ci: CI-orbital coupling blocks
!> - gradjerr: Gradient error estimation via block averaging
!>
!> @note Call allocate_m_gradhess() after set_gradhess_all_size() but before optimization.
!> @note Call deallocate_m_gradhess() to free all optimization memory at termination.
module m_gradhess
contains
!> @brief Allocate all gradient and Hessian data structures.
!> @details Calls allocation routines for all 9 gradient/Hessian sub-modules:
!> 1. allocate_gradhess_all(): Combined matrices for all parameters
!> 2. allocate_gradhessj(): Jastrow derivative components
!> 3. allocate_gradhessjo(): Old Jastrow values
!> 4. allocate_gradhess_ci(): CI matrices
!> 5. allocate_gradhess_jas(): Jastrow matrices
!> 6. allocate_gradhess_mix_jas_ci(): Jastrow-CI coupling
!> 7. allocate_gradhess_mix_jas_orb(): Jastrow-orbital coupling
!> 8. allocate_gradhess_mix_orb_ci(): CI-orbital coupling
!> 9. allocate_gradjerr(): Error estimation arrays
subroutine allocate_m_gradhess()
    use gradhess_all, only: allocate_gradhess_all
    use gradhess_ci, only: allocate_gradhess_ci
    use gradhess_jas, only: allocate_gradhess_jas
    use gradhess_mix_jas_ci, only: allocate_gradhess_mix_jas_ci
    use gradhess_mix_jas_orb, only: allocate_gradhess_mix_jas_orb
    use gradhess_mix_orb_ci, only: allocate_gradhess_mix_orb_ci
    use gradhessj, only: allocate_gradhessj
    use gradhessjo, only: allocate_gradhessjo
    use gradjerr, only: allocate_gradjerr

    implicit none

    call allocate_gradhess_all()
    call allocate_gradhessj()
    call allocate_gradhessjo()
    call allocate_gradhess_ci()
    call allocate_gradhess_jas()
    call allocate_gradhess_mix_jas_ci()
    call allocate_gradhess_mix_jas_orb()
    call allocate_gradhess_mix_orb_ci()
    call allocate_gradjerr()
end subroutine allocate_m_gradhess

!> @brief Deallocate all gradient and Hessian data structures.
!> @details Calls deallocation routines for all 9 gradient/Hessian sub-modules
!> in the same order as allocation.
subroutine deallocate_m_gradhess()
    use gradhess_all, only: deallocate_gradhess_all
    use gradhess_ci, only: deallocate_gradhess_ci
    use gradhess_jas, only: deallocate_gradhess_jas
    use gradhess_mix_jas_ci, only: deallocate_gradhess_mix_jas_ci
    use gradhess_mix_jas_orb, only: deallocate_gradhess_mix_jas_orb
    use gradhess_mix_orb_ci, only: deallocate_gradhess_mix_orb_ci
    use gradhessj, only: deallocate_gradhessj
    use gradhessjo, only: deallocate_gradhessjo
    use gradjerr, only: deallocate_gradjerr

    implicit none

    call deallocate_gradhess_all()
    call deallocate_gradhessj()
    call deallocate_gradhessjo()
    call deallocate_gradhess_ci()
    call deallocate_gradhess_jas()
    call deallocate_gradhess_mix_jas_ci()
    call deallocate_gradhess_mix_jas_orb()
    call deallocate_gradhess_mix_orb_ci()
    call deallocate_gradjerr()
end subroutine deallocate_m_gradhess
end module 
