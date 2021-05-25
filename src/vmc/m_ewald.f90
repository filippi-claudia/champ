module ewald_mod
     ! NP         is the number of factors of (r/cutr-1) we multiply polynomial by
     ! IVOL_RATIO is the ratio of the simulation to primitive cell volumes
     ! IBIG_RATIO is the ratio of the number of vectors or norms before the optimal Ewald separation
     !            to the number after the separation
     ! NSYM       is the ratio of the number of vectors to the number of norms
     !            and depends on the symmetry of the lattice.

     implicit none

     integer, parameter :: NCOEFX = 10, NPX = 4, IVOL_RATIO = 10, IBIG_RATIO = 15, NSYM = 8
     integer, parameter :: NGNORMX = 1000, NGVECX = NGNORMX*NSYM, NG1DX = 60
     integer, parameter :: NGNORM_SIMX = NGNORMX*IVOL_RATIO, NGVEC_SIMX = NGVECX*IVOL_RATIO
     integer, parameter :: NGNORM_BIGX = IBIG_RATIO*NGNORMX, NGVEC_BIGX = IBIG_RATIO*NGVECX
     integer, parameter :: NGNORM_SIM_BIGX = IBIG_RATIO*NGNORM_SIMX, NGVEC_SIM_BIGX = IBIG_RATIO*NGVEC_SIMX

     private
     public :: NCOEFX, NPX, IVOL_RATIO, IBIG_RATIO, NSYM
     public :: NGNORMX, NGVECX, NG1DX
     public :: NGNORM_SIMX, NGVEC_SIMX
     public :: NGNORM_BIGX, NGVEC_BIGX
     public :: NGNORM_SIM_BIGX, NGVEC_SIM_BIGX
     save
 end module ewald_mod

 module ewald
     !> Arguments: b_coul, b_coul_sim, y_coul, y_coul_sim
     use precision_kinds, only: dp
     use ewald_mod, only: NCOEFX, NGNORMX, NGNORM_SIMX

     implicit none

     real(dp), dimension(:), allocatable :: b_coul !(NCOEFX)
     real(dp), dimension(:), allocatable :: b_coul_sim !(NCOEFX)
     real(dp), dimension(:), allocatable :: y_coul !(NGNORMX)
     real(dp), dimension(:), allocatable :: y_coul_sim !(NGNORM_SIMX)

     private
     public   ::  b_coul, b_coul_sim, y_coul, y_coul_sim
     public :: allocate_ewald, deallocate_ewald
     save
 contains
     subroutine allocate_ewald()
         use precision_kinds, only: dp
         use ewald_mod, only: NCOEFX, NGNORMX, NGNORM_SIMX
         if (.not. allocated(b_coul)) allocate (b_coul(NCOEFX))
         if (.not. allocated(b_coul_sim)) allocate (b_coul_sim(NCOEFX))
         if (.not. allocated(y_coul)) allocate (y_coul(NGNORMX))
         if (.not. allocated(y_coul_sim)) allocate (y_coul_sim(NGNORM_SIMX))
     end subroutine allocate_ewald

     subroutine deallocate_ewald()
         if (allocated(y_coul_sim)) deallocate (y_coul_sim)
         if (allocated(y_coul)) deallocate (y_coul)
         if (allocated(b_coul_sim)) deallocate (b_coul_sim)
         if (allocated(b_coul)) deallocate (b_coul)
     end subroutine deallocate_ewald

 end module ewald

 module ewald_basis
     !> Arguments: vps_basis_fourier
     use precision_kinds, only: dp
     use ewald_mod, only: NGNORM_BIGX

     implicit none

     real(dp), dimension(:), allocatable :: vps_basis_fourier !(NGNORM_BIGX)

     private
     public   :: vps_basis_fourier
     public :: allocate_ewald_basis, deallocate_ewald_basis
     save
 contains
     subroutine allocate_ewald_basis()
         use precision_kinds, only: dp
         use ewald_mod, only: NGNORM_BIGX
         if (.not. allocated(vps_basis_fourier)) allocate (vps_basis_fourier(NGNORM_BIGX))
     end subroutine allocate_ewald_basis

     subroutine deallocate_ewald_basis()
         if (allocated(vps_basis_fourier)) deallocate (vps_basis_fourier)
     end subroutine deallocate_ewald_basis

 end module ewald_basis

 module periodic
     !> Arguments: cutg, cutg_big, cutg_sim, cutg_sim_big, cutr, cutr_sim, glatt, glatt_inv, glatt_sim, gnorm, gnorm_sim, gvec, gvec_sim, igmult, igmult_sim, igvec, igvec_sim, ireal_imag, isrange, k_inv, kvec, nband, ncoef_per, ng1d, ng1d_sim, ngnorm, ngnorm_big, ngnorm_orb, ngnorm_sim, ngnorm_sim_big, ngvec, ngvec_big, ngvec_orb, ngvec_sim, ngvec_sim_big, nkvec, np, npoly, rknorm, rkvec, rkvec_shift, rlatt, rlatt_inv, rlatt_sim, rlatt_sim_inv, vcell, vcell_sim, znuc2_sum, znuc_sum
     use ewald_mod, only: IVOL_RATIO
     use ewald_mod, only: NGNORM_BIGX, NGVEC_BIGX
     use ewald_mod, only: NGNORM_SIM_BIGX, NGVEC_SIM_BIGX
     use precision_kinds, only: dp
     use vmc_mod, only: MORB

     implicit none

     real(dp) :: cutg
     real(dp) :: cutg_big
     real(dp) :: cutg_sim
     real(dp) :: cutg_sim_big
     real(dp) :: cutr
     real(dp) :: cutr_sim
     real(dp), dimension(:, :), allocatable :: glatt !(3,3)
     real(dp), dimension(:, :), allocatable :: glatt_inv !(3,3)
     real(dp), dimension(:, :), allocatable :: glatt_sim !(3,3)
     real(dp), dimension(:), allocatable :: gnorm !(NGNORM_BIGX)
     real(dp), dimension(:), allocatable :: gnorm_sim !(NGNORM_SIM_BIGX)
     real(dp), dimension(:, :), allocatable :: gvec !(3,NGVEC_BIGX)
     real(dp), dimension(:, :), allocatable :: gvec_sim !(3,NGVEC_SIM_BIGX)
     integer, dimension(:), allocatable :: igmult !(NGNORM_BIGX)
     integer, dimension(:), allocatable :: igmult_sim !(NGNORM_SIM_BIGX)
     integer, dimension(:, :), allocatable :: igvec !(3,NGVEC_BIGX)
     integer, dimension(:, :), allocatable :: igvec_sim !(3,NGVEC_SIM_BIGX)
     integer, dimension(:), allocatable :: ireal_imag !(MORB)
     integer :: isrange
     integer, dimension(:), allocatable :: k_inv !(IVOL_RATIO)
     integer, dimension(:, :), allocatable :: kvec !(3,IVOL_RATIO)
     integer, dimension(:), allocatable :: nband !(IVOL_RATIO)
     integer :: ncoef_per
     integer, dimension(:), allocatable :: ng1d !(3)
     integer, dimension(:), allocatable :: ng1d_sim !(3)
     integer :: ngnorm
     integer :: ngnorm_big
     integer :: ngnorm_orb
     integer :: ngnorm_sim
     integer :: ngnorm_sim_big
     integer :: ngvec
     integer :: ngvec_big
     integer :: ngvec_orb
     integer :: ngvec_sim
     integer :: ngvec_sim_big
     integer :: nkvec
     integer :: np
     integer :: npoly
     real(dp), dimension(:), allocatable :: rknorm !(IVOL_RATIO)
     real(dp), dimension(:, :), allocatable :: rkvec !(3,IVOL_RATIO)
     real(dp), dimension(:), allocatable :: rkvec_shift !(3)
     real(dp), dimension(:, :), allocatable :: rlatt !(3,3)
     real(dp), dimension(:, :), allocatable :: rlatt_inv !(3,3)
     real(dp), dimension(:, :), allocatable :: rlatt_sim !(3,3)
     real(dp), dimension(:, :), allocatable :: rlatt_sim_inv !(3,3)
     real(dp) :: vcell
     real(dp) :: vcell_sim
     real(dp) :: znuc2_sum
     real(dp) :: znuc_sum

     private
     public :: cutg, cutg_big, cutg_sim, cutg_sim_big, cutr, cutr_sim, glatt, glatt_inv
     public :: glatt_sim, gnorm, gnorm_sim, gvec, gvec_sim, igmult, igmult_sim, igvec
     public :: igvec_sim, ireal_imag, isrange, k_inv, kvec, nband, ncoef_per, ng1d, ng1d_sim
     public :: ngnorm, ngnorm_big, ngnorm_orb, ngnorm_sim, ngnorm_sim_big, ngvec, ngvec_big
     public :: ngvec_orb, ngvec_sim, ngvec_sim_big, nkvec, np, npoly, rknorm, rkvec, rkvec_shift
     public :: rlatt, rlatt_inv, rlatt_sim, rlatt_sim_inv, vcell, vcell_sim, znuc2_sum, znuc_sum
     public :: allocate_periodic, deallocate_periodic
     save
 contains
     subroutine allocate_periodic()
         use ewald_mod, only: IVOL_RATIO
         use ewald_mod, only: NGNORM_BIGX, NGVEC_BIGX
         use ewald_mod, only: NGNORM_SIM_BIGX, NGVEC_SIM_BIGX
         use precision_kinds, only: dp
         use vmc_mod, only: MORB
         if (.not. allocated(glatt)) allocate (glatt(3, 3))
         if (.not. allocated(glatt_inv)) allocate (glatt_inv(3, 3))
         if (.not. allocated(glatt_sim)) allocate (glatt_sim(3, 3))
         if (.not. allocated(gnorm)) allocate (gnorm(NGNORM_BIGX))
         if (.not. allocated(gnorm_sim)) allocate (gnorm_sim(NGNORM_SIM_BIGX))
         if (.not. allocated(gvec)) allocate (gvec(3, NGVEC_BIGX))
         if (.not. allocated(gvec_sim)) allocate (gvec_sim(3, NGVEC_SIM_BIGX))
         if (.not. allocated(igmult)) allocate (igmult(NGNORM_BIGX))
         if (.not. allocated(igmult_sim)) allocate (igmult_sim(NGNORM_SIM_BIGX))
         if (.not. allocated(igvec)) allocate (igvec(3, NGVEC_BIGX))
         if (.not. allocated(igvec_sim)) allocate (igvec_sim(3, NGVEC_SIM_BIGX))
         if (.not. allocated(ireal_imag)) allocate (ireal_imag(MORB))
         if (.not. allocated(k_inv)) allocate (k_inv(IVOL_RATIO))
         if (.not. allocated(kvec)) allocate (kvec(3, IVOL_RATIO))
         if (.not. allocated(nband)) allocate (nband(IVOL_RATIO))
         if (.not. allocated(ng1d)) allocate (ng1d(3))
         if (.not. allocated(ng1d_sim)) allocate (ng1d_sim(3))
         if (.not. allocated(rknorm)) allocate (rknorm(IVOL_RATIO))
         if (.not. allocated(rkvec)) allocate (rkvec(3, IVOL_RATIO))
         if (.not. allocated(rkvec_shift)) allocate (rkvec_shift(3))
         if (.not. allocated(rlatt)) allocate (rlatt(3, 3))
         if (.not. allocated(rlatt_inv)) allocate (rlatt_inv(3, 3))
         if (.not. allocated(rlatt_sim)) allocate (rlatt_sim(3, 3))
         if (.not. allocated(rlatt_sim_inv)) allocate (rlatt_sim_inv(3, 3))
     end subroutine allocate_periodic

     subroutine deallocate_periodic()
         if (allocated(rlatt_sim_inv)) deallocate (rlatt_sim_inv)
         if (allocated(rlatt_sim)) deallocate (rlatt_sim)
         if (allocated(rlatt_inv)) deallocate (rlatt_inv)
         if (allocated(rlatt)) deallocate (rlatt)
         if (allocated(rkvec_shift)) deallocate (rkvec_shift)
         if (allocated(rkvec)) deallocate (rkvec)
         if (allocated(rknorm)) deallocate (rknorm)
         if (allocated(ng1d_sim)) deallocate (ng1d_sim)
         if (allocated(ng1d)) deallocate (ng1d)
         if (allocated(nband)) deallocate (nband)
         if (allocated(kvec)) deallocate (kvec)
         if (allocated(k_inv)) deallocate (k_inv)
         if (allocated(ireal_imag)) deallocate (ireal_imag)
         if (allocated(igvec_sim)) deallocate (igvec_sim)
         if (allocated(igvec)) deallocate (igvec)
         if (allocated(igmult_sim)) deallocate (igmult_sim)
         if (allocated(igmult)) deallocate (igmult)
         if (allocated(gvec_sim)) deallocate (gvec_sim)
         if (allocated(gvec)) deallocate (gvec)
         if (allocated(gnorm_sim)) deallocate (gnorm_sim)
         if (allocated(gnorm)) deallocate (gnorm)
         if (allocated(glatt_sim)) deallocate (glatt_sim)
         if (allocated(glatt_inv)) deallocate (glatt_inv)
         if (allocated(glatt)) deallocate (glatt)
     end subroutine deallocate_periodic

 end module periodic

 module pworbital
     !> Arguments: c_im, c_ip, c_rm, c_rp, icmplx, isortg, isortk, ngorb
     use ewald_mod, only: IVOL_RATIO
     use ewald_mod, only: NGVECX
     use precision_kinds, only: dp
     use vmc_mod, only: MORB

     implicit none

     real(dp), dimension(:, :), allocatable :: c_im !(NGVECX,MORB)
     real(dp), dimension(:, :), allocatable :: c_ip !(NGVECX,MORB)
     real(dp), dimension(:, :), allocatable :: c_rm !(NGVECX,MORB)
     real(dp), dimension(:, :), allocatable :: c_rp !(NGVECX,MORB)
     integer :: icmplx
     integer, dimension(:, :), allocatable :: isortg !(NGVECX,MORB)
     integer, dimension(:), allocatable :: isortk !(IVOL_RATIO)
     integer, dimension(:), allocatable :: ngorb !(MORB)

     private
     public :: c_im, c_ip, c_rm, c_rp, icmplx, isortg, isortk, ngorb
     public :: allocate_pworbital, deallocate_pworbital
     save
 contains
     subroutine allocate_pworbital()
         use ewald_mod, only: IVOL_RATIO
         use ewald_mod, only: NGVECX
         use precision_kinds, only: dp
         use vmc_mod, only: MORB
         if (.not. allocated(c_im)) allocate (c_im(NGVECX, MORB))
         if (.not. allocated(c_ip)) allocate (c_ip(NGVECX, MORB))
         if (.not. allocated(c_rm)) allocate (c_rm(NGVECX, MORB))
         if (.not. allocated(c_rp)) allocate (c_rp(NGVECX, MORB))
         if (.not. allocated(isortg)) allocate (isortg(NGVECX, MORB))
         if (.not. allocated(isortk)) allocate (isortk(IVOL_RATIO))
         if (.not. allocated(ngorb)) allocate (ngorb(MORB))
     end subroutine allocate_pworbital

     subroutine deallocate_pworbital()
         if (allocated(ngorb)) deallocate (ngorb)
         if (allocated(isortk)) deallocate (isortk)
         if (allocated(isortg)) deallocate (isortg)
         if (allocated(c_rp)) deallocate (c_rp)
         if (allocated(c_rm)) deallocate (c_rm)
         if (allocated(c_ip)) deallocate (c_ip)
         if (allocated(c_im)) deallocate (c_im)
     end subroutine deallocate_pworbital

 end module pworbital

 module test
     !> Arguments: f, vbare_coul, vbare_jas, vbare_psp
     use ewald_mod, only: NGNORM_BIGX
     use ewald_mod, only: NGNORM_SIM_BIGX
     use precision_kinds, only: dp

     implicit none

     real(dp) :: f
     real(dp), dimension(:), allocatable :: vbare_coul !(NGNORM_SIM_BIGX)
     real(dp), dimension(:), allocatable :: vbare_jas !(NGNORM_SIM_BIGX)
     real(dp), dimension(:), allocatable :: vbare_psp !(NGNORM_BIGX)

     private
     public :: f, vbare_coul, vbare_jas, vbare_psp
     public :: allocate_test, deallocate_test
     save
 contains
     subroutine allocate_test()
         use ewald_mod, only: NGNORM_BIGX
         use ewald_mod, only: NGNORM_SIM_BIGX
         use precision_kinds, only: dp
         if (.not. allocated(vbare_coul)) allocate (vbare_coul(NGNORM_SIM_BIGX))
         if (.not. allocated(vbare_jas)) allocate (vbare_jas(NGNORM_SIM_BIGX))
         if (.not. allocated(vbare_psp)) allocate (vbare_psp(NGNORM_BIGX))
     end subroutine allocate_test

     subroutine deallocate_test()
         if (allocated(vbare_psp)) deallocate (vbare_psp)
         if (allocated(vbare_jas)) deallocate (vbare_jas)
         if (allocated(vbare_coul)) deallocate (vbare_coul)
     end subroutine deallocate_test

 end module test

 module tempor
     !> Arguments: dist_nn
     use precision_kinds, only: dp

     implicit none

     real(dp) :: dist_nn

     private
     public :: dist_nn
     save
 end module tempor

 module tempor_test
     !> Arguments: c_imag, c_real, igvec_dft, iwgvec, ngg, ngvec_dft, rkvec_tmp, rkvec_tmp2
     use ewald_mod, only: IVOL_RATIO
     use ewald_mod, only: NGVEC_BIGX
     use precision_kinds, only: dp

     implicit none

     real(dp), dimension(:), allocatable :: c_imag !(NGVEC_BIGX)
     real(dp), dimension(:), allocatable :: c_real !(NGVEC_BIGX)
     integer, dimension(:, :), allocatable :: igvec_dft !(3,NGVEC_BIGX)
     integer, dimension(:), allocatable :: iwgvec !(NGVEC_BIGX)
     integer, dimension(:), allocatable :: ngg !(IVOL_RATIO)
     integer :: ngvec_dft
     real(dp), dimension(:), allocatable :: rkvec_tmp !(3)
     real(dp), dimension(:), allocatable :: rkvec_tmp2 !(3)

     private
     public :: c_imag, c_real, igvec_dft, iwgvec, ngg, ngvec_dft, rkvec_tmp, rkvec_tmp2
     public :: allocate_tempor_test, deallocate_tempor_test
     save
 contains
     subroutine allocate_tempor_test()
         use ewald_mod, only: IVOL_RATIO
         use ewald_mod, only: NGVEC_BIGX
         use precision_kinds, only: dp
         if (.not. allocated(c_imag)) allocate (c_imag(NGVEC_BIGX))
         if (.not. allocated(c_real)) allocate (c_real(NGVEC_BIGX))
         if (.not. allocated(igvec_dft)) allocate (igvec_dft(3, NGVEC_BIGX))
         if (.not. allocated(iwgvec)) allocate (iwgvec(NGVEC_BIGX))
         if (.not. allocated(ngg)) allocate (ngg(IVOL_RATIO))
         if (.not. allocated(rkvec_tmp)) allocate (rkvec_tmp(3))
         if (.not. allocated(rkvec_tmp2)) allocate (rkvec_tmp2(3))
     end subroutine allocate_tempor_test

     subroutine deallocate_tempor_test()
         if (allocated(rkvec_tmp2)) deallocate (rkvec_tmp2)
         if (allocated(rkvec_tmp)) deallocate (rkvec_tmp)
         if (allocated(ngg)) deallocate (ngg)
         if (allocated(iwgvec)) deallocate (iwgvec)
         if (allocated(igvec_dft)) deallocate (igvec_dft)
         if (allocated(c_real)) deallocate (c_real)
         if (allocated(c_imag)) deallocate (c_imag)
     end subroutine deallocate_tempor_test

 end module tempor_test

 subroutine allocate_m_ewald()
     use ewald, only: allocate_ewald
     use ewald_basis, only: allocate_ewald_basis
     use periodic, only: allocate_periodic
     use pworbital, only: allocate_pworbital
     use test, only: allocate_test
     use tempor_test, only: allocate_tempor_test

     implicit none

     call allocate_ewald()
     call allocate_ewald_basis()
     call allocate_periodic()
     call allocate_pworbital()
     call allocate_test()
     call allocate_tempor_test()
 end subroutine allocate_m_ewald

 subroutine deallocate_m_ewald()
     use ewald, only: deallocate_ewald
     use ewald_basis, only: deallocate_ewald_basis
     use periodic, only: deallocate_periodic
     use pworbital, only: deallocate_pworbital
     use test, only: deallocate_test
     use tempor_test, only: deallocate_tempor_test

     implicit none

     call deallocate_ewald()
     call deallocate_ewald_basis()
     call deallocate_periodic()
     call deallocate_pworbital()
     call deallocate_test()
     call deallocate_tempor_test()
 end subroutine deallocate_m_ewald
