module ewald_mod
     ! NP         is the number of factors of (r/cutr-1) we multiply polynomial by

     implicit none

!  Original parameters     
     integer, parameter :: NCOEFX = 10, NPX = 4
     integer, parameter :: NG1DX = 60, NGNORM_BIGX = 80000, NGVEC_BIGX = 80000

     private
     public :: NCOEFX, NPX, NG1DX
     public :: NGNORM_BIGX, NGVEC_BIGX
     save
 end module ewald_mod

 module periodic
     !> Arguments: cutg, cutg_big, cutr, glatt, glatt_inv,
     !> gnorm, gvec, igmult, igvec
     !> isrange, ncoef_per, ng1d, ngnorm, ngnorm_big,
     !> ngvec, ngvec_big
     !> np_coul, np_jas, npoly, rlatt, rlatt_inv, vcell,
     !> znuc2_sum, znuc_sum
      use ewald_mod, only: NGNORM_BIGX
      use ewald_mod, only: NGVEC_BIGX
      use precision_kinds, only: dp
      use vmc_mod, only: norb_tot
      use system, only : nctype
      
     implicit none

     real(dp) :: cutg
     real(dp) :: cutg_big
     real(dp) :: cutr
     real(dp), dimension(:, :), allocatable :: glatt !(3,3)
     real(dp), dimension(:, :), allocatable :: glatt_inv !(3,3)
     real(dp), dimension(:),    allocatable :: gnorm !(NGNORM_BIGX)
     real(dp), dimension(:, :), allocatable :: gvec !(3,NGVEC_BIGX)
     integer, dimension(:),     allocatable :: igmult !(NGNORM_BIGX)
     integer, dimension(:, :),  allocatable :: igvec !(3,NGVEC_BIGX)
     integer, dimension(:),     allocatable :: ng1d !(3)
     integer :: isrange, ncoef_per
     integer :: ngnorm, ngnorm_big, ngvec, ngvec_big
     integer :: np_coul, np_jas, npoly
     integer :: n_images
     real(dp), dimension(:, :), allocatable :: rlatt !(3,3)
     real(dp), dimension(:, :), allocatable :: rlatt_inv !(3,3)
     real(dp), dimension(:, :), allocatable :: ell !(3,n_images)
     real(dp) :: vcell
     real(dp) :: znuc2_sum
     real(dp) :: znuc_sum
     real(dp) :: alattice
     
     private
     public :: cutg, cutg_big, cutr, glatt, glatt_inv
     public :: gnorm, gvec, igmult, igvec
     public :: isrange, ncoef_per, ng1d
     public :: ngnorm, ngnorm_big, ngvec, ngvec_big
     public :: np_coul, np_jas, npoly
     public :: rlatt, rlatt_inv, vcell, znuc2_sum, znuc_sum
     public :: allocate_periodic, deallocate_periodic
     public :: alattice, n_images, ell
     save
 contains
     subroutine allocate_periodic()
      use ewald_mod, only: NGNORM_BIGX, NGVEC_BIGX
      use vmc_mod, only: norb_tot
         if (.not. allocated(glatt)) allocate (glatt(3, 3))
         if (.not. allocated(glatt_inv)) allocate (glatt_inv(3, 3))
         if (.not. allocated(gnorm)) allocate (gnorm(NGNORM_BIGX))
         if (.not. allocated(gvec)) allocate (gvec(3, NGVEC_BIGX))
         if (.not. allocated(igmult)) allocate (igmult(NGNORM_BIGX), source=0)
         if (.not. allocated(igvec)) allocate (igvec(3, NGVEC_BIGX), source=0)
         if (.not. allocated(ng1d)) allocate (ng1d(3), source=0)
         if (.not. allocated(rlatt)) allocate (rlatt(3, 3))
         if (.not. allocated(rlatt_inv)) allocate (rlatt_inv(3, 3))
         if (.not. allocated(ell)) allocate (ell(3, n_images))
     end subroutine allocate_periodic

     subroutine deallocate_periodic()
         if (allocated(rlatt_inv)) deallocate (rlatt_inv)
         if (allocated(rlatt)) deallocate (rlatt)
         if (allocated(ng1d)) deallocate (ng1d)
         if (allocated(igvec)) deallocate (igvec)
         if (allocated(igmult)) deallocate (igmult)
         if (allocated(gvec)) deallocate (gvec)
         if (allocated(gnorm)) deallocate (gnorm)
         if (allocated(glatt_inv)) deallocate (glatt_inv)
         if (allocated(glatt)) deallocate (glatt)
         if (allocated(ell)) deallocate (ell)
     end subroutine deallocate_periodic

 end module periodic

 module ewald
     !> Arguments: b_coul, y_coul
     use ewald_mod, only: NCOEFX
     use system, only: ncent_tot
     use precision_kinds, only: dp

     implicit none

     real(dp), dimension(:), allocatable :: b_coul !(NCOEFX)
     real(dp), dimension(:), allocatable :: b_jas  !(NCOEFX)
     real(dp), dimension(:), allocatable :: y_coul !(NGNORMX)
     real(dp), dimension(:), allocatable :: y_jas  !(NGNORMX)
     real(dp), dimension(:), allocatable :: cos_n_sum !(NGVECX)
     real(dp), dimension(:), allocatable :: sin_n_sum !(NGVECX)
     real(dp), dimension(:), allocatable :: cos_e_sum !(NGVECX)
     real(dp), dimension(:), allocatable :: sin_e_sum !(NGVECX)
     real(dp), dimension(:), allocatable :: cos_sum ! (NGVECX)
     real(dp), dimension(:), allocatable :: sin_sum ! (NGVECX)
     real(dp), dimension(:,:), allocatable :: cos_g ! (nelec,NGVECX)
     real(dp), dimension(:,:), allocatable :: sin_g ! (nelec,NGVECX)
     real(dp), dimension(:,:,:), allocatable :: dcos_g ! (3,nelec,NGVECX)
     real(dp), dimension(:,:,:), allocatable :: dsin_g ! (3,nelec,NGVECX)
     real(dp), dimension(:), allocatable :: cos_sum_new ! (NGVECX)
     real(dp), dimension(:), allocatable :: sin_sum_new ! (NGVECX)
     real(dp), dimension(:,:), allocatable :: cos_g_new ! (nelec,NGVECX)
     real(dp), dimension(:,:), allocatable :: sin_g_new ! (nelec,NGVECX)
     real(dp), dimension(:,:,:), allocatable :: dcos_g_new ! (3,nelec,NGVECX)
     real(dp), dimension(:,:,:), allocatable :: dsin_g_new ! (3,nelec,NGVECX)
     real(dp), dimension(:), allocatable :: sk !(NGNORMX number of shells)
     
     private
     public   ::  b_coul, b_jas, y_coul, y_jas
     public   ::  cos_n_sum, sin_n_sum, cos_e_sum, sin_e_sum
     public   ::  cos_sum, sin_sum, cos_g, sin_g, dcos_g, dsin_g
     public   ::  cos_sum_new, sin_sum_new, cos_g_new, sin_g_new, dcos_g_new, dsin_g_new
     public :: allocate_ewald, deallocate_ewald, sk
     save
 contains
     subroutine allocate_ewald()
      use ewald_mod, only: NCOEFX
      use periodic, only: ngnorm,ngvec
      use system, only: nelec
         if (.not. allocated(b_coul)) allocate (b_coul(NCOEFX))
         if (.not. allocated(b_jas))  allocate (b_jas(NCOEFX))
         if (.not. allocated(y_coul)) allocate (y_coul(ngnorm))
         if (.not. allocated(y_jas))  allocate (y_jas(ngnorm))
         if (.not. allocated(sin_n_sum)) allocate (sin_n_sum(ngvec))
         if (.not. allocated(cos_n_sum)) allocate (cos_n_sum(ngvec))
         if (.not. allocated(sin_e_sum)) allocate (sin_e_sum(ngvec))
         if (.not. allocated(cos_e_sum)) allocate (cos_e_sum(ngvec))
         if (.not. allocated(sin_sum)) allocate (sin_sum(ngvec))
         if (.not. allocated(cos_sum)) allocate (cos_sum(ngvec))
         if (.not. allocated(cos_g)) allocate (cos_g(nelec,ngvec))
         if (.not. allocated(sin_g)) allocate (sin_g(nelec,ngvec))
         if (.not. allocated(dcos_g)) allocate (dcos_g(3,nelec,ngvec))
         if (.not. allocated(dsin_g)) allocate (dsin_g(3,nelec,ngvec))
         if (.not. allocated(sin_sum_new)) allocate (sin_sum_new(ngvec))
         if (.not. allocated(cos_sum_new)) allocate (cos_sum_new(ngvec))
         if (.not. allocated(cos_g_new)) allocate (cos_g_new(nelec,ngvec))
         if (.not. allocated(sin_g_new)) allocate (sin_g_new(nelec,ngvec))
         if (.not. allocated(dcos_g_new)) allocate (dcos_g_new(3,nelec,ngvec))
         if (.not. allocated(dsin_g_new)) allocate (dsin_g_new(3,nelec,ngvec))
         if (.not. allocated(sk)) allocate (sk(ngnorm))
     end subroutine allocate_ewald

     subroutine deallocate_ewald()
         if (allocated(y_coul)) deallocate (y_coul)
         if (allocated(y_jas))  deallocate (y_jas)
         if (allocated(b_coul)) deallocate (b_coul)
         if (allocated(b_jas))  deallocate (b_jas)
         if (allocated(cos_n_sum)) deallocate (cos_n_sum)
         if (allocated(sin_n_sum)) deallocate (sin_n_sum)
         if (allocated(cos_e_sum)) deallocate (cos_e_sum)
         if (allocated(sin_e_sum)) deallocate (sin_e_sum)
         if (allocated(cos_sum)) deallocate (cos_sum)
         if (allocated(sin_sum)) deallocate (sin_sum)
         if (allocated(cos_g)) deallocate (cos_g)
         if (allocated(sin_g)) deallocate (sin_g)
         if (allocated(dcos_g)) deallocate (dcos_g)
         if (allocated(dsin_g)) deallocate (dsin_g)
         if (allocated(cos_sum_new)) deallocate (cos_sum_new)
         if (allocated(sin_sum_new)) deallocate (sin_sum_new)
         if (allocated(cos_g_new)) deallocate (cos_g_new)
         if (allocated(sin_g_new)) deallocate (sin_g_new)
         if (allocated(dcos_g_new)) deallocate (dcos_g_new)
         if (allocated(dsin_g_new)) deallocate (dsin_g_new)
         if (allocated(sk)) deallocate (sk)
     end subroutine deallocate_ewald

 end module ewald

 module ewald_test
     !> Arguments: f, vbare_coul, vbare_jas
      use ewald_mod, only: NGNORM_BIGX
      use precision_kinds, only: dp

     implicit none

     real(dp) :: f
     real(dp), dimension(:), allocatable :: vbare_coul !(NGNORM_BIGX)
     real(dp), dimension(:), allocatable :: vbare_jas !(NGNORM_BIGX)

     private
     public :: f, vbare_coul, vbare_jas
     public :: allocate_ewald_test, deallocate_ewald_test
     save
 contains
     subroutine allocate_ewald_test()
      use ewald_mod, only: NGNORM_BIGX
         if (.not. allocated(vbare_coul)) allocate (vbare_coul(NGNORM_BIGX))
         if (.not. allocated(vbare_jas)) allocate (vbare_jas(NGNORM_BIGX))
     end subroutine allocate_ewald_test

     subroutine deallocate_ewald_test()
         if (allocated(vbare_jas)) deallocate (vbare_jas)
         if (allocated(vbare_coul)) deallocate (vbare_coul)
     end subroutine deallocate_ewald_test
 end module ewald_test

module m_ewald
contains
! allocation of ewald done in set_ewald after definition of ngvec etc.
! subroutine allocate_m_ewald()
!     use ewald, only: allocate_ewald
!     use periodic, only: allocate_periodic
!     use ewald_test, only: allocate_ewald_test
!
!     implicit none
!
!     call allocate_periodic()
!     call allocate_ewald()
! end subroutine allocate_m_ewald

 subroutine deallocate_m_ewald()
     use ewald, only: deallocate_ewald
     use periodic, only: deallocate_periodic
     use ewald_test, only: deallocate_ewald_test

     implicit none

     call deallocate_periodic()
     call deallocate_ewald()
 end subroutine deallocate_m_ewald
end module

