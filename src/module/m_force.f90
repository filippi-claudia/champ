module multiple_geo

    use precision_kinds, only: dp
    use mstates_mod, only: MSTATES

     implicit none

     integer :: MFORCE
     integer, parameter :: MFORCE_WT_PRD = 1000
     integer, parameter :: MWF = 3
     integer :: nforce
     integer :: istrech
     real(dp) :: alfstr
     integer :: iwf
     integer, dimension(:), allocatable :: iwftype !(MFORCE)
     integer :: nwftype
     integer :: itausec
     integer :: nwprod
     real(dp), dimension(:, :, :), allocatable :: delc !(3,MCENT,MFORCE)

     real(dp), dimension(:, :), allocatable :: fcm2 !(MSTATES,MFORCE)
     real(dp), dimension(:, :), allocatable :: fcum !(MSTATES,MFORCE)
     ! DMC arrays:
     real(dp), dimension(:), allocatable :: fgcm2 !(MFORCE)
     real(dp), dimension(:), allocatable :: fgcum !(MFORCE)

     real(dp) :: pecent

     private
     public :: MFORCE, MFORCE_WT_PRD, MWF
     public :: nforce
     public   ::  istrech, alfstr
     public :: iwf, iwftype, nwftype
     !public :: allocate_wfsec
     public :: deallocate_wfsec
     public :: itausec, nwprod
     public   ::  delc
     public :: deallocate_forcestr
     public   ::  fcm2, fcum, fgcm2, fgcum
     public :: allocate_forcest, deallocate_forcest
     public :: pecent
     save

contains
    ! subroutine allocate_wfsec()
    !     use multiple_geo, only: MFORCE
    !     if (.not. allocated(iwftype)) allocate (iwftype(MFORCE))
    ! end subroutine allocate_wfsec

    subroutine deallocate_wfsec()
        if (allocated(iwftype)) deallocate (iwftype)
    end subroutine deallocate_wfsec

    subroutine deallocate_forcestr()
        if (allocated(delc)) deallocate (delc)
    end subroutine deallocate_forcestr

    subroutine allocate_forcest()
        use mstates_mod, only: MSTATES
        if (.not. allocated(fcm2)) allocate (fcm2(MSTATES, MFORCE))
        if (.not. allocated(fcum)) allocate (fcum(MSTATES, MFORCE))
        ! DMC arrays:
        if (.not. allocated(fgcm2)) allocate (fgcm2(MFORCE))
        if (.not. allocated(fgcum)) allocate (fgcum(MFORCE))
    end subroutine allocate_forcest

    subroutine deallocate_forcest()
        if (allocated(fcum)) deallocate (fcum)
        if (allocated(fcm2)) deallocate (fcm2)
        ! DMC arrays:
        if (allocated(fcm2)) deallocate (fgcm2)
        if (allocated(fcum)) deallocate (fgcum)
    end subroutine deallocate_forcest


end module multiple_geo


 module m_force_analytic
     !> Arguments: iforce_analy, iuse_zmat, alfgeo
     use precision_kinds, only: dp

     implicit none

     integer :: iforce_analy
     integer :: iuse_zmat
     real(dp) :: alfgeo

    ! from force_fin
     real(dp), dimension(:, :), allocatable :: da_energy_ave !(3,MCENT)
     real(dp), dimension(:), allocatable :: da_energy_err !(3)

     ! from force_mat_n
     real(dp), dimension(:, :), allocatable :: force_o !(6*MCENT,mconf)

     private
     public   :: iforce_analy, iuse_zmat, alfgeo
     public   :: da_energy_ave, da_energy_err
     public   :: allocate_force_fin, deallocate_force_fin
     public   ::  force_o
     public   :: allocate_force_mat_n, deallocate_force_mat_n
     save

 contains
     subroutine allocate_force_fin()
         use system, only: ncent_tot
         if (.not. allocated(da_energy_ave)) allocate (da_energy_ave(3, ncent_tot))
         if (.not. allocated(da_energy_err)) allocate (da_energy_err(3))
     end subroutine allocate_force_fin

     subroutine deallocate_force_fin()
         if (allocated(da_energy_err)) deallocate (da_energy_err)
         if (allocated(da_energy_ave)) deallocate (da_energy_ave)
     end subroutine deallocate_force_fin

     subroutine allocate_force_mat_n()
         use sr_mod, only: mconf
         use system, only: ncent_tot
         if (.not. allocated(force_o)) allocate (force_o(6*ncent_tot, mconf))
     end subroutine allocate_force_mat_n

     subroutine deallocate_force_mat_n()
         if (allocated(force_o)) deallocate (force_o)
     end subroutine deallocate_force_mat_n

 end module m_force_analytic

 module forcewt
     !> Arguments: wcum, wsum
     use multiple_geo, only: MFORCE
     use precision_kinds, only: dp
     use mstates_mod, only: MSTATES

     implicit none

     real(dp), dimension(:, :), allocatable :: wcum !(MSTATES,MFORCE)
     real(dp), dimension(:, :), allocatable :: wsum !(MSTATES,MFORCE)

     private
     public   ::  wcum, wsum
     public :: allocate_forcewt, deallocate_forcewt
     save
 contains
     subroutine allocate_forcewt()
         use multiple_geo, only: MFORCE
         use mstates_mod, only: MSTATES
         if (.not. allocated(wcum)) allocate (wcum(MSTATES, MFORCE))
         if (.not. allocated(wsum)) allocate (wsum(MSTATES, MFORCE))
     end subroutine allocate_forcewt

     subroutine deallocate_forcewt()
         if (allocated(wsum)) deallocate (wsum)
         if (allocated(wcum)) deallocate (wcum)
     end subroutine deallocate_forcewt

 end module forcewt

 module m_force
    contains
    subroutine allocate_m_force()
        use multiple_geo, only: allocate_forcest
        !  use multiple_geo, only: allocate_forcestr
        use forcewt, only: allocate_forcewt
        use m_force_analytic, only: allocate_force_fin
        use m_force_analytic, only: allocate_force_mat_n

        implicit none

        call allocate_forcest()
        !  call allocate_forcestr()
        call allocate_forcewt()
        call allocate_force_fin()
        call allocate_force_mat_n()
    end subroutine allocate_m_force

    subroutine deallocate_m_force()
        use multiple_geo, only: deallocate_forcest
        use multiple_geo, only: deallocate_forcestr
        use forcewt, only: deallocate_forcewt
        use m_force_analytic, only: deallocate_force_fin
        use m_force_analytic, only: deallocate_force_mat_n

        implicit none

        call deallocate_forcest()
        call deallocate_forcestr()
        call deallocate_forcewt()
        call deallocate_force_fin()
        call deallocate_force_mat_n()
    end subroutine deallocate_m_force
 end module 
