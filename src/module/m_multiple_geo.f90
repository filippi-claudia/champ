module multiple_geo

      use mstates_mod, only: MSTATES
      use precision_kinds, only: dp

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

