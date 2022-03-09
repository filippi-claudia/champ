module efield_mod
     !> Arguments: MCHARGES

     implicit none

     integer, parameter :: MCHARGES = 100
     private
     public :: MCHARGES
     save
 end module efield_mod

 module efield
     !> Arguments: iscreen, ncharges, iefield

     implicit none

     integer :: iefield
     integer :: iscreen
     integer :: ncharges

     private
     public :: iscreen, ncharges, iefield
     save
 end module efield

 module efield_blk
     !> Arguments: zcharge, bscreen, qcharge, ycharge, xcharge, ascreen
     use precision_kinds, only: dp
     use efield_mod, only: MCHARGES

     implicit none

     real(dp), dimension(:), allocatable :: ascreen !(MCHARGES)
     real(dp), dimension(:), allocatable :: bscreen !(MCHARGES)
     real(dp), dimension(:), allocatable :: qcharge !(MCHARGES)
     real(dp), dimension(:), allocatable :: xcharge !(MCHARGES)
     real(dp), dimension(:), allocatable :: ycharge !(MCHARGES)
     real(dp), dimension(:), allocatable :: zcharge !(MCHARGES)

     private
     public :: zcharge, bscreen, qcharge, ycharge, xcharge, ascreen
     public :: allocate_efield_blk, deallocate_efield_blk
     save
 contains
     subroutine allocate_efield_blk()
         use efield_mod, only: MCHARGES
         if (.not. allocated(ascreen)) allocate (ascreen(MCHARGES))
         if (.not. allocated(bscreen)) allocate (bscreen(MCHARGES))
         if (.not. allocated(qcharge)) allocate (qcharge(MCHARGES))
         if (.not. allocated(xcharge)) allocate (xcharge(MCHARGES))
         if (.not. allocated(ycharge)) allocate (ycharge(MCHARGES))
         if (.not. allocated(zcharge)) allocate (zcharge(MCHARGES))
     end subroutine allocate_efield_blk

     subroutine deallocate_efield_blk()
         if (allocated(zcharge)) deallocate (zcharge)
         if (allocated(ycharge)) deallocate (ycharge)
         if (allocated(xcharge)) deallocate (xcharge)
         if (allocated(qcharge)) deallocate (qcharge)
         if (allocated(bscreen)) deallocate (bscreen)
         if (allocated(ascreen)) deallocate (ascreen)
     end subroutine deallocate_efield_blk

 end module efield_blk

module m_efield
contains
 subroutine allocate_m_efield()
     use efield_blk, only: allocate_efield_blk

     implicit none

     call allocate_efield_blk()
 end subroutine allocate_m_efield

 subroutine deallocate_m_efield()
     use efield_blk, only: deallocate_efield_blk

     implicit none

     call deallocate_efield_blk()
 end subroutine deallocate_m_efield
end module 
