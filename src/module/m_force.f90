 module forcewt
     !> Arguments: wcum, wsum
      use mstates_mod, only: MSTATES
      use multiple_geo, only: MFORCE
      use precision_kinds, only: dp

     implicit none

     real(dp), dimension(:, :), allocatable :: wcum !(MSTATES,MFORCE)
     real(dp), dimension(:, :), allocatable :: wsum !(MSTATES,MFORCE)

     private
     public   ::  wcum, wsum
     public :: allocate_forcewt, deallocate_forcewt
     save
 contains
     subroutine allocate_forcewt()
      use mstates_mod, only: MSTATES
      use multiple_geo, only: MFORCE
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
      use forcewt, only: allocate_forcewt
      use m_force_analytic, only: allocate_force_fin
      use m_force_analytic, only: allocate_force_mat_n
      use multiple_geo, only: allocate_forcest
        !  use multiple_geo, only: allocate_forcestr

        implicit none

        call allocate_forcest()
        !  call allocate_forcestr()
        call allocate_forcewt()
        call allocate_force_fin()
        call allocate_force_mat_n()
    end subroutine allocate_m_force

    subroutine deallocate_m_force()
      use forcewt, only: deallocate_forcewt
      use m_force_analytic, only: deallocate_force_fin
      use m_force_analytic, only: deallocate_force_mat_n
      use multiple_geo, only: deallocate_forcest,deallocate_forcestr

        implicit none

        call deallocate_forcest()
        call deallocate_forcestr()
        call deallocate_forcewt()
        call deallocate_force_fin()
        call deallocate_force_mat_n()
    end subroutine deallocate_m_force
 end module 
