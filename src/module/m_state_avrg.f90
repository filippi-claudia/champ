module sa_check
     !> Arguments: energy_all, energy_err_all
     use precision_kinds, only: dp
     use mstates_mod, only: MSTATES

     real(dp), dimension(:), allocatable :: energy_all !(MSTATES)
     real(dp), dimension(:), allocatable :: energy_err_all !(MSTATES)

     private
     public :: energy_all, energy_err_all
     public :: allocate_sa_check, deallocate_sa_check
     save
 contains
     subroutine allocate_sa_check()
         use mstates_mod, only: MSTATES
         if (.not. allocated(energy_all)) allocate (energy_all(MSTATES))
         if (.not. allocated(energy_err_all)) allocate (energy_err_all(MSTATES))
     end subroutine allocate_sa_check

     subroutine deallocate_sa_check()
         if (allocated(energy_err_all)) deallocate (energy_err_all)
         if (allocated(energy_all)) deallocate (energy_all)
     end subroutine deallocate_sa_check

 end module sa_check

 module sa_weights
     !> Arguments: iweight, nweight, weights
     use precision_kinds, only: dp
     use mstates_mod, only: MSTATES

     integer, dimension(:), allocatable :: iweight !(MSTATES)
     integer :: nweight
     real(dp), dimension(:), allocatable :: weights !(MSTATES)

     private
     public :: iweight, nweight, weights
     public :: allocate_sa_weights, deallocate_sa_weights
     save
 contains
     subroutine allocate_sa_weights()
         use mstates_mod, only: MSTATES
         if (.not. allocated(iweight)) allocate (iweight(MSTATES), source=0)
         if (.not. allocated(weights)) allocate (weights(MSTATES))
     end subroutine allocate_sa_weights

     subroutine deallocate_sa_weights()
         if (allocated(weights)) deallocate (weights)
         if (allocated(iweight)) deallocate (iweight)
     end subroutine deallocate_sa_weights

 end module sa_weights

 subroutine allocate_m_state_avrg()
     use sa_check, only: allocate_sa_check
     use sa_weights, only: allocate_sa_weights

     call allocate_sa_check()
     call allocate_sa_weights()
 end subroutine allocate_m_state_avrg

 subroutine deallocate_m_state_avrg()
     use sa_check, only: deallocate_sa_check
     use sa_weights, only: deallocate_sa_weights

     call deallocate_sa_check()
     call deallocate_sa_weights()
 end subroutine deallocate_m_state_avrg
