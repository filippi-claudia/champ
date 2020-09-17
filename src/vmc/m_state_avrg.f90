 module sa_check
   !> Arguments: energy_all, energy_err_all
   use precision_kinds, only: dp
   use mstates_mod, only: MSTATES

   real(dp) :: energy_all(MSTATES)
   real(dp) :: energy_err_all(MSTATES)

   private
   public :: energy_all, energy_err_all
   save
 end module sa_check

 module sa_weights
   !> Arguments: iweight, nweight, weights
   use precision_kinds, only: dp
   use mstates_mod, only: MSTATES

   integer  :: iweight(MSTATES)
   integer  :: nweight
   real(dp) :: weights(MSTATES)

   private
   public :: iweight, nweight, weights
   save
 end module sa_weights