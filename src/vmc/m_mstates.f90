 module mstates_mod
  !> Arguments: MSTATES, MDETCSFX
  integer, parameter :: MSTATES=3
  integer, parameter :: MDETCSFX=20

  private 
  public :: MSTATES, MDETCSFX
  save 
 end module mstates_mod

module mstates_ctrl
   !> Arguments: iefficiency, nstates_psig, iguiding
 
   integer  :: iefficiency
   integer  :: iguiding
   integer  :: nstates_psig

   private
   public :: iefficiency, nstates_psig, iguiding
   save
 end module mstates_ctrl
 
 module mstates2
   !> Arguments: effcum, effcm2
   use precision_kinds, only: dp
   use mstates_mod, only: MSTATES

   real(dp) :: effcm2(MSTATES)
   real(dp) :: effcum(MSTATES)
   private

   public :: effcum, effcm2
   save
 end module mstates2

 module mstates3
   !> Arguments: weights_g, iweight_g
   use precision_kinds, only: dp
   use mstates_mod, only: MSTATES

   integer  :: iweight_g(MSTATES)
   real(dp) :: weights_g(MSTATES)
   private

   public :: weights_g, iweight_g
   save
 end module mstates3