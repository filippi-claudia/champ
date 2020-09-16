 module efield_mod
  !> Arguments: MCHARGES
  integer, parameter :: MCHARGES = 100
  private
  public :: MCHARGES
  save
 end module efield_mod

 module efield
   !> Arguments: iscreen, ncharges, iefield
 
   integer  :: iefield
   integer  :: iscreen
   integer  :: ncharges

   private
   public :: iscreen, ncharges, iefield
   save
 end module efield

 module efield_blk
   !> Arguments: zcharge, bscreen, qcharge, ycharge, xcharge, ascreen
   use precision_kinds, only: dp
   use efield_mod, only: MCHARGES
 
   real(dp) :: ascreen(MCHARGES)
   real(dp) :: bscreen(MCHARGES)
   real(dp) :: qcharge(MCHARGES)
   real(dp) :: xcharge(MCHARGES)
   real(dp) :: ycharge(MCHARGES)
   real(dp) :: zcharge(MCHARGES)

   private
   public :: zcharge, bscreen, qcharge, ycharge, xcharge, ascreen
   save
 end module efield_blk