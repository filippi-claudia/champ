module properties
  !> Arguments: MAXPROP
  integer, parameter :: MAXPROP=6
  private
  public :: MAXPROP
  save 
 end module properties

 module prp000
   !> Arguments: ipropprt, iprop, nprop
 
   integer  :: iprop
   integer  :: ipropprt
   integer  :: nprop

   private
   public :: ipropprt, iprop, nprop
   save
 end module prp000
 
 module prp001
   !> Arguments: vprop
   use properties, only: MAXPROP
   use precision_kinds, only: dp

   real(dp) :: vprop(MAXPROP)

   private
   public :: vprop
   save
 end module prp001

 module prp002
   !> Arguments: vprop_old
   use properties, only: MAXPROP
   use precision_kinds, only: dp
   include 'dmc.h'
 
   real(dp) :: vprop_old(MAXPROP,MWALK)
   real(dp) :: vprop_old2(MAXPROP)

   private
   public :: vprop_old, vprop_old2
   save
 end module prp002

 module prp003
   !> Arguments: vprop_cm2, vprop_cum, cc_nuc, vprop_sum
   use properties, only: MAXPROP
   use precision_kinds, only: dp

    real(dp) :: cc_nuc(3)
    real(dp) :: vprop_cm2(MAXPROP)
    real(dp) :: vprop_cum(MAXPROP)
    real(dp) :: vprop_sum(MAXPROP)

    private
    public :: vprop_cm2, vprop_cum, cc_nuc, vprop_sum
    save
 end module prp003