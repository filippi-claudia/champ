module const
   !> Arguments: pi, hb, etrial, delta, deltai, fbias, nelec, imetro, ipr
   use precision_kinds, only: dp

   real(dp) :: delta
   real(dp) :: deltai
   real(dp) :: etrial
   real(dp) :: fbias
   real(dp) :: hb
   integer  :: imetro
   integer  :: ipr
   integer  :: nelec
   real(dp) :: pi

   private
   public   :: pi, hb, etrial, delta, deltai, fbias, nelec, imetro, ipr
   save
 end module const

 module const2
  !> Arguments: deltar, deltat
  use precision_kinds, only: dp

   real(dp) :: deltar
   real(dp) :: deltat

   private
   public :: deltar, deltat
   save
 end module const2

 module constant
  !> Arguments: twopi
  use precision_kinds, only: dp

  real(dp) :: twopi

  private
  public :: twopi
  save
 end module constant