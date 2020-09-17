 module grdnthes
   !> Arguments: hessian_zmat
   use precision_kinds, only: dp
   use vmc, only: MCENT

   real(dp) :: hessian_zmat(3,MCENT)

   private
   public   ::  hessian_zmat 
   save
 end module grdnthes

 module grdntsmv
   !> Arguments: igrdaidx, igrdcidx, igrdmv
   use force_mod, only: MFORCE
   use vmc, only: MCENT

    integer  :: igrdaidx(MFORCE)
    integer  :: igrdcidx(MFORCE)
    integer  :: igrdmv(3,MCENT)

    private 
    public :: igrdaidx, igrdcidx, igrdmv 
    save
 end module grdntsmv

 module grdntspar
   !> Arguments: delgrdba, delgrdbl, delgrdda, delgrdxyz, igrdtype, ngradnts
   use precision_kinds, only: dp

   real(dp) :: delgrdba
   real(dp) :: delgrdbl
   real(dp) :: delgrdda
   real(dp) :: delgrdxyz
   integer  :: igrdtype
   integer  :: ngradnts

   private 
   public :: delgrdba, delgrdbl, delgrdda, delgrdxyz, igrdtype, ngradnts 
   save
 end module grdntspar