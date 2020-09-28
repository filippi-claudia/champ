module optci
  !> Arguments:
  !     flags and dimensions for generalized CI expectation values
  !     maximal number of terms, max dim of reduced matrices
  integer, parameter :: MXCITERM=5000
  integer, parameter :: MXCIREDUCED=1
  integer, parameter :: MXCIMATDIM=MXCITERM*(MXCIREDUCED+1)/2

  private
  public :: MXCITERM, MXCIREDUCED, MXCIMATDIM
  save 
 end module optci
  
 module ci000
   !> Arguments: iciprt, nciprim, nciterm
   !     iciprt : print flag 
   !     nciterm: number of terms (possibly after a transformation)
   !     nciprim: number of primitive terms (determinant ratios)
   
   integer  :: iciprt
   integer  :: nciprim
   integer  :: nciterm

   private
   public :: iciprt, nciprim, nciterm
   save
 end module ci000

 module ci001_blk
   !> Arguments: ci_oe, ci_o
   use optci, only: MXCITERM, MXCIREDUCED
   use precision_kinds, only: dp

   real(dp) :: ci_o(MXCITERM)
   real(dp) :: ci_oe(MXCITERM,MXCIREDUCED)

   private
   public :: ci_oe, ci_o
   save
 end module ci001_blk

 module ci002_blk
   !> Arguments: ci_o_old, ci_oe_old
   use optci, only: MXCITERM, MXCIREDUCED
   use precision_kinds, only: dp
 
    real(dp) :: ci_o_old(MXCITERM)
    real(dp) :: ci_oe_old(MXCITERM,MXCIREDUCED)

    private
    public :: ci_o_old, ci_oe_old
    save
 end module ci002_blk

 module ci003_blk
   !> Arguments: ci_e_old, ci_e
   use optci, only: MXCITERM
   use precision_kinds, only: dp

   real(dp) :: ci_e(MXCITERM)
   real(dp) :: ci_e_old(MXCITERM)

   private
   public :: ci_e_old, ci_e
   save
 end module ci003_blk

 module ci004_blk
   !> Arguments: ci_de, ci_de_old
   use optci, only: MXCITERM
   use precision_kinds, only: dp
 
   real(dp) :: ci_de(MXCITERM)
   real(dp) :: ci_de_old(MXCITERM)
 
   private
   public :: ci_de, ci_de_old
   save
 end module ci004_blk

 module ci005_blk
    !> Arguments: ci_o_cum, ci_o_sum
    use optci, only: MXCITERM
    use precision_kinds, only: dp
 
    real(dp) :: ci_o_cum(MXCITERM)
    real(dp) :: ci_o_sum(MXCITERM)

    private
    public :: ci_o_cum, ci_o_sum
    save
 end module ci005_blk
 
 module ci006_blk
   !> Arguments: ci_de_cum, ci_de_sum
   use optci, only: MXCITERM
   use precision_kinds, only: dp
 
    real(dp) :: ci_de_cum(MXCITERM)
    real(dp) :: ci_de_sum(MXCITERM)

    private
    public :: ci_de_cum, ci_de_sum
    save
 end module ci006_blk

 module ci008_blk
   !> Arguments: ci_oe_cm2, ci_oe_sum, ci_oe_cum
   use optci, only: MXCITERM, MXCIREDUCED
   use precision_kinds, only: dp
 
   real(dp) :: ci_oe_cm2(MXCITERM,MXCIREDUCED)
   real(dp) :: ci_oe_cum(MXCITERM,MXCIREDUCED)
   real(dp) :: ci_oe_sum(MXCITERM,MXCIREDUCED)

   private
   public :: ci_oe_cm2, ci_oe_sum, ci_oe_cum
   save
 end module ci008_blk

 module ci009_blk
   !> Arguments: ci_oo_sum, ci_oo_cm2, ci_oo_cum
   use optci, only: MXCIMATDIM
   use precision_kinds, only: dp
 
   real(dp) :: ci_oo_cm2(MXCIMATDIM)
   real(dp) :: ci_oo_cum(MXCIMATDIM)
   real(dp) :: ci_oo_sum(MXCIMATDIM)

   private
   public :: ci_oo_sum, ci_oo_cm2, ci_oo_cum
   save
 end module ci009_blk

 module ci010_blk
   !> Arguments: ci_ooe_cum, ci_ooe_sum
   use optci, only: MXCIMATDIM
   use precision_kinds, only: dp
 
   real(dp) :: ci_ooe_cum(MXCIMATDIM)
   real(dp) :: ci_ooe_sum(MXCIMATDIM)
   private
 
   public :: ci_ooe_cum, ci_ooe_sum
   save
 end module ci010_blk