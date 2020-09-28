module jasn
   !> Arguments: d2ijn, d2n, fijn, fjn, fsn, fsumn
   use precision_kinds, only: dp
   use vmc_mod, only: MELEC

   real(dp) :: d2ijn(MELEC,MELEC)
   real(dp) :: d2n
   real(dp) :: fijn(3,MELEC,MELEC)
   real(dp) :: fjn(3,MELEC)
   real(dp) :: fsn(MELEC,MELEC)
   real(dp) :: fsumn

   private 
   public :: d2ijn, d2n, fijn, fjn, fsn, fsumn 
   save
 end module jasn

 module jaso
   !> Arguments: d2ijo, d2o, fijo, fjo, fso, fsumo
   use precision_kinds, only: dp
   use vmc_mod, only: MELEC

    real(dp) :: d2ijo(MELEC,MELEC)
    real(dp) :: d2o
    real(dp) :: fijo(3,MELEC,MELEC)
    real(dp) :: fjo(3,MELEC)
    real(dp) :: fso(MELEC,MELEC)
    real(dp) :: fsumo

    private
    public :: d2ijo, d2o, fijo, fjo, fso, fsumo
    save
 end module jaso

 module jaspar
   !> Arguments: nspin1, nspin2, sspin, sspinn, is
   use precision_kinds, only: dp

   integer  :: is
   integer  :: nspin1
   integer  :: nspin2
   real(dp) :: sspin
   real(dp) :: sspinn

   private
   public   :: nspin1, nspin2, sspin, sspinn, is
   save
 end module jaspar

 module jaspar1
   !> Arguments: cjas1, cjas2
   use force_mod, only: MWF
   use precision_kinds, only: dp

   real(dp) :: cjas1(MWF)
   real(dp) :: cjas2(MWF)

   private
   public   ::  cjas1, cjas2
   save
 end module jaspar1

 module jaspar2
   !> Arguments: a1, a2
   use force_mod, only: MWF
   use precision_kinds, only: dp

   real(dp) :: a1(83,3,MWF)
   real(dp) :: a2(83,3,MWF)

   private 
   public :: a1, a2 
   save
 end module jaspar2

 module jaspar3
   !> Arguments: a, b, c, fck, nord, scalek
   use force_mod, only: MWF
   use precision_kinds, only: dp
   use vmc_mod, only: MCTYPE
   use vmc_mod, only: MORDJ1

   real(dp) :: a(MORDJ1,MWF)
   real(dp) :: b(MORDJ1,2,MWF)
   real(dp) :: c(83,MCTYPE,MWF)
   real(dp) :: fck(15,MCTYPE,MWF)
   integer  :: nord
   real(dp) :: scalek(MWF)

   private 
   public :: a, b, c, fck, nord, scalek 
   save
 end module jaspar3

 module jaspar4
   !> Arguments: a4, norda, nordb, nordc
   use force_mod, only: MWF
   use precision_kinds, only: dp
   use vmc_mod, only: MCTYPE
   use vmc_mod, only: MORDJ1

   real(dp) :: a4(MORDJ1,MCTYPE,MWF)
   integer  :: norda
   integer  :: nordb
   integer  :: nordc

   private 
   public :: a4, norda, nordb, nordc 
   save
 end module jaspar4

 module jaspar6
   !> Arguments: asymp_jasa, asymp_jasb, asymp_r, c1_jas6, c1_jas6i, c2_jas6, cutjas, cutjasi
   use precision_kinds, only: dp
   use vmc_mod, only: MCTYPE

   real(dp) :: asymp_jasa(MCTYPE)
   real(dp) :: asymp_jasb(2)
   real(dp) :: asymp_r
   real(dp) :: c1_jas6
   real(dp) :: c1_jas6i
   real(dp) :: c2_jas6
   real(dp) :: cutjas
   real(dp) :: cutjasi

   private 
   public :: asymp_jasa, asymp_jasb, asymp_r, c1_jas6, c1_jas6i, c2_jas6, cutjas, cutjasi 
   save
 end module jaspar6

 module jaspointer
   !> Arguments: npoint, npointa
   use vmc_mod, only: MCTYP3X

   integer  :: npoint(MCTYP3X)
   integer  :: npointa(3*MCTYP3X)

   private 
   public :: npoint, npointa 
   save
 end module jaspointer