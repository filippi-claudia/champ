module header
   !> Arguments: date, title

   character*20 title
   character*24 date

   private 
   public :: date, title 
   save
 end module header

 module inputflags
   !> Arguments: iznuc,igeometry,ibasis_num,ilcao,iexponents,
   !             ideterminants,ijastrow_parameter, ioptorb_def,ilattice,
   !             ici_def,iforces,icsfs,imstates,igradients,icharge_efield,
   !             imultideterminants,ioptorb_mixvirt,imodify_zmat,izmatrix_check,
   !             ihessian_zmat, node_cutoff, eps_node_cutoff, scalecoef, iqmmm
   use precision_kinds, only: dp

   integer :: iznuc
   integer :: igeometry
   integer :: ibasis_num
   integer :: ilcao
   integer :: iexponents
   integer :: ideterminants
   integer :: ijastrow_parameter
   integer :: ioptorb_def
   integer :: ilattice
   integer :: ici_def
   integer :: iforces
   integer :: icsfs
   integer :: imstates
   integer :: igradients
   integer :: icharge_efield
   integer :: imultideterminants
   integer :: ioptorb_mixvirt
   integer :: imodify_zmat
   integer :: izmatrix_check
   integer :: ihessian_zmat
   integer :: node_cutoff 
   real(dp) :: eps_node_cutoff
   real(dp) :: scalecoef
   integer :: iqmmm

   private
   public :: iznuc,igeometry,ibasis_num,ilcao, iexponents
   public :: ideterminants,ijastrow_parameter, ioptorb_def,ilattice
   public :: ici_def,iforces,icsfs,imstates,igradients,icharge_efield
   public :: imultideterminants,ioptorb_mixvirt,imodify_zmat,izmatrix_check
   public :: ihessian_zmat
   public :: node_cutoff, eps_node_cutoff, scalecoef
   public :: iqmmm
   save
 end module inputflags

 module general
    !> Arguments: pooldir, pp_id, bas_id, filename, filenames_bas_num,
    !>            atomtyp, atomsymbol, wforce
    character*256 :: pooldir
    character*256 :: pp_id
    character*256 :: bas_id 
    character*256 :: filename 
    character*256, allocatable, dimension(:) :: filenames_bas_num
    character*20  :: atomtyp 
    character*20  :: atomsymbol
    character*20  :: wforce 

    private
    public :: pooldir, pp_id, bas_id, atomtyp, filename, atomsymbol
    public :: filenames_bas_num, wforce
    save
 end module general

 module method_opt
   !> should be in the input somehow no ?
   !> Arguments: method
 
   character*20 :: method

   private
   public :: method
   save
 end module method_opt

 module pars
   !> only used in read input and only Z, a20, a21 are called in vmc
   !> more are called in dmc
   !> Arguments: Z, a00, a20, a21, c0000, c1110, c2000, eps_fock, xm1, xm12, xm2, xma, xms
   use precision_kinds, only: dp

   real(dp) :: Z
   real(dp) :: a00
   real(dp) :: a20
   real(dp) :: a21
   real(dp) :: c0000
   real(dp) :: c1110
   real(dp) :: c2000
   real(dp) :: eps_fock
   real(dp) :: xm1
   real(dp) :: xm12
   real(dp) :: xm2
   real(dp) :: xma
   real(dp) :: xms

   private
   public :: Z, a00, a20, a21, c0000, c1110, c2000
   public :: eps_fock, xm1, xm12, xm2, xma, xms
   save
 end module pars
