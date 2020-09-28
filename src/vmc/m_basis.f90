
 module basis
   !> Arguments: zex, betaq, n1s, n2s, n2p, n3s, n3p, n3dzr, n3dx2, n3dxy, n3dxz, n3dyz, n4s, n4p, n4fxxx, n4fyyy, n4fzzz, n4fxxy, n4fxxz, n4fyyx, n4fyyz, n4fzzx, n4fzzy, n4fxyz, nsa, npa, ndzra, ndz2a, ndxya, ndxza, ndyza
   use force_mod, only: MWF
   use precision_kinds, only: dp
   use vmc_mod, only: MBASIS, MCTYPE

   !  ncent  = number of centers
   !  zex    = screening constants for each basis function
   !  cent   = positions of each center
   !  pecent = potential energy of the centers
   !  znuc   = charge of the nuclei (centers)
   !  n1s    = number of 1s functions at each center
   !  n2s    = number of 2s functions at each center
   !  n2p    = number of 2p functions of each type at each center
   !  n3s    = number of 3s functions at each center
   !  n3p    = number of 3p functions of each type at each center
   !  n3dzr  = number of z**2-r**2 d functions at each center
   !  n3dx2  = number of x**2-y**2 d functions at each center
   !  n3dxy  = number of xy d functions at each center
   !  n3dxz  = number of xz d functions at each center
   !  n3dyz  = number of yz d functions at each center
   !  n4s    = number of 4s functions at each center
   !  n4p    = number of 4p functions of each type at each center

   real(dp) :: zex(MBASIS,MWF)
   real(dp) :: betaq
   integer  :: n1s(MCTYPE)
   integer  :: n2s(MCTYPE)
   integer  :: n2p(3,MCTYPE)
   integer  :: n3s(MCTYPE)
   integer  :: n3p(3,MCTYPE)
   integer  :: n3dzr(MCTYPE)
   integer  :: n3dx2(MCTYPE)
   integer  :: n3dxy(MCTYPE)
   integer  :: n3dxz(MCTYPE)
   integer  :: n3dyz(MCTYPE)
   integer  :: n4s(MCTYPE)
   integer  :: n4p(3,MCTYPE)
   integer  :: n4fxxx(MCTYPE)
   integer  :: n4fyyy(MCTYPE)
   integer  :: n4fzzz(MCTYPE)
   integer  :: n4fxxy(MCTYPE)
   integer  :: n4fxxz(MCTYPE)
   integer  :: n4fyyx(MCTYPE)
   integer  :: n4fyyz(MCTYPE)
   integer  :: n4fzzx(MCTYPE)
   integer  :: n4fzzy(MCTYPE)
   integer  :: n4fxyz(MCTYPE)
   integer  :: nsa(MCTYPE)
   integer  :: npa(3,MCTYPE)
   integer  :: ndzra(MCTYPE)
   integer  :: ndz2a(MCTYPE)
   integer  :: ndxya(MCTYPE)
   integer  :: ndxza(MCTYPE)
   integer  :: ndx2a(MCTYPE)
   integer  :: ndyza(MCTYPE)

   private
   public :: zex, betaq, n1s, n2s, n2p, n3s, n3p, n3dzr, n3dx2, n3dxy, n3dxz, n3dyz
   public :: n4s, n4p, n4fxxx, n4fyyy, n4fzzz, n4fxxy, n4fxxz, n4fyyx, n4fyyz
   public :: n4fzzx, n4fzzy, n4fxyz, nsa, npa, ndzra, ndz2a, ndxya, ndxza, ndyza
   public :: ndx2a
   save
 end module basis

 module numbas_mod
  !> Arguments: MRWF_PTS, MRWF
  integer, parameter :: MRWF_PTS=4000
  integer, parameter :: MRWF=200
  private 
  public :: MRWF, MRWF_PTS
  save 
 end module numbas_mod 

 module numexp
   !> Arguments: ae, ce
   use numbas_mod, only: MRWF
   use force_mod, only: MFORCE
   use precision_kinds, only: dp
   use vmc_mod, only: MCTYPE
   use vmc_mod, only: NCOEF
 
   real(dp) :: ae(2,MRWF,MCTYPE,MFORCE)
   real(dp) :: ce(NCOEF,MRWF,MCTYPE,MFORCE)
 
   private 
   public :: ae, ce 
   save
 end module numexp

 module numbas
   !> Arguments: arg, d2rwf, igrid, iwrwf, nr, nrbas, numr, r0, rwf
   use numbas_mod, only: MRWF, MRWF_PTS
   use force_mod, only: MWF
   use precision_kinds, only: dp
   use vmc_mod, only: MBASIS, MCTYPE

   real(dp) :: arg(MCTYPE)
   real(dp) :: d2rwf(MRWF_PTS,MRWF,MCTYPE,MWF)
   integer  :: igrid(MCTYPE)
   integer  :: iwrwf(MBASIS,MCTYPE)
   integer  :: nr(MCTYPE)
   integer  :: nrbas(MCTYPE)
   integer  :: numr
   real(dp) :: r0(MCTYPE)
   real(dp) :: rwf(MRWF_PTS,MRWF,MCTYPE,MWF)

   private
   public :: arg, d2rwf, igrid, iwrwf, nr, nrbas, numr, r0, rwf
   save
 end module numbas

 module numbas1
   !> Arguments: iwlbas, nbastyp
   use vmc_mod, only: MBASIS, MCTYPE

   integer  :: iwlbas(MBASIS,MCTYPE)
   integer  :: nbastyp(MCTYPE)

   private
   public :: iwlbas, nbastyp
   save
 end module numbas1

 module numbas2
   !> Arguments: ibas0, ibas1
   use vmc_mod, only: MCENT

   integer  :: ibas0(MCENT)
   integer  :: ibas1(MCENT)

   private
   public :: ibas0, ibas1
   save
 end module numbas2