module grdnthes
     !> Arguments: hessian_zmat
     use precision_kinds, only: dp
     use vmc_mod, only: MCENT

     real(dp), dimension(:, :), allocatable :: hessian_zmat !(3,MCENT)

     private
     public   ::  hessian_zmat
     public :: allocate_grdnthes, deallocate_grdnthes
     save
 contains
     subroutine allocate_grdnthes()
         use precision_kinds, only: dp
         use vmc_mod, only: MCENT
         if (.not. allocated(hessian_zmat)) allocate (hessian_zmat(3, MCENT))
     end subroutine allocate_grdnthes

     subroutine deallocate_grdnthes()
         if (allocated(hessian_zmat)) deallocate (hessian_zmat)
     end subroutine deallocate_grdnthes

 end module grdnthes

 module grdntsmv
     !> Arguments: igrdaidx, igrdcidx, igrdmv
     use force_mod, only: MFORCE
     use vmc_mod, only: MCENT

     integer, dimension(:), allocatable :: igrdaidx !(MFORCE)
     integer, dimension(:), allocatable :: igrdcidx !(MFORCE)
     integer, dimension(:, :), allocatable :: igrdmv !(3,MCENT)

     private
     public :: igrdaidx, igrdcidx, igrdmv
     public :: allocate_grdntsmv, deallocate_grdntsmv
     save
 contains
     subroutine allocate_grdntsmv()
         use force_mod, only: MFORCE
         use vmc_mod, only: MCENT
         if (.not. allocated(igrdaidx)) allocate (igrdaidx(MFORCE))
         if (.not. allocated(igrdcidx)) allocate (igrdcidx(MFORCE))
         if (.not. allocated(igrdmv)) allocate (igrdmv(3, MCENT))
     end subroutine allocate_grdntsmv

     subroutine deallocate_grdntsmv()
         if (allocated(igrdmv)) deallocate (igrdmv)
         if (allocated(igrdcidx)) deallocate (igrdcidx)
         if (allocated(igrdaidx)) deallocate (igrdaidx)
     end subroutine deallocate_grdntsmv

 end module grdntsmv

 module grdntspar
     !> Arguments: delgrdba, delgrdbl, delgrdda, delgrdxyz, igrdtype, ngradnts
     use precision_kinds, only: dp

     real(dp) :: delgrdba
     real(dp) :: delgrdbl
     real(dp) :: delgrdda
     real(dp) :: delgrdxyz
     integer :: igrdtype
     integer :: ngradnts

     private
     public :: delgrdba, delgrdbl, delgrdda, delgrdxyz, igrdtype, ngradnts
     save
 end module grdntspar

 subroutine allocate_m_grdnt()
     use grdnthes, only: allocate_grdnthes
     use grdntsmv, only: allocate_grdntsmv

     call allocate_grdnthes()
     call allocate_grdntsmv()
 end subroutine allocate_m_grdnt

 subroutine deallocate_m_grdnt()
     use grdnthes, only: deallocate_grdnthes
     use grdntsmv, only: deallocate_grdntsmv

     call deallocate_grdnthes()
     call deallocate_grdntsmv()
 end subroutine deallocate_m_grdnt
