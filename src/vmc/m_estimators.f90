module estcum
     !> Arguments: ecum, ecum1, iblk, pecum, r2cum, tjfcum, tpbcum, avcum
     use force_mod, only: MFORCE
     use precision_kinds, only: dp
     use mstates_mod, only: MSTATES

     real(dp), dimension(:, :), allocatable :: ecum !(MSTATES,MFORCE)
     real(dp), dimension(:), allocatable :: ecum1 !(MSTATES)
     integer :: iblk
     real(dp), dimension(:), allocatable :: pecum !(MSTATES)
     real(dp) :: r2cum
     real(dp), dimension(:), allocatable :: tjfcum !(MSTATES)
     real(dp), dimension(:), allocatable :: tpbcum !(MSTATES)
     real(dp), dimension(:), allocatable :: avcum !(MSTATES*3)

     private
     public   ::  ecum, ecum1, iblk, pecum, r2cum, tjfcum, tpbcum, avcum
     public :: allocate_estcum, deallocate_estcum
     save
 contains
     subroutine allocate_estcum()
         use force_mod, only: MFORCE
         use precision_kinds, only: dp
         use mstates_mod, only: MSTATES
         if (.not. allocated(ecum)) allocate (ecum(MSTATES, MFORCE))
         if (.not. allocated(ecum1)) allocate (ecum1(MSTATES))
         if (.not. allocated(pecum)) allocate (pecum(MSTATES))
         if (.not. allocated(tjfcum)) allocate (tjfcum(MSTATES))
         if (.not. allocated(tpbcum)) allocate (tpbcum(MSTATES))
         if (.not. allocated(avcum)) allocate (avcum(MSTATES*3))
     end subroutine allocate_estcum

     subroutine deallocate_estcum()
         if (allocated(avcum)) deallocate (avcum)
         if (allocated(tpbcum)) deallocate (tpbcum)
         if (allocated(tjfcum)) deallocate (tjfcum)
         if (allocated(pecum)) deallocate (pecum)
         if (allocated(ecum1)) deallocate (ecum1)
         if (allocated(ecum)) deallocate (ecum)
     end subroutine deallocate_estcum

 end module estcum

 module estsig
     !> Arguments: ecm21s, ecum1s
     use precision_kinds, only: dp
     use mstates_mod, only: MSTATES

     real(dp), dimension(:), allocatable :: ecm21s !(MSTATES)
     real(dp), dimension(:), allocatable :: ecum1s !(MSTATES)

     private
     public   ::  ecm21s, ecum1s
     public :: allocate_estsig, deallocate_estsig
     save
 contains
     subroutine allocate_estsig()
         use precision_kinds, only: dp
         use mstates_mod, only: MSTATES
         if (.not. allocated(ecm21s)) allocate (ecm21s(MSTATES))
         if (.not. allocated(ecum1s)) allocate (ecum1s(MSTATES))
     end subroutine allocate_estsig

     subroutine deallocate_estsig()
         if (allocated(ecum1s)) deallocate (ecum1s)
         if (allocated(ecm21s)) deallocate (ecm21s)
     end subroutine deallocate_estsig

 end module estsig

 module estsum
     !> Arguments: acc, esum, esum1, pesum, r2sum, tjfsum, tpbsum,
     !> efsum, efsum1, egsum, egsum1, ei1sum, ei2sum, ei3sum, esum1_dmc, esum_dmc,
     !> pesum_dmc, r2sum, risum, tausum, tjfsum_dmc, tpbsum_dmc, w_acc_sum, w_acc_sum1, wdsum,
     !> wdsum1, wfsum, wfsum1, wg_acc_sum, wg_acc_sum1, wgdsum, wgsum, wgsum1, wsum1, wsum_dmc

     use force_mod, only: MFORCE
     use precision_kinds, only: dp
     use mstates_mod, only: MSTATES

     real(dp) :: acc
     real(dp), dimension(:, :), allocatable :: esum !(MSTATES,MFORCE)
     real(dp), dimension(:), allocatable :: esum1 !(MSTATES)
     real(dp), dimension(:), allocatable :: pesum !(MSTATES)
     real(dp) :: r2sum
     real(dp), dimension(:), allocatable :: tjfsum !(MSTATES)
     real(dp), dimension(:), allocatable :: tpbsum !(MSTATES)
     !> DMC variables: 
     real(dp) :: efsum
     real(dp) :: efsum1
     real(dp), dimension(:), allocatable :: egsum  !(MFORCE)
     real(dp), dimension(:), allocatable :: egsum1 !(MFORCE)
     real(dp) :: ei1sum
     real(dp) :: ei2sum
     real(dp) :: ei3sum
     real(dp), dimension(:), allocatable :: esum1_dmc !(MFORCE)
     real(dp) :: esum_dmc
     real(dp), dimension(:), allocatable :: pesum_dmc !(MFORCE)
     real(dp) :: risum
     real(dp), dimension(:), allocatable :: tausum !(MFORCE)
     real(dp), dimension(:), allocatable :: tjfsum_dmc !(MFORCE)
     real(dp), dimension(:), allocatable :: tpbsum_dmc !(MFORCE)
     real(dp) :: w_acc_sum
     real(dp) :: w_acc_sum1
     real(dp) :: wdsum
     real(dp) :: wdsum1
     real(dp) :: wfsum
     real(dp) :: wfsum1
     real(dp) :: wg_acc_sum
     real(dp) :: wg_acc_sum1
     real(dp) :: wgdsum
     real(dp), dimension(:), allocatable :: wgsum !(MFORCE)
     real(dp), dimension(:), allocatable :: wgsum1 !(MFORCE)
     real(dp), dimension(:), allocatable :: wsum1 !(MFORCE)
     real(dp) :: wsum_dmc

     private
     public :: acc, esum, esum1, pesum, r2sum, tjfsum, tpbsum
     public :: allocate_estsum, deallocate_estsum
     public :: efsum, efsum1, egsum, egsum1, ei1sum, ei2sum, ei3sum, esum1_dmc, esum_dmc
     public :: pesum_dmc, risum, tausum, tjfsum_dmc, tpbsum_dmc, w_acc_sum
     public :: w_acc_sum1, wdsum, wdsum1, wfsum, wfsum1, wg_acc_sum, wg_acc_sum1
     public :: wgdsum, wgsum, wgsum1, wsum1, wsum_dmc
     public :: allocate_estsum_dmc, deallocate_estsum_dmc
     save

 contains
     subroutine allocate_estsum()
         use force_mod, only: MFORCE
         use mstates_mod, only: MSTATES
         if (.not. allocated(esum)) allocate (esum(MSTATES, MFORCE))
         if (.not. allocated(esum1)) allocate (esum1(MSTATES))
         if (.not. allocated(pesum)) allocate (pesum(MSTATES))
         if (.not. allocated(tjfsum)) allocate (tjfsum(MSTATES))
         if (.not. allocated(tpbsum)) allocate (tpbsum(MSTATES))
     end subroutine allocate_estsum

     subroutine deallocate_estsum()
         if (allocated(tpbsum)) deallocate (tpbsum)
         if (allocated(tjfsum)) deallocate (tjfsum)
         if (allocated(pesum)) deallocate (pesum)
         if (allocated(esum1)) deallocate (esum1)
         if (allocated(esum)) deallocate (esum)
     end subroutine deallocate_estsum

     subroutine allocate_estsum_dmc()
         use force_mod, only: MFORCE
         if (.not. allocated(egsum)) allocate(egsum(MFORCE))
         if (.not. allocated(egsum1)) allocate(egsum1(MFORCE))
         if (.not. allocated(esum1_dmc)) allocate(esum1_dmc(MFORCE))
         if (.not. allocated(pesum_dmc)) allocate(pesum_dmc(MFORCE))
         if (.not. allocated(tausum)) allocate(tausum(MFORCE))
         if (.not. allocated(tjfsum_dmc)) allocate(tjfsum_dmc(MFORCE))
         if (.not. allocated(tpbsum_dmc)) allocate(tpbsum_dmc(MFORCE))
         if (.not. allocated(wgsum)) allocate(wgsum(MFORCE))
         if (.not. allocated(wgsum1)) allocate(wgsum1(MFORCE))
         if (.not. allocated(wsum1)) allocate(wsum1(MFORCE))
     end subroutine allocate_estsum_dmc

     subroutine deallocate_estsum_dmc()
         if (allocated(egsum)) deallocate(egsum)
         if (allocated(egsum1)) deallocate(egsum1)
         if (allocated(esum1_dmc)) deallocate(esum1_dmc)
         if (allocated(pesum_dmc)) deallocate(pesum_dmc)
         if (allocated(tausum)) deallocate(tausum)
         if (allocated(tjfsum_dmc)) deallocate(tjfsum_dmc)
         if (allocated(tpbsum_dmc)) deallocate(tpbsum_dmc)
         if (allocated(wgsum)) deallocate(wgsum)
         if (allocated(wgsum1)) deallocate(wgsum1)
         if (allocated(wsum1)) deallocate(wsum1)
     end subroutine deallocate_estsum_dmc
 end module estsum

 module estpsi
     !> Arguments: apsi, aref, detref
     use precision_kinds, only: dp
     use mstates_mod, only: MSTATES

     real(dp), dimension(:), allocatable :: apsi !(MSTATES)
     real(dp) :: aref
     real(dp), dimension(:), allocatable :: detref !(2)

     private
     public   ::  apsi, aref, detref
     public :: allocate_estpsi, deallocate_estpsi
     save
 contains
     subroutine allocate_estpsi()
         use precision_kinds, only: dp
         use mstates_mod, only: MSTATES
         if (.not. allocated(apsi)) allocate (apsi(MSTATES))
         if (.not. allocated(detref)) allocate (detref(2))
     end subroutine allocate_estpsi

     subroutine deallocate_estpsi()
         if (allocated(detref)) deallocate (detref)
         if (allocated(apsi)) deallocate (apsi)
     end subroutine deallocate_estpsi

 end module estpsi

 module est2cm
     !> Arguments: ecm2, ecm21, pecm2, r2cm2, tjfcm2, tpbcm2, avcm2
     use force_mod, only: MFORCE
     use precision_kinds, only: dp
     use mstates_mod, only: MSTATES

     real(dp), dimension(:, :), allocatable :: ecm2 !(MSTATES,MFORCE)
     real(dp), dimension(:), allocatable :: ecm21 !(MSTATES)
     real(dp), dimension(:), allocatable :: pecm2 !(MSTATES)
     real(dp) :: r2cm2
     real(dp), dimension(:), allocatable :: tjfcm2 !(MSTATES)
     real(dp), dimension(:), allocatable :: tpbcm2 !(MSTATES)
     real(dp), dimension(:), allocatable :: avcm2 !(MSTATES*3)

     private
     public   :: ecm2, ecm21, pecm2, r2cm2, tjfcm2, tpbcm2, avcm2
     public :: allocate_est2cm, deallocate_est2cm
     save
 contains
     subroutine allocate_est2cm()
         use force_mod, only: MFORCE
         use precision_kinds, only: dp
         use mstates_mod, only: MSTATES
         if (.not. allocated(ecm2)) allocate (ecm2(MSTATES, MFORCE))
         if (.not. allocated(ecm21)) allocate (ecm21(MSTATES))
         if (.not. allocated(pecm2)) allocate (pecm2(MSTATES))
         if (.not. allocated(tjfcm2)) allocate (tjfcm2(MSTATES))
         if (.not. allocated(tpbcm2)) allocate (tpbcm2(MSTATES))
         if (.not. allocated(avcm2)) allocate (avcm2(MSTATES*3))
     end subroutine allocate_est2cm

     subroutine deallocate_est2cm()
         if (allocated(avcm2)) deallocate (avcm2)
         if (allocated(tpbcm2)) deallocate (tpbcm2)
         if (allocated(tjfcm2)) deallocate (tjfcm2)
         if (allocated(pecm2)) deallocate (pecm2)
         if (allocated(ecm21)) deallocate (ecm21)
         if (allocated(ecm2)) deallocate (ecm2)
     end subroutine deallocate_est2cm

 end module est2cm

 subroutine allocate_m_estimators()
     use estcum, only: allocate_estcum
     use estsig, only: allocate_estsig
     use estsum, only: allocate_estsum
     use estpsi, only: allocate_estpsi
     use est2cm, only: allocate_est2cm

     call allocate_estcum()
     call allocate_estsig()
     call allocate_estsum()
     call allocate_estpsi()
     call allocate_est2cm()
 end subroutine allocate_m_estimators
