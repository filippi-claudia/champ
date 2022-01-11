module estcum
     !> Arguments: ecum, ecum1, iblk, pecum, r2cum, tjfcum, tpbcum, avcum
     use force_mod, only: MFORCE
     use precision_kinds, only: dp
     use mstates_mod, only: MSTATES

     implicit none 

     real(dp), dimension(:, :), allocatable :: ecum !(MSTATES,MFORCE)
     real(dp), dimension(:), allocatable :: ecum1 !(MSTATES)
     integer :: iblk
     real(dp), dimension(:), allocatable :: pecum !(MSTATES)
     real(dp) :: r2cum
     real(dp), dimension(:), allocatable :: tjfcum !(MSTATES)
     real(dp), dimension(:), allocatable :: tpbcum !(MSTATES)
     real(dp), dimension(:), allocatable :: avcum !(MSTATES*3)
     !> DMC variables:
     real(dp) :: ecum1_dmc
     real(dp) :: ecum_dmc
     real(dp) :: efcum
     real(dp) :: efcum1
     real(dp), dimension(:), allocatable :: egcum  !(MFORCE)
     real(dp), dimension(:), allocatable :: egcum1 !(MFORCE)
     real(dp) :: ei1cum
     real(dp) :: ei2cum
     real(dp) :: ei3cum
     real(dp), dimension(:), allocatable :: pecum_dmc !(MFORCE)
     real(dp) :: r2cum_dmc
     real(dp) :: ricum
     real(dp), dimension(:), allocatable :: taucum !(MFORCE)
     real(dp), dimension(:), allocatable :: tjfcum_dmc !(MFORCE)
     real(dp), dimension(:), allocatable :: tpbcum_dmc !(MFORCE)
     real(dp) :: w_acc_cum
     real(dp) :: w_acc_cum1
     real(dp) :: wcum1
     real(dp) :: wcum_dmc
     real(dp) :: wdcum
     real(dp) :: wdcum1
     real(dp) :: wfcum
     real(dp) :: wfcum1
     real(dp) :: wg_acc_cum
     real(dp) :: wg_acc_cum1
     real(dp), dimension(:), allocatable :: wgcum !(MFORCE)
     real(dp), dimension(:), allocatable :: wgcum1 !(MFORCE)
     real(dp) :: wgdcum
     integer :: ipass

     private
     public :: ecum, ecum1, iblk, pecum, r2cum, tjfcum, tpbcum, avcum
     public :: allocate_estcum, deallocate_estcum
     !> DMC variables:
     public :: ecum1_dmc, ecum_dmc, efcum, efcum1, egcum, egcum1, ei1cum, ei2cum
     public :: ei3cum, pecum_dmc, r2cum_dmc, ricum, taucum, tjfcum_dmc, tpbcum_dmc
     public :: w_acc_cum, w_acc_cum1, wcum1, wcum_dmc, wdcum, wdcum1, wfcum, wfcum1
     public :: wg_acc_cum, wg_acc_cum1, wgcum, wgcum1
     public :: wgdcum, ipass
     public :: allocate_estcum_dmc, deallocate_estcum_dmc
     save
 contains
     subroutine allocate_estcum()
         use force_mod, only: MFORCE
         use mstates_mod, only: MSTATES
         if (.not. allocated(ecum)) allocate (ecum(MSTATES, MFORCE), source=0.0_dp)
         if (.not. allocated(ecum1)) allocate (ecum1(MSTATES), source=0.0_dp)
         if (.not. allocated(pecum)) allocate (pecum(MSTATES), source=0.0_dp)
         if (.not. allocated(tjfcum)) allocate (tjfcum(MSTATES), source=0.0_dp)
         if (.not. allocated(tpbcum)) allocate (tpbcum(MSTATES), source=0.0_dp)
         if (.not. allocated(avcum)) allocate (avcum(MSTATES*3), source=0.0_dp)
     end subroutine allocate_estcum

     subroutine deallocate_estcum()
         if (allocated(avcum)) deallocate (avcum)
         if (allocated(tpbcum)) deallocate (tpbcum)
         if (allocated(tjfcum)) deallocate (tjfcum)
         if (allocated(pecum)) deallocate (pecum)
         if (allocated(ecum1)) deallocate (ecum1)
         if (allocated(ecum)) deallocate (ecum)
     end subroutine deallocate_estcum

     subroutine allocate_estcum_dmc()
         if (.not. allocated(egcum)) allocate(egcum(MFORCE), source=0.0_dp)
         if (.not. allocated(egcum1)) allocate(egcum1(MFORCE), source=0.0_dp)
         if (.not. allocated(pecum_dmc)) allocate(pecum_dmc(MFORCE), source=0.0_dp)
         if (.not. allocated(taucum)) allocate(taucum(MFORCE), source=0.0_dp)
         if (.not. allocated(tjfcum_dmc)) allocate(tjfcum_dmc(MFORCE), source=0.0_dp)
         if (.not. allocated(tpbcum_dmc)) allocate(tpbcum_dmc(MFORCE), source=0.0_dp)
         if (.not. allocated(wgcum)) allocate(wgcum(MFORCE), source=0.0_dp)
         if (.not. allocated(wgcum1)) allocate(wgcum1(MFORCE), source=0.0_dp)
     end subroutine allocate_estcum_dmc

     subroutine deallocate_estcum_dmc()
         if (allocated(egcum)) deallocate(egcum)
         if (allocated(egcum1)) deallocate(egcum1)
         if (allocated(pecum_dmc)) deallocate(pecum_dmc)
         if (allocated(taucum)) deallocate(taucum)
         if (allocated(tjfcum_dmc)) deallocate(tjfcum_dmc)
         if (allocated(tpbcum_dmc)) deallocate(tpbcum_dmc)
         if (allocated(wgcum)) deallocate(wgcum)
         if (allocated(wgcum1)) deallocate(wgcum1)
     end subroutine deallocate_estcum_dmc
 end module estcum

 module estsig
     !> Arguments: ecm21s, ecum1s
     use precision_kinds, only: dp
     use mstates_mod, only: MSTATES

     implicit none 

     real(dp), dimension(:), allocatable :: ecm21s !(MSTATES)
     real(dp), dimension(:), allocatable :: ecum1s !(MSTATES)

     private
     public   ::  ecm21s, ecum1s
     public :: allocate_estsig, deallocate_estsig
     save
 contains
     subroutine allocate_estsig()
         use mstates_mod, only: MSTATES
         if (.not. allocated(ecm21s)) allocate (ecm21s(MSTATES), source=0.0_dp)
         if (.not. allocated(ecum1s)) allocate (ecum1s(MSTATES), source=0.0_dp)
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

     implicit none 

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
         if (.not. allocated(esum)) allocate (esum(MSTATES, MFORCE), source=0.0_dp)
         if (.not. allocated(esum1)) allocate (esum1(MSTATES), source=0.0_dp)
         if (.not. allocated(pesum)) allocate (pesum(MSTATES), source=0.0_dp)
         if (.not. allocated(tjfsum)) allocate (tjfsum(MSTATES), source=0.0_dp)
         if (.not. allocated(tpbsum)) allocate (tpbsum(MSTATES), source=0.0_dp)
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
         if (.not. allocated(egsum)) allocate(egsum(MFORCE), source=0.0_dp)
         if (.not. allocated(egsum1)) allocate(egsum1(MFORCE), source=0.0_dp)
         if (.not. allocated(esum1_dmc)) allocate(esum1_dmc(MFORCE), source=0.0_dp)
         if (.not. allocated(pesum_dmc)) allocate(pesum_dmc(MFORCE), source=0.0_dp)
         if (.not. allocated(tausum)) allocate(tausum(MFORCE), source=0.0_dp)
         if (.not. allocated(tjfsum_dmc)) allocate(tjfsum_dmc(MFORCE), source=0.0_dp)
         if (.not. allocated(tpbsum_dmc)) allocate(tpbsum_dmc(MFORCE), source=0.0_dp)
         if (.not. allocated(wgsum)) allocate(wgsum(MFORCE), source=0.0_dp)
         if (.not. allocated(wgsum1)) allocate(wgsum1(MFORCE), source=0.0_dp)
         if (.not. allocated(wsum1)) allocate(wsum1(MFORCE), source=0.0_dp)
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

     implicit none 

     real(dp), dimension(:), allocatable :: apsi !(MSTATES)
     real(dp) :: aref
     real(dp), dimension(:), allocatable :: detref !(2)

     private
     public   ::  apsi, aref, detref
     public :: allocate_estpsi, deallocate_estpsi
     save
 contains
     subroutine allocate_estpsi()
         use mstates_mod, only: MSTATES
         if (.not. allocated(apsi)) allocate (apsi(MSTATES), source=0.0_dp)
         if (.not. allocated(detref)) allocate (detref(2), source=0.0_dp)
     end subroutine allocate_estpsi

     subroutine deallocate_estpsi()
         if (allocated(detref)) deallocate (detref)
         if (allocated(apsi)) deallocate (apsi)
     end subroutine deallocate_estpsi

 end module estpsi

 module est2cm
     !> Arguments: ecm2, ecm21, pecm2, r2cm2, tjfcm2, tpbcm2, avcm2,
     !> ecm21_dmc, ecm2_dmc, efcm2, efcm21, egcm2, egcm21, ei1cm2, ei2cm2,
     !> ei3cm2, pecm2_dmc, r2cm2_dmc, ricm2, tjfcm_dmc, tpbcm2_dmc, wcm2, wcm21, wdcm2, wdcm21,
     !> wfcm2, wfcm21, wgcm2, wgcm21, wgdcm2

     use force_mod, only: MFORCE
     use precision_kinds, only: dp
     use mstates_mod, only: MSTATES

     implicit none 

     real(dp), dimension(:, :), allocatable :: ecm2 !(MSTATES,MFORCE)
     real(dp), dimension(:), allocatable :: ecm21 !(MSTATES)
     real(dp), dimension(:), allocatable :: pecm2 !(MSTATES)
     real(dp) :: r2cm2
     real(dp), dimension(:), allocatable :: tjfcm2 !(MSTATES)
     real(dp), dimension(:), allocatable :: tpbcm2 !(MSTATES)
     real(dp), dimension(:), allocatable :: avcm2 !(MSTATES*3)
     !> DMC variables:
     real(dp) :: ecm21_dmc
     real(dp) :: ecm2_dmc
     real(dp) :: efcm2
     real(dp) :: efcm21
     real(dp), dimension(:), allocatable :: egcm2  !(MFORCE)
     real(dp), dimension(:), allocatable :: egcm21 !(MFORCE)
     real(dp) :: ei1cm2
     real(dp) :: ei2cm2
     real(dp) :: ei3cm2
     real(dp), dimension(:), allocatable :: pecm2_dmc !(MFORCE)
     real(dp) :: r2cm2_dmc
     real(dp) :: ricm2
     real(dp), dimension(:), allocatable :: tjfcm_dmc  !(MFORCE)
     real(dp), dimension(:), allocatable :: tpbcm2_dmc !(MFORCE)
     real(dp) :: wcm2
     real(dp) :: wcm21
     real(dp) :: wdcm2
     real(dp) :: wdcm21
     real(dp) :: wfcm2
     real(dp) :: wfcm21
     real(dp), dimension(:), allocatable :: wgcm2 !(MFORCE)
     real(dp), dimension(:), allocatable :: wgcm21 !(MFORCE)
     real(dp) :: wgdcm2

     private
     public :: ecm2, ecm21, pecm2, r2cm2, tjfcm2, tpbcm2, avcm2
     public :: allocate_est2cm, deallocate_est2cm
     public :: ecm21_dmc, ecm2_dmc, efcm2, efcm21, egcm2, egcm21, ei1cm2, ei2cm2
     public :: ei3cm2, pecm2_dmc, r2cm2_dmc, ricm2, tjfcm_dmc, tpbcm2_dmc
     public :: wcm2, wcm21, wdcm2, wdcm21
     public :: wfcm2, wfcm21, wgcm2, wgcm21, wgdcm2
     public :: allocate_est2cm_dmc, deallocate_est2cm_dmc
     save

 contains
     subroutine allocate_est2cm()
         use force_mod, only: MFORCE
         use mstates_mod, only: MSTATES
         if (.not. allocated(ecm2))   allocate(ecm2(MSTATES, MFORCE), source=0.0_dp)
         if (.not. allocated(ecm21))  allocate(ecm21(MSTATES), source=0.0_dp)
         if (.not. allocated(pecm2))  allocate(pecm2(MSTATES), source=0.0_dp)
         if (.not. allocated(tjfcm2)) allocate(tjfcm2(MSTATES), source=0.0_dp)
         if (.not. allocated(tpbcm2)) allocate(tpbcm2(MSTATES), source=0.0_dp)
         if (.not. allocated(avcm2))  allocate(avcm2(MSTATES*3), source=0.0_dp)
     end subroutine allocate_est2cm

     subroutine deallocate_est2cm()
         if (allocated(avcm2))  deallocate(avcm2)
         if (allocated(tpbcm2)) deallocate(tpbcm2)
         if (allocated(tjfcm2)) deallocate(tjfcm2)
         if (allocated(pecm2))  deallocate(pecm2)
         if (allocated(ecm21))  deallocate(ecm21)
         if (allocated(ecm2))   deallocate(ecm2)
     end subroutine deallocate_est2cm

     subroutine allocate_est2cm_dmc()
         if (.not. allocated(egcm2))      allocate(egcm2(MFORCE), source=0.0_dp)
         if (.not. allocated(egcm21))     allocate(egcm21(MFORCE), source=0.0_dp)
         if (.not. allocated(pecm2_dmc))  allocate(pecm2_dmc(MFORCE), source=0.0_dp)
         if (.not. allocated(tjfcm_dmc))  allocate(tjfcm_dmc(MFORCE), source=0.0_dp)
         if (.not. allocated(tpbcm2_dmc)) allocate(tpbcm2_dmc(MFORCE), source=0.0_dp)
         if (.not. allocated(wgcm2))      allocate(wgcm2(MFORCE), source=0.0_dp)
         if (.not. allocated(wgcm21))     allocate(wgcm21(MFORCE), source=0.0_dp)
     end subroutine allocate_est2cm_dmc

     subroutine deallocate_est2cm_dmc()
         if (allocated(egcm2))      deallocate(egcm2)
         if (allocated(egcm21))     deallocate(egcm21)
         if (allocated(pecm2_dmc))  deallocate(pecm2_dmc)
         if (allocated(tjfcm_dmc))  deallocate(tjfcm_dmc)
         if (allocated(tpbcm2_dmc)) deallocate(tpbcm2_dmc)
         if (allocated(wgcm2))      deallocate(wgcm2)
         if (allocated(wgcm21))     deallocate(wgcm21)
     end subroutine deallocate_est2cm_dmc

 end module est2cm

 subroutine allocate_m_estimators()
     use estcum, only: allocate_estcum
     use estsig, only: allocate_estsig
     use estsum, only: allocate_estsum
     use estpsi, only: allocate_estpsi
     use est2cm, only: allocate_est2cm

     implicit none 

     call allocate_estcum()
     call allocate_estsig()
     call allocate_estsum()
     call allocate_estpsi()
     call allocate_est2cm()
 end subroutine allocate_m_estimators

 subroutine deallocate_m_estimators()
     use estcum, only: deallocate_estcum
     use estsig, only: deallocate_estsig
     use estsum, only: deallocate_estsum
     use estpsi, only: deallocate_estpsi
     use est2cm, only: deallocate_est2cm

     implicit none 

     call deallocate_estcum()
     call deallocate_estsig()
     call deallocate_estsum()
     call deallocate_estpsi()
     call deallocate_est2cm()
 end subroutine deallocate_m_estimators
