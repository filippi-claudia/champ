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
        use csfs, only: nstates
         use force_mod, only: MFORCE
         use precision_kinds, only: dp
         use mstates_mod, only: MSTATES
         if (.not. allocated(ecum)) allocate (ecum(nstates, MFORCE))
         if (.not. allocated(ecum1)) allocate (ecum1(nstates))
         if (.not. allocated(pecum)) allocate (pecum(nstates))
         if (.not. allocated(tjfcum)) allocate (tjfcum(nstates))
         if (.not. allocated(tpbcum)) allocate (tpbcum(nstates))
         if (.not. allocated(avcum)) allocate (avcum(nstates*3))
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
        use csfs, only: nstates
         use precision_kinds, only: dp
         use mstates_mod, only: MSTATES
         if (.not. allocated(ecm21s)) allocate (ecm21s(nstates))
         if (.not. allocated(ecum1s)) allocate (ecum1s(nstates))
     end subroutine allocate_estsig

     subroutine deallocate_estsig()
         if (allocated(ecum1s)) deallocate (ecum1s)
         if (allocated(ecm21s)) deallocate (ecm21s)
     end subroutine deallocate_estsig

 end module estsig

 module estsum
     !> Arguments: acc, esum, esum1, pesum, r2sum, tjfsum, tpbsum
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

     private
     public   :: acc, esum, esum1, pesum, r2sum, tjfsum, tpbsum
     public :: allocate_estsum, deallocate_estsum
     save
 contains
     subroutine allocate_estsum()
        use csfs, only: nstates
         use force_mod, only: MFORCE
         use precision_kinds, only: dp
         use mstates_mod, only: MSTATES
         if (.not. allocated(esum)) allocate (esum(nstates, MFORCE))
         if (.not. allocated(esum1)) allocate (esum1(nstates))
         if (.not. allocated(pesum)) allocate (pesum(nstates))
         if (.not. allocated(tjfsum)) allocate (tjfsum(nstates))
         if (.not. allocated(tpbsum)) allocate (tpbsum(nstates))
     end subroutine allocate_estsum

     subroutine deallocate_estsum()
         if (allocated(tpbsum)) deallocate (tpbsum)
         if (allocated(tjfsum)) deallocate (tjfsum)
         if (allocated(pesum)) deallocate (pesum)
         if (allocated(esum1)) deallocate (esum1)
         if (allocated(esum)) deallocate (esum)
     end subroutine deallocate_estsum

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
        use csfs, only: nstates
         use precision_kinds, only: dp
         use mstates_mod, only: MSTATES
         if (.not. allocated(apsi)) allocate (apsi(nstates))
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
        use csfs, only: nstates
         use force_mod, only: MFORCE
         use precision_kinds, only: dp
         use mstates_mod, only: MSTATES
         if (.not. allocated(ecm2)) allocate (ecm2(nstates, MFORCE))
         if (.not. allocated(ecm21)) allocate (ecm21(nstates))
         if (.not. allocated(pecm2)) allocate (pecm2(nstates))
         if (.not. allocated(tjfcm2)) allocate (tjfcm2(nstates))
         if (.not. allocated(tpbcm2)) allocate (tpbcm2(nstates))
         if (.not. allocated(avcm2)) allocate (avcm2(nstates*3))
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


