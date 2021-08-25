module mix_jas_ci
     !> Arguments: de_o_ci, dj_de_ci, dj_o_ci, dj_oe_ci
     use optjas, only: MPARMJ
     use precision_kinds, only: dp

     real(dp), dimension(:, :), allocatable :: de_o_ci !(MPARMJ,MDET)
     real(dp), dimension(:, :), allocatable :: dj_de_ci !(MPARMJ,MDET)
     real(dp), dimension(:, :), allocatable :: dj_o_ci !(MPARMJ,MDET)
     real(dp), dimension(:, :), allocatable :: dj_oe_ci !(MPARMJ,MDET)

     private
     public :: de_o_ci, dj_de_ci, dj_o_ci, dj_oe_ci
     public :: allocate_mix_jas_ci, deallocate_mix_jas_ci
     save
 contains
     subroutine allocate_mix_jas_ci()
         use optjas, only: MPARMJ
         use dets, only: ndet
         if (.not. allocated(de_o_ci)) allocate (de_o_ci(MPARMJ, ndet))
         if (.not. allocated(dj_de_ci)) allocate (dj_de_ci(MPARMJ, ndet))
         if (.not. allocated(dj_o_ci)) allocate (dj_o_ci(MPARMJ, ndet))
         if (.not. allocated(dj_oe_ci)) allocate (dj_oe_ci(MPARMJ, ndet))
     end subroutine allocate_mix_jas_ci

     subroutine deallocate_mix_jas_ci()
         if (allocated(dj_oe_ci)) deallocate (dj_oe_ci)
         if (allocated(dj_o_ci)) deallocate (dj_o_ci)
         if (allocated(dj_de_ci)) deallocate (dj_de_ci)
         if (allocated(de_o_ci)) deallocate (de_o_ci)
     end subroutine deallocate_mix_jas_ci

 end module mix_jas_ci

 module mix_jas_orb
     !> Arguments: de_o, dj_ho, dj_o, dj_oe
     use optorb_mod, only: MXREDUCED
     use optjas, only: MPARMJ
     use precision_kinds, only: dp
     use mstates_mod, only: MSTATES

     real(dp), dimension(:, :, :), allocatable :: de_o !(MPARMJ,MXREDUCED,MSTATES)
     real(dp), dimension(:, :, :), allocatable :: dj_ho !(MPARMJ,MXREDUCED,MSTATES)
     real(dp), dimension(:, :, :), allocatable :: dj_o !(MPARMJ,MXREDUCED,MSTATES)
     real(dp), dimension(:, :, :), allocatable :: dj_oe !(MPARMJ,MXREDUCED,MSTATES)

     private
     public :: de_o, dj_ho, dj_o, dj_oe
     public :: allocate_mix_jas_orb, deallocate_mix_jas_orb
     save
 contains
     subroutine allocate_mix_jas_orb()
         use optorb_mod, only: MXREDUCED
         use optjas, only: MPARMJ
         use mstates_mod, only: MSTATES
         if (.not. allocated(de_o)) allocate (de_o(MPARMJ, MXREDUCED, MSTATES))
         if (.not. allocated(dj_ho)) allocate (dj_ho(MPARMJ, MXREDUCED, MSTATES))
         if (.not. allocated(dj_o)) allocate (dj_o(MPARMJ, MXREDUCED, MSTATES))
         if (.not. allocated(dj_oe)) allocate (dj_oe(MPARMJ, MXREDUCED, MSTATES))
     end subroutine allocate_mix_jas_orb

     subroutine deallocate_mix_jas_orb()
         if (allocated(dj_oe)) deallocate (dj_oe)
         if (allocated(dj_o)) deallocate (dj_o)
         if (allocated(dj_ho)) deallocate (dj_ho)
         if (allocated(de_o)) deallocate (de_o)
     end subroutine deallocate_mix_jas_orb

 end module mix_jas_orb

 module mix_orb_ci
     !> Arguments: ci_de_o, ci_o_ho, ci_o_o, ci_o_oe
     use optorb_mod, only: MXREDUCED
     use optci, only: MXCITERM
     use precision_kinds, only: dp

     real(dp), dimension(:, :), allocatable :: ci_de_o !(MXCITERM,MXREDUCED)
     real(dp), dimension(:, :), allocatable :: ci_o_ho !(MXCITERM,MXREDUCED)
     real(dp), dimension(:, :), allocatable :: ci_o_o !(MXCITERM,MXREDUCED)
     real(dp), dimension(:, :), allocatable :: ci_o_oe !(MXCITERM,MXREDUCED)

     private
     public :: ci_de_o, ci_o_ho, ci_o_o, ci_o_oe
     public :: allocate_mix_orb_ci, deallocate_mix_orb_ci
     save
 contains
     subroutine allocate_mix_orb_ci()
         use optorb_mod, only: MXREDUCED
         use optci, only: MXCITERM
         if (.not. allocated(ci_de_o)) allocate (ci_de_o(MXCITERM, MXREDUCED))
         if (.not. allocated(ci_o_ho)) allocate (ci_o_ho(MXCITERM, MXREDUCED))
         if (.not. allocated(ci_o_o)) allocate (ci_o_o(MXCITERM, MXREDUCED))
         if (.not. allocated(ci_o_oe)) allocate (ci_o_oe(MXCITERM, MXREDUCED))
     end subroutine allocate_mix_orb_ci

     subroutine deallocate_mix_orb_ci()
         if (allocated(ci_o_oe)) deallocate (ci_o_oe)
         if (allocated(ci_o_o)) deallocate (ci_o_o)
         if (allocated(ci_o_ho)) deallocate (ci_o_ho)
         if (allocated(ci_de_o)) deallocate (ci_de_o)
     end subroutine deallocate_mix_orb_ci

 end module mix_orb_ci

 subroutine allocate_m_mixderiv()
     use mix_jas_ci, only: allocate_mix_jas_ci
     use mix_jas_orb, only: allocate_mix_jas_orb
     use mix_orb_ci, only: allocate_mix_orb_ci

     call allocate_mix_jas_ci()
     call allocate_mix_jas_orb()
     call allocate_mix_orb_ci()
 end subroutine allocate_m_mixderiv

 subroutine deallocate_m_mixderiv()
     use mix_jas_ci, only: deallocate_mix_jas_ci
     use mix_jas_orb, only: deallocate_mix_jas_orb
     use mix_orb_ci, only: deallocate_mix_orb_ci

     call deallocate_mix_jas_ci()
     call deallocate_mix_jas_orb()
     call deallocate_mix_orb_ci()
 end subroutine deallocate_m_mixderiv
