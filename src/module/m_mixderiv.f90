module mix_jas_ci
     !> Arguments: de_o_ci, dj_de_ci, dj_o_ci, dj_oe_ci
     use optwf_parms, only: nparmj
     use precision_kinds, only: dp

     real(dp), dimension(:, :, :), allocatable :: de_o_ci !(nparmj,MDET,nwftypemax)
     real(dp), dimension(:, :, :), allocatable :: dj_de_ci !(nparmj,MDET,nwftypemax)
     real(dp), dimension(:, :, :), allocatable :: dj_o_ci !(nparmj,MDET,nwftypemax)
     real(dp), dimension(:, :, :), allocatable :: dj_oe_ci !(nparmj,MDET,nwftypemax)

     private
     public :: de_o_ci, dj_de_ci, dj_o_ci, dj_oe_ci
     public :: allocate_mix_jas_ci, deallocate_mix_jas_ci
     save
 contains
     subroutine allocate_mix_jas_ci()
      use optwf_parms, only: nparmj
      use slater, only: ndet
      use vmc_mod, only: nwftypemax
         if (.not. allocated(de_o_ci)) allocate (de_o_ci(nparmj, ndet, nwftypemax))
         if (.not. allocated(dj_de_ci)) allocate (dj_de_ci(nparmj, ndet, nwftypemax))
         if (.not. allocated(dj_o_ci)) allocate (dj_o_ci(nparmj, ndet, nwftypemax))
         if (.not. allocated(dj_oe_ci)) allocate (dj_oe_ci(nparmj, ndet, nwftypemax))
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
     use mstates_mod, only: MSTATES
     use optorb_mod, only: mxreduced
     use optwf_parms, only: nparmj
     use precision_kinds, only: dp

     real(dp), dimension(:, :, :), allocatable :: de_o !(nparmj,mxreduced,MSTATES)
     real(dp), dimension(:, :, :), allocatable :: dj_ho !(nparmj,mxreduced,MSTATES)
     real(dp), dimension(:, :, :), allocatable :: dj_o !(nparmj,mxreduced,MSTATES)
     real(dp), dimension(:, :, :), allocatable :: dj_oe !(nparmj,mxreduced,MSTATES)

     private
     public :: de_o, dj_ho, dj_o, dj_oe
     public :: allocate_mix_jas_orb, deallocate_mix_jas_orb
     save
 contains
     subroutine allocate_mix_jas_orb()
      use mstates_mod, only: MSTATES
      use optorb_mod, only: mxreduced
      use optwf_parms, only: nparmj
         if (.not. allocated(de_o)) allocate (de_o(nparmj, mxreduced, MSTATES))
         if (.not. allocated(dj_ho)) allocate (dj_ho(nparmj, mxreduced, MSTATES))
         if (.not. allocated(dj_o)) allocate (dj_o(nparmj, mxreduced, MSTATES))
         if (.not. allocated(dj_oe)) allocate (dj_oe(nparmj, mxreduced, MSTATES))
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
     use optci, only: mxciterm
     use optorb_mod, only: mxreduced
     use precision_kinds, only: dp

     real(dp), dimension(:, :), allocatable :: ci_de_o !(mxciterm,mxreduced)
     real(dp), dimension(:, :), allocatable :: ci_o_ho !(mxciterm,mxreduced)
     real(dp), dimension(:, :), allocatable :: ci_o_o !(mxciterm,mxreduced)
     real(dp), dimension(:, :), allocatable :: ci_o_oe !(mxciterm,mxreduced)

     private
     public :: ci_de_o, ci_o_ho, ci_o_o, ci_o_oe
     public :: allocate_mix_orb_ci, deallocate_mix_orb_ci
     save
 contains
     subroutine allocate_mix_orb_ci()
      use optci, only: mxciterm
      use optorb_mod, only: mxreduced
         if (.not. allocated(ci_de_o)) allocate (ci_de_o(mxciterm, mxreduced))
         if (.not. allocated(ci_o_ho)) allocate (ci_o_ho(mxciterm, mxreduced))
         if (.not. allocated(ci_o_o)) allocate (ci_o_o(mxciterm, mxreduced))
         if (.not. allocated(ci_o_oe)) allocate (ci_o_oe(mxciterm, mxreduced))
     end subroutine allocate_mix_orb_ci

     subroutine deallocate_mix_orb_ci()
         if (allocated(ci_o_oe)) deallocate (ci_o_oe)
         if (allocated(ci_o_o)) deallocate (ci_o_o)
         if (allocated(ci_o_ho)) deallocate (ci_o_ho)
         if (allocated(ci_de_o)) deallocate (ci_de_o)
     end subroutine deallocate_mix_orb_ci

 end module mix_orb_ci

 module m_mixderiv
 contains
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
 end module 
