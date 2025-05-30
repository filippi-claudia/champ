module da_energy_sumcum
    !> Arguments: da_energy_cm2, da_energy_cum, da_energy_sum, da_psi_cum, da_psi_sum
      use precision_kinds, only: dp

    implicit none

    real(dp), dimension(:, :, :), allocatable :: da_energy_cm2 !(3, MCENT, PTH)
    real(dp), dimension(:, :, :), allocatable :: da_energy_cum !(3, MCENT, PTH)
    real(dp), dimension(:, :, :), allocatable :: da_energy_sum !(3, MCENT, PTH)
    real(dp), dimension(:, :, :), allocatable :: da_psi_cum !(3, MCENT, PTH)
    real(dp), dimension(:, :, :), allocatable :: da_psi_sum !(3, MCENT, PTH)

    private
    public :: da_energy_cm2, da_energy_cum, da_energy_sum, da_psi_cum, da_psi_sum
    public :: allocate_da_energy_sumcum, deallocate_da_energy_sumcum
    save
contains
    subroutine allocate_da_energy_sumcum()
      use system,  only: ncent_tot
      use force_pth, only: PTH
        if (.not. allocated(da_energy_cm2)) allocate (da_energy_cm2(3, ncent_tot, PTH))
        if (.not. allocated(da_energy_cum)) allocate (da_energy_cum(3, ncent_tot, PTH))
        if (.not. allocated(da_energy_sum)) allocate (da_energy_sum(3, ncent_tot, PTH))
        if (.not. allocated(da_psi_cum)) allocate (da_psi_cum(3, ncent_tot, PTH))
        if (.not. allocated(da_psi_sum)) allocate (da_psi_sum(3, ncent_tot, PTH))
    end subroutine allocate_da_energy_sumcum

    subroutine deallocate_da_energy_sumcum()
        if (allocated(da_psi_sum)) deallocate(da_psi_sum)
        if (allocated(da_psi_cum)) deallocate(da_psi_cum)
        if (allocated(da_energy_sum)) deallocate(da_energy_sum)
        if (allocated(da_energy_cum)) deallocate(da_energy_cum)
        if (allocated(da_energy_cm2)) deallocate(da_energy_cm2)
    end subroutine deallocate_da_energy_sumcum

end module da_energy_sumcum

module da_jastrow
    !> Arguments: da_d2j, da_j, da_vj
      use precision_kinds, only: dp

    implicit none

    real(dp), dimension(:, :), allocatable :: da_d2j !(3, MCENT)
    real(dp), dimension(:, :, :, :), allocatable :: da_j !(3, MELEC, MELEC, MCENT)
    real(dp), dimension(:, :, :, :), allocatable :: da_vj !(3, 3, MELEC, MCENT)

    private
    public   ::  da_d2j, da_j, da_vj
    public :: allocate_da_jastrow, deallocate_da_jastrow
    save
contains
    subroutine allocate_da_jastrow()
      use system,  only: ncent_tot,nelec
        if (.not. allocated(da_d2j)) allocate (da_d2j(3, ncent_tot))
        if (.not. allocated(da_j)) allocate (da_j(3, nelec, nelec, ncent_tot))
        if (.not. allocated(da_vj)) allocate (da_vj(3, 3, nelec, ncent_tot))
    end subroutine allocate_da_jastrow

    subroutine deallocate_da_jastrow
        if (allocated(da_vj)) deallocate(da_vj)
        if (allocated(da_j)) deallocate(da_j)
        if (allocated(da_d2j)) deallocate(da_d2j)
    end subroutine deallocate_da_jastrow

end module da_jastrow

module da_orbval
    !> Arguments: da_d2orb, da_dorb, da_orb
      use precision_kinds, only: dp
      use vmc_mod, only: norb_tot

    implicit none

    real(dp), dimension(:, :, :, :), allocatable :: da_d2orb !(3, MELEC, norb_tot, MCENT)
    real(dp), dimension(:, :, :, :, :), allocatable :: da_dorb !(3, 3, MELEC, norb_tot, MCENT)
    real(dp), dimension(:, :, :, :), allocatable :: da_orb !(3, MELEC, norb_tot, MCENT)

    private
    public   ::  da_d2orb, da_dorb, da_orb
    public :: allocate_da_orbval, deallocate_da_orbval
    save
contains
    subroutine allocate_da_orbval()
      use system, only: ncent, ncent_tot, nelec
      use vmc_mod, only: norb_tot
        if (.not. allocated(da_d2orb)) allocate (da_d2orb(3, nelec, norb_tot, ncent_tot))
        if (.not. allocated(da_dorb)) allocate (da_dorb(3, 3, nelec, norb_tot, ncent_tot))
        if (.not. allocated(da_orb)) allocate (da_orb(3, nelec, norb_tot, ncent_tot))
    end subroutine allocate_da_orbval

    subroutine deallocate_da_orbval()
        if (allocated(da_orb)) deallocate(da_orb)
        if (allocated(da_dorb)) deallocate(da_dorb)
        if (allocated(da_d2orb)) deallocate(da_d2orb)
    end subroutine deallocate_da_orbval

end module da_orbval

module da_pseudo
    !> Arguments: da_pecent, da_vps, da_nonloc, da_pe_en

    use precision_kinds, only: dp

    implicit none

    real(dp), dimension(:, :), allocatable :: da_pecent !(3, MCENT)
    real(dp), dimension(:, :, :, :), allocatable :: da_vps !(3, MELEC, MCENT, MPS_L)
    real(dp), dimension(:, :), allocatable :: da_nonloc !(3, MCENT)
    real(dp), dimension(:, :), allocatable :: da_pe_en  !(3, MCENT)

    private
    public   :: da_pecent, da_vps, da_nonloc, da_pe_en
    public :: allocate_da_pseudo, deallocate_da_pseudo
    save
contains
    subroutine allocate_da_pseudo()
      use pseudo_mod, only: MPS_L
      use system, only: ncent_tot, nelec
        if (.not. allocated(da_pecent)) allocate (da_pecent(3, ncent_tot))
        if (.not. allocated(da_vps)) allocate (da_vps(3, nelec, ncent_tot, MPS_L))
        if (.not. allocated(da_nonloc)) allocate (da_nonloc(3, ncent_tot))
        if (.not. allocated(da_pe_en)) allocate (da_pe_en(3, ncent_tot))

        da_nonloc = 0.0D0

    end subroutine allocate_da_pseudo

    subroutine deallocate_da_pseudo()
        if (allocated(da_nonloc)) deallocate(da_nonloc)
        if (allocated(da_vps)) deallocate(da_vps)
        if (allocated(da_pecent)) deallocate(da_pecent)
        if (allocated(da_pe_en)) deallocate(da_pe_en)
    end subroutine deallocate_da_pseudo

end module da_pseudo

module da_energy_now
    !> Arguments: da_energy, da_psi
      use precision_kinds, only: dp

    implicit none

    real(dp), dimension(:, :), allocatable :: da_energy !(3, MCENT)
    real(dp), dimension(:, :), allocatable :: da_psi !(3, MCENT)

    private
    public   ::  da_energy, da_psi
    public :: allocate_da_energy_now, deallocate_da_energy_now
    save
contains
    subroutine allocate_da_energy_now()
      use system, only: ncent_tot
        if (.not. allocated(da_energy)) allocate (da_energy(3, ncent_tot))
        if (.not. allocated(da_psi)) allocate (da_psi(3, ncent_tot))
    end subroutine allocate_da_energy_now

    subroutine deallocate_da_energy_now()
        if (allocated(da_psi)) deallocate(da_psi)
        if (allocated(da_energy)) deallocate(da_energy)
    end subroutine deallocate_da_energy_now

end module da_energy_now

module deloc_dj_m
    !> Arguments: denergy
    use mstates_mod, only: MSTATES
    use optwf_parms, only: nparmj
    use precision_kinds, only: dp

    implicit none

    real(dp), dimension(:, :), allocatable :: denergy !(nparmj, MSTATES)

    private
    public :: denergy
    public :: allocate_deloc_dj_m, deallocate_deloc_dj_m
    save
contains
    subroutine allocate_deloc_dj_m()
      use mstates_mod, only: MSTATES
      use optwf_parms, only: nparmj
        if (.not. allocated(denergy)) allocate (denergy(nparmj, MSTATES))
    end subroutine allocate_deloc_dj_m

    subroutine deallocate_deloc_dj_m()
        if (allocated(denergy)) deallocate(denergy)
    end subroutine deallocate_deloc_dj_m

end module deloc_dj_m

module denergy_det_m
    !> Arguments: denergy_det
      use precision_kinds, only: dp

    implicit none

    real(dp), dimension(:, :, :), allocatable :: denergy_det !(MDET, 2, nbjx)

    private
    public :: denergy_det
    public :: allocate_denergy_det_m, deallocate_denergy_det_m
    save
contains
    subroutine allocate_denergy_det_m()
      use slater, only: ndet
      use vmc_mod, only: nbjx
        if (.not. allocated(denergy_det)) allocate (denergy_det(ndet, 2, nbjx))
    end subroutine allocate_denergy_det_m

    subroutine deallocate_denergy_det_m()
        if (allocated(denergy_det)) deallocate(denergy_det)
    end subroutine deallocate_denergy_det_m

end module denergy_det_m

module derivjas
    !> Arguments: d2g, g, go, gvalue
    use precision_kinds, only: dp

    implicit none

    real(dp), dimension(:, :), allocatable :: d2g !(nparmj, nwftypejas)
    real(dp), dimension(:, :, :, :), allocatable :: g !(3, nelec, nparmj, nwftypejas)
    real(dp), dimension(:, :, :, :), allocatable :: go !(nelec, nelec, nparmj, nwftypejas)
    real(dp), dimension(:, :), allocatable :: gvalue !(nparmj, nwftypejas)

    private
    public   :: d2g, g, go, gvalue
    public :: allocate_derivjas, deallocate_derivjas
    save
contains
    subroutine allocate_derivjas()
      use optwf_parms, only: nparmj
      use system, only: nelec
      use vmc_mod, only: nwftypejas
        if (.not. allocated(d2g)) allocate (d2g(nparmj, nwftypejas))
        if (.not. allocated(g)) allocate (g(3, nelec, nparmj, nwftypejas))
        if (.not. allocated(go)) allocate (go(nelec, nelec, nparmj, nwftypejas))
        if (.not. allocated(gvalue)) allocate (gvalue(nparmj,nwftypejas))
    end subroutine allocate_derivjas

    subroutine deallocate_derivjas()
        if (allocated(gvalue)) deallocate(gvalue)
        if (allocated(go)) deallocate(go)
        if (allocated(g)) deallocate(g)
        if (allocated(d2g)) deallocate(d2g)
    end subroutine deallocate_derivjas

end module derivjas

module dorb_m
    !> Arguments: iworbd

    implicit none

    integer, dimension(:, :), allocatable :: iworbd !(MELEC, MDET)

    private
    public :: iworbd
    public :: allocate_dorb_m, deallocate_dorb_m
    save

contains

    subroutine allocate_dorb_m()
      use slater, only: ndet
      use system, only: nelec
        if (.not. allocated(iworbd)) allocate (iworbd(nelec, ndet), source=0)
    end subroutine allocate_dorb_m

    subroutine deallocate_dorb_m()
        if (allocated(iworbd)) deallocate(iworbd)
    end subroutine deallocate_dorb_m

end module dorb_m

module ijasnonlin
    !> Arguments: d1d2a, d1d2b, d2d2a, d2d2b
      use precision_kinds, only: dp

    implicit none

    real(dp), dimension(:), allocatable :: d1d2a !(MCTYPE)
    real(dp), dimension(:), allocatable :: d1d2b !(2)
    real(dp), dimension(:), allocatable :: d2d2a !(MCTYPE)
    real(dp), dimension(:), allocatable :: d2d2b !(2)

    private
    public :: d1d2a, d1d2b, d2d2a, d2d2b
    public :: allocate_ijasnonlin, deallocate_ijasnonlin
    save
contains
    subroutine allocate_ijasnonlin()
      use system, only: nctype_tot
        if (.not. allocated(d1d2a)) allocate (d1d2a(nctype_tot))
        if (.not. allocated(d1d2b)) allocate (d1d2b(2))
        if (.not. allocated(d2d2a)) allocate (d2d2a(nctype_tot))
        if (.not. allocated(d2d2b)) allocate (d2d2b(2))
    end subroutine allocate_ijasnonlin

    subroutine deallocate_ijasnonlin()
        if (allocated(d2d2b)) deallocate(d2d2b)
        if (allocated(d2d2a)) deallocate(d2d2a)
        if (allocated(d1d2b)) deallocate(d1d2b)
        if (allocated(d1d2a)) deallocate(d1d2a)
    end subroutine deallocate_ijasnonlin

end module ijasnonlin

module derivest
   !> DMC derivatives
   !> Arguments: derivcm2, derivcum, derivsum, derivtotave

   use precision_kinds, only: dp
   use system,  only: ncent_tot
   use force_pth, only: PTH

   implicit none


    real(dp), dimension(:,:,:), allocatable :: derivcm2 ! (3,NCENT,IPTH)
    real(dp), dimension(:,:,:,:), allocatable :: derivcum !(3,3,NCENT,IPTH)
    real(dp), dimension(:,:,:,:), allocatable :: derivsum !(3,3,NCENT,IPTH)
    real(dp), dimension(:,:,:), allocatable :: derivtotave ! (3,NCENT,IPTH)

    private
    public :: derivcm2, derivcum, derivsum, derivtotave
    public :: allocate_derivest, deallocate_derivest
    save

contains
    subroutine allocate_derivest()
        if (.not. allocated(derivcm2)) allocate(derivcm2(3,ncent_tot,PTH))
        if (.not. allocated(derivcum)) allocate(derivcum(3,3,ncent_tot, PTH))
        if (.not. allocated(derivsum)) allocate(derivsum(3,3,ncent_tot, PTH))
        if (.not. allocated(derivtotave)) allocate(derivtotave(3,ncent_tot, PTH))
    end subroutine allocate_derivest

    subroutine deallocate_derivest
        if (allocated(derivcm2)) deallocate(derivcm2)
        if (allocated(derivcum)) deallocate(derivcum)
        if (allocated(derivsum)) deallocate(derivsum)
        if (allocated(derivtotave)) deallocate(derivtotave)
    end subroutine deallocate_derivest
 end module derivest

module m_deriv
contains
subroutine allocate_m_deriv()
    use da_energy_now, only: allocate_da_energy_now
    use da_energy_sumcum, only: allocate_da_energy_sumcum
    use da_jastrow, only: allocate_da_jastrow
    use da_orbval, only: allocate_da_orbval
    use da_pseudo, only: allocate_da_pseudo
    use deloc_dj_m, only: allocate_deloc_dj_m
    use denergy_det_m, only: allocate_denergy_det_m
    use derivjas, only: allocate_derivjas
    use dorb_m, only: allocate_dorb_m
    use ijasnonlin, only: allocate_ijasnonlin

    implicit none

    call allocate_da_energy_sumcum()
    call allocate_da_jastrow()
    call allocate_da_orbval()
    call allocate_da_pseudo()
    call allocate_da_energy_now()
    call allocate_deloc_dj_m()
    call allocate_denergy_det_m()
    call allocate_derivjas()
    call allocate_dorb_m()
    call allocate_ijasnonlin()
end subroutine allocate_m_deriv

subroutine deallocate_m_deriv()
    use da_energy_now, only: deallocate_da_energy_now
    use da_energy_sumcum, only: deallocate_da_energy_sumcum
    use da_jastrow, only: deallocate_da_jastrow
    use da_orbval, only: deallocate_da_orbval
    use da_pseudo, only: deallocate_da_pseudo
    use deloc_dj_m, only: deallocate_deloc_dj_m
    use denergy_det_m, only: deallocate_denergy_det_m
    use derivjas, only: deallocate_derivjas
    use dorb_m, only: deallocate_dorb_m
    use ijasnonlin, only: deallocate_ijasnonlin

    implicit none

    call deallocate_da_energy_sumcum()
    call deallocate_da_jastrow()
    call deallocate_da_orbval()
    call deallocate_da_pseudo()
    call deallocate_da_energy_now()
    call deallocate_deloc_dj_m()
    call deallocate_denergy_det_m()
    call deallocate_derivjas()
    call deallocate_dorb_m()
    call deallocate_ijasnonlin()
end subroutine deallocate_m_deriv
end module
