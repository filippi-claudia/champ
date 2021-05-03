module da_energy_sumcum
    !> Arguments: da_energy_cm2, da_energy_cum, da_energy_sum, da_psi_cum, da_psi_sum
    use precision_kinds, only: dp
    use vmc_mod, only: MCENT
    use mstates_mod, only: MSTATES

    real(dp), dimension(:, :, :), allocatable :: da_energy_cm2 !(3, MCENT,MSTATES)
    real(dp), dimension(:, :, :), allocatable :: da_energy_cum !(3, MCENT,MSTATES)
    real(dp), dimension(:, :, :), allocatable :: da_energy_sum !(3, MCENT,MSTATES)
    real(dp), dimension(:, :, :), allocatable :: da_psi_cum !(3, MCENT,MSTATES)
    real(dp), dimension(:, :, :), allocatable :: da_psi_sum !(3, MCENT,MSTATES)

    private
    public :: da_energy_cm2, da_energy_cum, da_energy_sum, da_psi_cum, da_psi_sum
    public :: allocate_da_energy_sumcum, deallocate_da_energy_sumcum
    save
contains
    subroutine allocate_da_energy_sumcum()
        use precision_kinds, only: dp
        use vmc_mod, only: MCENT
        if (.not. allocated(da_energy_cm2)) allocate (da_energy_cm2(3, MCENT, MSTATES))
        if (.not. allocated(da_energy_cum)) allocate (da_energy_cum(3, MCENT, MSTATES))
        if (.not. allocated(da_energy_sum)) allocate (da_energy_sum(3, MCENT, MSTATES))
        if (.not. allocated(da_psi_cum)) allocate (da_psi_cum(3, MCENT, MSTATES))
        if (.not. allocated(da_psi_sum)) allocate (da_psi_sum(3, MCENT, MSTATES))
    end subroutine allocate_da_energy_sumcum

    subroutine deallocate_da_energy_sumcum()
        if (allocated(da_psi_sum)) deallocate(da_psi_sum)
        if (allocated(da_psi_cum)) deallocate(da_psi_cum)
        if (allocated(da_energy_sum)) deallocate(da_energy_sum)
        if (allocated(da_energy_cum)) deallocate(da_energy_cum)
        if (allocated(da_energy_cm2)) deallocate(da_energy_cm2)
    end subroutine deallocate_da_energy_sumcum

end module da_energy_sumcum

module da_jastrow4val
    !> Arguments: da_d2j, da_j, da_vj
    use precision_kinds, only: dp
    use mstates_mod, only: MSTATES
    use vmc_mod, only: MELEC, MCENT

    real(dp), dimension(:, :, :, :), allocatable :: da_d2j !(3, MELEC, MCENT, MSTATES)
    real(dp), dimension(:, :, :, :), allocatable :: da_j !(3, MELEC, MCENT, MSTATES)
    real(dp), dimension(:, :, :, :, :), allocatable :: da_vj !(3, 3, MELEC, MCENT, MSTATES)

    private
    public   ::  da_d2j, da_j, da_vj
    public :: allocate_da_jastrow4val, deallocate_da_jastrow4val
    save
contains
    subroutine allocate_da_jastrow4val()
        use precision_kinds, only: dp
        use vmc_mod, only: MELEC, MCENT
        if (.not. allocated(da_d2j)) allocate (da_d2j(3, MELEC, MCENT, MSTATES))
        if (.not. allocated(da_j)) allocate (da_j(3, MELEC, MCENT, MSTATES))
        if (.not. allocated(da_vj)) allocate (da_vj(3, 3, MELEC, MCENT, MSTATES))
    end subroutine allocate_da_jastrow4val

    subroutine deallocate_da_jastrow4val()
        if (allocated(da_vj)) deallocate(da_vj)
        if (allocated(da_j)) deallocate(da_j)
        if (allocated(da_d2j)) deallocate(da_d2j)
    end subroutine deallocate_da_jastrow4val

end module da_jastrow4val

module da_orbval
    !> Arguments: da_d2orb, da_dorb, da_orb
    use precision_kinds, only: dp
    use vmc_mod, only: MELEC, MORB, MCENT
    use mstates_mod, only: MSTATES

    real(dp), dimension(:, :, :, :, :), allocatable :: da_d2orb !(3, MELEC, MORB, MCENT, MSTATES)
    real(dp), dimension(:, :, :, :, :, :), allocatable :: da_dorb !(3, 3, MELEC, MORB, MCENT, MSTATES)
    real(dp), dimension(:, :, :, :, :), allocatable :: da_orb !(3, MELEC, MORB, MCENT, MSTATES)

    private
    public   ::  da_d2orb, da_dorb, da_orb
    public :: allocate_da_orbval, deallocate_da_orbval
    save
contains
    subroutine allocate_da_orbval()
        use precision_kinds, only: dp
        use vmc_mod, only: MELEC, MORB, MCENT
        if (.not. allocated(da_d2orb)) allocate (da_d2orb(3, MELEC, MORB, MCENT, MSTATES))
        if (.not. allocated(da_dorb)) allocate (da_dorb(3, 3, MELEC, MORB, MCENT, MSTATES))
        if (.not. allocated(da_orb)) allocate (da_orb(3, MELEC, MORB, MCENT, MSTATES))
    end subroutine allocate_da_orbval

    subroutine deallocate_da_orbval()
        if (allocated(da_orb)) deallocate(da_orb)
        if (allocated(da_dorb)) deallocate(da_dorb)
        if (allocated(da_d2orb)) deallocate(da_d2orb)
    end subroutine deallocate_da_orbval

end module da_orbval

module da_pseudo
    !> Arguments: da_pecent, da_vps, da_nonloc

    use pseudo_mod, only: MPS_L
    use precision_kinds, only: dp
    use vmc_mod, only: MELEC, MCENT

    real(dp), dimension(:, :), allocatable :: da_pecent !(3, MCENT)
    real(dp), dimension(:, :, :, :), allocatable :: da_vps !(3, MELEC, MCENT, MPS_L)
    real(dp), dimension(:, :), allocatable :: da_nonloc !(3, MCENT)

    private
    public   :: da_pecent, da_vps, da_nonloc
    public :: allocate_da_pseudo, deallocate_da_pseudo
    save
contains
    subroutine allocate_da_pseudo()
        use pseudo_mod, only: MPS_L
        use precision_kinds, only: dp
        use vmc_mod, only: MELEC, MCENT
        if (.not. allocated(da_pecent)) allocate (da_pecent(3, MCENT))
        if (.not. allocated(da_vps)) allocate (da_vps(3, MELEC, MCENT, MPS_L))
        if (.not. allocated(da_nonloc)) allocate (da_nonloc(3, MCENT))

        da_nonloc = 0.0D0

    end subroutine allocate_da_pseudo

    subroutine deallocate_da_pseudo()
        if (allocated(da_nonloc)) deallocate(da_nonloc)
        if (allocated(da_vps)) deallocate(da_vps)
        if (allocated(da_pecent)) deallocate(da_pecent)
    end subroutine deallocate_da_pseudo

end module da_pseudo

module da_energy_now
    !> Arguments: da_energy, da_psi
    use precision_kinds, only: dp
    use vmc_mod, only: MCENT
    use mstates_mod, only: MSTATES

    real(dp), dimension(:, :, :), allocatable :: da_energy !(3, MCENT, MSTATES)
    real(dp), dimension(:, :, :), allocatable :: da_psi !(3, MCENT, MSTATES)

    private
    public   ::  da_energy, da_psi
    public :: allocate_da_energy_now, deallocate_da_energy_now
    save
contains
    subroutine allocate_da_energy_now()
        use precision_kinds, only: dp
        use vmc_mod, only: MCENT
        if (.not. allocated(da_energy)) allocate (da_energy(3, MCENT, MSTATES))
        if (.not. allocated(da_psi)) allocate (da_psi(3, MCENT, MSTATES))
    end subroutine allocate_da_energy_now

    subroutine deallocate_da_energy_now()
        if (allocated(da_psi)) deallocate(da_psi)
        if (allocated(da_energy)) deallocate(da_energy)
    end subroutine deallocate_da_energy_now

end module da_energy_now

module deloc_dj_m
    !> Arguments: denergy
    use optjas, only: MPARMJ
    use precision_kinds, only: dp
    use mstates_mod, only: MSTATES

    real(dp), dimension(:, :), allocatable :: denergy !(MPARMJ, MSTATES)

    private
    public :: denergy
    public :: allocate_deloc_dj_m, deallocate_deloc_dj_m
    save
contains
    subroutine allocate_deloc_dj_m()
        use optjas, only: MPARMJ
        use precision_kinds, only: dp
        use mstates_mod, only: MSTATES
        if (.not. allocated(denergy)) allocate (denergy(MPARMJ, MSTATES))
    end subroutine allocate_deloc_dj_m

    subroutine deallocate_deloc_dj_m()
        if (allocated(denergy)) deallocate(denergy)
    end subroutine deallocate_deloc_dj_m

end module deloc_dj_m

module denergy_det_m
    !> Arguments: denergy_det
    use precision_kinds, only: dp
    use vmc_mod, only: MDET
    use mstates_mod, only: MSTATES

    real(dp), dimension(:, :, :), allocatable :: denergy_det !(MDET, 2, MSTATES)

    private
    public :: denergy_det
    public :: allocate_denergy_det_m, deallocate_denergy_det_m
    save
contains
    subroutine allocate_denergy_det_m()
        use precision_kinds, only: dp
        use vmc_mod, only: MDET
        if (.not. allocated(denergy_det)) allocate (denergy_det(MDET, 2, MSTATES))
    end subroutine allocate_denergy_det_m

    subroutine deallocate_denergy_det_m()
        if (allocated(denergy_det)) deallocate(denergy_det)
    end subroutine deallocate_denergy_det_m

end module denergy_det_m

module denupdn
    !> Arguments: rprobdn, rprobup
    use precision_kinds, only: dp
    use vmc_mod, only: nrad

    real(dp), dimension(:), allocatable :: rprobdn !(nrad)
    real(dp), dimension(:), allocatable :: rprobup !(nrad)

    private
    public   ::  rprobdn, rprobup
    public :: allocate_denupdn, deallocate_denupdn
    save
contains
    subroutine allocate_denupdn()
        use precision_kinds, only: dp
        use vmc_mod, only: nrad
        if (.not. allocated(rprobdn)) allocate (rprobdn(nrad))
        if (.not. allocated(rprobup)) allocate (rprobup(nrad))
    end subroutine allocate_denupdn

    subroutine deallocate_denupdn()
        if (allocated(rprobup)) deallocate(rprobup)
        if (allocated(rprobdn)) deallocate(rprobdn)
    end subroutine deallocate_denupdn

end module denupdn

module derivjas
    !> Arguments: d2g, g, go, gvalue
    use optjas, only: MPARMJ
    use precision_kinds, only: dp
    use mstates_mod, only: MSTATES
    use vmc_mod, only: MELEC

    real(dp), dimension(:, :), allocatable :: d2g !(MPARMJ, MSTATES)
    real(dp), dimension(:, :, :, :), allocatable :: g !(3, MELEC, MPARMJ, MSTATES)
    real(dp), dimension(:, :, :, :), allocatable :: go !(MELEC, MELEC, MPARMJ, MSTATES)
    real(dp), dimension(:, :), allocatable :: gvalue !(MPARMJ, MSTATES)

    private
    public   :: d2g, g, go, gvalue
    public :: allocate_derivjas, deallocate_derivjas
    save
contains
    subroutine allocate_derivjas()
        use optjas, only: MPARMJ
        use precision_kinds, only: dp
        use vmc_mod, only: MELEC
        if (.not. allocated(d2g)) allocate (d2g(MPARMJ, MSTATES))
        if (.not. allocated(g)) allocate (g(3, MELEC, MPARMJ, MSTATES))
        if (.not. allocated(go)) allocate (go(MELEC, MELEC, MPARMJ, MSTATES))
        if (.not. allocated(gvalue)) allocate (gvalue(MPARMJ, MSTATES))
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
    use vmc_mod, only: MELEC, MDET

    integer, dimension(:, :), allocatable :: iworbd !(MELEC, MDET)

    private
    public :: iworbd
    public :: allocate_dorb_m, deallocate_dorb_m
    save
contains
    subroutine allocate_dorb_m()
        use vmc_mod, only: MELEC, MDET
        if (.not. allocated(iworbd)) allocate(iworbd(MELEC, MDET))
    end subroutine allocate_dorb_m

    subroutine deallocate_dorb_m()
        if (allocated(iworbd)) deallocate(iworbd)
    end subroutine deallocate_dorb_m

end module dorb_m

module ijasnonlin
    !> Arguments: d1d2a, d1d2b, d2d2a, d2d2b
    use precision_kinds, only: dp
    use vmc_mod, only: MCTYPE

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
        use precision_kinds, only: dp
        use vmc_mod, only: MCTYPE
        if (.not. allocated(d1d2a)) allocate(d1d2a(MCTYPE))
        if (.not. allocated(d1d2b)) allocate(d1d2b(2))
        if (.not. allocated(d2d2a)) allocate(d2d2a(MCTYPE))
        if (.not. allocated(d2d2b)) allocate(d2d2b(2))
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
   !> Arguments: derivcm2, derivcum, derivsum, derivtotave_num_old

   use force_mod, only: MFORCE
   use precision_kinds, only: dp

    real(dp), dimension(:), allocatable :: derivcm2 !(MFORCE)
    real(dp), dimension(:,:), allocatable :: derivcum !(10,MFORCE)
    real(dp), dimension(:,:), allocatable :: derivsum !(10,MFORCE)
    real(dp), dimension(:), allocatable :: derivtotave_num_old !(MFORCE)

    private
    public :: derivcm2, derivcum, derivsum, derivtotave_num_old
    public :: allocate_derivest, deallocate_derivest 
    save

contains
    subroutine allocate_derivest()
        if (.not. allocated(derivcm2)) allocate(derivcm2(MFORCE))
        if (.not. allocated(derivcum)) allocate(derivcum(10,MFORCE))
        if (.not. allocated(derivsum)) allocate(derivsum(10,MFORCE))
        if (.not. allocated(derivtotave_num_old)) allocate(derivtotave_num_old(MFORCE))
    end subroutine allocate_derivest

    subroutine deallocate_derivest
        if (allocated(derivcm2)) deallocate(derivcm2)
        if (allocated(derivcum)) deallocate(derivcum)
        if (allocated(derivsum)) deallocate(derivsum)
        if (allocated(derivtotave_num_old)) deallocate(derivtotave_num_old)
    end subroutine deallocate_derivest
 end module derivest

subroutine allocate_m_deriv()
    use da_energy_sumcum, only: allocate_da_energy_sumcum
    use da_jastrow4val, only: allocate_da_jastrow4val
    use da_orbval, only: allocate_da_orbval
    use da_pseudo, only: allocate_da_pseudo
    use da_energy_now, only: allocate_da_energy_now
    use deloc_dj_m, only: allocate_deloc_dj_m
    use denergy_det_m, only: allocate_denergy_det_m
    use denupdn, only: allocate_denupdn
    use derivjas, only: allocate_derivjas
    use dorb_m, only: allocate_dorb_m
    use ijasnonlin, only: allocate_ijasnonlin

    call allocate_da_energy_sumcum()
    call allocate_da_jastrow4val()
    call allocate_da_orbval()
    call allocate_da_pseudo()
    call allocate_da_energy_now()
    call allocate_deloc_dj_m()
    call allocate_denergy_det_m()
    call allocate_denupdn()
    call allocate_derivjas()
    call allocate_dorb_m()
    call allocate_ijasnonlin()
end subroutine allocate_m_deriv

subroutine deallocate_m_deriv()
    use da_energy_sumcum, only: deallocate_da_energy_sumcum
    use da_jastrow4val, only: deallocate_da_jastrow4val
    use da_orbval, only: deallocate_da_orbval
    use da_pseudo, only: deallocate_da_pseudo
    use da_energy_now, only: deallocate_da_energy_now
    use deloc_dj_m, only: deallocate_deloc_dj_m
    use denergy_det_m, only: deallocate_denergy_det_m
    use denupdn, only: deallocate_denupdn
    use derivjas, only: deallocate_derivjas
    use dorb_m, only: deallocate_dorb_m
    use ijasnonlin, only: deallocate_ijasnonlin

    call deallocate_da_energy_sumcum()
    call deallocate_da_jastrow4val()
    call deallocate_da_orbval()
    call deallocate_da_pseudo()
    call deallocate_da_energy_now()
    call deallocate_deloc_dj_m()
    call deallocate_denergy_det_m()
    call deallocate_denupdn()
    call deallocate_derivjas()
    call deallocate_dorb_m()
    call deallocate_ijasnonlin()
end subroutine deallocate_m_deriv
