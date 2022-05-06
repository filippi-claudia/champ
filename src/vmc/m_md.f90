module md_mass
    use vmc_mod, only : mctype
    use precision_kinds, only: dp
    use atom, only: ncent

    integer, dimension(MCTYPE) :: ntype(MCTYPE)
    real(dp), dimension(:), allocatable :: amass
    character,dimension(MCTYPE):: asymb(MCTYPE)

    private
    public :: ntype, asymb, amass
    public :: allocate_mass,deallocate_mass

    save

contains
    subroutine allocate_mass()
        use precision_kinds, only: dp
        use atom, only: ncent
        if (.not. allocated(amass)) allocate (amass(ncent))
    end subroutine allocate_mass

    subroutine deallocate_mass()
        if (allocated(amass)) deallocate (amass)
    end subroutine deallocate_mass
end module md_mass

module md_var
   use atom, only:ncent
   use precision_kinds, only: dp

    integer ::  md_tau
    real(dp) :: md_dt
    real(dp), dimension(:,:), allocatable :: pos, forces_ave
    real(dp), dimension(:,:), allocatable :: sigma,sigma_md
    real(dp), dimension(:,:), allocatable :: acc, acc_new, vel

    private
    public :: md_dt, md_tau
    public :: pos, vel, acc, acc_new
    public :: sigma_md, forces_ave, sigma
    public :: allocate_md, deallocate_md
    save
contains
    subroutine allocate_md()
        use precision_kinds, only: dp
        use atom, only: ncent
        if (.not. allocated(pos)) allocate (pos(3, ncent))
        if (.not. allocated(vel)) allocate (vel(3, ncent))
        if (.not. allocated(acc)) allocate (acc(3, ncent))
        if (.not. allocated(acc_new)) allocate (acc_new(3, ncent))
        if (.not. allocated(forces_ave)) allocate (forces_ave(3,ncent))
        if (.not. allocated(sigma)) allocate (sigma(3,ncent))
        if (.not. allocated(sigma_md)) allocate (sigma_md(3,ncent))
    end subroutine allocate_md

    subroutine deallocate_md()
        if (allocated(pos)) deallocate (pos)
        if (allocated(vel)) deallocate (vel)
        if (allocated(acc)) deallocate (acc)
        if (allocated(acc_new)) deallocate (acc_new)
        if (allocated(forces_ave)) deallocate (forces_ave)
        if (allocated(sigma)) deallocate (sigma)
        if (allocated(sigma_md)) deallocate (sigma_md)
    end subroutine deallocate_md

end module md_var

module md_fit

   use atom, only:ncent
   use precision_kinds, only: dp

    integer :: nfit
    real(dp), dimension(:,:), allocatable:: x_save(:,:,:)
    real(dp), dimension(:,:), allocatable:: x_save_up(:,:,:)
    real(dp), dimension(:,:), allocatable:: wnoi(:,:)
    real(dp), dimension(:,:), allocatable:: aY(:,:)
    real(dp), dimension(:,:), allocatable:: a_save(:,:,:)
    real(dp), dimension(:,:), allocatable:: a_save_up(:,:,:)

    private
    public :: nfit, x_save, x_save_up, wnoi
    public :: aY, a_save, a_save_up
    public :: allocate_mdfit, deallocate_mdfit
    save

contains
    subroutine allocate_mdfit()
        use precision_kinds, only: dp
        use atom, only: ncent
        if (.not. allocated(x_save)) allocate (x_save(3, ncent, nfit))
        if (.not. allocated(x_save_up)) allocate (x_save_up(3, ncent,nfit))
        if (.not. allocated(wnoi)) allocate (wnoi(3, ncent))
        if (.not. allocated(aY)) allocate (aY(3, ncent))
        if (.not. allocated(a_save)) allocate (a_save(3, ncent,nfit))
        if (.not. allocated(a_save_up)) allocate (a_save_up(3,ncent,nfit))
    end subroutine allocate_mdfit

    subroutine deallocate_mdfit()
        if (allocated(x_save)) deallocate (x_save)
        if (allocated(x_save_up)) deallocate (x_save_up)
        if (allocated(wnoi)) deallocate (wnoi)
        if (allocated(aY)) deallocate (aY)
        if (allocated(a_save)) deallocate (a_save)
        if (allocated(a_save_up)) deallocate (a_save_up)
    end subroutine deallocate_mdfit
end module md_fit
