 module optwf_contrl
   !> Arguments: ioptci, ioptjas, ioptorb, nparm

    integer  :: ioptci
    integer  :: ioptjas
    integer  :: ioptorb
    integer  :: nparm

    private
    public :: ioptci, ioptjas, ioptorb, nparm
    save
 end module optwf_contrl

 module optwf_corsam
   !> Arguments: add_diag_tmp, energy, energy_err, force, force_err
   use force_mod, only: MFORCE
   use precision_kinds, only: dp

   real(dp) :: add_diag(MFORCE)
   real(dp) :: add_diag_tmp(MFORCE)
   real(dp) :: energy(MFORCE)
   real(dp) :: energy_err(MFORCE)
   real(dp) :: force(MFORCE)
   real(dp) :: force_err(MFORCE)

   private
   public :: add_diag, add_diag_tmp, energy, energy_err, force, force_err
   save
 end module optwf_corsam

 module optwf_func
   !> Arguments: ifunc_omega, omega, omega_hes
   use precision_kinds, only: dp

   integer  :: ifunc_omega
   real(dp) :: omega
   real(dp) :: omega_hes

   private
   public :: ifunc_omega, omega, omega_hes
   save
 end module optwf_func

 module optwf_nparmj
   !> Arguments: nparma, nparmb, nparmc, nparmf
   use vmc, only: MCTYPE, MCTYP3X

   integer  :: nparma(MCTYP3X)
   integer  :: nparmb(3)
   integer  :: nparmc(MCTYPE)
   integer  :: nparmf(MCTYPE)

   private
   public :: nparma, nparmb, nparmc, nparmf
   save
 end module optwf_nparmj

 module optwf_parms
   !> Arguments: nparmd, nparme, nparmg, nparmj, nparml, nparms

   integer  :: nparmd
   integer  :: nparme
   integer  :: nparmg
   integer  :: nparmj
   integer  :: nparml
   integer  :: nparms

   private
   public :: nparmd, nparme, nparmg, nparmj, nparml, nparms
   save
 end module optwf_parms

 module optwf_wjas
   !> Arguments: iwjasa, iwjasb, iwjasc, iwjasf
   use vmc, only: MCTYPE, MCTYP3X

   integer  :: iwjasa(83,MCTYP3X)
   integer  :: iwjasb(83,3)
   integer  :: iwjasc(83,MCTYPE)
   integer  :: iwjasf(15,MCTYPE)

   private
   public :: iwjasa, iwjasb, iwjasc, iwjasf
   save
 end module optwf_wjas