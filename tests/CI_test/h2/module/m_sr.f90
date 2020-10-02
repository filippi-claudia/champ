module sr_mod
     !> Arguments:
     integer, parameter :: MPARM = 15100
     integer, parameter :: MOBS = 10 + 6*MPARM
     integer, parameter :: MCONF = 10000
     integer, parameter :: MVEC = 160

     private
     public :: MPARM, MOBS, MCONF, MVEC
     save
 end module sr_mod

 module sr_index
     !> Arguments: jelo, jelo2, jelohfj

     integer  :: jelo
     integer  :: jelo2
     integer  :: jelohfj

     private
     public :: jelo, jelo2, jelohfj
     save
 end module sr_index

 module sr_mat_n
     !> Arguments: elocal, h_sr, jefj, jfj, jhfj, nconf_n, obs, s_diag, s_ii_inv, sr_ho, sr_o, wtg, obs_tot
     use sr_mod, only: MPARM, MOBS, MCONF
     use precision_kinds, only: dp
     use mstates_mod, only: MSTATES

     real(dp) :: elocal(MCONF, MSTATES)
     real(dp) :: h_sr(MPARM)
     integer  :: jefj
     integer  :: jfj
     integer  :: jhfj
     integer  :: nconf_n
     real(dp) :: obs(MOBS, MSTATES)
     real(dp) :: s_diag(MPARM, MSTATES)
     real(dp) :: s_ii_inv(MPARM)
     real(dp) :: sr_ho(MPARM, MCONF)
     real(dp) :: sr_o(MPARM, MCONF)
     real(dp) :: wtg(MCONF, MSTATES)
     real(dp) :: obs_tot(MOBS, MSTATES)

     private
     public :: elocal, h_sr, jefj, jfj, jhfj, nconf_n, obs, s_diag, s_ii_inv, sr_ho, sr_o, wtg, obs_tot
     save
 end module sr_mat_n
