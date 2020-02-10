!> \brief File collecting all modules that replace common blocks.  
!>
!> \author P. Lopez-Tarifa & F. Zapata NLeSC(2019)

 module precision_kinds
   ! named constants for 4, 2, and 1 byte integers:
   integer, parameter :: &
        i4b = selected_int_kind(9), &
        i2b = selected_int_kind(4), &
        i1b = selected_int_kind(2)
   ! single, double and quadruple precision reals:
   integer, parameter :: &
        sp = kind(1.0), &
        dp = selected_real_kind(2 * precision(1.0_sp)), &
        qp = selected_real_kind(2 * precision(1.0_dp))
 end module precision_kinds

 module atom
   !> Arguments: znuc, cent, pecent, iwctype, nctype, ncent
   use precision_kinds, only: dp  

   include 'vmc.h'

   real(dp) :: cent( 3, MCENT)
   real(dp) :: znuc( MCTYPE)
   real(dp) :: pecent
   integer  :: iwctype( MCENT), nctype, ncent

   private
   public   :: znuc, cent, pecent, iwctype, nctype, ncent 
   save
 end module atom
   
module const
    !> Arguments: pi, hb, etrial, delta, deltai, fbias, nelec, imetro, ipr
    use precision_kinds, only: dp
    include 'vmc.h'

    real(dp) :: delta
    real(dp) :: deltai
    real(dp) :: etrial
    real(dp) :: fbias
    real(dp) :: hb
    integer  :: imetro
    integer  :: ipr
    integer  :: nelec
    real(dp) :: pi

    save
end module const

 module contrl_per
   !> Arguments: iperiodic, ibasis 

   integer  :: iperiodic, ibasis

   private
   public   :: iperiodic, ibasis
   save
 end module contrl_per
 
 module da_pseudo
   !> Arguments: da_pecent, da_vps, da_nonloc  

   use precision_kinds, only: dp  

   include 'vmc.h'
   include 'pseudo.h'

   real(dp) :: da_pecent( 3, MCENT), da_vps( 3, MELEC, MCENT, MPS_L)
   real(dp) :: da_nonloc( 3, MCENT)= 0.0D0 

   private
   public   :: da_pecent, da_vps, da_nonloc 
   save
 end module da_pseudo 
 
 module force_analy 
   !> Arguments: iforce_analy 

   integer  :: iforce_analy 

   private
   public   :: iforce_analy 
   save
 end module force_analy 

module ghostatom
    !> Arguments: newghostype, nghostcent
    use precision_kinds, only: dp
    include 'vmc.h'

    integer  :: newghostype
    integer  :: nghostcent

    save
end module ghostatom

module jaspar
    !> Arguments: nspin1, nspin2, sspin, sspinn, is
    use precision_kinds, only: dp
    include 'vmc.h'

    integer  :: is
    integer  :: nspin1
    integer  :: nspin2
    real(dp) :: sspin
    real(dp) :: sspinn

    save
end module jaspar



module config
    !> Arguments: delttn, enew, eold, nearestn, nearesto, pen, peo, psi2n, psi2o, psido, psijo, rminn, rminno, rmino, rminon, rvminn, rvminno, rvmino, rvminon, tjfn, tjfo, tjfoo, vnew, vold, xnew, xold
    use precision_kinds, only: dp
    include 'vmc.h'
    include 'force.h'
    include 'mstates.h'

    real(dp) :: delttn(MELEC)
    real(dp) :: enew(MFORCE)
    real(dp) :: eold(MSTATES,MFORCE)
    integer  :: nearestn(MELEC)
    integer  :: nearesto(MELEC)
    real(dp) :: pen
    real(dp) :: peo(MSTATES)
    real(dp) :: psi2n(MFORCE)
    real(dp) :: psi2o(MSTATES,MFORCE)
    real(dp) :: psido(MSTATES)
    real(dp) :: psijo
    real(dp) :: rminn(MELEC)
    real(dp) :: rminno(MELEC)
    real(dp) :: rmino(MELEC)
    real(dp) :: rminon(MELEC)
    real(dp) :: rvminn(3,MELEC)
    real(dp) :: rvminno(3,MELEC)
    real(dp) :: rvmino(3,MELEC)
    real(dp) :: rvminon(3,MELEC)
    real(dp) :: tjfn
    real(dp) :: tjfo(MSTATES)
    integer  :: tjfoo
    real(dp) :: vnew(3,MELEC)
    real(dp) :: vold(3,MELEC)
    real(dp) :: xnew(3,MELEC)
    real(dp) :: xold(3,MELEC)

    save
end module config
