 module config
   !> looks a lot like sampling stuff
   !> Arguments: delttn, enew, eold, nearestn, nearesto, pen, peo,
   !> psi2n, psi2o, psido, psijo, rminn, rminno, rmino, rminon, rvminn,
   !> rvminno, rvmino, rvminon, tjfn, tjfo, tjfoo, vnew, vold, xnew, xold
   use force_mod, only: MFORCE
   use precision_kinds, only: dp
   use vmc_mod, only: MELEC
   use mstates_mod, only: MSTATES

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
   real(dp) :: tjfoo
   real(dp) :: vnew(3,MELEC)
   real(dp) :: vold(3,MELEC)
   real(dp) :: xnew(3,MELEC)
   real(dp) :: xold(3,MELEC)

   private
   public   :: delttn, enew, eold, nearestn, nearesto, pen, peo, psi2n
   public   :: psi2o, psido, psijo, rminn, rminno, rmino, rminon, rvminn
   public   :: rvminno, rvmino, rvminon, tjfn, tjfo, tjfoo, vnew, vold, xnew, xold
   save
 end module config
 
 module rnyucm
   !> I guess the random number generator
   !> used to move in the MC sampling
   !> Arguments: ll, lm

   integer  :: ll(4)
   integer  :: mm(4)
   data mm / 502,1521,4071,2107/
   data ll /   0,   0,   0,   1/

   private
   public :: ll, mm
   save
end module rnyucm

 module stats
   !> has to do with sampling
   !> Arguments: rejmax
   use precision_kinds, only: dp

   real(dp) :: rejmax

   private
   public :: rejmax
   save
 end module stats

 module step
   !> I guess has to do with the sampling
   !> Arguments: ekin, ekin2, rprob, suc, trunfb, try
   use precision_kinds, only: dp
   use vmc_mod, only: nrad

   real(dp) :: ekin(nrad)
   real(dp) :: ekin2(nrad)
   real(dp) :: rprob(nrad)
   real(dp) :: suc(nrad)
   real(dp) :: trunfb(nrad)
   real(dp) :: try(nrad)

   private
   public :: ekin, ekin2, rprob, suc, trunfb, try
   save
 end module step

 module tmpnode
   !> has to do with the sampling
   !> Arguments: distance_node_sum
   use precision_kinds, only: dp

   real(dp) :: distance_node_sum

   private
   public :: distance_node_sum
   save
 end module tmpnode

 module kinet
   !> kinetic energy ?
   !> only used in metropolis
   !> Arguments: dtdx2n, dtdx2o
   use precision_kinds, only: dp
   use vmc_mod, only: MELEC

   real(dp) :: dtdx2n(MELEC)
   real(dp) :: dtdx2o(MELEC)

   private 
   public :: dtdx2n, dtdx2o 
   save
 end module kinet