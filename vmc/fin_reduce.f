      subroutine fin_reduce
c for compatibility with MPI version, written by Claudia Filippi

      implicit real*8(a-h,o-z)
      include 'vmc.h'
      include 'force.h'
      include 'mstates.h'

      common /contrl/ nstep,nblk,nblkeq,nconf,nconf_new,isite,idump,irstar
      common /estcum/ ecum1(MSTATES),ecum(MSTATES,MFORCE),pecum(MSTATES),tpbcum(MSTATES),tjfcum(MSTATES),r2cum,iblk
      common /forcewt/ wsum(MSTATES,MFORCE),wcum(MSTATES,MFORCE)

      efin=ecum1(1)/wcum(1,1)

      call optjas_fin(wcum(1,1),ecum1)

      call optci_fin(iblk,wcum(1,1),efin)

      call optorb_fin(wcum(1,1),ecum1)

      call optx_jas_ci_fin(wcum(1,1),efin)

      call optx_jas_orb_fin(wcum(1,1),ecum1)

      call optx_orb_ci_fin(wcum(1,1),efin)

      call pcm_fin(wcum(1,1),iblk)

      call mmpol_fin(wcum(1,1),iblk)

      call force_analy_fin(wcum(1,1),iblk,efin)

      return
      end
