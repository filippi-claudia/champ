      subroutine acues1_reduce
      implicit real*8(a-h,o-z)
      include 'dmc.h'
      include 'vmc.h'
      include 'force.h'

      common /estcum/ wcum,w_acc_cum,wfcum,wgcum(MFORCE),wg_acc_cum,wdcum,
     &wgdcum, wcum1,w_acc_cum1,wfcum1,wgcum1(MFORCE),wg_acc_cum1,
     &wdcum1, ecum,efcum,egcum(MFORCE),ecum1,efcum1,egcum1(MFORCE),
     &ei1cum,ei2cum,ei3cum, pecum(MFORCE),tpbcum(MFORCE),tjfcum(MFORCE),r2cum,
     &ricum,taucum(MFORCE)

      efin=egcum1(1)/wgcum1(1)

      call optjas_fin(wgcum1(1),egcum1(1))
      call optci_fin(iblk,wgcum1(1),efin)
      call optorb_fin(wgcum1(1),egcum1(1))
      call optx_jas_ci_fin(wgcum1(1),efin)
      call optx_jas_orb_fin(wgcum1(1),egcum1(1))
      call optx_orb_ci_fin(wgcum1(1),efin)

      return
      end
