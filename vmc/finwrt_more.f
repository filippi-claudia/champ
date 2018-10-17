      subroutine finwrt_more

      implicit real*8(a-h,o-z)

      include 'vmc.h'
      include 'force.h'
      include 'mstates.h'

      common /csfs/ ccsf(MDET,MSTATES,MWF),cxdet(MDET*MDETCSFX)
     &,icxdet(MDET*MDETCSFX),iadet(MDET),ibdet(MDET),ncsf,nstates

      common /contrl/ nstep,nblk,nblkeq,nconf,nconf_new,isite,idump,irstar
      common /estsum/ esum1(MSTATES),esum(MSTATES,MFORCE),pesum(MSTATES),tpbsum(MSTATES),tjfsum(MSTATES),r2sum,acc
      common /estcum/ ecum1(MSTATES),ecum(MSTATES,MFORCE),pecum(MSTATES),tpbcum(MSTATES),tjfcum(MSTATES),r2cum,iblk
      common /est2cm/ ecm21(MSTATES),ecm2(MSTATES,MFORCE),pecm2(MSTATES),tpbcm2(MSTATES),tjfcm2(MSTATES),r2cm2
      common /estsig/ ecum1s(MSTATES),ecm21s(MSTATES)
      common /estpsi/ apsi(MSTATES),aref

      passes=dfloat(iblk*nstep)

      write(6,'(''average psid, det_ref '',2d12.5)') (apsi(istate)/passes,istate=1,nstates),aref/passes

      return
      end
