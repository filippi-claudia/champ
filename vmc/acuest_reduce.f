      subroutine acuest_reduce(enow)
c Written by Claudia Filippi

      implicit real*8(a-h,o-z)

      include 'vmc.h'
      include 'force.h'
      include 'mstates.h'

      common /estcum/ ecum1(MSTATES),ecum(MSTATES,MFORCE),pecum(MSTATES),tpbcum(MSTATES),tjfcum(MSTATES),r2cum,iblk

      dimension enow(MSTATES,MFORCE)

      iblk=iblk+1

      call acuest_write(enow,1)

      return

      entry acues1_reduce

      call qpcm_update_vol(iupdate)

      if(iupdate.eq.1) call pcm_compute_penupv
        
      return
      end
