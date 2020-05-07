      subroutine acuest_reduce(enow)
c Written by Claudia Filippi

      use csfs, only: ccsf, cxdet, iadet, ibdet, icxdet, ncsf, nstates

      use est2cm, only: ecm2, ecm21, pecm2, r2cm2, tjfcm2, tpbcm2, avcm2
      use estcum, only: ecum, ecum1, iblk, pecum, r2cum, tjfcum, tpbcum, avcum
      use estpsi, only: apsi, aref, detref
      use estsig, only: ecm21s, ecum1s
      use estsum, only: acc, esum, esum1, pesum, r2sum, tjfsum, tpbsum
      use forcepar, only: deltot, istrech, nforce
      use forcest, only: fcm2, fcum
      use forcewt, only: wcum, wsum
      use mpiconf, only: idtask, nproc, wid
      implicit real*8(a-h,o-z)











      include 'vmc.h'
      include 'force.h'
      include 'mstates.h'
      include 'optorb.h'
      include 'mpif.h'

      parameter (MOBS=MSTATES*(8+5*MFORCE)+10)








      character*20 filename

      dimension obs(MOBS)
      dimension collect(MOBS),enow(MSTATES,MFORCE)

      iblk=iblk+nproc

      jo=0
      do 10 istate=1,nstates
        jo=jo+1
        obs(jo)=enow(istate,1)

        jo=jo+1
   10   obs(jo)=apsi(istate)

      jo=jo+1
      obs(jo)=aref

      do iab=1,2
        jo=jo+1
        obs(jo)=detref(iab)
      enddo

      do 20 ifr=1,nforce
        do 20 istate=1,nstates
          jo=jo+1
          obs(jo)=ecum(istate,ifr)

          jo=jo+1
          obs(jo)=ecm2(istate,ifr)

          jo=jo+1
          obs(jo)=fcum(istate,ifr)

          jo=jo+1
          obs(jo)=fcm2(istate,ifr)

          jo=jo+1
   20     obs(jo)=wcum(istate,ifr)

      do 30 i=1,MSTATES*3
        jo=jo+1
        obs(jo)=avcum(i)

        jo=jo+1
   30   obs(jo)=avcm2(i)
      
      jo=jo+1
      obs(jo)=acc

      jo_tot=jo

      if(jo_tot.gt.MOBS)  call fatal_error('ACUEST_REDUCE: increase MOBS')
 
      call mpi_reduce(obs,collect,jo_tot
     &,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

      call mpi_bcast(collect,jo_tot
     &     ,mpi_double_precision,0,MPI_COMM_WORLD,ierr)

      jo=0
      do 110 istate=1,nstates
        jo=jo+1
        enow(istate,1)=collect(jo)/nproc

        jo=jo+1
  110   apsi(istate)=collect(jo)/nproc
      
      jo=jo+1
      aref=collect(jo)/nproc
      
      do iab=1,2
        jo=jo+1
        detref(iab)=collect(jo)/nproc
      enddo

      do 120 ifr=1,nforce
        do 120 istate=1,nstates
          jo=jo+1
          ecum(istate,ifr)=collect(jo)

          jo=jo+1
          ecm2(istate,ifr)=collect(jo)

          jo=jo+1
          fcum(istate,ifr)=collect(jo)

          jo=jo+1
          fcm2(istate,ifr)=collect(jo)

          jo=jo+1
  120     wcum(istate,ifr)=collect(jo)

      do 130 i=1,MSTATES*3
        jo=jo+1
        avcum(i)=collect(jo)

        jo=jo+1
  130   avcm2(i)=collect(jo)
       
      jo=jo+1
      acollect=collect(jo)
       
c reduce properties
      call prop_reduce
c optimization reduced at the end of the run: large vectors to pass
c pcm averages reduced at the end of the run

      call mpi_barrier(MPI_COMM_WORLD,ierr)

      if(wid) then

        acc=acollect

c optorb reduced at the end of the run: set printout to 0
        iorbprt_sav=iorbprt
        iorbprt=0
        call acuest_write(enow,nproc)
        iorbprt=iorbprt_sav

       else

        do 210 ifr=1,nforce
          do 210 istate=1,nstates
           ecum(istate,ifr)=0
           ecm2(istate,ifr)=0
           fcum(istate,ifr)=0
           fcm2(istate,ifr)=0
  210      wcum(istate,ifr)=0

        do 220 i=1,MSTATES*3
          avcum(i)=0
  220     avcm2(i)=0

        acc=0
      endif

      return

      entry acues1_reduce

      call qpcm_update_vol(iupdate)
      if(iupdate.eq.1) then
        call pcm_reduce_chvol
        call pcm_compute_penupv
      endif

      return
      end
