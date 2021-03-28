      subroutine acuest_reduce(enow)
c Written by Claudia Filippi
      use force_mod, only: MFORCE
      use csfs, only: nstates
      use mstates_mod, only: MSTATES
      use estcum, only: ecum, pecum, tpbcum, tjfcum, iblk
      use est2cm, only: ecm2, pecm2, tpbcm2, tjfcm2
      use estpsi, only: apsi, aref, detref
      use estsum, only: acc
      use forcepar, only: nforce
      use forcest, only: fcm2, fcum
      use forcewt, only: wcum
      use mpiconf, only: nproc, wid
      use mpi

      ! this in not even in the master as the line
      ! is commented in optorb.h !
      use optorb_cblock, only: iorbprt, iorbprt_sav

      implicit real*8(a-h,o-z)

      parameter (MOBS=MSTATES*(8+5*MFORCE)+10)
      character*20 filename
      dimension obs(MOBS)
      dimension collect(MOBS),enow(MSTATES,MFORCE)
      ! ipudate was not declared anywhere
      integer :: iupdate = 0

      iblk=iblk+nproc

      jo=0
      do istate=1,nstates
         jo=jo+1
         obs(jo)=enow(istate,1)

         jo=jo+1
         obs(jo)=apsi(istate)
      enddo

      jo=jo+1
      obs(jo)=aref

      do istate=1,nstates
         do iab=1,2
            jo=jo+1
            obs(jo)=detref(iab,istate)
         enddo
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

      do 30 i=1,nstates
        jo=jo+1
        obs(jo)=pecum(i)

        jo=jo+1
        obs(jo)=tpbcum(i)

        jo=jo+1
        obs(jo)=tjfcum(i)

        jo=jo+1
        obs(jo)=pecm2(i)

        jo=jo+1
        obs(jo)=tpbcm2(i)

        jo=jo+1
   30   obs(jo)=tjfcm2(i)
      
      jo=jo+1
      obs(jo)=acc

      jo_tot=jo

      if(jo_tot.gt.MOBS) then
         call fatal_error('ACUEST_REDUCE: increase MOBS')
      endif
 
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
      
      do istate=1,nstates
         do iab=1,2
           jo=jo+1
           detref(iab,istate)=collect(jo)/nproc
         enddo
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

      do 130 i=1,nstates
        jo=jo+1
        pecum(i)=collect(jo)

        jo=jo+1
        tpbcum(i)=collect(jo)

        jo=jo+1
        tjfcum(i)=collect(jo)

        jo=jo+1
        pecm2(i)=collect(jo)

        jo=jo+1
        tpbcm2(i)=collect(jo)

        jo=jo+1
  130   tjfcm2(i)=collect(jo)
       
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

        do 220 i=1,nstates
          pecum(i)=0
          tpbcum(i)=0
          tjfcum(i)=0
          pecm2(i)=0
          tpbcm2(i)=0
  220     tjfcm2(i)=0

        acc=0
      endif

      return

      entry acues1_reduce
      
      call qpcm_update_vol()
      if(iupdate.eq.1) then
        call pcm_reduce_chvol
        call pcm_compute_penupv
      endif

      return
      end
