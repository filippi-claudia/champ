      subroutine acuest_reduce(enow)
c Written by Claudia Filippi

      use mpi
      use precision_kinds, only: dp
      use force_mod, only: MFORCE
      use csfs, only: nstates
      use mstates_mod, only: MSTATES

      use est2cm, only: ecm2, avcm2
      use estcum, only: ecum, iblk, avcum
      use estpsi, only: apsi, aref, detref
      use estsum, only: acc
      use forcepar, only: nforce
      use forcest, only: fcm2, fcum
      use forcewt, only: wcum
      use mpiconf, only: nproc, wid

      ! this in not even in the master as the line
      ! is commented in optorb.h !
      use optorb_cblock, only: iorbprt, iorbprt_sav

      implicit real*8(a-h,o-z)


      integer MOBS
      integer iupdate
      character*20 filename

      dimension enow(MSTATES,MFORCE)
      !real, dimension(MSTATES, MFORCE), INTENT(INOUT) :: enow

      real(dp), dimension(:), allocatable  :: local_obs
      real(dp), dimension(:), allocatable  :: collect
      

      MOBS = MSTATES*(8+5*MFORCE)+10
      allocate(local_obs(MOBS))
      allocate(collect(MOBS))
      

      ! ipudate was not declared anywhere
      iupdate = 0

      iblk=iblk+nproc

      jo=0
      do 10 istate=1,nstates
        jo=jo+1
        local_obs(jo)=enow(istate,1)

        jo=jo+1
   10   local_obs(jo)=apsi(istate)

      jo=jo+1
      local_obs(jo)=aref

      do iab=1,2
        jo=jo+1
        local_obs(jo)=detref(iab)
      enddo

      do 20 ifr=1,nforce
        do 20 istate=1,nstates
          jo=jo+1
          local_obs(jo)=ecum(istate,ifr)

          jo=jo+1
          local_obs(jo)=ecm2(istate,ifr)

          jo=jo+1
          local_obs(jo)=fcum(istate,ifr)

          jo=jo+1
          local_obs(jo)=fcm2(istate,ifr)

          jo=jo+1
   20     local_obs(jo)=wcum(istate,ifr)

      do 30 i=1,nstates*3
        jo=jo+1
        local_obs(jo)=avcum(i)

        jo=jo+1
   30   local_obs(jo)=avcm2(i)
      
      jo=jo+1
      local_obs(jo)=acc

      jo_tot=jo

      if(jo_tot.gt.MOBS)  call fatal_error('ACUEST_REDUCE: increase MOBS')
 
      call mpi_reduce(local_obs,collect,jo_tot
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

      do 130 i=1,nstates*3
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

        do 220 i=1,nstates*3
          avcum(i)=0
  220     avcm2(i)=0

        acc=0
      endif

      if (allocated(local_obs)) deallocate(local_obs)
      if (allocated(collect)) deallocate(collect)

      return

      entry acues1_reduce
      
      call qpcm_update_vol()
      if(iupdate.eq.1) then
        call pcm_reduce_chvol
        call pcm_compute_penupv
      endif

      return

      end
