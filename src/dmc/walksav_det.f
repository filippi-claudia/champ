      subroutine walksav_det(iw)
c Written by Claudia Filippi

      use precision_kinds, only: dp
      use vmc_mod, only: MORB, MDET
      use vmc_mod, only: MMAT_DIM
      use vmc_mod, only: MEXCIT
      use dmc_mod, only: MWALK

      use const, only: nelec
      use vmc_mod, only: MORB, MDET, MCTYP3X
      use vmc_mod, only: NSPLIN, nrad, MORDJ, MORDJ1, MMAT_DIM, MMAT_DIM2, MMAT_DIM20
      use vmc_mod, only: radmax, delri, NEQSX, MTERMS, MCENT3, NCOEF, MEXCIT
      use dmc_mod, only: MWALK, MFPROD, MFPRD1, MPATH
      use const, only: delta, deltai, etrial, fbias, hb, imetro, ipr, nelec, pi
      use forcepar, only: deltot, istrech, nforce
      use force_mod, only: MFORCE, MFORCE_WT_PRD
      use mstates_mod, only: MSTATES, MDETCSFX
      use branch, only: eest, eigv, eold, ff, fprod, nwalk, pwt, wdsumo
      use branch, only: wgdsumo, wt, wtgen, wthist
      use slater, only: d2dx2, ddx, fp, fpp, slmi
      use dets, only: cdet, ndet
      use elec, only: ndn, nup
      use orbval, only: dorb, orb
      use coefs, only: norb
      use csfs, only: nstates
      use ycompact, only: ymat
      use multislater, only: detiab
      use multidet, only: ivirt, kref, numrep_det
      use multimat, only: aa, wfmat
      use mpi

      implicit none

      integer :: i, iab, ierr, iorb, irecv
      integer :: irequest, irequest_array, isend, istate
      integer :: istatus, itag, iw, iw2
      integer :: j, k, kk
      integer :: ndim, nel


      integer, allocatable, save :: krefw(:)
      real(dp), allocatable, save :: slmuiw(:, :)
      real(dp), allocatable, save :: slmdiw(:, :)
      real(dp), allocatable, save :: fpuw(:, :, :)
      real(dp), allocatable, save :: fpdw(:, :, :)
      real(dp), allocatable, save :: fppuw(:, :)
      real(dp), allocatable, save :: fppdw(:, :)
      real(dp), allocatable, save :: ddxw(:, :, :)
      real(dp), allocatable, save :: d2dx2w(:, :)
      real(dp), allocatable, save :: detuw(:, :)
      real(dp), allocatable, save :: detdw(:, :)
      real(dp), allocatable, save :: aaw(:,:,:,:)
      real(dp), allocatable, save :: wfmatw(:,:,:,:)
      real(dp), allocatable, save :: ymatw(:,:,:,:,:)
      real(dp), allocatable, save :: orbw(:,:,:)
      real(dp), allocatable, save :: dorbw(:,:,:,:)

      dimension istatus(MPI_STATUS_SIZE)
      dimension irequest_array(MPI_STATUS_SIZE)

      if(.not.allocated(aaw)) allocate(aaw(nelec,MORB,MWALK,2))
      if(.not.allocated(wfmatw)) allocate(wfmatw(MEXCIT**2,MDET,MWALK,2))
      if(.not.allocated(ymatw)) allocate(ymatw(MORB,nelec,MWALK,2,MSTATES))
      if(.not.allocated(orbw)) allocate(orbw(nelec,MORB,MWALK))
      if(.not.allocated(dorbw)) allocate(dorbw(3,nelec,MORB,MWALK))

      if(.not.allocated(krefw)) allocate(krefw(MWALK))
      if(.not.allocated(slmuiw)) allocate(slmuiw(MMAT_DIM,MWALK))
      if(.not.allocated(slmdiw)) allocate(slmdiw(MMAT_DIM,MWALK))
      if(.not.allocated(fpuw)) allocate(fpuw(3, MMAT_DIM,MWALK))
      if(.not.allocated(fpdw)) allocate(fpdw(3, MMAT_DIM,MWALK))
      if(.not.allocated(fppuw)) allocate(fppuw(MMAT_DIM,MWALK))
      if(.not.allocated(fppdw)) allocate(fppdw(MMAT_DIM,MWALK))
      if(.not.allocated(ddxw)) allocate(ddxw(3, nelec,MWALK))
      if(.not.allocated(d2dx2w)) allocate(d2dx2w(nelec,MWALK))
      if(.not.allocated(detuw)) allocate(detuw(MDET,MWALK))
      if(.not.allocated(detdw)) allocate(detdw(MDET,MWALK))



       do 20 k=1,ndet
         detuw(k,iw)=detiab(k,1)
   20    detdw(k,iw)=detiab(k,2)

       krefw(iw)=kref
       do 40 j=1,nup*nup
         slmuiw(j,iw)=slmi(j,1)
         fpuw(1,j,iw)=fp(1,j,1)
         fpuw(2,j,iw)=fp(2,j,1)
   40    fpuw(3,j,iw)=fp(3,j,1)
       do 50 j=1,ndn*ndn
         slmdiw(j,iw)=slmi(j,2)
         fpdw(1,j,iw)=fp(1,j,2)
         fpdw(2,j,iw)=fp(2,j,2)
   50    fpdw(3,j,iw)=fp(3,j,2)
       do 55 i=1,nelec
         ddxw(1,i,iw)=ddx(1,i)
         ddxw(2,i,iw)=ddx(2,i)
   55    ddxw(3,i,iw)=ddx(3,i)

       do 62 iab=1,2
         nel=nup
         if(iab.eq.2) nel=ndn
         do 58 j=ivirt(iab),norb
          do 58 i=1,nel
            do 57 istate=1,nstates
   57         ymatw(j,i,iw,iab,istate)=ymat(j,i,iab,istate)
   58       aaw(i,j,iw,iab)=aa(i,j,iab)
          do 62 k=1,ndet
            if(k.ne.kref) then
              ndim=numrep_det(k,iab)
              do 60 i=1,ndim*ndim
   60           wfmatw(i,k,iw,iab)=wfmat(i,k,iab)
            endif
   62  continue

       do 63 i=1,nelec
         do 63 iorb=1,norb
           orbw(i,iorb,iw)=orb(i,iorb)
           do 63 kk=1,3
   63        dorbw(kk,i,iorb,iw)=dorb(kk,i,iorb)

      return

      entry walkstrdet(iw)

      do 70 k=1,ndet
        detiab(k,1)=detuw(k,iw)
   70   detiab(k,2)=detdw(k,iw)

      kref=krefw(iw)
      do 80 j=1,nup*nup
        slmi(j,1)=slmuiw(j,iw)
        fp(1,j,1)=fpuw(1,j,iw)
        fp(2,j,1)=fpuw(2,j,iw)
   80   fp(3,j,1)=fpuw(3,j,iw)
      do 90 j=1,ndn*ndn
        slmi(j,2)=slmdiw(j,iw)
        fp(1,j,2)=fpdw(1,j,iw)
        fp(2,j,2)=fpdw(2,j,iw)
   90   fp(3,j,2)=fpdw(3,j,iw)
      do 95 i=1,nelec
        ddx(1,i)=ddxw(1,i,iw)
        ddx(2,i)=ddxw(2,i,iw)
   95   ddx(3,i)=ddxw(3,i,iw)

       do 102 iab=1,2
         nel=nup
         if(iab.eq.2) nel=ndn
         do 98 j=ivirt(iab),norb
          do 98 i=1,nel
            do 97 istate=1,nstates
   97         ymat(j,i,iab,istate)=ymatw(j,i,iw,iab,istate)
   98       aa(i,j,iab)=aaw(i,j,iw,iab)
          do 102 k=1,ndet
            if(k.ne.kref) then
              ndim=numrep_det(k,iab)
              do 100 i=1,ndim*ndim
  100           wfmat(i,k,iab)=wfmatw(i,k,iw,iab)
            endif
  102  continue

       do 103 i=1,nelec
         do 103 iorb=1,norb
           orb(i,iorb)=orbw(i,iorb,iw)
           do 103 kk=1,3
  103        dorb(kk,i,iorb)=dorbw(kk,i,iorb,iw)

      return

      entry splitjdet(iw,iw2)

      do 110 k=1,ndet
        detuw(k,iw2)=detuw(k,iw)
  110   detdw(k,iw2)=detdw(k,iw)

      krefw(iw2)=krefw(iw)
      do 120 j=1,nup*nup
        slmuiw(j,iw2)=slmuiw(j,iw)
        fpuw(1,j,iw2)=fpuw(1,j,iw)
        fpuw(2,j,iw2)=fpuw(2,j,iw)
  120   fpuw(3,j,iw2)=fpuw(3,j,iw)
      do 130 j=1,ndn*ndn
        slmdiw(j,iw2)=slmdiw(j,iw)
        fpdw(1,j,iw2)=fpdw(1,j,iw)
        fpdw(2,j,iw2)=fpdw(2,j,iw)
  130   fpdw(3,j,iw2)=fpdw(3,j,iw)
      do 135 i=1,nelec
        ddxw(1,i,iw2)=ddxw(1,i,iw)
        ddxw(2,i,iw2)=ddxw(2,i,iw)
  135   ddxw(3,i,iw2)=ddxw(3,i,iw)

       do 142 iab=1,2
         nel=nup
         if(iab.eq.2) nel=ndn
         do 138 j=ivirt(iab),norb
          do 138 i=1,nel
            do 137 istate=1,nstates
  137         ymatw(j,i,iw2,iab,istate)=ymatw(j,i,iw,iab,istate)
  138       aaw(i,j,iw2,iab)=aaw(i,j,iw,iab)
          do 142 k=1,ndet
            if(k.ne.krefw(iw)) then
              ndim=numrep_det(k,iab)
              do 140 i=1,ndim*ndim
  140           wfmatw(i,k,iw2,iab)=wfmatw(i,k,iw,iab)
            endif
  142  continue

       do 143 i=1,nelec
         do 143 iorb=1,norb
           orbw(i,iorb,iw2)=orbw(i,iorb,iw)
           do 143 kk=1,3
  143        dorbw(kk,i,iorb,iw2)=dorbw(kk,i,iorb,iw)

      return

      entry send_det(irecv)

      itag=0
      call mpi_isend(detuw(1,nwalk),ndet,mpi_double_precision,irecv
     &,itag+1,MPI_COMM_WORLD,irequest,ierr)
      call mpi_isend(detdw(1,nwalk),ndet,mpi_double_precision,irecv
     &,itag+2,MPI_COMM_WORLD,irequest,ierr)
      itag=itag+2

      call mpi_isend(slmuiw(1,nwalk),nup*nup,mpi_double_precision
     &  ,irecv,itag+1,MPI_COMM_WORLD,irequest,ierr)
      call mpi_isend(fpuw(1,1,nwalk),3*nup*nup,mpi_double_precision
     &  ,irecv,itag+2,MPI_COMM_WORLD,irequest,ierr)
      call mpi_isend(slmdiw(1,nwalk),ndn*ndn,mpi_double_precision
     &  ,irecv,itag+3,MPI_COMM_WORLD,irequest,ierr)
      call mpi_isend(fpdw(1,1,nwalk),3*ndn*ndn,mpi_double_precision
     &  ,irecv,itag+4,MPI_COMM_WORLD,irequest,ierr)
      call mpi_isend(ddxw(1,1,nwalk),3*nelec,mpi_double_precision
     &  ,irecv,itag+5,MPI_COMM_WORLD,irequest,ierr)
      call mpi_isend(krefw(nwalk),1,mpi_integer
     &  ,irecv,itag+6,MPI_COMM_WORLD,irequest,ierr)
      itag=itag+6

      call mpi_isend(aaw(1,1,nwalk,1),nelec*norb,mpi_double_precision
     & ,irecv,itag+1,MPI_COMM_WORLD,irequest,ierr)
      call mpi_isend(aaw(1,1,nwalk,2),nelec*norb,mpi_double_precision
     & ,irecv,itag+2,MPI_COMM_WORLD,irequest,ierr)
      itag=itag+2

      do 150 istate=1,nstates
        do 150 iab=1,2
        itag=itag+1
 150    call mpi_isend(ymatw(1,1,nwalk,iab,istate),MORB*nelec,mpi_double_precision
     &   ,irecv,itag,MPI_COMM_WORLD,irequest,ierr)

      do 160 iab=1,2
        do 160 k=1,ndet
          ndim=numrep_det(k,iab)
          if(k.ne.krefw(nwalk).and.ndim.gt.0) then
            itag=itag+1
            call mpi_isend(wfmatw(1,k,nwalk,iab),ndim*ndim,mpi_double_precision
     &     ,irecv,itag,MPI_COMM_WORLD,irequest,ierr)
          endif
 160  continue

      call mpi_isend(orbw(1,1,nwalk),nelec*norb,mpi_double_precision
     &  ,irecv,itag+1,MPI_COMM_WORLD,irequest,ierr)
      call mpi_isend(dorbw(1,1,1,nwalk),3*nelec*norb,mpi_double_precision
     &  ,irecv,itag+2,MPI_COMM_WORLD,irequest,ierr)
      itag=itag+2

      return

      entry recv_det(isend)

      itag=0
      call mpi_recv(detuw(1,nwalk),ndet,mpi_double_precision,isend
     &,itag+1,MPI_COMM_WORLD,istatus,ierr)
      call mpi_recv(detdw(1,nwalk),ndet,mpi_double_precision,isend
     &,itag+2,MPI_COMM_WORLD,istatus,ierr)
      itag=itag+2

      call mpi_recv(slmuiw(1,nwalk),nup*nup,mpi_double_precision
     &  ,isend,itag+1,MPI_COMM_WORLD,istatus,ierr)
      call mpi_recv(fpuw(1,1,nwalk),3*nup*nup,mpi_double_precision
     &  ,isend,itag+2,MPI_COMM_WORLD,istatus,ierr)
      call mpi_recv(slmdiw(1,nwalk),ndn*ndn,mpi_double_precision
     &  ,isend,itag+3,MPI_COMM_WORLD,istatus,ierr)
      call mpi_recv(fpdw(1,1,nwalk),3*ndn*ndn,mpi_double_precision
     &  ,isend,itag+4,MPI_COMM_WORLD,istatus,ierr)
      call mpi_recv(ddxw(1,1,nwalk),3*nelec,mpi_double_precision
     &  ,isend,itag+5,MPI_COMM_WORLD,istatus,ierr)
      call mpi_recv(krefw(nwalk),1,mpi_integer
     &  ,isend,itag+6,MPI_COMM_WORLD,irequest_array,ierr)
      itag=itag+6

      call mpi_recv(aaw(1,1,nwalk,1),nelec*norb,mpi_double_precision
     & ,isend,itag+1,MPI_COMM_WORLD,istatus,ierr)
      call mpi_recv(aaw(1,1,nwalk,2),nelec*norb,mpi_double_precision
     & ,isend,itag+2,MPI_COMM_WORLD,istatus,ierr)
      itag=itag+2

      do 250 istate=1,nstates
        do 250 iab=1,2
        itag=itag+1
 250    call mpi_recv(ymatw(1,1,nwalk,iab,istate),MORB*nelec,mpi_double_precision
     &   ,isend,itag,MPI_COMM_WORLD,istatus,ierr)

      do 260 iab=1,2
        do 260 k=1,ndet
          ndim=numrep_det(k,iab)
          if(k.ne.krefw(nwalk).and.ndim.gt.0) then
            itag=itag+1
            call mpi_recv(wfmatw(1,k,nwalk,iab),ndim*ndim,mpi_double_precision
     &     ,isend,itag,MPI_COMM_WORLD,istatus,ierr)
        endif
 260  continue

      call mpi_recv(orbw(1,1,nwalk),nelec*morb,mpi_double_precision
     &  ,isend,itag+1,MPI_COMM_WORLD,istatus,ierr)
      call mpi_recv(dorbw(1,1,1,nwalk),3*nelec*morb,mpi_double_precision
     &  ,isend,itag+2,MPI_COMM_WORLD,istatus,ierr)
      itag=itag+2

      return
      end
