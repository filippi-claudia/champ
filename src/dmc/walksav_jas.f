      subroutine walksav_jas(iw)
c Written by Claudia Filippi

      use vmc_mod, only: MELEC, MORB, MBASIS, MDET, MCENT, MCTYPE, MCTYP3X,
     &NSPLIN, nrad, MORDJ, MORDJ1, MMAT_DIM, MMAT_DIM2, MMAT_DIM20,
     &radmax, delri, NEQSX, MTERMS, MCENT3, NCOEF, MEXCIT
      use dmc_mod, only: MWALK, MFPROD, MFPRD1, MPATH
      use const, only: delta, deltai, etrial, fbias, hb, imetro, ipr, nelec, pi
      use force_mod, only: MFORCE, MFORCE_WT_PRD, MWF

      implicit real*8(a-h,o-z)

      include 'mpif.h'


      common /branch/ wtgen(0:MFPRD1),ff(0:MFPRD1),eold(MWALK,MFORCE),
     &pwt(MWALK,MFORCE),wthist(MWALK,0:MFORCE_WT_PRD,MFORCE),
     &wt(MWALK),eigv,eest,wdsumo,wgdsumo,fprod,nwalk

      common /jaso/ fso(MELEC,MELEC),fijo(3,MELEC,MELEC)
     &,d2ijo(MELEC,MELEC),d2o,fsumo,fjo(3,MELEC)

      common /velocity_jastrow/vj(3,MELEC),vjn(3,MELEC)

      dimension fsow(MELEC,MELEC,MWALK),fijow(3,MELEC,MELEC,MWALK)
     &,fsumow(MWALK),fjow(3,MELEC,MWALK),d2ow(MWALK),d2ijow(MELEC,MELEC,MWALK)

      dimension vjw(3,MELEC,MWALK)

      dimension istatus(MPI_STATUS_SIZE)

      save fsow,fijow,fsumow,fjow,d2ow,d2ijow

      save vjw

      fsumow(iw)=fsumo

      do 10 i=1,nelec
        fjow(1,i,iw)=fjo(1,i)
        fjow(2,i,iw)=fjo(2,i)
  10    fjow(3,i,iw)=fjo(3,i)

      do 20 i=2,nelec
        do 20 j=1,i-1
        fsow(i,j,iw)=fso(i,j)
        fijow(1,i,j,iw)=fijo(1,i,j)
        fijow(2,i,j,iw)=fijo(2,i,j)
        fijow(3,i,j,iw)=fijo(3,i,j)
        fijow(1,j,i,iw)=fijo(1,j,i)
        fijow(2,j,i,iw)=fijo(2,j,i)
  20    fijow(3,j,i,iw)=fijo(3,j,i)

      do 25 i=1,nelec
        fsow(i,i,iw)=fso(i,i)
        fijow(1,i,i,iw)=fijo(1,i,i)
        fijow(2,i,i,iw)=fijo(2,i,i)
  25    fijow(3,i,i,iw)=fijo(3,i,i)

      do 26 i=1,nelec
        do 26 kk=1,3
  26      vjw(kk,i,iw)=vj(kk,i)
    
      return

      entry walkstrjas(iw)

      fsumo=fsumow(iw)

      do 30 i=1,nelec
        fjo(1,i)=fjow(1,i,iw)
        fjo(2,i)=fjow(2,i,iw)
  30    fjo(3,i)=fjow(3,i,iw)

      do 40 i=2,nelec
        do 40 j=1,i-1
        fso(i,j)=fsow(i,j,iw)
        fijo(1,i,j)=fijow(1,i,j,iw)
        fijo(2,i,j)=fijow(2,i,j,iw)
        fijo(3,i,j)=fijow(3,i,j,iw)
        fijo(1,j,i)=fijow(1,j,i,iw)
        fijo(2,j,i)=fijow(2,j,i,iw)
  40    fijo(3,j,i)=fijow(3,j,i,iw)

      do 45 i=1,nelec
        fso(i,i)=fsow(i,i,iw)
        fijo(1,i,i)=fijow(1,i,i,iw)
        fijo(2,i,i)=fijow(2,i,i,iw)
  45    fijo(3,i,i)=fijow(3,i,i,iw)

      do 46 i=1,nelec
        do 46 kk=1,3
  46      vj(kk,i)=vjw(kk,i,iw)
    
      return

      entry splitjjas(iw,iw2)

      fsumow(iw2)=fsumow(iw)

      do 50 i=1,nelec
        fjow(1,i,iw2)=fjow(1,i,iw)
        fjow(2,i,iw2)=fjow(2,i,iw)
  50    fjow(3,i,iw2)=fjow(3,i,iw)

      do 60 i=2,nelec
        do 60 j=1,i-1
        fsow(i,j,iw2)=fsow(i,j,iw)
        fijow(1,i,j,iw2)=fijow(1,i,j,iw)
        fijow(2,i,j,iw2)=fijow(2,i,j,iw)
        fijow(3,i,j,iw2)=fijow(3,i,j,iw)
        fijow(1,j,i,iw2)=fijow(1,j,i,iw)
        fijow(2,j,i,iw2)=fijow(2,j,i,iw)
  60    fijow(3,j,i,iw2)=fijow(3,j,i,iw)

      do 65 i=1,nelec
        fsow(i,i,iw2)=fsow(i,i,iw)
        fijow(1,i,i,iw2)=fijow(1,i,i,iw)
        fijow(2,i,i,iw2)=fijow(2,i,i,iw)
  65    fijow(3,i,i,iw2)=fijow(3,i,i,iw)

      do 66 i=1,nelec
        do 66 kk=1,3
  66      vjw(kk,i,iw2)=vjw(kk,i,iw)
    
      return

      entry send_jas(irecv)

      itag=0
      call mpi_isend(fsumow(nwalk),1,mpi_double_precision,irecv
     &,itag+1,MPI_COMM_WORLD,irequest,ierr)
      call mpi_isend(fjow(1,1,nwalk),3*nelec,mpi_double_precision,irecv
     &,itag+2,MPI_COMM_WORLD,irequest,ierr)
      call mpi_isend(fsow(1,1,nwalk),MELEC*nelec,mpi_double_precision
     &,irecv,itag+3,MPI_COMM_WORLD,irequest,ierr)
      call mpi_isend(fijow(1,1,1,nwalk),3*MELEC*nelec
     &,mpi_double_precision,irecv,itag+4,MPI_COMM_WORLD,irequest,ierr)

      call mpi_isend(vjw(1,1,nwalk),3*nelec,mpi_double_precision,irecv
     &,itag+5,MPI_COMM_WORLD,irequest,ierr)

      return

      entry recv_jas(isend)

      itag=0
      call mpi_recv(fsumow(nwalk),1,mpi_double_precision,isend
     &,itag+1,MPI_COMM_WORLD,istatus,ierr)
      call mpi_recv(fjow(1,1,nwalk),3*nelec,mpi_double_precision,isend
     &,itag+2,MPI_COMM_WORLD,istatus,ierr)
      call mpi_recv(fsow(1,1,nwalk),MELEC*nelec,mpi_double_precision
     &,isend,itag+3,MPI_COMM_WORLD,istatus,ierr)
      call mpi_recv(fijow(1,1,1,nwalk),3*MELEC*nelec
     &,mpi_double_precision,isend,itag+4,MPI_COMM_WORLD,istatus,ierr)

      call mpi_recv(vjw(1,1,nwalk),3*nelec,mpi_double_precision,isend
     &,itag+5,MPI_COMM_WORLD,istatus,ierr)

      return
      end
