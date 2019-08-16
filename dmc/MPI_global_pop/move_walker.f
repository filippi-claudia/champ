      subroutine send_walker(irecv)
c Written by Claudia Filippi
      implicit real*8(a-h,o-z)
      include 'vmc.h'
      include 'dmc.h'
      include 'force.h'
      include 'mpif.h'

      parameter (zero=0.d0,one=1.d0,two=2.d0,three=3.d0,four=4.d0)
      parameter (ten=10.d0,tenth=.1d0,half=.5d0,third=1.d0/3.d0)

      common /forcepar/ deltot(MFORCE),nforce,istrech
      common /forcest/ fgcum(MFORCE),fgcm2(MFORCE)
      common /force_dmc/ itausec,nwprod

      common /const/ pi,hb,etrial,delta,deltai,fbias,nelec,imetro,ipr
      common /config/ xold(3,MELEC,MWALK,MFORCE),vold(3,MELEC,MWALK,MFORCE),
     &psido(MWALK,MFORCE),psijo(MWALK,MFORCE),peo(MWALK,MFORCE),d2o(MWALK,MFORCE)
      common /velratio/ fratio(MWALK,MFORCE)
      common /age/ iage(MWALK),ioldest,ioldestmx
      common /stats/ dfus2ac,dfus2un,dr2ac,dr2un,acc,trymove,nacc,
     &nbrnch,nodecr
      common /branch/ wtgen(0:MFPRD1),ff(0:MFPRD1),eold(MWALK,MFORCE),
     &pwt(MWALK,MFORCE),wthist(MWALK,0:MFORCE_WT_PRD,MFORCE),
     &wt(MWALK),eigv,eest,wdsumo,wgdsumo,fprod,nwalk
      common /jacobsave/ ajacob,ajacold(MWALK,MFORCE)

      common /mpiconf/ idtask,nproc,wid

      dimension istatus(MPI_STATUS_SIZE)

      call mpi_isend(wt(nwalk),1,mpi_double_precision,irecv,1
     &,MPI_COMM_WORLD,irequest,ierr)
      call mpi_isend(iage(nwalk),1,mpi_integer,irecv,2
     &,MPI_COMM_WORLD,irequest,ierr)

      itag=2
      do 15 ifr=1,nforce
        call mpi_isend(ajacold(nwalk,ifr),1,mpi_double_precision,irecv
     &  ,itag+1,MPI_COMM_WORLD,irequest,ierr)
        call mpi_isend(eold(nwalk,ifr),1,mpi_double_precision,irecv
     &  ,itag+2,MPI_COMM_WORLD,irequest,ierr)
        call mpi_isend(psido(nwalk,ifr),1,mpi_double_precision,irecv
     &  ,itag+3,MPI_COMM_WORLD,irequest,ierr)
        call mpi_isend(psijo(nwalk,ifr),1,mpi_double_precision,irecv
     &  ,itag+4,MPI_COMM_WORLD,irequest,ierr)
        call mpi_isend(peo(nwalk,ifr),1,mpi_double_precision,irecv
     &  ,itag+5,MPI_COMM_WORLD,irequest,ierr)
        call mpi_isend(d2o(nwalk,ifr),1,mpi_double_precision,irecv
     &  ,itag+6,MPI_COMM_WORLD,irequest,ierr)
        call mpi_isend(pwt(nwalk,ifr),1,mpi_double_precision,irecv
     &  ,itag+7,MPI_COMM_WORLD,irequest,ierr)
        call mpi_isend(fratio(nwalk,ifr),1,mpi_double_precision,irecv
     &  ,itag+8,MPI_COMM_WORLD,irequest,ierr)
        call mpi_isend(vold(1,1,nwalk,ifr),3*nelec,mpi_double_precision
     &  ,irecv,itag+9,MPI_COMM_WORLD,irequest,ierr)
        call mpi_isend(xold(1,1,nwalk,ifr),3*nelec,mpi_double_precision
     &  ,irecv,itag+10,MPI_COMM_WORLD,irequest,ierr)
        itag=itag+10
        do 15 ip=0,nwprod-1
        itag=itag+1
        call mpi_isend(wthist(nwalk,ip,ifr),1,mpi_double_precision,irecv
     &  ,itag,MPI_COMM_WORLD,irequest,ierr)
   15 continue

c     call send_det(itag,irecv)
c     call send_jas(itag,irecv)

c     nwalk=nwalk-1

      call prop_send(irecv,itag)
      call pcm_send(irecv,itag)
      call mmpol_send(irecv,itag)

      return

      entry recv_walker(isend)

c     nwalk=nwalk+1

      call mpi_recv(wt(nwalk),1,mpi_double_precision,isend,1
     &,MPI_COMM_WORLD,istatus,ierr)
      call mpi_recv(iage(nwalk),1,mpi_integer,isend,2
     &,MPI_COMM_WORLD,istatus,ierr)

      itag=2
      do 25 ifr=1,nforce
        call mpi_recv(ajacold(nwalk,ifr),1,mpi_double_precision,isend
     &  ,itag+1,MPI_COMM_WORLD,istatus,ierr)
        call mpi_recv(eold(nwalk,ifr),1,mpi_double_precision,isend
     &  ,itag+2,MPI_COMM_WORLD,istatus,ierr)
        call mpi_recv(psido(nwalk,ifr),1,mpi_double_precision,isend
     &  ,itag+3,MPI_COMM_WORLD,istatus,ierr)
        call mpi_recv(psijo(nwalk,ifr),1,mpi_double_precision,isend
     &  ,itag+4,MPI_COMM_WORLD,istatus,ierr)
        call mpi_recv(peo(nwalk,ifr),1,mpi_double_precision,isend
     &  ,itag+5,MPI_COMM_WORLD,istatus,ierr)
        call mpi_recv(d2o(nwalk,ifr),1,mpi_double_precision,isend
     &  ,itag+6,MPI_COMM_WORLD,istatus,ierr)
        call mpi_recv(pwt(nwalk,ifr),1,mpi_double_precision,isend
     &  ,itag+7,MPI_COMM_WORLD,istatus,ierr)
        call mpi_recv(fratio(nwalk,ifr),1,mpi_double_precision,isend
     &  ,itag+8,MPI_COMM_WORLD,istatus,ierr)
        call mpi_recv(vold(1,1,nwalk,ifr),3*nelec,mpi_double_precision
     &  ,isend,itag+9,MPI_COMM_WORLD,istatus,ierr)
        call mpi_recv(xold(1,1,nwalk,ifr),3*nelec,mpi_double_precision
     &  ,isend,itag+10,MPI_COMM_WORLD,istatus,ierr)
        itag=itag+10
        do 25 ip=0,nwprod-1
        itag=itag+1
        call mpi_recv(wthist(nwalk,ip,ifr),1,mpi_double_precision,isend
     &  ,itag,MPI_COMM_WORLD,istatus,ierr)
   25 continue

c     call recv_det(itag,isend)
c     call recv_jas(itag,isend)

      call prop_recv(isend,itag)
      call pcm_recv(isend,itag)
      call mmpol_recv(isend,itag)

      return
      end
