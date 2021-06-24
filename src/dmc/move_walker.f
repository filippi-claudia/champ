      subroutine send_walker(irecv)
c Written by Claudia Filippi
      implicit real*8(a-h,o-z)
      include 'vmc.h'
      include 'dmc.h'
      include 'force.h'
      include 'mpif.h'

      common /const/ pi,hb,etrial,delta,deltai,fbias,nelec,imetro,ipr
      common /config/ xold(3,MELEC,MWALK,MFORCE),vold(3,MELEC,MWALK,MFORCE),
     &psido(MWALK,MFORCE),psijo(MWALK,MFORCE),peo(MWALK,MFORCE),d2o(MWALK,MFORCE)
      common /velratio/ fratio(MWALK,MFORCE)
      common /age/ iage(MWALK),ioldest,ioldestmx
      common /branch/ wtgen(0:MFPRD1),ff(0:MFPRD1),eold(MWALK,MFORCE),
     &pwt(MWALK,MFORCE),wthist(MWALK,0:MFORCE_WT_PRD,MFORCE),
     &wt(MWALK),eigv,eest,wdsumo,wgdsumo,fprod,nwalk
      common /jacobsave/ ajacob,ajacold(MWALK,MFORCE)
      common /forcepar/ deltot(MFORCE),nforce,istrech
      common /forcest/ fgcum(MFORCE),fgcm2(MFORCE)
      common /force_dmc/ itausec,nwprod
      common /derivanaly/ deriv_energy_sum(10,3,MCENT,PTH),deriv_energy_cum(10,3,MCENT,PTH),
     &energy_snake(3,MCENT,MWALK,PTH),energy_hist(3,MCENT,MWALK,0:MFORCE_WT_PRD,PTH),
     &deriv_energy_old(3,MCENT,MWALK),pathak_old(MWALK,PTH),eps_pathak(PTH),ipathak
      common /atom/ znuc(MCTYPE),cent(3,MCENT),pecent
     &,iwctype(MCENT),nctype,ncent
      common /force_analy/ iforce_analy

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

      if(iforce_analy.eq.1) then
        itag=itag+1
        call mpi_isend(pathak_old(nwalk,1),1,mpi_double_precision,irecv
     &  ,itag,MPI_COMM_WORLD,irequest,ierr)
        itag=itag+1
        call mpi_isend(deriv_energy_old(1,1,nwalk),3*ncent,mpi_double_precision,irecv
     &  ,itag,MPI_COMM_WORLD,irequest,ierr)
        itag=itag+1
        call mpi_isend(energy_snake(1,1,nwalk,1),3*ncent*int(PHT),mpi_double_precision,irecv
     &  ,itag,MPI_COMM_WORLD,irequest,ierr)
        do 20 ip=0,nwprod-1
        itag=itag+1
   20   call mpi_isend(energy_hist(1,1,nwalk,ip,1),3*ncent*int(PHT),mpi_double_precision,irecv
     &  ,itag,MPI_COMM_WORLD,irequest,ierr)
      endif

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

      if(iforce_analy.eq.1) then
        itag=itag+1
        call mpi_recv(pathak_old(nwalk,1),1,mpi_double_precision,isend
     &  ,itag,MPI_COMM_WORLD,istatus,ierr)
        itag=itag+1
        call mpi_recv(deriv_energy_old(1,1,nwalk),3*ncent,mpi_double_precision,isend
     &  ,itag,MPI_COMM_WORLD,istatus,ierr)
        itag=itag+1
        call mpi_recv(energy_snake(1,1,nwalk,1),3*ncent*int(PHT),mpi_double_precision,isend
     &  ,itag,MPI_COMM_WORLD,istatus,ierr)
        do 30 ip=0,nwprod-1
        itag=itag+1
   30   call mpi_recv(energy_hist(1,1,nwalk,ip,1),3*ncent*int(PHT),mpi_double_precision,isend
     &  ,itag,MPI_COMM_WORLD,istatus,ierr)
      endif

      call prop_recv(isend,itag)
      call pcm_recv(isend,itag)
      call mmpol_recv(isend,itag)

      return
      end
