      subroutine send_walker(irecv)
c Written by Claudia Filippi

      use const, only: nelec
      use forcepar, only: nforce
      use age, only: iage
      use config, only: d2o, peo_dmc, psido_dmc, psijo_dmc, vold_dmc, xold_dmc
      use force_dmc, only: nwprod
      use branch, only: eold, nwalk, pwt, wt
      use branch, only: wthist
      use jacobsave, only: ajacold
      use velratio, only: fratio
      use mpi

      implicit none

      integer :: ierr, ifr, ip, irecv, irequest
      integer :: isend, itag
      integer, dimension(MPI_STATUS_SIZE) :: istatus




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
        call mpi_isend(psido_dmc(nwalk,ifr),1,mpi_double_precision,irecv
     &  ,itag+3,MPI_COMM_WORLD,irequest,ierr)
        call mpi_isend(psijo_dmc(nwalk,ifr),1,mpi_double_precision,irecv
     &  ,itag+4,MPI_COMM_WORLD,irequest,ierr)
        call mpi_isend(peo_dmc(nwalk,ifr),1,mpi_double_precision,irecv
     &  ,itag+5,MPI_COMM_WORLD,irequest,ierr)
        call mpi_isend(d2o(nwalk,ifr),1,mpi_double_precision,irecv
     &  ,itag+6,MPI_COMM_WORLD,irequest,ierr)
        call mpi_isend(pwt(nwalk,ifr),1,mpi_double_precision,irecv
     &  ,itag+7,MPI_COMM_WORLD,irequest,ierr)
        call mpi_isend(fratio(nwalk,ifr),1,mpi_double_precision,irecv
     &  ,itag+8,MPI_COMM_WORLD,irequest,ierr)
        call mpi_isend(vold_dmc(1,1,nwalk,ifr),3*nelec,mpi_double_precision
     &  ,irecv,itag+9,MPI_COMM_WORLD,irequest,ierr)
        call mpi_isend(xold_dmc(1,1,nwalk,ifr),3*nelec,mpi_double_precision
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
        call mpi_recv(psido_dmc(nwalk,ifr),1,mpi_double_precision,isend
     &  ,itag+3,MPI_COMM_WORLD,istatus,ierr)
        call mpi_recv(psijo_dmc(nwalk,ifr),1,mpi_double_precision,isend
     &  ,itag+4,MPI_COMM_WORLD,istatus,ierr)
        call mpi_recv(peo_dmc(nwalk,ifr),1,mpi_double_precision,isend
     &  ,itag+5,MPI_COMM_WORLD,istatus,ierr)
        call mpi_recv(d2o(nwalk,ifr),1,mpi_double_precision,isend
     &  ,itag+6,MPI_COMM_WORLD,istatus,ierr)
        call mpi_recv(pwt(nwalk,ifr),1,mpi_double_precision,isend
     &  ,itag+7,MPI_COMM_WORLD,istatus,ierr)
        call mpi_recv(fratio(nwalk,ifr),1,mpi_double_precision,isend
     &  ,itag+8,MPI_COMM_WORLD,istatus,ierr)
        call mpi_recv(vold_dmc(1,1,nwalk,ifr),3*nelec,mpi_double_precision
     &  ,isend,itag+9,MPI_COMM_WORLD,istatus,ierr)
        call mpi_recv(xold_dmc(1,1,nwalk,ifr),3*nelec,mpi_double_precision
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
