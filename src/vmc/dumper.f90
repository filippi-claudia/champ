module dumper_mod
! MPI version created by Claudia Filippi starting from serial version
! routine to pick up and dump everything needed to restart
! job where it left off

      use config, only: xold
      use contrl_file,    only: ounit
      use csfs, only: nstates
      use est2cm, only: ecm2, ecm21, pecm2, tjfcm2, tpbcm2
      use estcum, only: ecum, ecum1, pecum, tjfcum, tpbcum
      use estsig, only: ecm21s, ecum1s
      use estsum, only: acc
      use forcewt, only: wcum
      use multiple_geo, only: nforce, fcm2, fcum
      use multiple_states, only: efficiency_init
      use mpi
      use mpiconf, only: idtask, nproc, wid
      use precision_kinds, only: dp
      use pseudo, only: nloc
      use qua, only: nquad, wq, xq, yq, zq
      use system, only: nelec

      use random_mod,      only: savern
      use dumper_more_mod, only: dumper_more, startr_more
      use error, only: fatal_error
      use random_mod, only: setrn, random_dp
      use optci_mod, only: optci_init
      use optorb_f_mod, only: optorb_init
      use optjas_mod, only: optjas_init
      use properties_mod, only: prop_init
      use pcm_mod, only: pcm_init
      use mmpol, only: mmpol_init
      use force_analytic, only: force_analy_init
      use forcewt, only: wcum
      use mmpol,   only: mmpol_init
      use mpiconf, only: idtask,nproc,wid
      use multiple_geo, only: fcm2,fcum,nforce
      use multiple_states, only: efficiency_init
      use optci_mod, only: optci_init
      use optjas_mod, only: optjas_init
      use optorb_f_mod, only: optorb_init
      use pcm_mod, only: pcm_init
      use precision_kinds, only: dp
      use properties_mod, only: prop_init
      use random_mod, only: random_dp,savern,setrn

      implicit none

contains
      subroutine dumper
      implicit none
      integer :: i, id, idfrom, idget, ierr
      integer :: ifr, istate, j, k
      integer :: nelecx, nforcex, nlocx, nproco
      integer :: nq_id, nqd_id, nqx
      integer, dimension(8,0:nproc) :: irn
      integer, dimension(MPI_STATUS_SIZE) :: istatus
      real(dp) :: rnd, wq_id, x_id, xq_id, yq_id
      real(dp) :: zq_id

      rewind 10


      call savern(irn(1,idtask))


      if(idtask.ne.0) then
        call mpi_send(irn(:,idtask), 8, mpi_integer, 0, 1, MPI_COMM_WORLD, ierr)
        call mpi_send(xold,3*nelec,mpi_double_precision,0,1,MPI_COMM_WORLD,ierr)
        call mpi_send(xq,nquad,mpi_double_precision,0,2,MPI_COMM_WORLD,ierr)
        call mpi_send(yq,nquad,mpi_double_precision,0,3,MPI_COMM_WORLD,ierr)
        call mpi_send(zq,nquad,mpi_double_precision,0,4,MPI_COMM_WORLD,ierr)
       else
        do id=1, nproc-1
          call mpi_recv(irn(:, id), 8, mpi_integer, id,1,MPI_COMM_WORLD,istatus,ierr)
        enddo
        write(10) nproc
        write(10) ((irn(i,j),i=1,8),j=0,nproc-1)
        write(10) nelec,nforce,nloc
        write(10) ((xold(k,i),k=1,3),i=1,nelec)
        if(nloc.gt.0) write(10) nquad,(xq(i),yq(i),zq(i),wq(i),i=1,nquad)
        do id=1,nproc-1
          call mpi_recv(xold,3*nelec,mpi_double_precision,id &
          ,1,MPI_COMM_WORLD,istatus,ierr)
          call mpi_recv(xq,nquad,mpi_double_precision,id &
          ,2,MPI_COMM_WORLD,istatus,ierr)
          call mpi_recv(yq,nquad,mpi_double_precision,id &
          ,3,MPI_COMM_WORLD,istatus,ierr)
          call mpi_recv(zq,nquad,mpi_double_precision,id &
          ,4,MPI_COMM_WORLD,istatus,ierr)
          write(10) ((xold(k,i),k=1,3),i=1,nelec)
          if(nloc.gt.0) write(10) nquad,(xq(i),yq(i),zq(i),wq(i),i=1,nquad)
        enddo
      endif



      if(.not.wid) return

      call dumper_more

      end subroutine

!-----------------------------------------------------------------------
      subroutine startr
      implicit none
      integer :: i, id, idfrom, idget, ierr
      integer :: ifr, istate, j, k
      integer :: nelecx, nforcex, nlocx, nproco
      integer :: nq_id, nqd_id, nqx, nscounts
      integer, dimension(8,0:nproc) :: irn
      integer, dimension(MPI_STATUS_SIZE) :: istatus
      real(dp) :: rnd, wq_id, x_id, xq_id, yq_id
      real(dp) :: zq_id

      write(ounit,'(1x,''attempting restart from unit 10'')')

      rewind 10
      read(10) nproco
      if(nproco.ne.nproc) write(ounit,'(''Warning: different number of processors'',/ &
       ,9x,''old number of processors'',i3,/,9x,''continuing with'',i3,'' processors'')') &
       nproco,nproc
      read(10) ((irn(i,j),i=1,8),j=0,nproco-1)
      if(idtask.le.nproco-1) call setrn(irn(1,idtask))
      read(10) nelecx,nforcex,nlocx
      if (nelecx.ne.nelec) call fatal_error('STARTR: nelec')
      if (nforcex.ne.nforce) call fatal_error('STARTR: nforce')
      if (nlocx.ne.nloc) call fatal_error('STARTR: nloc')
      if(idtask.le.nproco-1) then
        do id=0,idtask
          read(10) ((xold(k,i),k=1,3),i=1,nelec)
          if(nloc.gt.0) read(10) nqx,(xq(i),yq(i),zq(i),wq(i),i=1,nqx)
        enddo
        if(nqx.ne.nquad) call fatal_error('STARTR: nquad')
        do id=idtask+1,nproco-1
          read(10) (x_id,i=1,3*nelec)
          if(nloc.gt.0) read(10) nq_id,(xq_id,yq_id,zq_id,wq_id,i=1,nqd_id)
        enddo
       else
        do id=0,nproco-1
          read(10) (x_id,i=1,3*nelec)
          if(nloc.gt.0) read(10) nq_id,(xq_id,yq_id,zq_id,wq_id,i=1,nqd_id)
        enddo
      endif

      if(nproc.gt.nproco) then
        if(idtask.le.nproco-1) then
          do idget=nproco,nproc-1
! xold from idtask to idget
            if(idtask.eq.mod(idget,nproco)) &
            call mpi_send(xold,3*nelec,mpi_double_precision,idget &
            ,idget,MPI_COMM_WORLD,ierr)
!    &      ,idget,MPI_COMM_WORLD,irequest,ierr)
          enddo
         else
          do id=1,(3*nelec)*idtask
            rnd=random_dp()
          enddo
          idfrom=mod(idtask,nproco)
! xold from idfrom to idtask
          call mpi_recv(xold,3*nelec,mpi_double_precision,idfrom &
          ,idtask,MPI_COMM_WORLD,istatus,ierr)
        endif
      endif

      call startr_more

      if (wid) return

      acc=0

      do istate=1,nstates

      pecum(istate)=0
      tpbcum(istate)=0

      pecm2(istate)=0
      tpbcm2(istate)=0

      ecum1(istate)=0
      ecum1s(istate)=0
      ecm21(istate)=0
      ecm21s(istate)=0

      do ifr=1,nforce
        ecum(istate,ifr)=0
        ecm2(istate,ifr)=0
        wcum(istate,ifr)=0
        fcum(istate,ifr)=0
        fcm2(istate,ifr)=0
      enddo
      enddo

      call optjas_init
      call optorb_init(0)
      call optci_init(0)
      call prop_init(0)
      call pcm_init(0)
      call mmpol_init(0)
      call force_analy_init(0)
      call efficiency_init

      end
end module
