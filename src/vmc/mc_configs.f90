module mc_configs

      use mpi
      use config, only: xnew, xold
      use control_vmc, only: vmc_irstar, vmc_isite, vmc_nconf_new, vmc_icharged_atom
      use contrl_file, only: ounit
      use error, only: fatal_error
      use fin_reduce_mod, only: fin_reduce
      use mpiconf, only: idtask,nproc
      use pcm_mod, only: pcm_qvol
      use sites_mod, only: sites
      use system,  only: iwctype,ncent,ncent_tot,nelec,znuc

      implicit none

contains
      subroutine mc_configs_start
      implicit none

      integer :: i, icharge_system, id
      integer :: ierr
      integer :: k, l, ntotal_sites
      integer, dimension(ncent_tot) :: nsite
      character(len=64) filename
      character(len=32) :: cnum
      logical :: configs_file_open
      logical :: use_sites

! The parser already sets a different random seed on each processor.
      if(vmc_irstar.ne.1) then

        use_sites = vmc_isite.eq.1
        if (.not.use_sites) then
          open(unit=9,iostat=ierr,file='mc_configs_start')
          use_sites = ierr.ne.0
        endif

        if (.not.use_sites) then
          rewind 9
          do id=0,idtask
            read(9,*,iostat=ierr) ((xold(k,i),k=1,3),i=1,nelec)
            if (ierr.ne.0) then
              use_sites = .true.
              close(9)
              exit
            endif
          enddo
        endif

        if (.not.use_sites) then
          write(ounit,'(/,''initial configuration from unit 9'')')
        else
! Generate the initial configuration from atomic sites.
          ntotal_sites=0
          do i=1,ncent
            ntotal_sites=ntotal_sites+int(znuc(iwctype(i))+0.5d0)
          enddo
          icharge_system=ntotal_sites-nelec

          l=0
          do i=1,ncent
            nsite(i)=int(znuc(iwctype(i))+0.5d0)
            if (vmc_icharged_atom.eq.i) then
              nsite(i)=int(znuc(iwctype(i))+0.5d0)-icharge_system
              if (nsite(i).lt.0) call fatal_error('MC_CONFIG: error in icharged_atom')
            endif
            l=l+nsite(i)
            if (l.gt.nelec) then
              nsite(i)=nsite(i)-(l-nelec)
              l=nelec
            endif
          enddo
          if (l.lt.nelec) nsite(1)=nsite(1)+(nelec-l)

          call sites(xold,nelec,nsite)
          open(unit=9,file='mc_configs_start')
          rewind 9

          write(ounit,'(/,''initial configuration from sites'')')
        endif

! If we are moving one electron at a time, then we need to initialize
! xnew, since only the first electron gets initialized in metrop
        do i=1,nelec
          do k=1,3
            xnew(k,i)=xold(k,i)
          enddo
        enddo
      endif

      inquire(unit=9, opened=configs_file_open)
      if (.not.configs_file_open) open(unit=9,file='mc_configs_start')

! If nconf_new > 0 then we want to dump configurations for a future
! optimization or dmc calculation.
      if (vmc_nconf_new.gt.0) then
        write(cnum,*) idtask
        filename='mc_configs_new'//trim(adjustl(cnum))
        open(unit=7,form='formatted',file=filename)
        rewind 7
      endif
      call pcm_qvol(nproc)
      end subroutine

!-----------------------------------------------------------------------
      subroutine mc_configs_write
      implicit none
      integer :: i, ic, id, ierr
      integer, dimension(MPI_STATUS_SIZE) :: istatus

      if(idtask.ne.0) then
        call mpi_send(xold,3*nelec,mpi_double_precision,0 &
        ,1,MPI_COMM_WORLD,ierr)
!    &  ,1,MPI_COMM_WORLD,irequest,ierr)
       else
        rewind 9
        write(9,*) ((xold(ic,i),ic=1,3),i=1,nelec)
        do id=1,nproc-1
          call mpi_recv(xnew,3*nelec,mpi_double_precision,id &
          ,1,MPI_COMM_WORLD,istatus,ierr)
          write(9,*) ((xnew(ic,i),ic=1,3),i=1,nelec)
        enddo
      endif
      close(9)

! reduce cum1 estimates, density and related quantities
      call fin_reduce

      return
      end
end module
