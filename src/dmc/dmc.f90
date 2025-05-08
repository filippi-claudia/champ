module dmc_f_mod
contains
      subroutine dmc
! Written by Cyrus Umrigar with major contributions by Claudia Filippi.
! Uses the diffusion Monte Carlo algorithm described in:
! A Diffusion Monte Carlo Algorithm with Very Small Time-Step Errors,
! C.J. Umrigar, M.P. Nightingale and K.J. Runge, J. Chem. Phys., 99, 2865 (1993).

      use acues1_mod, only: acues1
      use acues1_reduce_mod, only: acues1_reduce
      use acuest_mod, only: acuest
      use branch, only: eest
      use const, only: etrial
      ! use averages, only: average,average_write,init_averages_index
      use constants, only: pi
      use contrl_file, only: ounit
      use contrldmc, only: idmc
      use control_dmc, only: dmc_idump,dmc_irstar,dmc_nblk,dmc_nblkeq
      use control_dmc, only: dmc_nconf,dmc_nstep
      use dmc_ps_mov1, only: dmc_ps
      use dumper_mod, only: dumper
      use error,   only: fatal_error
      use estcum,  only: ipass
      use finwrt_mod, only: finwrt
      use fragments, only: etrialfrag, eestfrag, nfrag
      use init_mod, only: init
      use mc_configs_mod, only: mc_configs,mc_configs_write
      use mpitimer, only: elapsed_time
      use mpiconf, only: nproc
      use mpi
      use multiple_geo, only: iwftype,nforce,nwftype,nwprod
      use m_force_analytic, only: iforce_analy
      use precision_kinds, only: dp
      use pseudo,  only: nloc
      use rotqua_mod, only: rotqua
      use strech_mod, only: setup_force
      use vd_mod,         only: dmc_ivd
      use zerest_mod, only: zerest


      implicit none

      integer :: i, j, irun, lpass, ifrag, ierr
      real(dp), parameter :: one = 1.d0
      real(dp), parameter :: four = 4.d0
      real(dp) :: etrialcollect
      real(dp), dimension(:), allocatable :: etrialfragcollect

      if(nforce.gt.1) then
        call setup_force
       else
        ! debug line. ravindra
        if (.not. allocated(iwftype)) allocate (iwftype(nforce), source=0)
        if(iforce_analy.eq.1.and.dmc_ivd.eq.0) nwprod=1
        nwftype=1
        iwftype(1)=1
      endif

! read walker configurations
      call mc_configs

! get initial value of cpu time

      call elapsed_time("DMC : Reading initial walker configuration : ")

! initialize sums and averages
      ! call init_averages_index
      if(dmc_irstar.ne.1) call init
      ! if(dmc_irstar.ne.1) call average(0)

! forces implemented only for certain dmc control options
      if(nforce.gt.1) write(ounit,'(''Possible Warning: force implemented for certain dmc control options'')')

!     call flush(6)
! loops for dmc calculation
      lpass=0
      irun=0
      do i=1,dmc_nblk+2*dmc_nblkeq
        if((i.eq.dmc_nblkeq+1.or.i.eq.2*dmc_nblkeq+1).and.dmc_irstar.ne.1) then

          call elapsed_time("DMC : equilibrium CP : ")

          call zerest
          ! call average(0)
          call elapsed_time("DMC : zero out estimators and averages : ")
        endif
        do j=1,dmc_nstep
            if(i.gt.2*dmc_nblkeq.or.dmc_irstar.eq.1) then
                  lpass=lpass+1
                  irun=1
               endif
          ipass=ipass+1
          if (nloc.gt.0) call rotqua
          if(iabs(idmc).eq.1) then
!           call dmc_brock
           elseif(iabs(idmc).eq.2) then
            if (nloc.gt.0) then
              call dmc_ps(lpass,irun)
             else
!             call dmc_good
            endif
           else
            call fatal_error('DMC: iabs(idmc) must be 1 or 2')
          endif
          call mc_configs_write(i,ipass)
          call acues1
        enddo
        ! call average(2)
        ! call average_write
        call acuest
      enddo

      call acues1_reduce

      call finwrt
      call elapsed_time("DMC : all CP : ")

      if (dmc_idump.eq.1) call dumper
      close (unit=9)
      if (dmc_nconf.ne.0) close (unit=7)
      call elapsed_time("dumping restart files : ")

      end
end module
