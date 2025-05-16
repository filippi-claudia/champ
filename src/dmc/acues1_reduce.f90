module acues1_reduce_mod
contains
      subroutine acues1_reduce

      use control, only: mode
      use contrl_per, only: iperiodic
      use contrl_file, only: ounit
      use control, only: mode
      use derivest, only: derivcm2, derivcum, derivtotave
      use estcum, only: iblk
      use stats, only: acc, nacc, nodecr, trymove
      use estcum, only: ecum1_dmc, efcum1, egcum, egcum1
      use estcum, only: wcum1, wfcum1, wgcum, wgcum1
      use est2cm, only: ecm21_dmc, efcm21, egcm21
      use est2cm, only: wcm21
      use est2cm, only: wfcm21, wgcm21, wgcm2
      use force_analytic, only: force_analy_fin
      use force_analy_reduce_mod, only: force_analy_reduce
      use force_pth, only: PTH
      use m_force_analytic, only: iforce_analy
      use fragments, only: egcm21frag, egcum1frag, nfrag
      use mpiconf, only: nproc, wid
      use mpi
      use multiple_geo, only: MFORCE,nforce
      use optci_mod, only: optci_fin
      use optci_reduce_mod, only: optci_reduce
      use optjas_mod, only: optjas_fin
      use optjas_reduce_mod, only: optjas_reduce
      use optorb_f_mod, only: optorb_fin
      use optorb_reduce_mod, only: optorb_reduce
      use optx_jas_ci, only: optx_jas_ci_fin
      use optx_jas_ci_reduce_mod, only: optx_jas_ci_reduce
      use optx_jas_orb, only: optx_jas_orb_fin
      use optx_jas_orb_reduce_mod, only: optx_jas_orb_reduce
      use optx_orb_ci, only: optx_orb_ci_fin
      use optx_orb_ci_reduce_mod, only: optx_orb_ci_reduce
      use precision_kinds, only: dp
      use pathak_mod, only: ipathak      
      use stats,   only: acc,nacc,nodecr,trymove
      use system,    only: ncent


      implicit none

      integer :: ierr, ifr, i, nodecr_collect, nacc_collect, id, irequest, ic, k, iph
      integer, dimension(MPI_STATUS_SIZE) :: istatus
      real(dp) :: e1collect, e21collect, w1collect, w21collect, ef1collect
      real(dp) :: ef21collect, wf1collect, wf21collect, trymove_collect
      real(dp) :: acc_collect, efin
      real(dp) :: rn_eff


      real(dp), dimension(MFORCE) :: eg1collect
      real(dp), dimension(MFORCE) :: eg21collect
      real(dp), dimension(nfrag) :: eg1collectfrag
      real(dp), dimension(nfrag) :: eg21collectfrag
      real(dp), dimension(MFORCE) :: wg1collect
      real(dp), dimension(MFORCE) :: wg21collect
      real(dp), dimension(3,ncent,PTH) :: derivgerr

      if(mode.eq.'dmc_one_mpi2') return

      call mpi_reduce(ecum1_dmc,e1collect,1,mpi_double_precision &
      ,mpi_sum,0,MPI_COMM_WORLD,ierr)
      call mpi_reduce(ecm21_dmc,e21collect,1,mpi_double_precision &
      ,mpi_sum,0,MPI_COMM_WORLD,ierr)
      call mpi_reduce(wcum1,w1collect,1,mpi_double_precision &
      ,mpi_sum,0,MPI_COMM_WORLD,ierr)
      call mpi_reduce(wcm21,w21collect,1,mpi_double_precision &
      ,mpi_sum,0,MPI_COMM_WORLD,ierr)

      ecum1_dmc=e1collect
      wcum1=w1collect
      ecm21_dmc=e21collect
      wcm21=w21collect

      call mpi_reduce(efcum1,ef1collect,1,mpi_double_precision &
      ,mpi_sum,0,MPI_COMM_WORLD,ierr)
      call mpi_reduce(efcm21,ef21collect,1,mpi_double_precision &
      ,mpi_sum,0,MPI_COMM_WORLD,ierr)
      call mpi_reduce(wfcum1,wf1collect,1,mpi_double_precision &
      ,mpi_sum,0,MPI_COMM_WORLD,ierr)
      call mpi_reduce(wfcm21,wf21collect,1,mpi_double_precision &
      ,mpi_sum,0,MPI_COMM_WORLD,ierr)

      efcum1=ef1collect
      wfcum1=wf1collect
      efcm21=ef21collect
      wfcm21=wf21collect

      call mpi_reduce(egcum1,eg1collect,MFORCE,mpi_double_precision &
      ,mpi_sum,0,MPI_COMM_WORLD,ierr)
      call mpi_reduce(egcm21,eg21collect,MFORCE,mpi_double_precision &
      ,mpi_sum,0,MPI_COMM_WORLD,ierr)
      call mpi_reduce(wgcum1,wg1collect,MFORCE,mpi_double_precision &
      ,mpi_sum,0,MPI_COMM_WORLD,ierr)
      call mpi_reduce(wgcm21,wg21collect,MFORCE,mpi_double_precision &
      ,mpi_sum,0,MPI_COMM_WORLD,ierr)
      
      if (nfrag.gt.1) then
        call mpi_reduce(egcum1frag,eg1collectfrag,nfrag,mpi_double_precision &
        ,mpi_sum,0,MPI_COMM_WORLD,ierr)
        call mpi_reduce(egcm21frag,eg21collectfrag,nfrag,mpi_double_precision &
        ,mpi_sum,0,MPI_COMM_WORLD,ierr)

        egcum1frag(:)=eg1collectfrag(:)
        egcm21frag(:)=eg21collectfrag(:)
      endif

      do ifr=1,nforce
        egcum1(ifr)=eg1collect(ifr)
        wgcum1(ifr)=wg1collect(ifr)
        egcm21(ifr)=eg21collect(ifr)
        wgcm21(ifr)=wg21collect(ifr)
      enddo

      call mpi_reduce(nodecr,nodecr_collect,1,mpi_integer,mpi_sum,0,MPI_COMM_WORLD,ierr)
      call mpi_reduce(trymove,trymove_collect,1,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)
      call mpi_reduce(acc,acc_collect,1,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)
      call mpi_reduce(nacc,nacc_collect,1,mpi_integer,mpi_sum,0,MPI_COMM_WORLD,ierr)
      nodecr=nodecr_collect
      trymove=trymove_collect
      acc=acc_collect
      nacc=nacc_collect

      call optjas_reduce
      call optorb_reduce
      call optci_reduce
      call optx_jas_orb_reduce
      call optx_jas_ci_reduce
      call optx_orb_ci_reduce
      call force_analy_reduce

      if(wid) then
        do id=1,nproc-1
          call mpi_isend(egcum1,1,mpi_double_precision,id &
          ,1,MPI_COMM_WORLD,irequest,ierr)
          call mpi_isend(wgcum1,1,mpi_double_precision,id &
          ,2,MPI_COMM_WORLD,irequest,ierr)
        enddo
       else
        call mpi_recv(egcum1,1,mpi_double_precision,0 &
        ,1,MPI_COMM_WORLD,istatus,ierr)
        call mpi_recv(wgcum1,1,mpi_double_precision,0 &
        ,2,MPI_COMM_WORLD,istatus,ierr)
      endif

      efin=egcum(1)/wgcum(1)
      rn_eff=(wgcum(1)**2/wgcm2(1))-1.d0

      call optjas_fin(wgcum1(1),egcum1(1))
      call optci_fin(iblk,wgcum1(1),efin)
      call optorb_fin(wgcum1(1),egcum1(1))
      call optx_jas_ci_fin(wgcum1(1),efin)
      call optx_jas_orb_fin(wgcum1(1),egcum1(1))
      call optx_orb_ci_fin(wgcum1(1),efin)
      if(wid) call force_analy_fin(wgcum(1),rn_eff,efin)
      
      if(wid) then
         if (iforce_analy.gt.0) then
            do iph=1,PTH
               do ic=1,ncent
                  do k=1,3
                     derivtotave(k,ic,iph)=(derivcum(1,k,ic,iph)+2.d0*derivcum(2,k,ic,iph)-2.d0*efin*derivcum(3,k,ic,iph))/wgcum(1)
                     derivgerr(k,ic,iph)=errg(derivtotave(k,ic,iph)*wgcum(1),derivcm2(k,ic,iph),wgcum(1),rn_eff)
                  enddo
               enddo
            enddo
            do iph=1,PTH
               do ic=1,ncent
                  if (ipathak.gt.0) then
                     write(ounit,'(i5,i5,1p6e14.5)')iph,ic,(derivtotave(k,ic,iph),k=1,3),(derivgerr(k,ic,iph),k=1,3)
                  else
                     write(ounit,'(i5,1p6e14.5)') ic,(derivtotave(k,ic,iph),k=1,3),(derivgerr(k,ic,iph),k=1,3)
                  endif
               enddo
            enddo
         endif
      endif

      return
contains
        elemental pure function error(x,x2,w,rn)
          implicit none
          real(dp), intent(in) :: x, x2,w,rn
          real(dp)             :: error
          error=dsqrt(max((x2/w-(x/w)**2)/rn,0.d0))
        end function
        elemental pure function errg(x,x2,w,rn)
          implicit none
          real(dp), intent(in) :: x, x2, w, rn
          real(dp)             :: errg
          errg=error(x,x2,w,rn)  
        end function
      end
end module
