module acuest_mod
contains
      subroutine acuest
! MPI version created by Claudia Filippi starting from serial version
! routine to accumulate estimators for energy etc.

      use acuest_gpop_mod, only: acuest_gpop
      use age, only: ioldest
      use contrldmc, only: idmc
      use contrl_file, only: ounit
      use contrldmc, only: icut_e
      use control, only: mode
      use control_dmc, only: dmc_nstep
      use derivest, only: derivcm2, derivcum, derivsum, derivtotave
      use estcum, only: iblk
      use estsum, only: efsum, egsum, esum_dmc
      use estsum, only: pesum_dmc, tausum, tpbsum_dmc
      use estsum, only: wfsum, wgsum, wsum_dmc
      use estcum, only: ecum_dmc, efcum, egcum
      use estcum, only: pecum_dmc, taucum, tpbcum_dmc
      use estcum, only: wcum_dmc, wfcum, wgcum
      use est2cm, only: ecm2_dmc, efcm2, egcm2
      use est2cm, only: pecm2_dmc, tpbcm2_dmc, wcm2
      use est2cm, only: wfcm2, wgcm2
      use force_analytic,   only: force_analy_init, force_analy_cum
      use force_analy_reduce_mod, only: force_analy_reduce
      use force_pth, only: PTH
      use fragments, only: esum_i, ecum_i
      use fragments, only: esumfrag, ecumfrag, nfrag
      use fragments, only: egsumfrag, egcumfrag, egcm2frag
      use m_force_analytic, only: iforce_analy
      use mmpol,   only: mmpol_init
      use mmpol_dmc, only: mmpol_prt
      use mmpol_reduce_mod, only: mmpol_reduce
      use mpiconf, only: nproc, wid
      use mpiblk, only: iblk_proc
      use mpi
      use multiple_geo, only: fgcm2, fgcum, nforce, MFORCE
      use optci_mod, only: optci_cum,optci_init
      use optjas_mod, only: optjas_cum
      use optorb_f_mod, only: optorb_cum,optorb_init
      use pathak_mod, only: ipathak
      use pcm_dmc, only: pcm_prt
      use pcm_mod, only: pcm_init
      use pcm_reduce_mod, only: pcm_reduce
      use precision_kinds, only: dp
      use prop_dmc, only: prop_prt_dmc
      use prop_reduce_mod, only: prop_reduce
      use properties_mod, only: prop_init
      use system,    only: ncent, nelec

      implicit none

      integer :: i, iegerr, ierr, ifgerr
      integer :: ifr, ipeerr, itpber, ic, iph
      integer :: k, npass
      real(dp) :: e2collect
      real(dp) :: e2sum, ecollect, ef2collect, ef2sum
      real(dp) :: efcollect, efnow, egave, egave1
      real(dp) :: egerr, egnow
      real(dp) :: enow, fgave
      real(dp) :: fgerr, peave, peerr, penow
      real(dp) :: tpbave, tpberr
      real(dp) :: tpbnow, w, w2, w2collect
      real(dp) :: w2sum, wcollect, wf2collect, wf2sum
      real(dp) :: wfcollect, wfnow, wgnow, wnow
      real(dp) :: x, x2
      real(dp), dimension(3,ncent,PTH) :: derivgerr
      integer, dimension(3,ncent,PTH) :: iderivgerr
      real(dp), dimension(MFORCE) :: egcollect
      real(dp), dimension(MFORCE) :: wgcollect
      real(dp), dimension(MFORCE) :: pecollect
      real(dp), dimension(MFORCE) :: tpbcollect
      real(dp), dimension(MFORCE) :: eg2collect
      real(dp), dimension(MFORCE) :: wg2collect
      real(dp), dimension(MFORCE) :: pe2collect
      real(dp), dimension(MFORCE) :: tpb2collect
      real(dp), dimension(MFORCE) :: fsum
      real(dp), dimension(MFORCE) :: f2sum
      real(dp), dimension(MFORCE) :: eg2sum
      real(dp), dimension(MFORCE) :: wg2sum
      real(dp), dimension(MFORCE) :: pe2sum
      real(dp), dimension(MFORCE) :: tpb2sum
      real(dp), dimension(MFORCE) :: taucollect
      real(dp), dimension(MFORCE) :: fcollect
      real(dp), dimension(MFORCE) :: f2collect
      real(dp), dimension(nelec) :: ecollect_i
      real(dp), dimension(nfrag) :: ecollectfrag
      real(dp), dimension(nfrag) :: egcollectfrag
      real(dp), dimension(nfrag) :: eg2collectfrag
      real(dp), dimension(nfrag) :: egnowfrag
      real(dp), dimension(nfrag) :: eg2sumfrag
      real(dp), dimension(3,3,ncent,PTH) :: derivcollect
      real(dp), dimension(3,ncent,PTH) :: derivsum2
      real(dp), dimension(3,ncent,PTH) :: derivcollect2
      real(dp), parameter :: zero = 0.d0
      real(dp), parameter :: one = 1.d0

      if(mode.eq.'dmc_one_mpi2') then
        call acuest_gpop
        return
      endif

! wt   = weight of configurations
! xsum = sum of values of x from dmc
! xnow = average of values of x from dmc
! xcum = accumulated sums of xnow
! xcm2 = accumulated sums of xnow**2
! xave = current average value of x
! xerr = current error of x

      iblk=iblk+1
      iblk_proc=iblk_proc+nproc

      npass=iblk_proc*dmc_nstep

      if(idmc.gt.0) then
         wnow=wsum_dmc/dmc_nstep
         wfnow=wfsum/dmc_nstep
         enow=esum_dmc/wsum_dmc
         efnow=efsum/wfsum

         w2sum=wsum_dmc**2
         wf2sum=wfsum**2
         e2sum=esum_dmc*enow
         ef2sum=efsum*efnow
      endif

      do ifr=1,nforce
        wgnow=wgsum(ifr)/dmc_nstep
        egnow=egsum(ifr)/wgsum(ifr)
        penow=pesum_dmc(ifr)/wgsum(ifr)
        tpbnow=tpbsum_dmc(ifr)/wgsum(ifr)

        wg2sum(ifr)=wgsum(ifr)**2
        eg2sum(ifr)=egsum(ifr)*egnow
        pe2sum(ifr)=pesum_dmc(ifr)*penow
        tpb2sum(ifr)=tpbsum_dmc(ifr)*tpbnow
        if(ifr.gt.1) then
          fsum(ifr)=wgsum(1)*(egnow-egsum(1)/wgsum(1))
          f2sum(ifr)=wgsum(1)*(egnow-egsum(1)/wgsum(1))**2
        endif
      enddo

      if (icut_e.lt.0) then
        call mpi_allreduce(esum_i,ecollect_i,nelec &
        ,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)
        ecum_i(:) = ecum_i(:) + ecollect_i(:)
      end if
      if (nfrag.gt.1) then
        egnowfrag(:)=egsumfrag(:)/wgsum(1)
        eg2sumfrag(:)=egsumfrag(:)*egnowfrag(:)
        
        call mpi_allreduce(egsumfrag,egcollectfrag,nfrag &
        ,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)
        call mpi_allreduce(esumfrag,ecollectfrag,nfrag &
        ,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)

        egcumfrag(:)=egcumfrag(:)+egcollectfrag(:)
        ecumfrag(:) = ecumfrag(:) + ecollectfrag(:)
        call mpi_reduce(eg2sumfrag,eg2collectfrag,nfrag,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)
        egcm2frag(:)=egcm2frag(:)+eg2collectfrag(:)
      endif

      call mpi_allreduce(wgsum,wgcollect,MFORCE &
      ,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(egsum,egcollect,MFORCE &
      ,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(tausum,taucollect,MFORCE &
      ,mpi_double_precision,mpi_sum,MPI_COMM_WORLD,ierr)

      do ifr=1,nforce
        wgcum(ifr)=wgcum(ifr)+wgcollect(ifr)
        egcum(ifr)=egcum(ifr)+egcollect(ifr)
        taucum(ifr)=taucum(ifr)+taucollect(ifr)
      enddo

      egave=egcum(1)/wgcum(1)
      if (iforce_analy.gt.0) then
        do iph=1,PTH
          do ic=1,ncent
            do k=1,3
              derivsum2(k,ic,iph) = (derivsum(1,k,ic,iph)+2.d0*derivsum(2,k,ic,iph)-2.d0*egave*derivsum(3,k,ic,iph))**2/wgsum(1)
            enddo
          enddo
        enddo
      endif

      call mpi_reduce(pesum_dmc,pecollect,MFORCE,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)
      call mpi_reduce(tpbsum_dmc,tpbcollect,MFORCE,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

      call mpi_reduce(wg2sum,wg2collect,MFORCE,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)
      call mpi_reduce(eg2sum,eg2collect,MFORCE,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)
      call mpi_reduce(pe2sum,pe2collect,MFORCE,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)
      call mpi_reduce(tpb2sum,tpb2collect,MFORCE,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

      call mpi_reduce(fsum,fcollect,MFORCE,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)
      call mpi_reduce(f2sum,f2collect,MFORCE,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

      call mpi_reduce(derivsum,derivcollect,3*3*ncent*PTH,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

      call mpi_reduce(derivsum2,derivcollect2,3*ncent*PTH,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

      call mpi_reduce(esum_dmc,ecollect,1,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)
      call mpi_reduce(wsum_dmc,wcollect,1,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)
      call mpi_reduce(efsum,efcollect,1,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)
      call mpi_reduce(wfsum,wfcollect,1,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)

      call mpi_reduce(e2sum,e2collect,1,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)
      call mpi_reduce(w2sum,w2collect,1,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)
      call mpi_reduce(ef2sum,ef2collect,1,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)
      call mpi_reduce(wf2sum,wf2collect,1,mpi_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)


      call optjas_cum(wgsum(1),egnow)
      call optorb_cum(wgsum(1),egsum(1))
      call optci_cum(wgsum(1))
      call force_analy_cum(wgsum(1),egcum(1)/wgcum(1))

      call prop_reduce(wgsum(1))
      call pcm_reduce(wgsum(1))
      call mmpol_reduce(wgsum(1))

      if(.not.wid) goto 17

      wcm2=wcm2+w2collect
      wfcm2=wfcm2+wf2collect
      ecm2_dmc=ecm2_dmc+e2collect
      efcm2=efcm2+ef2collect

      wcum_dmc=wcum_dmc+wcollect
      wfcum=wfcum+wfcollect
      ecum_dmc=ecum_dmc+ecollect
      efcum=efcum+efcollect

      do ifr=1,nforce
        wgcm2(ifr)=wgcm2(ifr)+wg2collect(ifr)
        egcm2(ifr)=egcm2(ifr)+eg2collect(ifr)
        pecm2_dmc(ifr)=pecm2_dmc(ifr)+pe2collect(ifr)
        tpbcm2_dmc(ifr)=tpbcm2_dmc(ifr)+tpb2collect(ifr)

        pecum_dmc(ifr)=pecum_dmc(ifr)+pecollect(ifr)
        tpbcum_dmc(ifr)=tpbcum_dmc(ifr)+tpbcollect(ifr)

        if(iblk.eq.1) then
          egerr=0
          peerr=0
          tpberr=0
         else
          egerr=errg(egcum(ifr),egcm2(ifr),ifr)
          peerr=errg(pecum_dmc(ifr),pecm2_dmc(ifr),ifr)
          tpberr=errg(tpbcum_dmc(ifr),tpbcm2_dmc(ifr),ifr)
        endif

        egave=egcum(ifr)/wgcum(ifr)
        peave=pecum_dmc(ifr)/wgcum(ifr)
        tpbave=tpbcum_dmc(ifr)/wgcum(ifr)

        if(ifr.gt.1) then
          fgcum(ifr)=fgcum(ifr)+fcollect(ifr)
          fgcm2(ifr)=fgcm2(ifr)+f2collect(ifr)
          fgave=egcum(1)/wgcum(1)-egcum(ifr)/wgcum(ifr)
          if(iblk.eq.1) then
            fgerr=0
            ifgerr=0
           else
            fgerr=errg(fgcum(ifr),fgcm2(ifr),1)
            ifgerr=nint(1e12* fgerr)
          endif

        else

          egave=egcum(1)/wgcum(1)
          if (iforce_analy.gt.0) then
            do iph=1,PTH
              do ic=1,ncent
                do k=1,3
                  derivcum(1,k,ic,iph)=derivcum(1,k,ic,iph)+derivcollect(1,k,ic,iph)
                  derivcum(2,k,ic,iph)=derivcum(2,k,ic,iph)+derivcollect(2,k,ic,iph)
                  derivcum(3,k,ic,iph)=derivcum(3,k,ic,iph)+derivcollect(3,k,ic,iph)
                  derivcm2(k,ic,iph)=derivcm2(k,ic,iph)+derivcollect2(k,ic,iph)
                enddo
              enddo
            enddo
          end if
          call prop_prt_dmc(iblk,0,wgcum,wgcm2)
          call pcm_prt(iblk,wgcum,wgcm2)
          call mmpol_prt(iblk,wgcum,wgcm2)
        endif

! write out header first time

        if (iblk.eq.1.and.ifr.eq.1) then
          if(nforce.gt.1) then
            write(ounit,'(t5,''egnow'',t15,''egave'',t21,''(egerr)'' ,t32 &
              ,''peave'',t38,''(peerr)'',t49,''tpbave'',t55,''(tpberr)'',t66 &
              ,''fgave'',t79,''(fgerr)'',t93,''npass'',t102,''wgsum'',t112   &
              ,''ioldest'')')
          else
            write(ounit,'(t5,''egnow'',t15,''egave'',t21,''(egerr)'' ,t32&
              ,''peave'',t38,''(peerr)'',t49,''tpbave'',t55,''(tpberr)'',t67&
              ,''npass'',t77,''wgsum'',t85,''ioldest'')')
          endif
        endif

! write out current values of averages etc.

        iegerr=nint(100000* egerr)
        ipeerr=nint(100000* peerr)
        itpber=nint(100000*tpberr)

        if(ifr.eq.1) then
          if(nforce.gt.1) then
            write(ounit,'(f10.5,3(f10.5,''('',i5,'')''),62x,3i10)') &
            egcollect(ifr)/wgcollect(ifr), &
            egave,iegerr,peave,ipeerr,tpbave,itpber,npass, &
            nint(wgcollect(ifr)/nproc),ioldest
           else
            write(ounit,'(f10.5,3(f10.5,''('',i5,'')''),3i10)') &
            egcollect(ifr)/wgcollect(ifr), &
            egave,iegerr,peave,ipeerr,tpbave,itpber,npass, &
            nint(wgcollect(ifr)/nproc),ioldest
          endif
         else
          write(ounit,'(f10.5,3(f10.5,''('',i5,'')''),f17.12, ''('',i12,'')'',10x,i10)') &
            egcollect(ifr)/wgcollect(ifr), &
            egave,iegerr,peave,ipeerr,tpbave,itpber, &
            fgave,ifgerr,nint(wgcollect(ifr)/nproc)
        endif
      enddo

! zero out xsum variables for metrop

      17 wsum_dmc=zero
      wfsum=zero
      esum_dmc=zero
      efsum=zero
      if (icut_e.lt.0) then
        esum_i = zero
      endif

      if (nfrag.gt.1) then
        esumfrag = zero
        egsumfrag(:)=zero
      endif

      do ifr=1,nforce
        egsum(ifr)=zero
        wgsum(ifr)=zero
        pesum_dmc(ifr)=zero
        tpbsum_dmc(ifr)=zero
        tausum(ifr)=zero
      enddo
      do iph=1,PTH
        do k=1,3
          do ic=1,ncent
            derivsum(1,k,ic,iph)=zero
            derivsum(2,k,ic,iph)=zero
            derivsum(3,k,ic,iph)=zero
          enddo
        enddo
      enddo

      call optorb_init(1)
      call optci_init(1)

      call prop_init(1)
      call pcm_init(1)
      call mmpol_init(1)
      call force_analy_init(1)

      return
contains
        elemental pure function rn_eff(w,w2)
          implicit none
          real(dp), intent(in) :: w, w2
          real(dp)             :: rn_eff
          rn_eff=w**2/w2
        end function
        elemental pure function error(x,x2,w,w2)
          implicit none
          real(dp), intent(in) :: x, x2,w,w2
          real(dp)             :: error
          error=dsqrt(max((x2/w-(x/w)**2)/(rn_eff(w,w2)-1),0.d0))
        end function
        elemental pure function errg(x,x2,i)
          implicit none
          real(dp), intent(in) :: x, x2
          integer, intent(in)  :: i
          real(dp)             :: errg
          errg=error(x,x2,wgcum(i),wgcm2(i))
        end function
      end
end module
