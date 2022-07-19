      module pcm_dmc
      contains
      subroutine pcm_prt(iblk,wgcum,wgcm2)

      use force_mod, only: MFORCE
!      use contrl, only: nconf, nstep
      use control_dmc, only: dmc_nconf, dmc_nstep
      use pcm_cntrl, only: ipcm, ipcmprt
      use pcm_averages, only: spcmcum, spcmcm2, vpcmcum, vpcmcm2
      use pcm_averages, only: qopcm_cum, qopcm_cm2
      use precision_kinds, only: dp
      use contrl_file, only: ounit
      implicit none

      integer :: i, iblk, iqopcm_err, ispcmerr, ivpcmerr
      real(dp) :: errg, error, evalg_eff, hatokc, qopcm_ave
      real(dp) :: qopcm_err, rn_eff, rtevalg_eff1, sepcmkcal
      real(dp) :: spcmave, spcmerr, spcmkcal, vepcmkcal
      real(dp) :: vpcmave, vpcmerr, vpcmkcal, w
      real(dp) :: w2, x, x2
      real(dp), dimension(MFORCE) :: wgcum
      real(dp), dimension(MFORCE) :: wgcm2

      data hatokc/627.509541d0/

c Statement functions for error calculation, it might be reaplaced in the near future:
      rn_eff(w,w2)=w**2/w2
      error(x,x2,w,w2)=dsqrt(max((x2/w-(x/w)**2)/(rn_eff(w,w2)-1),0.d0))
      errg(x,x2,i)=error(x,x2,wgcum(i),wgcm2(i))

      if(ipcm.eq.0.or.ipcmprt.eq.0) return

      spcmave=spcmcum/wgcum(1)
      vpcmave=vpcmcum/wgcum(1)
      qopcm_ave=qopcm_cum/wgcum(1)

      if(iblk.eq.1) then
        spcmerr=0
        vpcmerr=0
        qopcm_err=0
       else
        spcmerr=errg(spcmcum,spcmcm2,1)
        vpcmerr=errg(vpcmcum,vpcmcm2,1)
        qopcm_err=errg(qopcm_cum,qopcm_cm2,1)

        ispcmerr=nint(100000*spcmerr)
        ivpcmerr=nint(100000*vpcmerr)
        iqopcm_err=nint(100000*qopcm_err)
      endif

      evalg_eff=dmc_nconf*dmc_nstep*rn_eff(wgcum(1),wgcm2(1))
      rtevalg_eff1=dsqrt(evalg_eff-1)

      spcmkcal=spcmave*hatokc
      vpcmkcal=vpcmave*hatokc
      sepcmkcal=spcmerr*hatokc
      vepcmkcal=vpcmerr*hatokc
      write(ounit,'(''pcm dG(surf) ='',t17,f12.7,'' +-'',f11.7,f9.5,2x,f12.7,'' +-'',f11.7)')
     & spcmave,spcmerr,spcmerr*rtevalg_eff1,spcmkcal,sepcmkcal
      write(ounit,'(''pcm dG(vol)  ='',t17,f12.7,'' +-'',f11.7,f9.5,2x,f12.7,''+-'',f11.7)')
     & vpcmave,vpcmerr,vpcmerr*rtevalg_eff1,vpcmkcal,vepcmkcal
c     write(ounit,'(''pcm qout     ='',t17,f12.7,'' +-'',f11.7,f9.5)')
c    & qopcm_ave,qopcm_err,qopcm_err*rtevalg_eff1

c     gpcmkcal=spcmkcal+vpcmkcal

      return
      end
c-----------------------------------------------------------------------
      subroutine pcm_fin(iblk,wgcum,wgcm2)

      use force_mod, only: MFORCE
      use pcm_cntrl, only: ipcm, ipcmprt

      use precision_kinds, only: dp
      implicit none

      integer :: iblk, ipcmprt_sav

      real(dp), dimension(MFORCE) :: wgcum
      real(dp), dimension(MFORCE) :: wgcm2



      if(ipcm.eq.0) return

      ipcmprt_sav=ipcmprt
      ipcmprt=1
      call pcm_prt(iblk,wgcum(1),wgcm2(1))
      ipcmprt=ipcmprt_sav

      return
      end
c-----------------------------------------------------------------------
      subroutine pcm_save(iw)

      use pcm_hpsi, only: pcms, pcmv, qopcm, enfpcm
      use pcmo, only: spcmo_dmc, vpcmo_dmc, qopcmo_dmc, enfpcmo_dmc
      use pcm_cntrl, only: ipcm
      use pcm_parms, only: nchs

      implicit none

      integer :: i, iw


      if(ipcm.eq.0) return

      spcmo_dmc(iw)=pcms
      vpcmo_dmc(iw)=pcmv
      qopcmo_dmc(iw)=qopcm

c     write(ounit,*) 'CIAO',qopcm,qopcmo_dmc(iw),iw,spcmo_dmc(iw),vpcmo_dmc(iw)
      do i=1,nchs
      enfpcmo_dmc(iw,i)=enfpcm(i)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine pcm_sum(p,q,iw)

      use pcm_hpsi, only: pcms, pcmv, qopcm, enfpcm
      use pcmo, only: spcmo_dmc, vpcmo_dmc, qopcmo_dmc, enfpcmo_dmc
      use pcm_cntrl, only: ipcm
      use pcm_parms, only: nchs
      use pcm_averages, only: spcmsum, vpcmsum
      use pcm_averages, only: qopcm_sum, enfpcm_sum

      use precision_kinds, only: dp
      implicit none

      integer :: i, iw
      real(dp) :: p, q

      if(ipcm.eq.0) return

      spcmsum=spcmsum+p*pcms+q*spcmo_dmc(iw)
      vpcmsum=vpcmsum+p*pcmv+q*vpcmo_dmc(iw)
      qopcm_sum=qopcm_sum+p*qopcm+q*qopcmo_dmc(iw)

c     write(ounit,*) 'HELLO',qopcm,qopcmo_dmc(iw),iw

      do i=1,nchs
      enfpcm_sum(i)= enfpcm_sum(i)+p*enfpcm(i)+q*enfpcmo_dmc(iw,i)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine pcm_cum(wsum_dmc)

      use pcm_cntrl, only: ipcm
      use pcm_parms, only: nchs
      use pcm_averages, only: spcmsum, spcmcum, spcmcm2, vpcmsum, vpcmcum, vpcmcm2
      use pcm_averages, only: qopcm_sum, qopcm_cum, qopcm_cm2, enfpcm_sum, enfpcm_cum, enfpcm_cm2

      use precision_kinds, only: dp
      implicit none

      integer :: i
      real(dp) :: enfpcm_now, qopcm_now, spcmnow, vpcmnow, wsum_dmc

      if(ipcm.eq.0) return

      spcmnow=spcmsum/wsum_dmc
      vpcmnow=vpcmsum/wsum_dmc
      qopcm_now=qopcm_sum/wsum_dmc

      spcmcm2=spcmcm2+spcmsum*spcmnow
      vpcmcm2=vpcmcm2+vpcmsum*vpcmnow
      qopcm_cm2=qopcm_cm2+qopcm_sum*qopcm_now

      spcmcum=spcmcum+spcmsum
      vpcmcum=vpcmcum+vpcmsum
      qopcm_cum=qopcm_cum+qopcm_sum
c     write (6,*) 'HELLO-CIAO', qopcm_cum

      do i=1,nchs
      enfpcm_now=enfpcm_sum(i)/wsum_dmc
      enfpcm_cm2(i)=enfpcm_cm2(i)+enfpcm_sum(i)*enfpcm_now
      enfpcm_cum(i)=enfpcm_cum(i)+enfpcm_sum(i)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine pcm_splitj(iw,iw2)

      use pcmo, only: spcmo_dmc, vpcmo_dmc, qopcmo_dmc, enfpcmo_dmc
      use pcm_parms, only: nchs

      implicit none

      integer :: i, iw, iw2


      spcmo_dmc(iw2)=spcmo_dmc(iw)
      vpcmo_dmc(iw2)=vpcmo_dmc(iw)
      qopcmo_dmc(iw2)=qopcmo_dmc(iw)

      do i=1,nchs
      enfpcmo_dmc(iw2,i)=enfpcmo_dmc(iw,i)
      enddo

      return
      end
      end module
