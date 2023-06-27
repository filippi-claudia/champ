module pcm_dmc
contains
      subroutine pcm_prt(iblk,wgcum,wgcm2)

      use contrl_file, only: ounit
      use control_dmc, only: dmc_nconf,dmc_nstep
      use multiple_geo, only: MFORCE
      use pcm_averages, only: qopcm_cm2,qopcm_cum,spcmcm2,spcmcum
      use pcm_averages, only: vpcmcm2,vpcmcum
      use pcm_cntrl, only: ipcm,ipcmprt
      use precision_kinds, only: dp
!      use contrl, only: nconf, nstep
      implicit none

      integer :: i, iblk, iqopcm_err, ispcmerr, ivpcmerr
      real(dp) :: evalg_eff, hatokc, qopcm_ave
      real(dp) :: qopcm_err, rtevalg_eff1, sepcmkcal
      real(dp) :: spcmave, spcmerr, spcmkcal, vepcmkcal
      real(dp) :: vpcmave, vpcmerr, vpcmkcal, w
      real(dp), dimension(MFORCE) :: wgcum
      real(dp), dimension(MFORCE) :: wgcm2

      data hatokc/627.509541d0/

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
      write(ounit,'(''pcm dG(surf) ='',t17,f12.7,'' +-'',f11.7,f9.5,2x,f12.7,'' +-'',f11.7)') &
       spcmave,spcmerr,spcmerr*rtevalg_eff1,spcmkcal,sepcmkcal
      write(ounit,'(''pcm dG(vol)  ='',t17,f12.7,'' +-'',f11.7,f9.5,2x,f12.7,''+-'',f11.7)') &
       vpcmave,vpcmerr,vpcmerr*rtevalg_eff1,vpcmkcal,vepcmkcal
!     write(ounit,'(''pcm qout     ='',t17,f12.7,'' +-'',f11.7,f9.5)')
!    & qopcm_ave,qopcm_err,qopcm_err*rtevalg_eff1

!     gpcmkcal=spcmkcal+vpcmkcal

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
!-----------------------------------------------------------------------
      subroutine pcm_fin(iblk,wgcum,wgcm2)

      use multiple_geo, only: MFORCE
      use pcm_cntrl, only: ipcm,ipcmprt
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
!-----------------------------------------------------------------------
      subroutine pcm_save(iw)

      use pcm_cntrl, only: ipcm
      use pcm_hpsi, only: enfpcm,pcms,pcmv,qopcm
      use pcm_parms, only: nchs
      use pcmo,    only: enfpcmo_dmc,qopcmo_dmc,spcmo_dmc,vpcmo_dmc

      implicit none

      integer :: i, iw


      if(ipcm.eq.0) return

      spcmo_dmc(iw)=pcms
      vpcmo_dmc(iw)=pcmv
      qopcmo_dmc(iw)=qopcm

!     write(ounit,*) 'CIAO',qopcm,qopcmo_dmc(iw),iw,spcmo_dmc(iw),vpcmo_dmc(iw)
      do i=1,nchs
      enfpcmo_dmc(iw,i)=enfpcm(i)
      enddo

      return
      end
!-----------------------------------------------------------------------
      subroutine pcm_sum(p,q,iw)

      use pcm_averages, only: enfpcm_sum,qopcm_sum,spcmsum,vpcmsum
      use pcm_cntrl, only: ipcm
      use pcm_hpsi, only: enfpcm,pcms,pcmv,qopcm
      use pcm_parms, only: nchs
      use pcmo,    only: enfpcmo_dmc,qopcmo_dmc,spcmo_dmc,vpcmo_dmc
      use precision_kinds, only: dp

      implicit none

      integer :: i, iw
      real(dp) :: p, q

      if(ipcm.eq.0) return

      spcmsum=spcmsum+p*pcms+q*spcmo_dmc(iw)
      vpcmsum=vpcmsum+p*pcmv+q*vpcmo_dmc(iw)
      qopcm_sum=qopcm_sum+p*qopcm+q*qopcmo_dmc(iw)

!     write(ounit,*) 'HELLO',qopcm,qopcmo_dmc(iw),iw

      do i=1,nchs
      enfpcm_sum(i)= enfpcm_sum(i)+p*enfpcm(i)+q*enfpcmo_dmc(iw,i)
      enddo

      return
      end
!-----------------------------------------------------------------------
      subroutine pcm_cum(wsum_dmc)

      use pcm_averages, only: enfpcm_cm2,enfpcm_cum,enfpcm_sum,qopcm_cm2
      use pcm_averages, only: qopcm_cum,qopcm_sum,spcmcm2,spcmcum
      use pcm_averages, only: spcmsum,vpcmcm2,vpcmcum,vpcmsum
      use pcm_cntrl, only: ipcm
      use pcm_parms, only: nchs
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
!     write (6,*) 'HELLO-CIAO', qopcm_cum

      do i=1,nchs
      enfpcm_now=enfpcm_sum(i)/wsum_dmc
      enfpcm_cm2(i)=enfpcm_cm2(i)+enfpcm_sum(i)*enfpcm_now
      enfpcm_cum(i)=enfpcm_cum(i)+enfpcm_sum(i)
      enddo

      return
      end
!-----------------------------------------------------------------------
      subroutine pcm_splitj(iw,iw2)

      use pcm_parms, only: nchs
      use pcmo,    only: enfpcmo_dmc,qopcmo_dmc,spcmo_dmc,vpcmo_dmc

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
