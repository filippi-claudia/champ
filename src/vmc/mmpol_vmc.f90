module mmpol_vmc
contains
      subroutine mmpol_prt(wcum,iblk)

      use contrl_file, only: ounit
      use mmpol_averages, only: cmmpol_cm2,cmmpol_cum,dmmpol_cm2
      use mmpol_averages, only: dmmpol_cum
      use mmpol_cntrl, only: immpol,immpolprt
      use precision_kinds, only: dp
      implicit none

      integer :: iblk, icmmpol_err, idmmpol_err
      real(dp) :: cmmpol_ave, cmmpol_err, cmmpol_merr, dmmpol_ave, dmmpol_err
      real(dp) :: qmmpol_ave, qmmpol_err, rtpass
      real(dp) :: wcum



      if(immpol.eq.0.or.immpolprt.eq.0) return

      dmmpol_ave=dmmpol_cum/wcum
      cmmpol_ave=cmmpol_cum/wcum

      if(iblk.eq.1) then
        dmmpol_err=0
        cmmpol_err=0

       else
        dmmpol_err=err(dmmpol_cum,dmmpol_cm2)
        cmmpol_err=err(cmmpol_cum,cmmpol_cm2)

        idmmpol_err=nint(100000*dmmpol_err)
        icmmpol_err=nint(100000*cmmpol_err)
      endif

      rtpass=dsqrt(wcum)
      if(immpol.eq.2)then
        write(ounit,'(''mmpol dE (induced dipoles) ='',t17,f12.7,'' +-'',f11.7,f9.5)') dmmpol_ave,dmmpol_err,dmmpol_err*rtpass
        write(ounit,'(''mmpol dE  (fixed charges)  ='',t17,f12.7,'' +-'',f11.7,f9.5)') cmmpol_ave,cmmpol_err,cmmpol_merr*rtpass
        write(ounit,'(''<H(QM/MM)>'',t17,f12.7,'' +-'',f11.7,f9.5)') dmmpol_ave+cmmpol_ave
      endif
      if(immpol.eq.3)then
        qmmpol_ave=dmmpol_ave+cmmpol_ave
        qmmpol_err=dsqrt(dmmpol_err**2.0d0+cmmpol_err**2.0d0)
        write(ounit,'(''mmpol <H(QM/MM)> ='',t17,f12.7,'' +-'',f11.7,f9.5)') qmmpol_ave,qmmpol_err,qmmpol_err*rtpass
      endif

      return
contains
        elemental pure function err(x,x2)
          implicit none
          real(dp), intent(in) :: x, x2
          real(dp)             :: err
          err=dsqrt(abs(x2/wcum-(x/wcum)**2)/iblk)
        end function
      end
!-----------------------------------------------------------------------
      subroutine mmpol_fin(wcum,iblk)

      use contrl_file, only: ounit
      use mmpol,   only: mmpol_dipoles
      use mmpol_averages, only: cmmpol_cm2,cmmpol_cum,dmmpol_cm2
      use mmpol_averages, only: dmmpol_cum,eek1_cm2,eek1_cum,eek2_cm2
      use mmpol_averages, only: eek2_cum,eek3_cm2,eek3_cum
      use mmpol_cntrl, only: immpol
      use mmpol_mod, only: MCHMM
      use mmpol_parms, only: nchmm
      use precision_kinds, only: dp
      implicit none

      integer :: i, iblk, icmmpol_err, idmmpol_err
      real(dp) :: cekcal, cmmpol_ave, cmmpol_err, cmmpol_kcal, cmmpol_merr
      real(dp) :: dckcal, dekcal, dmmpol_ave, dmmpol_err
      real(dp) :: dmmpol_kcal, hatokc, qmmpol_ave
      real(dp) :: qmmpol_err, rtpass, wcum
      real(dp), dimension(MCHMM) :: eek1_ave
      real(dp), dimension(MCHMM) :: eek1_err
      real(dp), dimension(MCHMM) :: eek2_ave
      real(dp), dimension(MCHMM) :: eek2_err
      real(dp), dimension(MCHMM) :: eek3_ave
      real(dp), dimension(MCHMM) :: eek3_err
      real(dp), dimension(3, MCHMM) :: eek_ave
      real(dp), dimension(3, MCHMM) :: eek_err



      data hatokc/627.509541d0/

      if(immpol.eq.0) return

      rtpass=dsqrt(wcum)

      dmmpol_ave=dmmpol_cum/wcum
      cmmpol_ave=cmmpol_cum/wcum

      dmmpol_err=err(dmmpol_cum,dmmpol_cm2)
      cmmpol_err=err(cmmpol_cum,cmmpol_cm2)

      idmmpol_err=nint(100000*dmmpol_err)
      icmmpol_err=nint(100000*cmmpol_err)

      if(immpol.eq.2)then
        write(ounit,*)
        write(ounit,'(''mmpol dE (induced dipoles) ='',t17,f12.7,'' +-'',f11.7,f9.5)') dmmpol_ave,dmmpol_err,dmmpol_err*rtpass
        write(ounit,'(''mmpol dE  (fixed charges)  ='',t17,f12.7,'' +-'',f11.7,f9.5)') cmmpol_ave,cmmpol_err,cmmpol_merr*rtpass
      endif

      if(immpol.eq.3)then
        qmmpol_ave=dmmpol_ave+cmmpol_ave
        qmmpol_err=dsqrt(dmmpol_err**2.0d0+cmmpol_err**2.0d0)
        write(ounit,'(''mmpol <H(QM/MM)> ='',t17,f12.7,'' +-'',f11.7,f9.5)') qmmpol_ave,qmmpol_err,qmmpol_err*rtpass
      endif

      if(immpol.ne.3)then
        do  i=1,nchmm
          eek1_ave(i)=eek1_cum(i)/wcum
          eek1_err(i)=err(eek1_cum(i),eek1_cm2(i))
          eek_ave(1,i)=eek1_cum(i)/wcum
          eek_err(1,i)=err(eek1_cum(i),eek1_cm2(i))

          eek2_ave(i)=eek2_cum(i)/wcum
          eek2_err(i)=err(eek2_cum(i),eek2_cm2(i))
          eek_ave(2,i)=eek2_cum(i)/wcum
          eek_err(2,i)=err(eek2_cum(i),eek2_cm2(i))

          eek3_ave(i)=eek3_cum(i)/wcum
          eek3_err(i)=err(eek3_cum(i),eek3_cm2(i))
          eek_ave(3,i)=eek3_cum(i)/wcum
          eek_err(3,i)=err(eek3_cum(i),eek3_cm2(i))
        enddo
         call mmpol_dipoles(eek_ave,eek_err)

        dmmpol_kcal=dmmpol_ave*hatokc
        cmmpol_kcal=cmmpol_ave*hatokc
        dekcal=dmmpol_err*hatokc
        cekcal=cmmpol_err*hatokc
        dckcal=dmmpol_kcal+cmmpol_kcal
        write(ounit,*)'    <H(QM/MM)/dipoles/charges/tot +- err (kcal/mol) '
        write(ounit,1000)dmmpol_kcal,dekcal,cmmpol_kcal,cekcal,dckcal
        write(ounit,*)'    <H(QM/MM)/dipoles/charges/tot +- err (hartree) '
        write(ounit,1000)dmmpol_ave,dmmpol_err,cmmpol_ave,cmmpol_err,dckcal/hatokc
        write(ounit,*)
      endif


      1000 format(9F12.5)

      return
contains
        elemental pure function err(x,x2)
          implicit none
          real(dp), intent(in) :: x, x2
          real(dp)             :: err
          err=dsqrt(abs(x2/wcum-(x/wcum)**2)/iblk)
        end function
      end
!-----------------------------------------------------------------------
      subroutine mmpol_save
      use mmpol_cntrl, only: immpol
      use mmpol_hpsi, only: eek_pol
      use mmpol_parms, only: nchmm
      use mmpolo,  only: cmmpolo,dmmpolo,eeko
      use precision_kinds, only: dp
      implicit none

      integer :: i
      real(dp) :: QMdp, QMq






      if(immpol.eq.0) return
      dmmpolo=QMdp
      cmmpolo=QMq

      do i=1,nchmm
        eeko(1,i)=eek_pol(1,i)
        eeko(2,i)=eek_pol(2,i)
        eeko(3,i)=eek_pol(3,i)
      enddo

      return
      end
!-----------------------------------------------------------------------
      subroutine mmpol_sum(p,q)

      use mmpol_averages, only: cmmpol_sum,dmmpol_sum,eek_sum
      use mmpol_cntrl, only: immpol
      use mmpol_hpsi, only: eek_pol
      use mmpol_parms, only: nchmm
      use mmpolo,  only: cmmpolo,dmmpolo,eeko
      use precision_kinds, only: dp

      implicit none

      integer :: i
      real(dp) :: QMdp, QMq, p, q


      if(immpol.eq.0) return
      dmmpol_sum=dmmpol_sum+p*QMdp+q*dmmpolo
      cmmpol_sum=cmmpol_sum+p*QMq+q*cmmpolo

      do i=1,nchmm
        eek_sum(1,i)= eek_sum(1,i)+p*eek_pol(1,i)+q*eeko(1,i)
        eek_sum(2,i)= eek_sum(2,i)+p*eek_pol(2,i)+q*eeko(2,i)
        eek_sum(3,i)= eek_sum(3,i)+p*eek_pol(3,i)+q*eeko(3,i)
      enddo

      return
      end
!-----------------------------------------------------------------------
      subroutine mmpol_cum(wsum)

      use mmpol_averages, only: cmmpol_cm2,cmmpol_cum,cmmpol_sum
      use mmpol_averages, only: dmmpol_cm2,dmmpol_cum,dmmpol_sum
      use mmpol_averages, only: eek1_cm2,eek1_cum,eek2_cm2,eek2_cum
      use mmpol_averages, only: eek3_cm2,eek3_cum,eek_sum
      use mmpol_cntrl, only: immpol
      use mmpol_parms, only: nchmm
      use precision_kinds, only: dp

      implicit none

      integer :: i
      real(dp) :: cmmpolnow, dmmpolnow, eek_now1, eek_now2, eek_now3
      real(dp) :: wsum


      if(immpol.eq.0) return
      dmmpolnow=dmmpol_sum/wsum
      cmmpolnow=cmmpol_sum/wsum

      dmmpol_cm2=dmmpol_cm2+dmmpol_sum*dmmpolnow
      cmmpol_cm2=cmmpol_cm2+cmmpol_sum*cmmpolnow

      dmmpol_cum=dmmpol_cum+dmmpol_sum
      cmmpol_cum=cmmpol_cum+cmmpol_sum

      do i=1,nchmm
        eek_now1=eek_sum(1,i)/wsum
        eek1_cm2(i)=eek1_cm2(i)+eek_sum(1,i)*eek_now1
        eek1_cum(i)=eek1_cum(i)+eek_sum(1,i)

        eek_now2=eek_sum(2,i)/wsum
        eek2_cm2(i)=eek2_cm2(i)+eek_sum(2,i)*eek_now2
        eek2_cum(i)=eek2_cum(i)+eek_sum(2,i)

        eek_now3=eek_sum(3,i)/wsum
        eek3_cm2(i)=eek3_cm2(i)+eek_sum(3,i)*eek_now3
        eek3_cum(i)=eek3_cum(i)+eek_sum(3,i)
      enddo

      return
      end
end module
