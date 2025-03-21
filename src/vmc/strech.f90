module strech_mod
      use precision_kinds, only: dp
      use contrl_file, only: ounit
      use control, only: mode
      use ewald, only: cos_n_sum, sin_n_sum
      use error,   only: fatal_error
      use matinv_mod, only: matinv
      use multiple_geo, only: MFORCE,MFORCE_WT_PRD,alfstr,delc,istrech
      use multiple_geo, only: itausec,iwftype,nforce,nwprod,pecent
      use optwf_control, only: ioptwf
      use pcm,     only: MCHS,MCHV
      use pcm_ameta, only: eta
      use pcm_cntrl, only: ipcm
      use pcm_force, only: sch_s
      use pcm_inda, only: inda
      use pcm_mod, only: sigma_R
      use pcm_parms, only: ch,nch,nchs,nesph,xpol
      use pcm_pot, only: penups,penupv
      use pot,     only: pot_nn
      use system,  only: cent,iwctype,ncent,ncent_tot,nelec,znuc

      real(dp), ALLOCATABLE, save :: centsav(:,:)
      real(dp), ALLOCATABLE, save :: pecentn(:)
      real(dp), ALLOCATABLE, save :: xpolsav(:,:)

contains
      subroutine strech(x,xstrech,ajacob,ifr,istrech_el)
! Written by Cyrus Umrigar and Claudia Filippi
! Modified by A. Amovilli for forces in PCM
! Uses the coordinate transform described in:
! 1) Two Aspects of Quantum Monte Carlo: Determination of Accurate Wavefunctions and
!    Determination of Potential Energy Surfaces of Molecules, C.J. Umrigar,
!    Int. J. Quant. Chem. Symp., 23, 217 (1989).
! 2) Correlated sampling in quantum Monte Carlo: A route to forces,
!    Claudia Filippi and C. J. Umrigar, Phys. Rev. B., 61, R16291, (2000).

! stretch space so that electrons close to a nucleus move almost
! rigidly with that nucleus

      implicit none

      integer :: i, ic, icent, ifl, ifr
      integer :: index, is, istrech_el, j
      integer :: jc, js, k, l
      real(dp) :: ajacob, cc, cc1, cc2
      real(dp) :: cc3, delta_gpol_fc, delta_qs
      real(dp) :: det, dist, dist2, enk
      real(dp) :: env, penups_fc, penupv_fc, rcm
      real(dp) :: rnp, rnp2, rr2, rr3
      real(dp) :: rsq, rsq1, wtsm, wtsmi
      real(dp) :: xi, xx, yi, yy
      real(dp) :: zi, zz
      real(dp), dimension(3,nelec) :: x
      real(dp), dimension(3,nelec) :: xstrech
      real(dp), dimension(ncent_tot) :: wt
      real(dp), dimension(3,3) :: dvol
      real(dp), dimension(3,ncent_tot) :: dwt
      real(dp), dimension(3) :: dwtsm
      real(dp), dimension(3,ncent_tot) :: cent_str
      real(dp), dimension(MCHS) :: q_strech
      real(dp), dimension(MCHS) :: efsol
      real(dp), dimension(ncent_tot) :: wt_pcm
      real(dp), parameter :: zero = 0.d0
      real(dp), parameter :: one = 1.d0

      if(.not.allocated(centsav)) allocate(centsav(3, ncent_tot))
      if(.not.allocated(pecentn)) allocate(pecentn(MFORCE))
      if(.not.allocated(xpolsav)) allocate(xpolsav(3,MCHV))

! set center and n-n potential for secondary geometries
      pecent=pecentn(ifr)
      do icent=1,ncent
        do k=1,3
          cent(k,icent)=centsav(k,icent)+delc(k,icent,ifr)
        enddo
      enddo

      ajacob=one

! PCM
      if(ipcm.eq.3.and.ioptwf.eq.0) then

! positions of surface charges on spheres rigidly displaced
          do j=1,nchs
            is=inda(j)
            do k=1,3
             xpol(k,j)=xpolsav(k,j)+delc(k,is,ifr)
            enddo
          enddo
! positions of volume charges space warped
          do j=nchs+1,nch
            wtsm=zero
            do icent=1,ncent
              dist2=zero
              do k=1,3
                dist2=dist2+(xpolsav(k,j)-centsav(k,icent))**2
              enddo
                dist=dsqrt(dist2)
                if(istrech.eq.1) wt_pcm(icent)=dexp(-alfstr*dist)
                if(istrech.eq.2) wt_pcm(icent)=one/dist**alfstr
                if(istrech.eq.3) wt_pcm(icent)=dexp(alfstr/dist)
                wtsm=wtsm+wt_pcm(icent)
            enddo
            wtsmi=one/wtsm
            do icent=1,ncent
              wt_pcm(icent)=wt_pcm(icent)*wtsmi
              do k=1,3
                xpol(k,j)=xpolsav(k,j)+wt_pcm(icent)*delc(k,icent,ifr)
              enddo
            enddo
          enddo

          do i=1,nchs
             ch(i)=sch_s(i,ifr)
          enddo

! endif PCM
      endif

      if(istrech_el.eq.0) then
        return
      endif

      do i=1,nelec
        do k=1,3
          xstrech(k,i)=x(k,i)
        enddo
      enddo

      if(istrech.eq.0) then
        return
      endif

      do i=1,nelec

        wtsm=zero
! initialize volume change matrix
        do k=1,3
          dwtsm(k)=zero
          do j=1,3
            dvol(j,k)=zero
            if(j.eq.k) dvol(j,k)=one
          enddo
        enddo

        do icent=1,ncent
          dist2=zero
          do k=1,3
            dist2=dist2+(x(k,i)-centsav(k,icent))**2
          enddo
            dist=dsqrt(dist2)
            if(istrech.eq.1) wt(icent)=dexp(-alfstr*dist)
            if(istrech.eq.2) wt(icent)=one/dist**alfstr
            if(istrech.eq.3) wt(icent)=dexp(alfstr/dist)
            wtsm=wtsm+wt(icent)
            do k=1,3
              if(istrech.eq.1) dwt(k,icent)=-alfstr*dexp(-alfstr*dist)* &
              (x(k,i)-centsav(k,icent))/dist
              if(istrech.eq.2) dwt(k,icent)=-alfstr* &
              (x(k,i)-centsav(k,icent))/dist**(alfstr+2)
              if(istrech.eq.3) dwt(k,icent)=-alfstr*dexp(alfstr/dist)* &
              (x(k,i)-centsav(k,icent))/(dist2*dist)
              dwtsm(k)=dwtsm(k)+dwt(k,icent)
            enddo
        enddo
          wtsmi=one/wtsm

          do icent=1,ncent
            do k=1,3
              dwt(k,icent)=(wtsm*dwt(k,icent)-wt(icent)*dwtsm(k))/ &
              wtsm**2
              dvol(1,k)=dvol(1,k)+dwt(k,icent)*delc(1,icent,ifr)
              dvol(2,k)=dvol(2,k)+dwt(k,icent)*delc(2,icent,ifr)
              dvol(3,k)=dvol(3,k)+dwt(k,icent)*delc(3,icent,ifr)
            enddo
          enddo
          call matinv(dvol,3,det)
          ajacob=ajacob*det
          do icent=1,ncent
            wt(icent)=wt(icent)*wtsmi
            do k=1,3
              xstrech(k,i)=xstrech(k,i)+wt(icent)*delc(k,icent,ifr)
! end loop over electrons
            enddo
          enddo
      enddo

      return
      end subroutine

! Set up n-n potential energy (and PCM related quantities) at displaced positions
      subroutine setup_force
      implicit none

      integer :: i, ic, icent, ifl, ifr
      integer :: index, is, istrech_el, j
      integer :: jc, js, k, l
      real(dp) :: ajacob, cc, cc1, cc2
      real(dp) :: cc3, delta_gpol_fc, delta_qs
      real(dp) :: det, dist, dist2, enk
      real(dp) :: env, penups_fc, penupv_fc, rcm
      real(dp) :: rnp, rnp2, rr2, rr3
      real(dp) :: rsq, rsq1, wtsm, wtsmi
      real(dp) :: xi, xx, yi, yy
      real(dp) :: zi, zz
      real(dp), dimension(3,nelec) :: x
      real(dp), dimension(3,nelec) :: xstrech
      real(dp), dimension(ncent_tot) :: wt
      real(dp), dimension(3,3) :: dvol
      real(dp), dimension(3,ncent_tot) :: dwt
      real(dp), dimension(3) :: dwtsm
      real(dp), dimension(3,ncent_tot) :: cent_str
      real(dp), dimension(MCHS) :: q_strech
      real(dp), dimension(MCHS) :: efsol
      real(dp), dimension(ncent_tot) :: wt_pcm
      real(dp), parameter :: zero = 0.d0
      real(dp), parameter :: one = 1.d0

      if(.not.allocated(centsav)) allocate(centsav(3, ncent_tot))
      if(.not.allocated(pecentn)) allocate(pecentn(MFORCE))
      if(.not.allocated(xpolsav)) allocate(xpolsav(3,MCHV))

      write(ounit,'(''istrech,alfstr ='',i4,2f10.5)') istrech,alfstr

      do i=1,nforce
        write(ounit,*) '--------------'
        do ic=1,ncent
          write(ounit,'(''center '',i2,'' conf '',i2,'' displace '',3f15.7)') ic,i,(delc(k,ic,i),k=1,3)
        enddo
      enddo
      write(ounit,'(''iwftypes'',20i2)') (iwftype(i),i=1,nforce)

      if(index(mode,'dmc').ne.0) then
        if(nwprod.gt.MFORCE_WT_PRD) call fatal_error('STRETCH: nwprod gt MFORCE_WT_PRD')
        write(ounit,'(''nwprod,itausec='',2i4)') nwprod,itausec
      endif

      do icent=1,ncent
        do k=1,3
          centsav(k,icent)=cent(k,icent)
        enddo
      enddo

!' PCM
      if(ipcm.eq.3) then
        do j=1,nch
          do k=1,3
            xpolsav(k,j)=xpol(k,j)
          enddo
        enddo

! Interatomic forces (and not wave function optimization)
        if(ioptwf.eq.0) then

          open(54,file='field',status='old',form='formatted')
          rewind(54)
          do i=1,nchs
            do j=1,nesph
             read(54,*)
            enddo
          enddo
          do i=1,nchs
            read(54,*) xi,yi,zi,efsol(i)
          enddo
          close(54)
          do k=1,nchs
            enk=0.0d0
            do l=1,ncent
              xx=xpolsav(1,k)-cent(1,l)
              yy=xpolsav(2,k)-cent(2,l)
              zz=xpolsav(3,k)-cent(3,l)
              rr2=xx**2+yy**2+zz**2
              rr3=rr2**1.5d0
              cc1=xx*eta(1,k)
              cc2=yy*eta(2,k)
              cc3=zz*eta(3,k)
              cc=cc1+cc2+cc3
              enk=enk+znuc(iwctype(l))*cc/rr3
            enddo
            env=0.0d0
!           do 70 l=nchs+1,nch
!             xx=xpolsav(1,k)-xpolsav(1,l)
!             yy=xpolsav(2,k)-xpolsav(2,l)
!             zz=xpolsav(3,k)-xpolsav(3,l)
!             rr2=xx**2+yy**2+zz**2
!             rr3=rr2**1.5d0
!             cc1=xx*eta(1,k)
!             cc2=yy*eta(2,k)
!             cc3=zz*eta(3,k)
!             cc=cc1+cc2+cc3
!  70         env=env+ch(l)*cc/rr3
            efsol(k)=efsol(k)+enk+env
          enddo
! endif interatomic forces are being computed
        endif
! endif PCM
      endif

! loop over geometries (if wf optimization, geometries for different adiag are equal)
      do ifl=1,nforce

        do i=1,ncent
            do k=1,3
              cent_str(k,i)=cent(k,i)+delc(k,i,ifl)
            enddo
          call pot_nn(cent_str,znuc,iwctype,ncent,pecentn(ifl),cos_n_sum(1,ifl),sin_n_sum(1,ifl))
        enddo

! PCM
! for wave function optimization or energy calculation, positions/charges unchanged
         if(ipcm.eq.3.and.ioptwf.eq.0) then

            do j=1,nchs
              is=inda(j)
              do k=1,3
               xpol(k,j)=xpolsav(k,j)+delc(k,is,ifl)
              enddo
            enddo

            call sigma_R(efsol,q_strech)

            do i=1,nchs
              sch_s(i,ifl)=2*q_strech(i)
            enddo

            delta_qs=0.d0
            do i=1,nchs
              delta_qs=delta_qs+dabs(q_strech(i)-ch(i))
            enddo
! check deviation of surface charges from charges of primary geometry
            write (6,'(''Geometry'',i4,'' : Deviation of surface charges from primary charges'',1p1d14.5)') ifl,delta_qs
! printout charges if deviation is  big
            if (delta_qs.gt.1.d-3) then
              do i=1,nchs
                write(ounit,'(''Warning: Large deviation in surface charges'',2f16.8)') q_strech(i),ch(i)
              enddo
            endif

            do j=nchs+1,nch
              wtsm=zero
              do icent=1,ncent
                dist2=zero
                do k=1,3
                  dist2=dist2+(xpolsav(k,j)-centsav(k,icent))**2
                enddo
                  dist=dsqrt(dist2)
                  if(istrech.eq.1) wt_pcm(icent)=dexp(-alfstr*dist)
                  if(istrech.eq.2) wt_pcm(icent)=one/dist**alfstr
                  if(istrech.eq.3) wt_pcm(icent)=dexp(alfstr/dist)
                  wtsm=wtsm+wt_pcm(icent)
              enddo
              wtsmi=one/wtsm
              do icent=1,ncent
                wt_pcm(icent)=wt_pcm(icent)*wtsmi
                do k=1,3
                  xpol(k,j)=xpolsav(k,j)+wt_pcm(icent)*delc(k,icent,ifl)
                enddo
              enddo
            enddo

              penups_fc=0.d0
              penupv_fc=0.d0
                do i=1,ncent
                  do j=1,nchs
                    js=inda(j)
                    rnp2=0.d0
                    do k=1,3
                      rnp2=rnp2+(xpol(k,j)-cent_str(k,i))**2.0d0
                    enddo
                    rnp=dsqrt(rnp2)
                    penups_fc=penups_fc+0.5d0*znuc(iwctype(i))*sch_s(j,ifl)/rnp
                  enddo

                do j=nchs+1,nch
                  rnp2=0.d0
                  do k=1,3
                    rnp2=rnp2+(xpol(k,j)-cent_str(k,i))**2.0d0
                  enddo
                  rnp=dsqrt(rnp2)
                  penupv_fc=penupv_fc+0.5d0*znuc(iwctype(i))*ch(j)/rnp
                enddo
                enddo
              delta_gpol_fc=penups_fc-penups+penupv_fc-penupv
              write(ounit,'(''nuclear delta_gpol contribution to force'',i4,'' ='',1p1d14.5)') ifl,delta_gpol_fc
              pecentn(ifl)=pecentn(ifl)+delta_gpol_fc
! endif PCM
          endif

! end loop forces
      enddo

      write(ounit,'(''n-n potential energies '',10f10.5)') (pecentn(ifl),ifl=1,nforce)

      return
      end
end module
