      module mmpol
      use error, only: fatal_error
      interface !LAPACK interface
        SUBROUTINE dgetrf( M, N, A, LDA, IPIV, INFO )
!*  -- LAPACK computational routine --
!*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
          INTEGER            INFO, LDA, M, N
          INTEGER            IPIV( * )
          DOUBLE PRECISION   A( LDA, * )
        END SUBROUTINE
        SUBROUTINE DGETRI( N, A, LDA, IPIV, WORK, LWORK, INFO )
!*  -- LAPACK routine (version 3.1) --
!*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!*     November 2006
          INTEGER            INFO, LDA, LWORK, N
          INTEGER            IPIV( * )
          DOUBLE PRECISION   A( LDA, * ), WORK( * )
        END SUBROUTINE
      end interface
      contains

      subroutine mmpol_extpot_read

c Written by Amovilli-Floris
c Modified by Riccardo and Claudia
c...........................................................
c     read data for MM-POL calculations
c     computes nuclei-charges  interactions (penu_q)
c     computes nuclei-dipoles  interactions (penu_dp)
c     computes charges-charges interactions (peqq)
c     computes charges-dipoles interactions (peq_dp)
c...........................................................
      use mmpol_mod, only: mmpolfile_sites
      use mmpol_cntrl, only: icall_mm, ich_mmpol, immpol, isites_mmpol
      use mmpol_parms, only: chmm, nchmm, x_mmpol
      use mmpol_dipol, only: alfa, dipo
      use mmpol_pot, only: penu_dp, penu_q, peq_dp, peqq, u_self
      use mmpol_inds, only: inds_pol
      use precision_kinds, only: dp
      use contrl_file,    only: ounit

      implicit none

      integer :: i, k
      real(dp) :: PI, penu_mmpol

      data PI/3.1415927D0/

      if(immpol.eq.0) return

c isites_mmpol is always set to zero in read_input
      if(isites_mmpol.eq.0) then

c Read site positions, fixed charges, and polarizabilities
        write(ounit,'(''mmpol open file'',a20)') mmpolfile_sites
        open (55,file=mmpolfile_sites,form='formatted',status='unknown')
        rewind 55
        read(55,*)
        read(55,*) nchmm
        read(55,*)

        do i=1,nchmm
          read(55,*) inds_pol(i),(x_mmpol(k,i),k=1,3),chmm(i),alfa(i),(dipo(k,i),k=1,3)
        enddo
        close(55)
c endif isites_mmpol
      endif

      icall_mm=0
      write(ounit,'(''mmpol nchmm ='',i10)') nchmm

c Compute site-site distances and screening functions, and charge-charge interaction
      call mmpol_compute_peq_q

c Compute nuclei-charge interaction
      call mmpol_compute_penu_q

c Compute dipole-dipole interaction and self-energy
      call mmpol_compute_pedp_dp

c Compute nuclei-dipole interaction
      call mmpol_compute_penu_dp

c Compute nuclei electric field and nuclei-dipole interaction
      call mmpol_field_nu

c Compute charge electric field and charge-dipole interaction
      call mmpol_field_q

c Compute Ainv of mu= Ainv E and electric field due to nuclei and MM charges
      if (ich_mmpol.eq.1) then
         call mmpol_Ainv
      endif

      penu_mmpol=penu_dp+penu_q
      write(ounit,*)
      write(ounit,'(''mmpol epot nuclei-sites dipoles ='',f12.6)') penu_dp
      write(ounit,'(''mmpol epot nuclei-sites charges ='',f12.6)') penu_q
      write(ounit,'(''mmpol epot nuclei-sites ='',f12.6)') penu_mmpol
      write(ounit,'(''mmpol epot charges-charges ='',f12.6)') peqq
      write(ounit,'(''mmpol epot charges-dipoles ='',f12.6)') peq_dp
      write(ounit,'(''mmpol self energy          ='',f12.6)') u_self
      write(ounit,*)

 1000 format(I4,2x,3F12.5,2x,F12.5,2x,F12.5)
      return
      end
c-----------------------------------------------------------------------
      subroutine mmpol_compute_peq_q
c............................................................
c     compute distances and screening between sites and
c     peqq  (charges-charges interaction)
c............................................................
      use mmpol_cntrl, only: immpol
      use mmpol_parms, only: chmm, nchmm, rqq, x_mmpol
      use mmpol_dipol, only: alfa
      use mmpol_pot, only: peqq
      use mmpol_fdc, only: a_cutoff, screen1, screen2
      use mmpol_inds, only: inds_pol
      use precision_kinds, only: dp

      implicit none

      integer :: i, j, k, l
      real(dp) :: dist, r_cutoff, ratio, sixth
      real(dp) :: xx, yy, zz

      if(immpol.eq.0) return

      sixth=1.d0/6.d0
c Compute site-site distances (also charge-charge)
      do k=1,nchmm-1
        do l=k+1,nchmm
          xx=(x_mmpol(1,k)-x_mmpol(1,l))**2.0d0
          yy=(x_mmpol(2,k)-x_mmpol(2,l))**2.0d0
          zz=(x_mmpol(3,k)-x_mmpol(3,l))**2.0d0
          dist=dsqrt(xx+yy+zz)
          rqq(k,l)=dist

          if(inds_pol(k).ne.inds_pol(l)) then
            screen1(k,l)=1.d0
            screen2(k,l)=1.d0
            r_cutoff=a_cutoff*(alfa(k)*alfa(l))**sixth
            ratio=rqq(k,l)/r_cutoff
            if(ratio.lt.1.d0) then
              screen2(k,l)=ratio**4
              screen1(k,l)=4.d0*ratio**3-3.d0*screen2(k,l)
            endif
           else
            screen1(k,l)=0.d0
            screen2(k,l)=0.d0
          endif

          rqq(l,k)=rqq(k,l)
          screen1(l,k)=screen1(k,l)
          screen2(l,k)=screen2(k,l)
        enddo
      enddo

      do k=1,nchmm
       rqq(k,k)=0.d0
      enddo

c Compute charge-charge interaction
      peqq=0.0d0
      do i=1,nchmm-1
        do j=i+1,nchmm
          if(inds_pol(i).ne.inds_pol(j)) then
            peqq=peqq+chmm(i)*chmm(j)/rqq(i,j)
          endif
         enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine mmpol_compute_pedp_dp
c............................................................
c     compute u_self and u_dd (self and dipole-dipole interaction)
c............................................................
      use mmpol_cntrl, only: immpol
      use mmpol_parms, only: nchmm, rqq, x_mmpol
      use mmpol_dipol, only: alfa, dipo
      use mmpol_pot, only: u_dd, u_self
      use mmpol_fdc, only: screen1, screen2
      use mmpol_inds, only: inds_pol
      use precision_kinds, only: dp

      implicit none

      integer :: i, j, k
      real(dp) :: dipo_mod, qiki, qikk, riki, riki3
      real(dp) :: riki5, sik

      if(immpol.eq.0) return

c Riccardo
c Compute self-energy and dipole-dipole interaction
      u_self=0.d0
      do i=1,nchmm
        dipo_mod=0.d0
        do j=1,3
          dipo_mod=dipo_mod+dipo(j,i)**2
        enddo
        u_self=u_self+dipo_mod/alfa(i)
      enddo
      u_self=0.5*u_self

      u_dd=0.d0
      do i=1,nchmm
        do k=1,nchmm
          if(inds_pol(k).ne.inds_pol(i)) then
            qiki=(x_mmpol(1,i)-x_mmpol(1,k))*dipo(1,i)+
     +           (x_mmpol(2,i)-x_mmpol(2,k))*dipo(2,i)+
     +           (x_mmpol(3,i)-x_mmpol(3,k))*dipo(3,i)
            qikk=(x_mmpol(1,i)-x_mmpol(1,k))*dipo(1,k)+
     +           (x_mmpol(2,i)-x_mmpol(2,k))*dipo(2,k)+
     +           (x_mmpol(3,i)-x_mmpol(3,k))*dipo(3,k)
            sik=dipo(1,i)*dipo(1,k)+
     +          dipo(2,i)*dipo(2,k)+
     +          dipo(3,i)*dipo(3,k)
              riki=1.d0/rqq(i,k)
              riki3=riki**3
              riki5=riki3*riki*riki
              u_dd=u_dd+0.5*(screen1(i,k)*sik*riki3-3.d0*screen2(i,k)*qiki*qikk*riki5)
          endif
        enddo
      enddo
c Riccardo

      return
      end
c-----------------------------------------------------------------------
      subroutine mmpol_compute_penu_dp
c............................................................
c     compute penu_dp (nuclei-induced dipoles on MM sites interaction)
c............................................................
      use atom, only: znuc, cent, iwctype, ncent
      use mmpol_cntrl, only: immpol
      use mmpol_parms, only: nchmm, x_mmpol
      use mmpol_dipol, only: dipo
      use mmpol_pot, only: penu_dp
      use precision_kinds, only: dp
      use contrl_file,    only: ounit
      implicit none

      integer :: i, j
      real(dp) :: dpdr, rnp, rnp2, rnp3, xx
      real(dp) :: yy, zz

      if(immpol.eq.0) return

      penu_dp=0.0d0
      do i=1,ncent
        do j=1,nchmm
          xx=(x_mmpol(1,j)-cent(1,i))
          yy=(x_mmpol(2,j)-cent(2,i))
          zz=(x_mmpol(3,j)-cent(3,i))
          rnp2=xx**2+yy**2+zz**2
          rnp=dsqrt(rnp2)
          rnp3=rnp**3.0d0
          dpdr=xx*dipo(1,j)+ yy*dipo(2,j)+ zz*dipo(3,j)
          penu_dp=penu_dp-znuc(iwctype(i))*dpdr/rnp3
        enddo
      enddo

      write(ounit,*) 'Compute dipoles penu_dp =' ,penu_dp

      return
      end
c-----------------------------------------------------------------------
      subroutine mmpol_compute_penu_q
c............................................................
c     compute penu_q  (nuclei-charges interaction)
c............................................................
      use atom, only: znuc, cent, iwctype, ncent
      use mmpol_cntrl, only: immpol
      use mmpol_parms, only: chmm, nchmm, x_mmpol
      use mmpol_pot, only: penu_q
      use precision_kinds, only: dp

      implicit none

      integer :: i, j
      real(dp) :: rnp, rnp2, xx, yy, zz

      if(immpol.eq.0) return

      penu_q=0.0d0
      do i=1,ncent
        do j=1,nchmm
          xx=(x_mmpol(1,j)-cent(1,i))**2.0d0
          yy=(x_mmpol(2,j)-cent(2,i))**2.0d0
          zz=(x_mmpol(3,j)-cent(3,i))**2.0d0
          rnp2=xx+yy+zz
          rnp=dsqrt(rnp2)
          penu_q=penu_q+znuc(iwctype(i))*chmm(j)/rnp
         enddo
       enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine mmpol_Ainv
C     ***************************************************************
C     contribution inverse A matrix to compute dipoles as mu=Ainv E
C     ***************************************************************
      use mmpol_mod, only: MCHMM
      use mmpol_parms, only: nchmm, rqq, x_mmpol
      use mmpol_dipol, only: alfa
      use mmpol_ahpol, only: ah_pol
      use mmpol_fdc, only: screen1, screen2
      use mmpol_inds, only: inds_pol
      use precision_kinds, only: dp

      implicit none

      integer :: i, ij, j, j1, j2
      integer :: j3, k, ki, l
      integer :: lj, nch3
      real(dp) :: delta_ij
      real(dp) :: det, dxx_kilj, tmat_kilj
      real(dp), dimension(3*MCHMM*3*MCHMM) :: ahpol_vec

C
c............................................................
c     The matrix ah_pol is computed and inverted
c     (3*nchmm)*(3*nchmm)
c............................................................

      do k=1,nchmm
        do i=1,3
          ki=(k-1)*3+i
          do l=1,nchmm
            do j=1,3
              lj=(l-1)*3+j
              ah_pol(ki,lj)=0.d0
              if(inds_pol(k).ne.inds_pol(l)) then
                 delta_ij=0.d0
                 if (i.eq.j) delta_ij=1.d0
                 dxx_kilj=(x_mmpol(i,k)-x_mmpol(i,l))*(x_mmpol(j,k)-x_mmpol(j,l))
                 tmat_kilj=3.d0*screen2(k,l)*dxx_kilj/rqq(k,l)**5-delta_ij*screen1(k,l)/rqq(k,l)**3
                 ah_pol(ki,lj)=-alfa(k)*tmat_kilj
              endif
            enddo
          enddo
          ah_pol(ki,ki)=1.d0
        enddo
      enddo

      nch3=3*nchmm
      do i=1,nch3
        do j=1,nch3
          ij=(j-1)*nch3+i
          ahpol_vec(ij)=ah_pol(i,j)
        enddo
      enddo

      call mmpol_matinv(ahpol_vec,nch3,det)

      do i=1,nch3
        do j=1,nch3
          ij=(j-1)*nch3+i
          ah_pol(i,j)=ahpol_vec(ij)
        enddo
      enddo

      do i=1,nch3
        do j=1,nchmm
          j3=3*J
          J2=J3-1
          J1=J3-2
          ah_pol(i,j1)=ah_pol(i,j1)*alfa(j)
          ah_pol(i,j2)=ah_pol(i,j2)*alfa(j)
          ah_pol(i,j3)=ah_pol(i,j3)*alfa(j)
        enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine mmpol_field_nu
C     ***************************************************************
C     contribution from nuclei to polarization charghes
C     ***************************************************************
      use atom, only: znuc, cent, iwctype, ncent
      use mmpol_cntrl, only: immpol
      use mmpol_parms, only: nchmm, x_mmpol
      use mmpol_dipol, only: dipo
      use mmpol_field, only: enk_pol
      use mmpol_pot, only: penu_dp
      use precision_kinds, only: dp
      implicit none

      integer :: i, k, l
      real(dp) :: rr2, rr3, xx, yy, zz





      real*8 ah_vec,det
C

C
c
      if(immpol.eq.0) return
c..............................................................
c     enk_pol=nuclear field on MM sites
c..............................................................
      do k=1,nchmm
        enk_pol(1,k)=0.0d0
        enk_pol(2,k)=0.0d0
        enk_pol(3,k)=0.0d0
        do l=1,ncent
          xx=x_mmpol(1,k)-cent(1,l)
          yy=x_mmpol(2,k)-cent(2,l)
          zz=x_mmpol(3,k)-cent(3,l)
          rr2=xx**2+yy**2+zz**2
          rr3=rr2**1.5d0
          enk_pol(1,k)=enk_pol(1,k)+znuc(iwctype(l))*xx/rr3
          enk_pol(2,k)=enk_pol(2,k)+znuc(iwctype(l))*yy/rr3
          enk_pol(3,k)=enk_pol(3,k)+znuc(iwctype(l))*zz/rr3
        enddo
      enddo

c Compute nuclei-dipole interaction
      penu_dp=0.0d0
      do k=1,nchmm
        do i=1,3
          penu_dp=penu_dp-dipo(i,k)*enk_pol(i,k)
        enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine mmpol_field_q
c............................................................
c    computes electric field due to fixed MM charges on MM sites
c    and potential interaction MM charges-dipoles (peq_dp)
c............................................................
      use mmpol_cntrl, only: immpol
      use mmpol_parms, only: chmm, nchmm, rqq, x_mmpol
      use mmpol_dipol, only: dipo
      use mmpol_field, only: eqk_pol
      use mmpol_pot, only: peq_dp
      use mmpol_fdc, only: screen1
      use mmpol_inds, only: inds_pol
      use contrl_file,    only: ounit
      use precision_kinds, only: dp
      implicit none

      integer :: i, k, l
      real(dp) :: rr3, xx, yy, zz









      if(immpol.eq.0) return

c Compute MM charge electric field
      do k=1,nchmm
        eqk_pol(1,k)=0.0d0
        eqk_pol(2,k)=0.0d0
        eqk_pol(3,k)=0.0d0
        do l=1,nchmm
          if(inds_pol(k).ne.inds_pol(l)) then
            xx=x_mmpol(1,k)-x_mmpol(1,l)
            yy=x_mmpol(2,k)-x_mmpol(2,l)
            zz=x_mmpol(3,k)-x_mmpol(3,l)
            rr3=rqq(k,l)**3
            eqk_pol(1,k)=eqk_pol(1,k)+screen1(k,l)*chmm(l)*xx/rr3
            eqk_pol(2,k)=eqk_pol(2,k)+screen1(k,l)*chmm(l)*yy/rr3
            eqk_pol(3,k)=eqk_pol(3,k)+screen1(k,l)*chmm(l)*zz/rr3
          endif
        enddo
      enddo

c Compute MM charge-dipole interaction
      peq_dp=0.0d0
      do k=1,nchmm
        do i=1,3
          peq_dp=peq_dp-dipo(i,k)*eqk_pol(i,k)
        enddo
      enddo
      write(ounit,*) '(charges-dipoles)=',peq_dp

      return
      end
c-----------------------------------------------------------------------
      subroutine mmpol_efield(nelec,coord)
C     ***************************************************************
c     For the accepted configuration, the electronic field is computed
c     on MM sites  eek_pol(1,k),eek_pol(2,k),eek_pol(3,k)
C     ***************************************************************
      use mmpol_hpsi, only: eek_pol
      use mmpol_cntrl, only: immpol
      use mmpol_parms, only: nchmm, x_mmpol
      use mmpol_fdc, only: rcolm
      use precision_kinds, only: dp
      implicit none

      integer :: i, k, nelec
      real(dp) :: AV, GC, PI, hatokc, rcolm3
      real(dp) :: rr1, rr2, rr3, rr3i
      real(dp) :: xx, yy, zz
      real(dp), dimension(3,*) :: coord







      DATA PI/3.1415927D0/,GC/1.9872159D0/,AV/0.60228D0/
      data hatokc/627.509541d0/

      if (immpol.eq.3) return

      rcolm3=rcolm**3
      do k=1,nchmm
        eek_pol(1,k)=0
        eek_pol(2,k)=0
        eek_pol(3,k)=0
        do i=1,nelec
          xx=x_mmpol(1,k)-coord(1,i)
          yy=x_mmpol(2,k)-coord(2,i)
          zz=x_mmpol(3,k)-coord(3,i)
          rr2=xx**2+yy**2+zz**2
          rr1=dsqrt(rr2)
          rr3=rr1*rr2
c..............................................................
c  correction for collision between e- and fixed charges
c..............................................................
          if (rr1.lt.rcolm) rr3=rcolm3
          rr3i=1.d0/rr3
          eek_pol(1,k)=eek_pol(1,k)-xx*rr3i
          eek_pol(2,k)=eek_pol(2,k)-yy*rr3i
          eek_pol(3,k)=eek_pol(3,k)-zz*rr3i
        enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine mmpol_extpot_ene(coord,nelec,peQMdp,peQMq)

c Written by Amovilli-Floris
c Modified by Riccardo
c......................................................
c       Calculate e-/dipoles interactions (MM_POL)
c                 e-/charges interactions (MM_POL)
c       and adds nuclei-dipoles interactions and
c                nuclei-charges interactions and
c......................................................
c     H(environment) is defined.
c     It includes both H(QM/MM) and H(MM).
c     peQMdp and peQMq are therefore mixed terms
c......................................................
      use mmpol_cntrl, only: icall_mm
      use mmpol_parms, only: nchmm
      use mmpol_pot, only: penu_dp, penu_q, pepol_dp, pepol_q, peq_dp, peqq, u_dd, u_self
      use precision_kinds, only: dp
      use contrl_file,    only: ounit
      implicit none

      integer :: i, nelec
      real(dp) :: AV, GC, PI, hatokc, peQMdp
      real(dp) :: peQMq
      real(dp), dimension(3,*) :: coord




      DATA PI/3.1415927D0/,GC/1.9872159D0/,AV/0.60228D0/
      data hatokc/627.509541d0/
      icall_mm=icall_mm+1

      peQMdp=penu_dp+peq_dp+u_self+u_dd
      peQMq =penu_q+peqq
      do i=1,nelec
        call mmpol_extpot_ene_elec(coord(1,i))
        peQMq=peQMq+pepol_q
        peQMdp=peQMdp+pepol_dp
      enddo

      if (icall_mm.eq.1)then
        write(ounit,*)
        write(ounit,*)'nchmm=',nchmm
        write(ounit,*)'u_self= ',u_self
        write(ounit,*)'u_dd = ',u_dd
        write(ounit,'(''QM-MM epot  ='',f12.6)') peQMdp+peQMq
        write(ounit,'(''QM-MM e-/dipoles  ='',f12.6)') pepol_dp
        write(ounit,'(''QM-MM e-/charges  ='',f12.6)') pepol_q
        write(ounit,*)
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine mmpol_extpot_ene_elec(x)

c Written by Amovilli-Floris
c......................................................
c       Calculate e-/dipoles interactions
c......................................................
      use mmpol_parms, only: chmm, nchmm, x_mmpol
      use mmpol_dipol, only: dipo
      use mmpol_pot, only: pepol_dp, pepol_q
      use mmpol_fdc, only: rcolm
      use precision_kinds, only: dp
      implicit none

      integer :: j
      real(dp) :: AV, GC, PI, dpdr, rcolm3
      real(dp) :: repol, repol2, repol3, xx
      real(dp) :: yy, zz
      real(dp), dimension(3) :: x





      DATA PI/3.1415927D0/,GC/1.9872159D0/,AV/0.60228D0/

      rcolm3=rcolm**3

      pepol_dp=0.0d0
      pepol_q=0.0d0
      do j=1,nchmm
        xx=(x_mmpol(1,j)-x(1))
        yy=(x_mmpol(2,j)-x(2))
        zz=(x_mmpol(3,j)-x(3))
        repol2=xx**2+yy**2+zz**2
        repol=dsqrt(repol2)
        repol3=repol*repol2
        dpdr=xx*dipo(1,j)+ yy*dipo(2,j)+ zz*dipo(3,j)
c......................................................
c    corrections for collisions electrons-dipoles
c    corrections for collisions electrons-charges
c......................................................
        if(repol.lt.rcolm) then
          repol=rcolm
          repol3=rcolm3
        endif

c......................................................
c     interaction with dipoles
c......................................................
        pepol_dp=pepol_dp+dpdr/repol3

c......................................................
c     interaction with fixed charges
c......................................................
        pepol_q=pepol_q-chmm(j)/repol
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine mmpol_dipoles(eek_ave,eek_err)
      use mmpol_mod, only: MCHMM
      use mmpol_parms, only: chmm, nchmm, x_mmpol
      use mmpol_dipol, only: alfa, dipo
      use mmpol_ahpol, only: ah_pol, bh_pol
      use mmpol_field, only: enk_pol, eqk_pol
      use mmpol_pot, only: u_dd, u_self
      use mmpol_inds, only: inds_pol

      use precision_kinds, only: dp
      implicit none

      integer :: i, j, k, k1, k2
      integer :: k3, nch3
      real(dp) :: efree1, efree2, efree3, efree_k1, efree_k2
      real(dp) :: efree_k3, u_ef, u_pol
      real(dp), dimension(3,MCHMM) :: eek_ave
      real(dp), dimension(3,MCHMM) :: eek_err
      real(dp), dimension(3*MCHMM) :: dp_array

      open (56,file='dipo_new',form='formatted',status='unknown')
      open (57,file='efield_electron',form='formatted',status='unknown')
      rewind 56
      rewind 57
c Riccardo
      open (58,file='E0',form='formatted',status='unknown')
      rewind 58
c Riccardo

      u_ef=0.d0
      do k=1,nchmm
        k3=3*k
        k2=k3-1
        k1=k3-2
        efree_k1=enk_pol(1,k)+eqk_pol(1,k)+eek_ave(1,k)
        efree_k2=enk_pol(2,k)+eqk_pol(2,k)+eek_ave(2,k)
        efree_k3=enk_pol(3,k)+eqk_pol(3,k)+eek_ave(3,k)

        bh_pol(k1)=efree_k1
        bh_pol(k2)=efree_k2
        bh_pol(k3)=efree_k3

        u_ef=u_ef-efree_k1*dipo(1,k)-efree_k2*dipo(2,k)-efree_k3*dipo(3,k)
      enddo

      write (6,*)
      write (6,*) 'INITIAL DIPOLES'
      write (6,*)
c u_ef = -mu*E + U_self + U_dd
      u_pol=u_dd+u_self+u_ef
      write (6,*) 'MMpol:  -mu*E    = ',u_ef
      write (6,*) 'MMpol:  U_self   = ',u_self
      write (6,*) 'MMpol:  U_dd     = ',u_dd
      write (6,*) 'MMpol:  U_pol    = ',u_pol
      write (6,*) 'MMpol: -0.5mu*E  = ',u_ef/2.d0
      write (6,*)

      nch3=3*nchmm
      do i=1,nch3
        dp_array(i)=0.0d0
        do j=1,nch3
          dp_array(i)=dp_array(i)+ah_pol(i,j)*bh_pol(j)
        enddo
      enddo
c...................................................................
c     write in 'efield' the average of solute electrons field
c     write in 'E0' the total static field
c     overwrite dipoles
c...................................................................
      write(56,*)'nchmm'
      write(56,*) nchmm
      write(56,*)'inds, coordinates,chmm,alfa,dipoles'

      do k=1,nchmm
        k3=3*k
        k2=k3-1
        k1=k3-2
        dipo(1,k)=dp_array(k1)
        dipo(2,k)=dp_array(k2)
        dipo(3,k)=dp_array(k3)
        write(56,1001) inds_pol(k),(x_mmpol(i,k),i=1,3),chmm(k),alfa(k),(dipo(i,k),i=1,3)
        write(57,1000) (x_mmpol(i,k),i=1,3),(eek_ave(i,k),i=1,3),(eek_err(i,k),i=1,3)
c Riccardo write E0
        write(58,1000) (x_mmpol(i,k),i=1,3),(enk_pol(i,k)+eqk_pol(i,k)+eek_ave(i,k),i=1,3)
c Riccardo
      enddo

      write (6,*) 'FINAL DIPOLES'
      write (6,*)

      call mmpol_compute_pedp_dp

      u_ef=0.d0
      do i=1,nchmm
        k3=3*i
        k2=k3-1
        k1=k3-2
        efree1=bh_pol(k1)
        efree2=bh_pol(k2)
        efree3=bh_pol(k3)
        u_ef=u_ef-efree1*dipo(1,i)-efree2*dipo(2,i)-efree3*dipo(3,i)
      enddo

      u_pol=u_dd+u_self+u_ef
      write (6,*) 'MMpol:  -mu*E    = ',u_ef
      write (6,*) 'MMpol:  U_self   = ',u_self
      write (6,*) 'MMpol:  U_dd     = ',u_dd
      write (6,*) 'MMpol:  U_pol    = ',u_pol
      write (6,*) 'MMpol: -0.5mu*E  = ',u_ef/2.d0

 1000 format(9F12.5)
 1001 format(I7,9F12.5)
 1200 format(I4,9F12.5)

      return
      end
c-----------------------------------------------------------------------
c    AVERAGES   subroutines
c-----------------------------------------------------------------------
      subroutine mmpol_init(iflag)

      use mmpol_cntrl, only: immpol
      use mmpol_parms, only: nchmm
      use mmpol_averages, only: cmmpol_cm2, cmmpol_cum, cmmpol_sum, dmmpol_cm2, dmmpol_cum, dmmpol_sum
      use mmpol_averages, only: eek1_cm2, eek1_cum, eek2_cm2, eek2_cum, eek3_cm2, eek3_cum, eek_sum

      implicit none

      integer :: i, iflag



      if(immpol.eq.0) return

      dmmpol_sum=0
      cmmpol_sum=0

      do i=1,nchmm
        eek_sum(1,i)=0.0d0
        eek_sum(2,i)=0.0d0
        eek_sum(3,i)=0.0d0
      enddo

      if(iflag.gt.0) return
      dmmpol_cum=0
      cmmpol_cum=0

      dmmpol_cm2=0
      cmmpol_cm2=0

      do i=1,nchmm
        eek1_cum(i)=0.0d0
        eek2_cum(i)=0.0d0
        eek3_cum(i)=0.0d0
        eek1_cm2(i)=0.0d0
        eek2_cm2(i)=0.0d0
        eek3_cm2(i)=0.0d0
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine mmpol_dump(iu)

      use mmpol_cntrl, only: immpol
      use mmpol_parms, only: nchmm
      use mmpol_averages, only: cmmpol_cm2, cmmpol_cum, dmmpol_cm2, dmmpol_cum
      use mmpol_averages, only: eek1_cm2, eek1_cum, eek2_cm2, eek2_cum, eek3_cm2, eek3_cum

      implicit none

      integer :: i, iu



      if(immpol.eq.0) return
      write(iu) dmmpol_cum,dmmpol_cm2
      write(iu) cmmpol_cum,cmmpol_cm2

      do i=1,nchmm
        write(iu) eek1_cum(i),eek1_cm2(i)
        write(iu) eek2_cum(i),eek2_cm2(i)
        write(iu) eek3_cum(i),eek3_cm2(i)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine mmpol_rstrt(iu)

      use mmpol_cntrl, only: immpol
      use mmpol_parms, only: nchmm
      use mmpol_averages, only: cmmpol_cm2, cmmpol_cum, dmmpol_cm2, dmmpol_cum
      use mmpol_averages, only: eek1_cm2, eek1_cum, eek2_cm2, eek2_cum, eek3_cm2, eek3_cum

      implicit none

      integer :: i, iu



      if(immpol.eq.0) return
      read(iu) dmmpol_cum,dmmpol_cm2
      read(iu) cmmpol_cum,cmmpol_cm2

      do i=1,nchmm
        read(iu) eek1_cum(i),eek1_cm2(i)
        read(iu) eek2_cum(i),eek2_cm2(i)
        read(iu) eek3_cum(i),eek3_cm2(i)
      enddo

      return
      end
c......................................................
c     end AVERAGES subroutines
c......................................................
      subroutine mmpol_matinv(a,nsub,determinant)
      use mmpol_mod, only: MCHMM
      use precision_kinds, only: dp
      use contrl_file,    only: ounit
      implicit none

      integer :: i, info, nsub
      integer, dimension(3*MCHMM) :: ipvt
      real(dp) :: determinant, ten
      real(dp), dimension(nsub,nsub) :: a
      real(dp), dimension(3*MCHMM) :: work
      real(dp), dimension(2) :: det

c routine to calculate inverse and determinant of matrix a
c assumed to be dimensioned a(nsub,nsub).
c the matrix a is replaced by its inverse.


      call dgetrf(nsub,nsub,a,nsub,ipvt,info)
      if(info.gt.0) then
            write(ounit,'(''mmpol.f'')')
        write(ounit,'(''MATINV: u(k,k)=0 with k= '',i5)') info
        call fatal_error('MATINV: info ne 0 in dgetrf')
      endif

      det(1) = 1.0d0
      det(2) = 0.0d0
      ten = 10.0d0
      do i = 1, nsub
        if (ipvt(i) .ne. i) det(1) = -det(1)
        det(1) = a(i,i)*det(1)
c        ...exit
        if (det(1) .eq. 0.0d0) go to 60
   10   if (dabs(det(1)) .ge. 1.0d0) go to 20
        det(1) = ten*det(1)
        det(2) = det(2) - 1.0d0
        go to 10
   20   continue
   30   if (dabs(det(1)) .lt. ten) go to 40
          det(1) = det(1)/ten
          det(2) = det(2) + 1.0d0
        go to 30
   40   continue
      enddo
   60 continue

      determinant = det(1)*10.0**det(2)

      call dgetri(nsub,a,nsub,ipvt,work,3*MCHMM,info)

      return

      end
      end module
