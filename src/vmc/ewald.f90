      module ewald_breakup
      use error, only: fatal_error
      interface                 !LAPACK interface
      SUBROUTINE DPOSV( UPLO, N, NRHS, A, LDA, B, LDB, INFO )
!     *  -- LAPACK computational routine --
!     *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!     *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      CHARACTER          UPLO
      INTEGER            INFO, LDA, LDB, N, NRHS
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * )

      END SUBROUTINE

      end interface

      contains
      subroutine set_ewald
! Written by Cyrus Umrigar
! Modified to use periodic GTO's by Edgar Josue Landinez Borda

      use system, only: znuc, iwctype, ncent
      use multiple_geo, only: pecent
      use control, only: ipr
      use contrl_file, only: ounit
      use contrl_per, only: iperiodic
      use ewald, only: b_coul, b_jas, y_coul, y_jas
      use ewald, only: allocate_ewald, deallocate_ewald
      use ewald_test, only: f, vbare_coul, vbare_jas
      use ewald_mod, only: NG1DX, NGVEC_BIGX
      use periodic, only: cutg, cutg_big, cutr, glatt
      use periodic, only: glatt_inv, gnorm, gvec, igmult, igvec
      use periodic, only: ncoef_per, ng1d, ngnorm, ngnorm_big
      use periodic, only: ngvec, ngvec_big
      use periodic, only: np_coul, np_jas, npoly, vcell
      use periodic, only: rlatt, rlatt_inv
      use periodic, only: znuc2_sum, znuc_sum
      use periodic, only : n_images, ell
      use error, only: fatal_error
      use find_pimage, only: check_lattice
      use matinv_mod, only: matinv

      use precision_kinds, only: dp

      implicit none

      integer :: i, ict, ifcon, ig, in
      integer :: ir, j, k, lowest_pow
      integer :: npts
      real(dp) :: pi, twopi
      real(dp) :: b0, chisq
      real(dp) :: datan, det, det1
      real(dp) :: dist_min, dx, gdistmin
      real(dp) :: rms, rr, testv, test_s
      real(dp) :: true, true_s, vgcell
      real(dp), dimension(3) :: rdist
      real(dp), dimension(3) :: gdist
      real(dp), dimension(3) :: r_tmp
      real(dp), parameter :: eps = 1.d-12

! for images evaluation
      integer :: ix,iy,iz
      integer :: nix,niy,niz
      integer :: imcount
      integer :: nisum, ngmx

      pi=4.d0*datan(1.d0)
      twopi=2*pi

! Number of coefficients polynomial expansion in the breakup
      ncoef_per=npoly+1

! Check that the lattice vectors are the smallest possible ones and return the smallest
! which is used to set the range of the real-space Ewald sums so that only one image
! of a nucleus or an electron is present within cutr and cutr_sim respectively.
      call check_lattice(rlatt,cutr,0)

      write(ounit, *)  "================================================================================="
      write(ounit, *)  "==  Setting-up Ewald-Breakup                                                   =="
      write(ounit, *)  "================================================================================="
      write(ounit, *)  " "

      write(ounit,'(''cutr ='',9f9.5)') cutr

! Calculate inverse transformations (from lattice coordinates to real coordinates)
! and cell volumes
      do i=1,3
        do k=1,3
          rlatt_inv(k,i)=rlatt(k,i)
        enddo
      enddo
      call matinv(rlatt_inv,3,det)

      det1=twopi/det
      glatt(1,1)=det1*(rlatt(2,2)*rlatt(3,3)-rlatt(2,3)*rlatt(3,2))
      glatt(2,1)=det1*(rlatt(3,2)*rlatt(1,3)-rlatt(3,3)*rlatt(1,2))
      glatt(3,1)=det1*(rlatt(1,2)*rlatt(2,3)-rlatt(1,3)*rlatt(2,2))
      glatt(1,2)=det1*(rlatt(2,3)*rlatt(3,1)-rlatt(2,1)*rlatt(3,3))
      glatt(2,2)=det1*(rlatt(3,3)*rlatt(1,1)-rlatt(3,1)*rlatt(1,3))
      glatt(3,2)=det1*(rlatt(1,3)*rlatt(2,1)-rlatt(1,1)*rlatt(2,3))
      glatt(1,3)=det1*(rlatt(2,1)*rlatt(3,2)-rlatt(2,2)*rlatt(3,1))
      glatt(2,3)=det1*(rlatt(3,1)*rlatt(1,2)-rlatt(3,2)*rlatt(1,1))
      glatt(3,3)=det1*(rlatt(1,1)*rlatt(2,2)-rlatt(1,2)*rlatt(2,1))

      write(ounit,'(/,''Reciprocal lattice basis vectors'',3(/,3f10.6))') &
       ((glatt(k,j),k=1,3),j=1,3)

      vcell=dabs(det)
      write(ounit,'(/,''Cell volume'',f14.8)') vcell

! Calculate inverse transformation for reciprocal lattice (from lattice coordinates to real coordinates)
! Needed to transform k-vectors
      do i=1,3
        do k=1,3
          glatt_inv(k,i)=glatt(k,i)
        enddo
      enddo
      call matinv(glatt_inv,3,det)

! real-space distances
! primitive cell
      call short_distance(rlatt,vcell,dist_min,rdist)

! reciprocal-space distances
! primitive cell
      vgcell=twopi**3/vcell
      call short_distance(glatt,vgcell,gdistmin,gdist)

! Estimate maximum shells number before build it

      ngvec=1
      do k=1,3
         ng1d(k)=int(cutg/gdist(k))
         ngvec=ngvec*(2*ng1d(k)+1)
         write(ounit,*)"Estimated ng1d given cutg",cutg, k, 2*ng1d(k)+1
      enddo
      write(ounit,*) "Estimated max ngvec", ngvec

      if(ngvec.gt.NGVEC_BIGX) then
         write(ounit,*) "Compute a new cutg to lower number of g-vectors"
         ngmx=7
         cutg=gdistmin*ngmx
         write(ounit,*) "New cutg value", cutg
      endif

! generate shells of primitive cell g-vectors
      call shells(cutg_big,glatt,gdist,igvec,gvec,gnorm,igmult,ngvec_big, &
      ngnorm_big,ng1d)

      ngnorm=ngnorm_big
      ngvec=0
      do k=1,ngnorm_big
        if(gnorm(k).gt.cutg+eps) then
          ngnorm=k-1
          goto 20
        endif
        ngvec=ngvec+igmult(k)
      enddo
   20 write(ounit,*) "Recomputed after shells generation ngvec", ngvec

      call allocate_ewald()

      if(ipr.ge.4) then
        write(ounit,'(/,''Shells within cutg_big,cutg'',2i8)') ngnorm_big,ngnorm
        write(ounit,'(/,''Vects. within cutg_big,cutg'',2i8)') ngvec_big,ngvec
        write(ounit,'(/,''ng1d for primitive cell'',3i4)') (ng1d(k),k=1,3)
      endif

      do k=1,3
        if(ng1d(k).gt.NG1DX) then
          if(ipr.ge.4) write(ounit,'(''k,ng1d(k),NG1DX='',i1,2i8)') k,ng1d(k),NG1DX
          call fatal_error ('ng1d(k)>NG1DX in set_ewald')
        endif
      enddo

      !open(1,file='gvectors_qmc')
      !write(1,'(i5,'' ngvec (half of them only)'')') ngvec
      !write(1,'(3i5,3f8.4)') ((igvec(k,i),k=1,3),(gvec(k,i),k=1,3),i=1,ngvec)
      !close(1)

! Coulomb interactions in primitive and simulation cells
      lowest_pow=-1
      b0=1.d0
      ifcon=1

! n-n, e-n interactions (primitive cell)
! put in uniform background by setting k=0 term to zero
      vbare_coul(1)=0.d0
      do k=2,ngnorm_big
!     Fourier transfom of 1/r
         vbare_coul(k)=2*twopi/(vcell*gnorm(k)**2)
      enddo

      if(ipr.ge.4) then
         write(ounit,'(''Compare Coulomb reconstructed'')')
         write(ounit,'(''1/r  vs vps_rest  diff'')')
         write(ounit,'(''start vbare_coul separate primitive'')')
         write(ounit,'(''check vbare_coul = '',20d12.4)') (vbare_coul(k),k=1,ngnorm)
      endif

      ict=1
      call separate(vbare_coul,b0,lowest_pow,ngnorm_big,igmult,gnorm, &
       ngnorm,cutr,vcell,ncoef_per,np_coul,b_coul,y_coul,chisq,ifcon)
      if(ipr.ge.4) write(ounit,'(''finish vbare_coul separate primitive'')')

      if(chisq.gt.0) then
        write(ounit,'(''Rms error in 1/r separation in primitive cell'',d12.5)') dsqrt(chisq)
       else
        write(ounit,'(''Warning: Rms error missing, chisq negative in 1/r primitive separate'',d12.4)') chisq
        if(chisq.lt.0.d0) call fatal_error ('chisq<0 in separate')
      endif

      if(ipr.ge.4) then
         write(ounit,'(/,''Separation of Coulomb interaction in primitive cell'')')
         write(ounit,'(''vbare_coul = '',20d12.4)') (vbare_coul(k),k=1,ngnorm)
         write(ounit,'(''vbare_coul = '',20d12.4)') (vbare_coul(k),k=1,ngnorm_big)
         write(ounit,'(''y_coul = '',20d12.4)') (y_coul(k),k=1,ngnorm)
         write(ounit,'(''b_coul = '',20d12.4)') (b_coul(k),k=1,ncoef_per)
      endif

! e-e interactions (simulation cell) (we can reuse vbare_coul)
! put in uniform background by setting k=0 term to zero

! TMP
      f=1.d0
! TMP
      vbare_jas(1)=0.d0
      do k=2,ngnorm_big
! Fourier transform of -1/r*(1-exp(-r/f)) for Jastrow
        vbare_jas(k)=-vbare_coul(k)/(1+(f*gnorm(k))**2)
      enddo

!     e-e Jastrow
      lowest_pow=0
      b0=0.5d0/f**2
      ifcon=1

      if(ipr.ge.4) write(ounit,'(''Start vbare_jas separate'')')
      call separate(vbare_jas,b0,lowest_pow,ngnorm_big,igmult, &
           gnorm, ngnorm, cutr,vcell,ncoef_per,np_jas, &
           b_jas,y_jas,chisq,ifcon)
      if(ipr.ge.4) write(ounit,'(''End vbare_jas separate'')')

      if(chisq.gt.0) then
        write(ounit,'(''Rms error in Jastrow separation'',d12.5)') dsqrt(chisq)
       else
        write(ounit,'(''Warning: Rms error missing, chisq negative in Jastrow separate'',d12.4)') chisq
        if(chisq.lt.0.d0) call fatal_error ('chisq<0 in separate')
      endif

      if(ipr.ge.4) then
         write(ounit,'(/,''Separation of Jastrow in simulation cell'')')
         write(ounit,'(''vbare_jas = '',20d12.4)') (vbare_jas(k),k=1,ngnorm)
         write(ounit,'(''vbare_jas = '',20d12.4)') (vbare_jas(k),k=1,ngnorm_big)
         write(ounit,'(''y_jas = '',20d12.4)') (y_jas(k),k=1,ngnorm)
      !  write(ounit,'(''b_jas = '',20d12.4)') (b_jas(k),k=1,ncoef_per)
      endif

      write(ounit,'(''b_jas = '',20d12.4)') (b_jas(k),k=1,ncoef_per)

! debug e-e Jastrow
! Since Jastrow has singlularity at 0, cannot match there, so evaluate
! rms error only for latter 3/4 of interval
      if(ipr.ge.4) then
         write(ounit,'(/,''Debugging e-e Jastrow separation'')')
         write(ounit,'(''      r       "true"       test      test-true -1/r*(1-exp(-r/f) d_true d_test'')')
         lowest_pow=0
         npts=101
         dx=cutr/(npts-1)
         rms=0.d0
         do i=1,npts
            r_tmp(1)=(i-1)*dx+1.d-20
            r_tmp(2)=0.d0
            r_tmp(3)=0.d0
            rr=sqrt(r_tmp(1)**2+r_tmp(2)**2+r_tmp(3)**2)
            testv=vsrange(rr,cutr,lowest_pow,ncoef_per,np_jas,b_jas)
            testv=testv+vlrange_r(r_tmp,gvec,ngnorm,igmult,y_jas)
            true=vlrange_r(r_tmp,gvec,ngnorm_big,igmult,vbare_jas)
            if(4*i.ge.npts) rms=rms+(true-testv)**2
            if(i.eq.1) then
               write(ounit,'(''jas'',f8.4,4f12.6,2f8.3)') rr,true,testv,testv-true
            else
               write(ounit,'(''jas'',f8.4,4f12.6,2f8.3)') rr,true,testv,testv-true &
                    ,-1/rr*(1-exp(-rr/f)),(true-true_s)/dx,(testv-test_s)/dx
            endif
            true_s=true
            test_s=testv
         enddo
         rms=sqrt(4*rms/(3*npts))
         write(ounit,'(''Rms error of jas fit on larger 3/4 interval='',d12.4)') rms
         if(rms.gt.1.d-3) write(ounit,'(''Warning: rms error of jas fit too large'',d12.4)') rms
         write(ounit,'(/,''End debugging e-e Jastrow separation'')')
      endif

      znuc_sum=0
      znuc2_sum=0
      do i=1,ncent
        znuc_sum=znuc_sum+znuc(iwctype(i))
        znuc2_sum=znuc2_sum+znuc(iwctype(i))**2
      enddo

      call pot_nn_ewald
      write(ounit,'(''pecent='',f12.6)') pecent

! images for periodic basis functions
      if(n_images.ge.1)then

         if (ipr.ge.4 ) write(ounit,*) "Initialization periodic images distances"
!     assuming same number of images in each direction
         nix=n_images
         niy=n_images
         niz=n_images

         n_images=(2*nix+1)
         n_images=n_images*n_images*n_images
!     decreasing number of images by 1 temporally reagarding the ao's implementation at basis_fns
         n_images=n_images-1
         if (ipr.ge.4 ) write(ounit,*) "Total number of periodic images set", n_images

         deallocate (ell)
         allocate (ell(3, n_images))
         ell=0.d0

! assuming cubic boxes
! set images counter to zero
         imcount=0
         do iz=-niz,niz,1
            do iy=-niy,niy,1
               do ix=-nix,nix,1
                  nisum=abs(ix)+abs(iy)+abs(iz)
                  if (nisum.gt.0) then
!                     print*,ix,iy,iz
                     imcount=imcount+1
                     if(ix.ne.0) ell(1,imcount)=1.d0*ix*rlatt(1,1)
                     if(iy.ne.0) ell(2,imcount)=1.d0*iy*rlatt(2,2)
                     if(iz.ne.0) ell(3,imcount)=1.d0*iz*rlatt(3,3)
                  endif
               enddo
            enddo
         enddo
         if (ipr.ge.4) write(ounit,*) "Total number of peridic images counted: ",imcount

      endif

      return
      end
!-----------------------------------------------------------------------

!> Written by Cyrus Umrigar
!> distcell(i) is the perpendicular distance between cell faces parallel
!> to the other 2 directions from i.
!> dist_min is the shortest of these three.
!> By choosing the range of the short-range part of the Ewald sums to be
!> <= half the shortest perpendicular distance we ensure that the short-range
!> part has zero or one terms.
      subroutine short_distance(vector,volume,dist_min,distcell)

      use contrl_file,    only: ounit
      use precision_kinds, only: dp
      implicit none

      real(dp) :: dist_min, vlen, volume
      real(dp), dimension(3,3) :: vector
      real(dp), dimension(3) :: v1
      real(dp), dimension(3) :: v2
      real(dp), dimension(3) :: v3
      real(dp), dimension(3) :: distcell

      v1(1)=vector(1,2)
      v1(2)=vector(2,2)
      v1(3)=vector(3,2)

      v2(1)=vector(1,3)
      v2(2)=vector(2,3)
      v2(3)=vector(3,3)

      call cross(v1,v2,v3)
      vlen=sqrt(v3(1)**2+v3(2)**2+v3(3)**2)
      distcell(1)=volume/vlen
      dist_min=distcell(1)

      v1(1)=vector(1,3)
      v1(2)=vector(2,3)
      v1(3)=vector(3,3)

      v2(1)=vector(1,1)
      v2(2)=vector(2,1)
      v2(3)=vector(3,1)

      call cross(v1,v2,v3)
      vlen=sqrt(v3(1)**2+v3(2)**2+v3(3)**2)
      distcell(2)=volume/vlen
      dist_min=min(dist_min,distcell(2))

      v1(1)=vector(1,1)
      v1(2)=vector(2,1)
      v1(3)=vector(3,1)

      v2(1)=vector(1,2)
      v2(2)=vector(2,2)
      v2(3)=vector(3,2)

      call cross(v1,v2,v3)
      vlen=sqrt(v3(1)**2+v3(2)**2+v3(3)**2)
      distcell(3)=volume/vlen
      dist_min=min(dist_min,distcell(3))

      return
      end
!-----------------------------------------------------------------------

!> evaluates the cross-product of v1 and v2 and puts it in v3
!> @author Edgar Josue Landinez Borda
      subroutine cross(v1,v2,v3)

      use precision_kinds, only: dp
      implicit none

      real(dp), dimension(3) :: v1
      real(dp), dimension(3) :: v2
      real(dp), dimension(3) :: v3

      v3(1) = v1(2) * v2(3) - v1(3) * v2(2)
      v3(2) = v1(3) * v2(1) - v1(1) * v2(3)
      v3(3) = v1(1) * v2(2) - v1(2) * v2(1)

      return
      end
!-----------------------------------------------------------------------

      subroutine shells(cutg,glatt,gdist,igvec,gvec,gnorm,igmult,ngvec_big,ngnorm_big,ng1d)
      use ewald_mod, only: NGVEC_BIGX
! Written by Cyrus Umrigar

      use precision_kinds, only: dp
      use error, only: fatal_error
      implicit none

      integer :: i1, i2, i2min, i3, i3min
      integer :: k, ngnorm_big, ngvec_big
      integer, dimension(3,*) :: igvec
      integer, dimension(*) :: igmult
      integer, dimension(*) :: ng1d
      real(dp) :: cutg, cutg2, glen2, gx, gy
      real(dp) :: gz
      real(dp), dimension(3,*) :: glatt
      real(dp), dimension(3) :: gdist
      real(dp), dimension(3,*) :: gvec
      real(dp), dimension(*) :: gnorm
      real(dp), dimension(NGVEC_BIGX) :: gnorm_tmp

      do k=1,3
         ng1d(k)=int(cutg/gdist(k))
      enddo

      cutg2=cutg**2
      ngvec_big=0
!     do 10 i1=-ng1d(1),ng1d(1)
      do i1=0,ng1d(1)
        if(i1.ne.0) then
          i2min=-ng1d(2)
         else
          i2min=0
        endif
!       do 10 i2=-ng1d(2),ng1d(2)
        do i2=i2min,ng1d(2)
          if(i2.ne.0.or.i1.ne.0) then
            i3min=-ng1d(3)
           else
            i3min=0
          endif
!     do 10 i3=-ng1d(3),ng1d(3)
          do i3=i3min,ng1d(3)

            gx=i1*glatt(1,1)+i2*glatt(1,2)+i3*glatt(1,3)
            gy=i1*glatt(2,1)+i2*glatt(2,2)+i3*glatt(2,3)
            gz=i1*glatt(3,1)+i2*glatt(3,2)+i3*glatt(3,3)

            glen2=gx*gx+gy*gy+gz*gz

            if(glen2.le.cutg2) then
              ngvec_big=ngvec_big+1
              if(ngvec_big.gt.NGVEC_BIGX) then
                call fatal_error ('ngvec_big > NGVEC_BIGX in shells')
              endif

              igvec(1,ngvec_big)=i1
              igvec(2,ngvec_big)=i2
              igvec(3,ngvec_big)=i3

              gvec(1,ngvec_big)=gx
              gvec(2,ngvec_big)=gy
              gvec(3,ngvec_big)=gz

              gnorm_tmp(ngvec_big)=dsqrt(glen2)
            endif
          enddo
        enddo
      enddo

      call sort(igvec,gvec,gnorm_tmp,gnorm,igmult,ngvec_big,ngnorm_big)

      return
      end
!-----------------------------------------------------------------------

      subroutine sort(igvec,gvec,gnorm_tmp,gnorm,igmult,ngvec_big,ngnorm_big)
      use ewald_mod, only: NGNORM_BIGX
! Written by Cyrus Umrigar
      use contrl_file,    only: ounit
      use control, only: ipr
      use error, only: fatal_error
      use precision_kinds, only: dp
      implicit none

      integer :: i, icheck, icount, it
      integer :: j, k, l, lognb2
      integer :: m, ngnorm_big, ngvec_big, nn
      integer, dimension(3,*) :: igvec
      integer, dimension(*) :: igmult
      real(dp) :: t
      real(dp), dimension(3,*) :: gvec
      real(dp), dimension(*) :: gnorm_tmp
      real(dp), dimension(*) :: gnorm
      real(dp), parameter :: eps = 1.d-12

      lognb2=int(dlog(dfloat(ngvec_big))/dlog(2.d0)+1.d-14)
      m=ngvec_big
      do nn=1,lognb2
        m=m/2
        k=ngvec_big-m
        do j=1,k
          do i=j,1,-m
            l=i+m
            if (gnorm_tmp(l).gt.gnorm_tmp(i)-eps) goto 20
            t=gnorm_tmp(i)
            gnorm_tmp(i)=gnorm_tmp(l)
            gnorm_tmp(l)=t
            do k=1,3
              it=igvec(k,i)
              igvec(k,i)=igvec(k,l)
              igvec(k,l)=it
              t=gvec(k,i)
              gvec(k,i)=gvec(k,l)
              gvec(k,l)=t
            enddo
          enddo
      20     continue
        enddo
      enddo

! figure out the multiplicities and convert gnorm from being ngvec_big long to being ngnorm_big long
      ngnorm_big=1
      icount=0
      do i=2,ngvec_big
        icount=icount+1
        if(gnorm_tmp(i)-gnorm_tmp(i-1).gt.eps) then
          igmult(ngnorm_big)=icount
          gnorm(ngnorm_big)=gnorm_tmp(i-1)
          ngnorm_big=ngnorm_big+1
          if(ngnorm_big.gt.NGNORM_BIGX) then
            !if(ipr.ge.4) write(ounit,'(''ngnorm_big='',i8)') ngnorm_big
            write(ounit,'(''ngnorm_big='',i8)') ngnorm_big
            call fatal_error ('ngnorm_big > NGNORM_BIGX in sort')
          endif
          icount=0
        endif
      enddo
      igmult(ngnorm_big)=icount+1
      gnorm(ngnorm_big)=gnorm_tmp(ngvec_big)

      icheck=0
      do i=1,ngnorm_big
        icheck=icheck+igmult(i)
      enddo
      if(icheck.ne.ngvec_big) call fatal_error ('problem in sort')

      return
      end
!-----------------------------------------------------------------------

      subroutine separate(v,b0,lowest_pow,ngnorm_big,igmult,gnorm,ngnorm, &
       cutr,vcell,ncoef_per,np,b,y,chisq,ifcon)
! Written by Cyrus Umrigar and Claudia Filippi

      use ewald_mod, only: NCOEFX, NPX
!      use constant, only: twopi
      use precision_kinds, only: dp
      use contrl_file,    only: ounit
      use control,    only: ipr
      use error, only: fatal_error
      implicit none

      integer :: i, i0, ifcon, ig, info
      integer :: j, k, lowest_pow
      integer :: ncoef_per, nfree, ngnorm, ngnorm_big
      integer :: np
      integer, dimension(*) :: igmult
      real(dp) :: anorm, b0, beta1, beta2, chisq
      real(dp) :: cutr, gr, rcond, vcell
      real(dp) :: vk, pi,twopi
      real(dp), dimension(NCOEFX,NCOEFX) :: a
      real(dp), dimension(NCOEFX+NPX) :: c
      real(dp), dimension(NCOEFX) :: work
      real(dp), dimension(*) :: v
      real(dp), dimension(*) :: b
      real(dp), dimension(*) :: y
      real(dp), dimension(*) :: gnorm

!     parameter(NPX=6)

      pi=4.d0*datan(1.d0)
      twopi=2*pi

      if(ncoef_per+np.gt.NCOEFX+NPX) call fatal_error ('ncoef_per+np > NCOEFX+NPX in separate')

      anorm=2*twopi*cutr**3/vcell

! check for imposed conditions
      if(ifcon.ne.1) then
        nfree=ncoef_per
        i0=1
        beta1=0.d0
        beta2=0.d0
       else
! one less equation because of cusp conditions
        nfree=ncoef_per-1
        i0=2
! setting up cusp condition constraints
        if(lowest_pow.eq.-1) then
! 1/r behavior
          beta1=1/cutr
          beta2=0.d0
         elseif(lowest_pow.eq.0) then
! e-e cusp conditions
          beta1=-b0*cutr/np
          beta2=1.d0/np
         else
          call fatal_error ('lowest_pow must be -1 or 0')
        endif
      endif

      if(ipr.ge.4) then
         write(ounit,'(/,''Ncoef ='',i5)') ncoef_per
         write(ounit,'(''Beta1,beta2 ='',2f15.8)') beta1,beta2
      endif

! zero right and left hand side of fitting equation
      do i=1,ncoef_per
        b(i)=0.d0
        do j=1,ncoef_per
          a(j,i)=0.d0
        enddo
      enddo

      chisq=0.d0
      if(ipr.ge.4) then
         write(ounit,*) "chisq", chisq
!     go over k values larger than those explicitly used
         write(ounit,'(''do ewald breakup compute ckn'')')
         write(ounit,'(''ngnorm, ngnorm_big'',2i9)') ngnorm, ngnorm_big
      endif

      do k=ngnorm+1,ngnorm_big
        gr=gnorm(k)*cutr
        ig=k
        call integral_sin_poly(gr,lowest_pow,ncoef_per,np,anorm,c)

!     Constraints.  That for c(2) is for cusp constraint only.

        vk=v(k)-beta1*c(1)
        c(2)=c(2)+beta2*c(1)

        chisq=chisq+igmult(k)*vk**2

!        write(ounit,*) "k",k,"chisq(k)", chisq

!        write(ounit,'(''vk='',d14.5)') vk

! add to right hand side
        do i=i0,ncoef_per
          b(i)=b(i)+igmult(k)*vk*c(i)
! add to left hand side
          do j=i0,ncoef_per
      20       a(j,i)=a(j,i)+igmult(k)*c(i)*c(j)
          enddo
        enddo
      enddo

      if(ipr.ge.4) then
         write(ounit,'(''c='',10d14.5)') (c(i),i=i0,ncoef_per)
         write(ounit,'(''a='',10d14.5)') ((a(i,j),i=i0,ncoef_per),j=i0,ncoef_per)
         write(ounit,'(''b='',10d14.5)') (b(i),i=i0,ncoef_per)
      endif

! to be reinserted
! invert right hand side
!      if(nfree.gt.0) then
!        call dpoco(a(i0,i0),NCOEFX,nfree,rcond,work,info)
!        write(ounit,'(''condition #, rcond, after return from dpoco'',d12.4)') rcond
!        if(rcond.lt.1.d-14) call fatal_error ('rcond too small in dpoco')
!        if(info.ne.0) call fatal_error ('info in dpoco.ne.0 when called from separate')
!      endif

! make a spare copy of right hand side
      do i=i0,ncoef_per
        work(i)=b(i)
      enddo

! solve linear equations
!     call dposl(a(i0,i0),NCOEFX,nfree,b(i0))
      call dposv('U',nfree,1,a(i0,i0),NCOEFX,b(i0),nfree,info)

!     write(ounit,*) (b(i),i=i0,ncoef_per)

! b is now the solution (t in Ceperley's paper)
      do i=i0,ncoef_per
        chisq=chisq-work(i)*b(i)
      enddo

! this is cusp constraint
      if(ifcon.eq.1) b(1)=beta1+beta2*b(2)

! subtract effect of short range potential on fourier components

      do k=1,ngnorm
        gr=gnorm(k)*cutr
        ig=k
        call integral_sin_poly(gr,lowest_pow,ncoef_per,np,anorm,c)
        y(k)=v(k)
        do i=1,ncoef_per
      50     y(k)=y(k)-c(i)*b(i)
        enddo
      enddo

!     write(ounit,'(''Poly coefs (t) = '',5d14.6)') (b(i),i=1,ncoef_per)
!     write(ounit,'(''Yk = '',20d12.4)') (y(k),k=1,ngnorm)

      return
      end
!-----------------------------------------------------------------------

      subroutine integral_sin_poly(g,lowest_pow,n,np,anorm,c)
! Written by Cyrus Umrigar and Claudia Filippi
! anorm = 4*pi*cutr^3/volume
! g = g*cutr
! x = r/cutr
! output coefficients c

      use precision_kinds, only: dp
      implicit none

      integer :: i, j, k, lowest_pow, n
!     integer :: np, npts
      integer :: np
      real(dp) :: anorm,  dcmplx, dx, g
!     real(dp) :: gi, sin, ti,et,em, x
      real(dp) :: gi, sin, x
      real(dp), dimension(*) :: c
      integer, parameter :: NPTS = 2001
      real(dp), dimension(NPTS) :: y

      complex*16 ti,et,em

! integrates sin(g*x)*x**i for i=lowest_pow+1 to n+np+lowest_pow and x from 0 to 1
      if(dabs(g).gt.1.d-10) then
        gi=1.d0/g
        ti=dcmplx(0.d0,-gi)
        et=dcmplx(dsin(g)*gi,-dcos(g)*gi)
        em=ti*(et-ti)
        do i=1,n+np+lowest_pow+1
          if(i.gt.lowest_pow+1) c(i-lowest_pow-1)=dreal(em)
          em=ti*(et-i*em)
        enddo
       else
        do i=1,n+np+lowest_pow+1
          c(i)=1.d0/(i+2+lowest_pow)
        enddo
      endif

! take care that expansion functions are h_i(x) = x**i*(1-x)**np
! Warning check if we need to go one more.
      do k=1,np
        do i=1,n+np-k
          c(i)=c(i)-c(i+1)
        enddo
      enddo

!     write(ounit,'(''g,c1='',f5.1,9f9.5)') g,(c(i),i=1,n)

! Calculate c from numerical integral rather than recursion for small non-zero g's
       if(g.ne.0.d0 .and. g.lt.10.d0) then
       dx=1.d0/(NPTS-1)
       do i=1,n
         do j=1,NPTS
           x=(j-1)*dx
           if(g.gt.1.d-9) then
             if(i+lowest_pow.ne.0) then
               y(j)=x**(i+lowest_pow)*(1-x)**np*sin(g*x)
              else
               y(j)=(1-x)**np*sin(g*x)
             endif
            else
             y(j)=x**(i+1+lowest_pow)*(1-x)**np
           endif
         enddo
         c(i)=bode(y,dx,NPTS)
         if(g.gt.1.d-6) then
           c(i)=c(i)/g
         endif
       enddo

!     write(ounit,'(''g,c2='',f5.1,9f9.5)') g,(c(i),i=1,n)
      endif

! multiply by anorm
      do i=1,n
        c(i)=anorm*c(i)
      enddo

      return
      end
!-----------------------------------------------------------------------

      real*8 function vsrange(r,cutr,lowest_pow,ncoef_per,np,b)
! Written by Cyrus Umrigar and Claudia Filippi
! h(x)= \sum_{i=1}^ncoef_per b_i x^{i-1} (1-x)^np, x=r/cutr

      use precision_kinds, only: dp
      implicit none

      integer :: i, lowest_pow, ncoef_per, np
      real(dp) :: cutr, r, x
      real(dp), dimension(*) :: b

      x=r/cutr
      vsrange=0.d0
      if(x.gt.1.d0) return
      do i=1,ncoef_per
        vsrange=b(ncoef_per-i+1)+x*vsrange
      enddo
      vsrange=vsrange*(1-x)**np
      if(lowest_pow.eq.-1) vsrange=vsrange/x

      return
      end
!-----------------------------------------------------------------------

      real*8 function ewald_pot(rvec,rr,gvec,gnorm,ngnorm,igmult,y,cutr,vcell)
! Written by Cyrus Umrigar

!      use constants, only: pi
      use precision_kinds, only: dp
      implicit none

      integer :: im, ivec, k, ngnorm
      integer, dimension(*) :: igmult
!     real(dp) :: cos, cutr, derfc, expon, pi
      real(dp) :: cutr, expon, pi, g2expinv
      real(dp) :: gauss_exp, product, rr, vcell
      real(dp), dimension(3) :: rvec
      real(dp), dimension(3,*) :: gvec
      real(dp), dimension(*) :: gnorm
      real(dp), dimension(*) :: y

      pi=4.d0*datan(1.d0)

      gauss_exp=5.0/cutr

      g2expinv=1.d0/(4.0*gauss_exp**2)

      ivec=1
! The factor of 2 in the next line is just to compensate for the 2 in the
! last line, which is there because we keep only half the vectors in the star.
!     ewald_pot=-pi/(2*vcell*gaus_exp**2)
      ewald_pot=-(2.d0*pi/vcell)*g2expinv
      do k=2,ngnorm
!     expon=exp(-(gnorm(k)/(2*gauss_exp))**2)
        expon=exp(-g2expinv*gnorm(k))
        do im=1,igmult(k)
          ivec=ivec+1
          product=rvec(1)*gvec(1,ivec)+ &
                  rvec(2)*gvec(2,ivec)+ &
                  rvec(3)*gvec(3,ivec)
          ewald_pot=ewald_pot+cos(product)*y(k)*expon
        enddo
      enddo
      ewald_pot=2*ewald_pot+y(1)+erfc(gauss_exp*rr)/rr

      return
      end

!-----------------------------------------------------------------------

      real*8 function vlrange_r(rvec,gvec,ngnorm,igmult,y)
! Written by Cyrus Umrigar

      use precision_kinds, only: dp
      implicit none

      integer :: im, ivec, k, ngnorm
      integer, dimension(*) :: igmult
      real(dp) :: cos, product, vlrange
      real(dp), dimension(3) :: rvec
      real(dp), dimension(3,*) :: gvec
      real(dp), dimension(*) :: y

      ivec=1
      vlrange=0.d0
      do k=2,ngnorm
        do im=1,igmult(k)
          ivec=ivec+1
          product=rvec(1)*gvec(1,ivec)+ &
                  rvec(2)*gvec(2,ivec)+ &
                  rvec(3)*gvec(3,ivec)
          vlrange=vlrange+cos(product)*y(k)
        enddo
      enddo
      vlrange_r=2*vlrange+y(1)

      return
      end
!-----------------------------------------------------------------------

      real*8 function vlrange(ngnorm,igmult,cos1_sum,cos2_sum,sin1_sum,sin2_sum,y)
! Written by Cyrus Umrigar

      use precision_kinds, only: dp
      use contrl_file,    only: ounit

      implicit none

      integer :: im, ivec, k, ngnorm
      integer, dimension(*) :: igmult
      real(dp) :: vl
      real(dp), dimension(*) :: cos1_sum
      real(dp), dimension(*) :: cos2_sum
      real(dp), dimension(*) :: sin1_sum
      real(dp), dimension(*) :: sin2_sum
      real(dp), dimension(*) :: y

      ivec=1
      vl=0.5d0*y(1)*(cos1_sum(1)*cos2_sum(1)+sin1_sum(1)*sin2_sum(1))
      do k=2,ngnorm
        do im=1,igmult(k)
          ivec=ivec+1
          vl=vl+y(k)*(cos1_sum(ivec)*cos2_sum(ivec)+sin1_sum(ivec)*sin2_sum(ivec))
        enddo
      enddo
      vlrange=vl

      return
      end
!-----------------------------------------------------------------------

      real*8 function bode(y,h,n)

      use precision_kinds, only: dp
      implicit none

!  need to recover pw_ewald original Cyrus Umrigar implementation

      real(dp), parameter :: oneb45=1.d0/45.d0
      real(dp), parameter :: c0=28.d0
      real(dp), parameter :: c1=64.d0
      real(dp), parameter :: c2=24.d0
      real(dp), parameter :: half=.5d0
      real(dp), dimension (*) :: y
      real(dp) :: s, h

      integer :: n, n1, i

      if(mod(n-1,4).ne.0) stop 'n must be 4*n+1 in bode'

      n1=n-1
      s=half*c0*(y(n)-y(1))
      do i=1,n1,4
         s=s+c0*y(i)+c1*(y(i+1)+y(i+3))+c2*y(i+2)
      enddo
      bode=s*h*oneb45
      return
      end

!-----------------------------------------------------------------------

      subroutine pot_nn_ewald
! Written by Cyrus Umrigar

      use contrl_file,    only: ounit
      use control,    only: ipr
      use system, only: znuc, cent, iwctype, ncent
      use multiple_geo, only: pecent

      use ewald, only: b_coul, y_coul, cos_n_sum, sin_n_sum

      use periodic, only: cutr, glatt
      use periodic, only: igmult, igvec
      use periodic, only: ncoef_per, ng1d, ngnorm
      use periodic, only: ngvec
      use periodic, only: np_coul, vcell
      use periodic, only: vcell, znuc2_sum, znuc_sum
      use precision_kinds, only: dp
      use find_pimage, only: find_image_pbc
      implicit none

      integer :: i, j, k, lowest_pow
      real(dp) :: c0, rnorm, vl
      real(dp) :: vs, zprod
      real(dp), dimension(3) :: r

! short-range sum
      lowest_pow=-1
      c0=(b_coul(2)-np_coul*b_coul(1))/2
      vs=c0*znuc2_sum
      do i=1,ncent
        do j=1,i-1
          zprod=znuc(iwctype(i))*znuc(iwctype(j))
          do k=1,3
            r(k)=cent(k,j)-cent(k,i)
          enddo
          call find_image_pbc(r,rnorm)
          vs=vs+zprod*vsrange(rnorm,cutr,lowest_pow,ncoef_per,np_coul,b_coul)
        enddo
      enddo

! long-range sum
      call cossin_n(znuc,iwctype,glatt,igvec,ngvec,cent,ncent,ng1d,cos_n_sum,sin_n_sum)

      vl=vlrange(ngnorm,igmult,cos_n_sum,cos_n_sum,sin_n_sum,sin_n_sum,y_coul)

      pecent=vs+vl
      vs=vs*2/ncent
      vl=vl*2/ncent
      if(ipr.ge.4) then
         write(ounit,'(''v_nn,vs,vl,vs1,vl1='',9f12.8)') pecent*2/ncent,vs,vl,znuc2_sum*c0*2/ncent,y_coul(1)*znuc_sum**2/ncent
      endif

      return
      end

!-----------------------------------------------------------------------

      subroutine pot_en_ewald(x,pe_en)
! Written by Cyrus Umrigar

      use contrl_file,    only: ounit
      use vmc_mod, only: nmat_dim2
      use system, only: znuc, cent, iwctype, ncent, ncent_tot
      use control, only: ipr
      use system, only: nelec
      use ewald, only: b_coul, y_coul, y_jas, b_jas
      use ewald, only: cos_n_sum, sin_n_sum, cos_e_sum, sin_e_sum
      use periodic, only: cutr, glatt
      use periodic, only: igmult, igvec
      use periodic, only: ncoef_per, ng1d, ngnorm
      use periodic, only: ngvec
      use periodic, only: np_coul
      use periodic, only: znuc_sum
      use pseudo, only: lpot, nloc
      use distance_mod, only: r_en, rvec_en

      use find_pimage, only: find_image_pbc

      use precision_kinds, only: dp
      implicit none

      integer :: i, ict, j, k, lowest_pow
      real(dp) :: pe_en, vl
      real(dp) :: vs, vs_aux
      real(dp), dimension(3,*) :: x

! short-range sum
! Warning: I need to call the appropriate vsrange
      vs=0.d0
      do i=1,ncent
        ict=iwctype(i)
        do j=1,nelec
          do k=1,3
             rvec_en(k,j,i)=x(k,j)-cent(k,i)
          enddo
          call find_image_pbc(rvec_en(1,j,i),r_en(j,i))
      lowest_pow=-1
          vs=vs-znuc(ict)*vsrange(r_en(j,i),cutr,lowest_pow,ncoef_per,np_coul,b_coul)
        enddo
      enddo

! long-range sum
      call cossin_e(glatt,igvec,ngvec,x,nelec,ng1d,cos_e_sum,sin_e_sum)

      vl=-2*vlrange(ngnorm,igmult,cos_n_sum,cos_e_sum,sin_n_sum,sin_e_sum,y_coul)

      pe_en=vs+vl

      vs=vs/nelec
      vl=vl/nelec
      if(ipr.ge.2) then
        if(nloc.eq.0) write(ounit,'(''v_en,vs,vl,vl1='',9f12.8)') pe_en/nelec,vs,vl,-y_coul(1)*znuc_sum
      endif

      return
      end

!-----------------------------------------------------------------------

      subroutine pot_ee_ewald(x,pe_ee)
! Written by Cyrus Umrigar

      use contrl_file,    only: ounit
      use vmc_mod, only: nmat_dim2
      use control, only: ipr
      use system, only: nelec
      use ewald, only: b_coul, y_coul
      use ewald, only: cos_n_sum, sin_n_sum, cos_e_sum, sin_e_sum

      use periodic, only: cutr
      use periodic, only: glatt, igmult, igvec
      use periodic, only: ncoef_per, ng1d
      use periodic, only: ngnorm, ngvec
      use periodic, only: np_coul
      use distance_mod, only: r_en, rvec_en, r_ee, rvec_ee
      use find_pimage, only: find_image_pbc
      use pseudo, only: lpot, nloc

      use precision_kinds, only: dp
      implicit none

      integer :: i, ij, j, k, lowest_pow
      real(dp) :: c0, pe_ee, vl
      real(dp) :: vs
      real(dp), dimension(3,*) :: x

!      write(ounit,*) "inside pot_ee_ewald nloc", nloc

! short-range sum
      lowest_pow=-1
      c0=(b_coul(2)-np_coul*b_coul(1))/2
      vs=c0*nelec
      ij=0
      do i=1,nelec
        do j=1,i-1
          ij=ij+1
          do k=1,3
            rvec_ee(k,ij)=x(k,i)-x(k,j)
          enddo
          call find_image_pbc(rvec_ee(1,ij),r_ee(ij))
          vs=vs+vsrange(r_ee(ij),cutr,lowest_pow,ncoef_per,np_coul,b_coul)
        enddo
      enddo

! long-range sum
      call cossin_e(glatt,igvec,ngvec,x,nelec,ng1d,cos_e_sum,sin_e_sum)

      vl=vlrange(ngnorm,igmult,cos_e_sum,cos_e_sum,sin_e_sum,sin_e_sum,y_coul)

      pe_ee=vs+vl

      vs=vs*2/nelec
      vl=vl*2/nelec
      if(ipr.ge.2) write(ounit,'(''v_ee,vs,vl,vs1,vl1='',9f12.8)') pe_ee*2/nelec,vs,vl,c0*2,y_coul(1)*nelec

      return
      end

!-----------------------------------------------------------------------

      subroutine cossin_psi(iel,glatt,gnorm,gvec,igvec,ngvec,r,nr,ng1d,cos_g,sin_g &
      ,dcos_g,dsin_g,ddcos_g,ddsin_g)

! Written by Cyrus Umrigar
! iflag = 0 Calculate cos(gr) and sin(gr) and first 2 derivs at electron positions.
!       = 1 Calculate cos(kr) and sin(kr) and first 2 derivs at electron positions.
! Needed for orbitals and their Laplacian.
! Presently using cossin_psi_g and cossin_psi_k instead.

      use contrl_file,    only: ounit
      use ewald_mod, only: NG1DX
      use system, only: nelec
      use periodic, only: igmult,ngnorm
      use precision_kinds, only: dp
      implicit none

      integer :: i, iel, im, ir, i1, i2, k, n
      integer :: ngvec, nr
      integer, dimension(3,*) :: igvec
      integer, dimension(3) :: ng1d
      real(dp) :: cos, cos_tmp0, cos_tmp1, cos_tmp2, dot
      real(dp) :: sin, sin_tmp0, sin_tmp1, sin_tmp2
      real(dp), dimension(3,3) :: glatt
      real(dp), dimension(*) :: gnorm
      real(dp), dimension(3,*) :: gvec
      real(dp), dimension(3,*) :: r
      real(dp), dimension(nelec,*) :: cos_g
      real(dp), dimension(nelec,*) :: sin_g
      real(dp), dimension(3,nelec,*) :: dcos_g
      real(dp), dimension(3,nelec,*) :: dsin_g
      real(dp), dimension(nelec,*) :: ddcos_g
      real(dp), dimension(nelec,*) :: ddsin_g
      real(dp), dimension(-NG1DX:NG1DX,3) :: cos_gr
      real(dp), dimension(-NG1DX:NG1DX,3) :: sin_gr

! Calculate cosines and sines for recip. lattice vectors along axes first.
      if(iel.eq.0)then
        i1=1
        i2=nr
       else
        i1=iel
        i2=iel
      endif
      do ir=i1,i2

      do i=1,3
        dot=0
        do k=1,3
          dot=dot+glatt(k,i)*r(k,ir)
        enddo
        cos_gr(1,i)=cos(dot)
        sin_gr(1,i)=sin(dot)
        cos_gr(-1,i)=cos_gr(1,i)
        sin_gr(-1,i)=-sin_gr(1,i)
        cos_gr(0,i)=1.d0
        sin_gr(0,i)=0.d0
        do n=2,ng1d(i)
          cos_gr(n,i)=cos_gr(n-1,i)*cos_gr(1,i)-sin_gr(n-1,i)*sin_gr(1,i)
          sin_gr(n,i)=sin_gr(n-1,i)*cos_gr(1,i)+cos_gr(n-1,i)*sin_gr(1,i)
          cos_gr(-n,i)=cos_gr(n,i)
          sin_gr(-n,i)=-sin_gr(n,i)
        enddo
      enddo

! If the calculation is for g-vectors then no shift; if for k-vectors there could be one.
      cos_tmp0=1.d0
      sin_tmp0=0.d0

      do i=1,ngvec
        cos_tmp1=cos_tmp0*cos_gr(igvec(1,i),1) &
                -sin_tmp0*sin_gr(igvec(1,i),1)
        sin_tmp1=sin_tmp0*cos_gr(igvec(1,i),1) &
                +cos_tmp0*sin_gr(igvec(1,i),1)
        cos_tmp2=cos_tmp1*cos_gr(igvec(2,i),2) &
                -sin_tmp1*sin_gr(igvec(2,i),2)
        sin_tmp2=sin_tmp1*cos_gr(igvec(2,i),2) &
                +cos_tmp1*sin_gr(igvec(2,i),2)
        cos_g(ir,i)=cos_tmp2*cos_gr(igvec(3,i),3) &
                   -sin_tmp2*sin_gr(igvec(3,i),3)
        sin_g(ir,i)=sin_tmp2*cos_gr(igvec(3,i),3) &
                   +cos_tmp2*sin_gr(igvec(3,i),3)
        do k=1,3
          dcos_g(k,ir,i)=-gvec(k,i)*sin_g(ir,i)
          dsin_g(k,ir,i)= gvec(k,i)*cos_g(ir,i)
        enddo
!       if(i.lt.5) write(ounit,'(''ir,i,gnorm(i),cos_g(ir,i),sin_g(ir,i),dcos_g(k,ir,i),dsin_g(k,ir,i)'',2i5,9d12.4)')
!    & ir,i,gnorm(i),cos_g(ir,i),sin_g(ir,i),(dcos_g(k,ir,i),dsin_g(k,ir,i),k=1,3)
      enddo

!     if(iel.eq.0) then
        i=0
        do k=1,ngnorm
          do im=1,igmult(k)
            i=i+1
            ddcos_g(ir,i)=-gnorm(k)*gnorm(k)*cos_g(ir,i)
            ddsin_g(ir,i)=-gnorm(k)*gnorm(k)*sin_g(ir,i)
          enddo
        enddo
!     endif

      enddo
      return
      end

!-----------------------------------------------------------------------
! The only diff. between cossin_n, cossin_p and cossin_e is whether the
! nuclear charge, pseudopotential or electronic charge is used, so they
! could be merged, but I did not do that to get a small gain in efficiency.
!-----------------------------------------------------------------------

      subroutine cossin_n(znuc,iwctype,glatt,igvec,ngvec,r,nr,ng1d,cos_sum,sin_sum)
      use ewald_mod, only: NG1DX
      use system, only: ncent_tot
! Written by Cyrus Umrigar
! Calculate cos_sum and sin_sum for nuclei

      use precision_kinds, only: dp
      implicit none

      integer :: i, ir, k, n, ngvec
      integer :: nr
      integer, dimension(*) :: iwctype
      integer, dimension(3,*) :: igvec
      integer, dimension(3) :: ng1d
      real(dp) :: cos, cos_tmp, dot, sin, sin_tmp
      real(dp), dimension(*) :: znuc
      real(dp), dimension(3,3) :: glatt
      real(dp), dimension(3,*) :: r
      real(dp), dimension(*) :: cos_sum
      real(dp), dimension(*) :: sin_sum
      real(dp), dimension(-NG1DX:NG1DX,3,ncent_tot) :: cos_gr
      real(dp), dimension(-NG1DX:NG1DX,3,ncent_tot) :: sin_gr

! Calculate cosines and sines for all positions and reciprocal lattice vectors
      do ir=1,nr
        do i=1,3
          dot=0
          do k=1,3
            dot=dot+glatt(k,i)*r(k,ir)
          enddo
          cos_gr(1,i,ir)=cos(dot)
          sin_gr(1,i,ir)=sin(dot)
          cos_gr(-1,i,ir)=cos_gr(1,i,ir)
          sin_gr(-1,i,ir)=-sin_gr(1,i,ir)
          cos_gr(0,i,ir)=1.d0
          sin_gr(0,i,ir)=0.d0
          do n=2,ng1d(i)
            cos_gr(n,i,ir)=cos_gr(n-1,i,ir)*cos_gr(1,i,ir)-sin_gr(n-1,i,ir)*sin_gr(1,i,ir)
            sin_gr(n,i,ir)=sin_gr(n-1,i,ir)*cos_gr(1,i,ir)+cos_gr(n-1,i,ir)*sin_gr(1,i,ir)
            cos_gr(-n,i,ir)=cos_gr(n,i,ir)
            sin_gr(-n,i,ir)=-sin_gr(n,i,ir)
          enddo
        enddo
      enddo

      do i=1,ngvec
        cos_sum(i)=0
        sin_sum(i)=0
        do ir=1,nr
          cos_tmp=cos_gr(igvec(1,i),1,ir)*cos_gr(igvec(2,i),2,ir) &
                 -sin_gr(igvec(1,i),1,ir)*sin_gr(igvec(2,i),2,ir)
          sin_tmp=sin_gr(igvec(1,i),1,ir)*cos_gr(igvec(2,i),2,ir) &
                 +cos_gr(igvec(1,i),1,ir)*sin_gr(igvec(2,i),2,ir)
          cos_sum(i)=cos_sum(i)+znuc(iwctype(ir))* &
                     (cos_tmp*cos_gr(igvec(3,i),3,ir) &
                     -sin_tmp*sin_gr(igvec(3,i),3,ir))
          sin_sum(i)=sin_sum(i)+znuc(iwctype(ir))* &
                     (sin_tmp*cos_gr(igvec(3,i),3,ir) &
                     +cos_tmp*sin_gr(igvec(3,i),3,ir))
        enddo
      enddo

      return
      end

!-----------------------------------------------------------------------

      subroutine cossin_e(glatt,igvec,ngvec,r,nr,ng1d,cos_sum,sin_sum)
      use ewald_mod, only: NG1DX
      use system, only: nelec
! Written by Cyrus Umrigar
! Calculate cos_sum and sin_sum for electrons

      use precision_kinds, only: dp
      implicit none

      integer :: i, ir, k, n, ngvec
      integer :: nr
      integer, dimension(3,*) :: igvec
      integer, dimension(3) :: ng1d
      real(dp) :: cos, cos_tmp, dot, sin, sin_tmp
      real(dp), dimension(3,3) :: glatt
      real(dp), dimension(3,*) :: r
      real(dp), dimension(*) :: cos_sum
      real(dp), dimension(*) :: sin_sum
      real(dp), dimension(-NG1DX:NG1DX,3,nelec) :: cos_gr
      real(dp), dimension(-NG1DX:NG1DX,3,nelec) :: sin_gr

! Calculate cosines and sines for all positions and reciprocal lattice vectors
      do ir=1,nr
        do i=1,3
          dot=0
          do k=1,3
            dot=dot+glatt(k,i)*r(k,ir)
          enddo
          cos_gr(1,i,ir)=cos(dot)
          sin_gr(1,i,ir)=sin(dot)
          cos_gr(-1,i,ir)=cos_gr(1,i,ir)
          sin_gr(-1,i,ir)=-sin_gr(1,i,ir)
          cos_gr(0,i,ir)=1.d0
          sin_gr(0,i,ir)=0.d0
          do n=2,ng1d(i)
            cos_gr(n,i,ir)=cos_gr(n-1,i,ir)*cos_gr(1,i,ir)-sin_gr(n-1,i,ir)*sin_gr(1,i,ir)
            sin_gr(n,i,ir)=sin_gr(n-1,i,ir)*cos_gr(1,i,ir)+cos_gr(n-1,i,ir)*sin_gr(1,i,ir)
            cos_gr(-n,i,ir)=cos_gr(n,i,ir)
            sin_gr(-n,i,ir)=-sin_gr(n,i,ir)
          enddo
        enddo
      enddo

      do i=1,ngvec
        cos_sum(i)=0
        sin_sum(i)=0
        do ir=1,nr
          cos_tmp=cos_gr(igvec(1,i),1,ir)*cos_gr(igvec(2,i),2,ir) &
                 -sin_gr(igvec(1,i),1,ir)*sin_gr(igvec(2,i),2,ir)
          sin_tmp=sin_gr(igvec(1,i),1,ir)*cos_gr(igvec(2,i),2,ir) &
                 +cos_gr(igvec(1,i),1,ir)*sin_gr(igvec(2,i),2,ir)
          cos_sum(i)=cos_sum(i)+ &
                     (cos_tmp*cos_gr(igvec(3,i),3,ir) &
                     -sin_tmp*sin_gr(igvec(3,i),3,ir))
          sin_sum(i)=sin_sum(i)+ &
                     (sin_tmp*cos_gr(igvec(3,i),3,ir) &
                     +cos_tmp*sin_gr(igvec(3,i),3,ir))
        enddo
      enddo

      return
      end

!-----------------------------------------------------------------------
      subroutine jastrow_longrange(iel,x,psij_per,d2_per,v_per,iflag)

      use contrl_file,    only: ounit
      use control, only: ipr
      use system, only: nelec
      use ewald, only: y_jas, b_jas
      use ewald, only: cos_sum, cos_g , dcos_g
      use ewald, only: sin_sum, sin_g , dsin_g
      use ewald, only: cos_sum_new, cos_g_new, dcos_g_new
      use ewald, only: sin_sum_new, sin_g_new, dsin_g_new
      use periodic, only: glatt, gnorm, gvec
      use periodic, only: igmult, igvec
      use periodic, only: ng1d, ngnorm
      use periodic, only: ngvec

      use ewald_mod, only: NG1DX

      use precision_kinds, only: dp
      implicit none

      integer :: i, iel, iflag, ir
      real(dp) :: psij_per, cos_tmp, sin_tmp, d2_per
      real(dp), save :: psij_per_old

      real(dp), dimension(3,*) :: x
      real(dp), dimension(3,*) :: v_per
      !real(dp), dimension(ngvec) :: cos_sum
      !real(dp), dimension(ngvec) :: sin_sum
      !real(dp), dimension(nelec,ngvec) :: cos_g
      !real(dp), dimension(nelec,ngvec) :: sin_g
      !real(dp), dimension(3,nelec,ngvec) :: dcos_g
      !real(dp), dimension(3,nelec,ngvec) :: dsin_g
      real(dp), dimension(nelec,ngvec) :: ddcos_g
      real(dp), dimension(nelec,ngvec) :: ddsin_g

      if(iel.eq.0) then

      call cossin_psi(iel,glatt,gnorm,gvec,igvec,ngvec,x,nelec,ng1d,cos_g,sin_g &
      ,dcos_g,dsin_g,ddcos_g,ddsin_g)

      do i=1,ngvec
        cos_sum(i)=0
        sin_sum(i)=0
        do ir=1,nelec
          cos_sum(i)=cos_sum(i)+ cos_g(ir,i)
          sin_sum(i)=sin_sum(i)+ sin_g(ir,i)
        enddo
      enddo

      psij_per=vlrange(ngnorm,igmult,cos_sum,cos_sum,sin_sum,sin_sum,y_jas)
      call dvlrange_ee(ngnorm,gnorm,igmult,cos_sum,sin_sum,dcos_g,dsin_g,ddcos_g,ddsin_g,y_jas,v_per,d2_per)

      psij_per_old=psij_per

      else

      cos_g_new(:,1:ngvec)=cos_g(:,1:ngvec)
      sin_g_new(:,1:ngvec)=sin_g(:,1:ngvec)
      dcos_g_new(:,:,1:ngvec)=dcos_g(:,:,1:ngvec)
      dsin_g_new(:,:,1:ngvec)=dsin_g(:,:,1:ngvec)

      call cossin_psi(iel,glatt,gnorm,gvec,igvec,ngvec,x,nelec,ng1d,cos_g_new,sin_g_new &
      ,dcos_g_new,dsin_g_new,ddcos_g,ddsin_g)

      do i=1,ngvec
        cos_sum_new(i)=cos_sum(i)+cos_g_new(iel,i)-cos_g(iel,i)
        sin_sum_new(i)=sin_sum(i)+sin_g_new(iel,i)-sin_g(iel,i)
      enddo
      psij_per=vlrange(ngnorm,igmult,cos_sum_new,cos_sum_new,sin_sum_new,sin_sum_new,y_jas)

      if(iflag.eq.0) then
        call dvlrange_ee(ngnorm,gnorm,igmult,cos_sum_new,sin_sum_new,dcos_g_new,dsin_g_new,ddcos_g,ddsin_g,y_jas,v_per,d2_per)
       else
! only compute ratio for non-loc potential
! WARNING ICASULA 3 might not work
         psij_per=psij_per_old-psij_per
      endif

      endif

      return
      end

      subroutine dvlrange_ee(ngnorm,gnorm,igmult,cos_sum,sin_sum,dcos_g,dsin_g,ddcos_g,ddsin_g,y,dvl,ddvl)
! Written by Cyrus Umrigar

      use precision_kinds, only: dp
      use system, only: nelec
      implicit none

      integer :: im, ir, ivec, j, k, ngnorm
      integer, dimension(*) :: igmult
      real(dp) :: ddvl, vl
      real(dp), dimension(3,*) :: dvl
      real(dp), dimension(*) :: cos_sum
      real(dp), dimension(*) :: sin_sum
      real(dp), dimension(*) :: y
      real(dp), dimension(*) :: gnorm
      real(dp), dimension(3,nelec,*) :: dcos_g, dsin_g
      real(dp), dimension(nelec,*) :: ddcos_g,ddsin_g

      ivec=1
      ddvl=0.d0
      do ir=1,nelec
       do j=1,3
        dvl(j,ir)=2*y(1)*(dcos_g(j,ir,ivec)*cos_sum(1)+dsin_g(j,ir,ivec)*sin_sum(1))
       enddo
       ddvl=ddvl+2*y(1)*(ddcos_g(ir,ivec)*cos_sum(1)+ddsin_g(ir,ivec)*sin_sum(1)+gnorm(1)*gnorm(1))
      enddo

      do k=2,ngnorm
        do im=1,igmult(k)
          ivec=ivec+1
          do ir=1,nelec
             do j=1,3
               dvl(j,ir)=dvl(j,ir)+2*y(k)*(dcos_g(j,ir,ivec)*cos_sum(ivec)+dsin_g(j,ir,ivec)*sin_sum(ivec))
             enddo
             ddvl=ddvl+2*y(k)*(ddcos_g(ir,ivec)*cos_sum(ivec)+ddsin_g(ir,ivec)*sin_sum(ivec)+gnorm(k)*gnorm(k))
          enddo
        enddo
      enddo

      return
      end
      end module
