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

      use pseudo_mod, only: MPS_L, MPS_GRID
      use ewald_mod, only: NGNORMX, NGVECX, NG1DX
      use ewald_mod, only: NGNORM_SIMX, NGVEC_SIMX, NCOEFX
      use system, only: znuc, iwctype, nctype, ncent, nctype_tot, ncent_tot
      use multiple_geo, only: pecent
      use control, only: ipr
      use ewald, only: b_coul, b_coul_sim, b_psp, b_jas, y_coul, y_coul_sim, y_psp, y_jas
      use ewald_basis, only: vps_basis_fourier
      use periodic, only: cutg, cutg_big, cutg_sim, cutg_sim_big, cutr, cutr_sim, glatt
      use periodic, only: glatt_inv, glatt_sim, gnorm, gnorm_sim, gvec, gvec_sim, igmult, igmult_sim, igvec, igvec_sim
      use periodic, only: isrange, ncoef_per, ng1d, ng1d_sim, ngnorm, ngnorm_big
      use periodic, only: ngnorm_sim, ngnorm_sim_big, ngvec, ngvec_big, ngvec_sim, ngvec_sim_big
      use periodic, only: np, npoly, rkvec_shift, rlatt_sim, rlatt_sim_inv, vcell
      use periodic, only: rlatt, rlatt_inv
      use periodic, only: vcell_sim, znuc2_sum, znuc_sum
      use tempor, only: dist_nn
      use test, only: f, vbare_coul, vbare_jas, vbare_psp
      use contrl_per, only: iperiodic
      use pseudo, only: lpot, nloc, vps
      use contrl_file,    only: ounit
      use grid3d_param, only: origin
      use error, only: fatal_error
      use find_pimage, only: check_lattice
      use matinv_mod, only: matinv
      use spline2_mod, only: spline2
      use readps_gauss, only: gauss_pot
      use periodic, only : n_images, ell


      use precision_kinds, only: dp

      implicit none

      integer :: i, ict, ifcon, ig, in
      integer :: ir, j, k, lowest_pow
      integer :: npts
      real(dp) :: pi, twopi
      real(dp) :: alpha, as
      real(dp) :: b0, because
      real(dp) :: chisq
      real(dp) :: coefs, components, constv, short_ewald
      real(dp) :: datan, derf, det, det1
      real(dp) :: det_sim, discontinuity, dist_min, dpot1
      real(dp) :: dpotn, dx, ewa
      real(dp) :: g2, g2a
      real(dp) :: gdistmin, gdistmin_sim
      real(dp) :: psp, rms
      real(dp) :: rr, sum, testv, test_s
      real(dp) :: those, true, true_s, vgcell
      real(dp) :: vgcell_sim, vl, vpot, dvpot, rg
      real(dp) :: vs, wt
      real(dp) :: deltar
      real(dp) :: vps_test
!      real(dp), dimension(MPS_GRID) :: r
!      real(dp), dimension(MPS_GRID) :: vps_short
!     real(dp), dimension(MPS_GRID) :: work
      real(dp), allocatable :: r(:)
      real(dp), allocatable :: vps_short(:)
      real(dp), dimension(3) :: rdist
      real(dp), dimension(3) :: gdist
      real(dp), dimension(3) :: rdist_sim
      real(dp), dimension(3) :: gdist_sim
      real(dp), dimension(3) :: rkvec_shift_tmp
      real(dp), dimension(3) :: r_tmp
      real(dp), dimension(3) :: r_fk
      real(dp), parameter :: eps = 1.d-12
      real(dp), dimension(nctype_tot) :: r0
      real(dp), dimension(nctype_tot) :: arg
      integer, dimension(nctype_tot) :: nr_ps
! this array is just for testing purposes
      real(dp), dimension(101) :: w_gauss
      real(dp), dimension(3) :: cutgv

! for images evaluation
      integer :: ix,iy,iz
      integer :: nix,niy,niz
      integer :: imcount
      integer :: nisum, ngmx

      pi=4.d0*datan(1.d0)
      twopi=2*pi

!     These are parameters defined by default in Champ's converter to get the grid
!     these parameters are the same for any element to get the basis grid file
!     the only different parameter is nrbas which define the number of angular shells
!     temporal parameters grid short range potential fft log-shifted grid
      r0=20.d0
      nr_ps=2000

!      r0=cutr/(arg**(nr_ps-1)-1)
!      arg(1)=cutr/((nr_ps(1)-1)*r0(1)))

!     allocate grid and vps_short / a.k.a fft  compontents
! asuming nr
      allocate(r(nr_ps(1)))
      allocate(vps_short(nr_ps(1)))

!     Number of coefficients polynomial expansion the breakup
      ncoef_per=npoly+1

! Check that the lattice vectors are the smallest possible ones and return the smallest
! which is used to set the range of the real-space Ewald sums so that only one image
! of a nucleus or an electron is present within cutr and cutr_sim respectively.
      call check_lattice(rlatt,cutr,0)
      call check_lattice(rlatt_sim,cutr_sim,1)

      write(ounit, *)  "================================================================================="
      write(ounit, *)  "=^..^=                                                                     =^..>="
      write(ounit, *)  "==                                                                             =="
      write(ounit, *)  "==  Setting-up Ewald-Breakup                                                   =="
      write(ounit, *)  "==                                                                             =="
      write(ounit, *)  "=^..^=                                                                     =<..^="
      write(ounit, *)  "================================================================================="

      write(ounit, *)  " "

      write(ounit,'(''cutr,cutr_sim ='',9f9.5)') cutr,cutr_sim

! this parameters setiing are necessary for FFT on the pseudopotential integral exponential grid
! same for any element since depends on r0,cutr and nr_ps=nr
!     computing arg depending on how to compute r (upper limit) it should be cutr in principle for vps(short)
!      print*,"cutr,deltar",cutr,deltar
      deltar=1.0*cutr/(r0(1)*(nr_ps(1)-1))
!      deltar=10.2*cutr/(r0(1)*(nr_ps(1)-1))
!      write(ounit,*) "cutr",cutr,deltar
      arg=(1.0d0+deltar)
!      write(ounit,*) "arg",arg(1)




! Calculate inverse transformations (from lattice coordinates to real coordinates)
! and cell volumes
      do i=1,3
        do k=1,3
          rlatt_inv(k,i)=rlatt(k,i)
          rlatt_sim_inv(k,i)=rlatt_sim(k,i)
        enddo
      enddo
      call matinv(rlatt_inv,3,det)
      call matinv(rlatt_sim_inv,3,det_sim)


! Primitive cell volume and reciprocal lattice
!     det=rlatt(1,1)*rlatt(2,2)*rlatt(3,3)
!    &   +rlatt(2,1)*rlatt(3,2)*rlatt(1,3)
!    &   +rlatt(3,1)*rlatt(1,2)*rlatt(2,3)
!    &   -rlatt(3,1)*rlatt(2,2)*rlatt(1,3)
!    &   -rlatt(1,1)*rlatt(3,2)*rlatt(2,3)
!    &   -rlatt(2,1)*rlatt(1,2)*rlatt(3,3)
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

! Simulation cell volume and reciprocal lattice
!     det=rlatt_sim(1,1)*rlatt_sim(2,2)*rlatt_sim(3,3)
!    &   +rlatt_sim(2,1)*rlatt_sim(3,2)*rlatt_sim(1,3)
!    &   +rlatt_sim(3,1)*rlatt_sim(1,2)*rlatt_sim(2,3)
!    &   -rlatt_sim(3,1)*rlatt_sim(2,2)*rlatt_sim(1,3)
!    &   -rlatt_sim(1,1)*rlatt_sim(3,2)*rlatt_sim(2,3)
!    &   -rlatt_sim(2,1)*rlatt_sim(1,2)*rlatt_sim(3,3)
      det1=twopi/det_sim
      glatt_sim(1,1)=det1*(rlatt_sim(2,2)*rlatt_sim(3,3)-rlatt_sim(2,3)*rlatt_sim(3,2))
      glatt_sim(2,1)=det1*(rlatt_sim(3,2)*rlatt_sim(1,3)-rlatt_sim(3,3)*rlatt_sim(1,2))
      glatt_sim(3,1)=det1*(rlatt_sim(1,2)*rlatt_sim(2,3)-rlatt_sim(1,3)*rlatt_sim(2,2))
      glatt_sim(1,2)=det1*(rlatt_sim(2,3)*rlatt_sim(3,1)-rlatt_sim(2,1)*rlatt_sim(3,3))
      glatt_sim(2,2)=det1*(rlatt_sim(3,3)*rlatt_sim(1,1)-rlatt_sim(3,1)*rlatt_sim(1,3))
      glatt_sim(3,2)=det1*(rlatt_sim(1,3)*rlatt_sim(2,1)-rlatt_sim(1,1)*rlatt_sim(2,3))
      glatt_sim(1,3)=det1*(rlatt_sim(2,1)*rlatt_sim(3,2)-rlatt_sim(2,2)*rlatt_sim(3,1))
      glatt_sim(2,3)=det1*(rlatt_sim(3,1)*rlatt_sim(1,2)-rlatt_sim(3,2)*rlatt_sim(1,1))
      glatt_sim(3,3)=det1*(rlatt_sim(1,1)*rlatt_sim(2,2)-rlatt_sim(1,2)*rlatt_sim(2,1))

      write(ounit,'(/,''Simulation cell reciprocal lattice basis vectors'',3(/,3f10.6))') &
       ((glatt_sim(k,j),k=1,3),j=1,3)

      vcell_sim=dabs(det_sim)
      write(ounit,'(/,''Simulation cell volume'',f14.8)') vcell_sim
      if((vcell_sim/vcell)-nint(vcell_sim/vcell).gt.1.d-9) then
        write(ounit,'(''Warning: vcell_sim/vcell='',f9.5, '' not an integer'')') vcell_sim/vcell
        call fatal_error ('Simulation cell volume is not a multiple of the primitive cell volume')
      endif

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
!     cutr=0.5d0*dist_min
!     cutr=min(rmax(ict),cutr)
      write(ounit,'(/,''Shortest distance to cell boundary'',f14.8)') dist_min/2

! simulation cell
      call short_distance(rlatt_sim,vcell_sim,dist_min,rdist_sim)
!     cutr_sim=0.5d0*dist_min
!      write(ounit,'(/,''Shortest distance to sim. cell boundary'',f14.8)') dist_min/2

! reciprocal-space distances
! primitive cell
      vgcell=twopi**3/vcell
      call short_distance(glatt,vgcell,gdistmin,gdist)
      write(ounit,'(/,''Shortest distance to recip. cell boundary'',f14.8)') gdistmin

! simulation cell
      vgcell_sim=twopi**3/vcell_sim
      call short_distance(glatt_sim,vgcell_sim,gdistmin_sim,gdist_sim)
!      write(ounit,'(/,''Shortest distance to sim. recip. cell boundary'',f14.8)') gdistmin_sim


! Estimate maximum shells number before build it

      ngvec=1
      do k=1,3
         ng1d(k)=int(cutg/gdist(k))
         ngvec=ngvec*(2*ng1d(k)+1)
         write(ounit,*)"Estimated ng1d given cutg",cutg, k, 2*ng1d(k)+1
      enddo
      write(ounit,*) "Estimated max ngvec", ngvec

      if(ngvec.gt.NGVECX) then
         write(ounit,*) "Compute a new cutg to lower  number og g-vectors"
         ngmx=7
         cutgv=gdistmin*ngmx
         cutg=MAXVAL(cutgv)
         write(ounit,*) "New cutg value", cutg

      endif


! generate shells of primitive cell g-vectors
      call shells(cutg_big,glatt,gdist,igvec,gvec,gnorm,igmult,ngvec_big, &
      ngnorm_big,ng1d,0)

      ngnorm=ngnorm_big
      ngvec=0
      do k=1,ngnorm_big
        if(gnorm(k).gt.cutg+eps) then
          ngnorm=k-1
          goto 20
        endif
        ngvec=ngvec+igmult(k)
      enddo
      write(ounit,*) "Recomputed after shells generation ngvec", ngvec

   20 if(ipr.ge.4) write(ounit,'(/,''Shells within cutg_big,cutg'',2i8)') ngnorm_big,ngnorm
      if(ipr.ge.4) write(ounit,'(/,''Vects. within cutg_big,cutg'',2i8)') ngvec_big,ngvec
      if(ipr.ge.4) write(ounit,'(/,''ng1d for primitive cell'',3i4)') (ng1d(k),k=1,3)
      if(ngvec.gt.NGVECX) then
        if(ipr.ge.4) write(ounit,'(''ngvec,NGVECX='',2i8)') ngvec,NGVECX
        call fatal_error ('ngvec>NGVECX in set_ewald')
      endif
      if(ngnorm.gt.NGNORMX) then
        if(ipr.ge.4) write(ounit,'(''ngnorm,NGNORMX='',2i8)') ngnorm,NGNORMX
        call fatal_error ('ngnorm>NGNORMX in set_ewald')
      endif
      do k=1,3
        if(ng1d(k).gt.NG1DX) then
          if(ipr.ge.4) write(ounit,'(''k,ng1d(k),NG1DX='',i1,2i8)') k,ng1d(k),NG1DX
          call fatal_error ('ng1d(k)>NG1DX in set_ewald')
        endif
      enddo

      open(1,file='gvectors_qmc')
      write(1,'(i5,'' ngvec (half of them only)'')') ngvec
      write(1,'(3i5,3f8.4)') ((igvec(k,i),k=1,3),(gvec(k,i),k=1,3),i=1,ngvec)
      close(1)

! generate shells of simulation cell g-vectors
      call shells(cutg_sim_big,glatt_sim,gdist_sim,igvec_sim,gvec_sim,gnorm_sim,igmult_sim,ngvec_sim_big, &
       ngnorm_sim_big,ng1d_sim,1)

      write(ounit,*) "Number of gvectors simulation cell: ", ngvec_sim, ngvec_sim_big
      ngnorm_sim=ngnorm_sim_big
      ngvec_sim=0
      do k=1,ngnorm_sim_big
        if(gnorm_sim(k).gt.cutg_sim+eps) then
          ngnorm_sim=k-1
          goto 50
        endif
        ngvec_sim=ngvec_sim+igmult_sim(k)
      enddo

   50 if(ipr.ge.4) then
         write(ounit,'(/,''Shells within cutg_sim_big,cutg_sim'',2i8)') ngnorm_sim_big,ngnorm_sim
         write(ounit,'(/,''Vects. within cutg_sim_big,cutg_sim'',2i8)') ngvec_sim_big,ngvec_sim
         write(ounit,'(/,''ng1d for simulation cell'',3i4)') (ng1d_sim(k),k=1,3)
      endif
      if(ngvec_sim.gt.NGVEC_SIMX) call fatal_error ('ngvec_sim>NGVEC_SIMX in set_ewald')
      if(ngnorm_sim.gt.NGNORM_SIMX) call fatal_error ('ngnorm_sim>NGNORM_SIMX in set_ewald')
      do k=1,3
        if(ng1d_sim(k).gt.NG1DX) call fatal_error ('ng1d_sim(k)>NG1DX in shells')
      enddo

! Convert k-vector shift from simulation-cell recip. lattice vector units to cartesian coordinates
      do k=1,3
         rkvec_shift_tmp(k)=rkvec_shift(k)
      enddo
      do k=1,3
         rkvec_shift(k)=0
         do i=1,3
            rkvec_shift(k)=rkvec_shift(k)+rkvec_shift_tmp(i)*glatt_sim(k,i)
         enddo
      enddo
      if(ipr.ge.4) then
         write(ounit,'(/,''rkvec_shift in sim-cell recip. lat. vec. units'',9f9.4)') (rkvec_shift_tmp(k),k=1,3)
         write(ounit,'(''rkvec_shift in cartesian coodinates'',9f9.4)') (rkvec_shift(k),k=1,3)
      endif
!     Generate k-vectors, i.e. simulation-cell recip. lattice vectors shiftedby rkvec_shift that are
! not related by a primitive-cell recip. lattice vector.
      call k_vectors



! Coulomb interactions in primitive and simulation cells
      lowest_pow=-1
      b0=1.d0
      ifcon=1
      isrange=0



! n-n, e-n interactions (primitive cell)
! put in uniform background by setting k=0 term to zero
      vbare_coul(1)=0.d0
      do k=2,ngnorm_big
!     Fourier transfom of 1/r
         vbare_coul(k)=2*twopi/(vcell*gnorm(k)**2)
      enddo



      if(ipr.ge.4) write(ounit,'(''Compare Coulomb reconstructed'')')
      if(ipr.ge.4) write(ounit,'(''1/r  vs vps_rest  diff'')')
      r_fk=0.d0
      ict=1
      write(ounit,*) "nr_ps(ict)", nr_ps(ict)
      do ir=2,nr_ps(ict),100
         r_fk(1)=r0(ict)*(arg(ict)**(ir-1)-1.d0)
!         write(ounit,*) "ir",ir,"r_fk(ir)",r_fk(ir)
         vps_test=0.d0
         r_fk(1)=max(1.0d-10,r_fk(1))
         call fourier_it(r_fk,gvec,ngnorm_big,igmult,vbare_coul, vps_test)
!     write(ounit, *) ir, 1/r_fk(1), vps_test
         if(ipr.ge.4)       write(ounit,*) ir, r_fk(1), 1.d0/r_fk(1), vps_test, (1.d0/r_fk(1))-vps_test
      enddo
      if(ipr.ge.4) write(ounit,'(''Finnish check Coulomb reconstruted'')')



!      write(ounit,*) "parameters for separate"
!      write(ounit,*) 'ngnorm_big',ngnorm_big
!      write(ounit,*) 'igmult',igmult
!      write(ounit,*) "gnorm",gnorm
!      write(ounit,*) "ngnorm",ngnorm
!      write(ounit,*) "cutr",cutr
!      write(ounit,*) "ncoef_per",ncoef_per


      if(ipr.ge.4) then
         write(ounit,'(''start vbare_coul separate primitive'')')
         write(ounit,'(''check vbare_coul = '',20d12.4)') (vbare_coul(k),k=1,ngnorm)
      endif
      call separate(vbare_coul,b0,lowest_pow,ngnorm_big,igmult,gnorm, &
       ngnorm,cutr,vcell,ncoef_per,np,b_coul,y_coul,chisq,ifcon,isrange)
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



! debug n-n interaction (primitive cell)
      if(ipr.ge.0) then
        write(ounit, &
       '(''      r       "true"      ewald       test        test-true    1/r     d_true  d_test   vsrange   vlrange'')')
        lowest_pow=-1
        npts=101
        dx=cutr/(npts-1)

        sum=0
        do i=1,npts
          rr=(i-1)*dx+1.d-20
          if(i.eq.1.or.i.eq.npts) then
            wt=0.5d0
           else
            wt=1
          endif
          sum=sum+wt*rr**2*vsrange(rr,cutr,lowest_pow,ncoef_per,np,b_coul)
        enddo
        constv=2*twopi*sum*dx/vcell
        write(ounit,'(''const='',9f12.8)') constv

        rms=0.d0
        do i=1,npts
          r_tmp(1)=(i-1)*dx+1.d-20
          r_tmp(2)=0.d0
          r_tmp(3)=0.d0
          rr=sqrt(r_tmp(1)**2+r_tmp(2)**2+r_tmp(3)**2)
          vs=vsrange(rr,cutr,lowest_pow,ncoef_per,np,b_coul)
          vl=vlrange_old(r_tmp,gvec,ngnorm,igmult,y_coul)
          testv=vs+vl
!         true=vlrange_old(r_tmp,gvec,ngnorm_big,igmult,vbare_coul)
          true=ewald_pot(r_tmp,rr,gvec,gnorm,ngnorm_big,igmult,vbare_coul,cutr,vcell)
          ewa=ewald_pot(r_tmp,rr,gvec,gnorm,ngnorm,igmult,vbare_coul,cutr,vcell)
!          short_ewald=derfc((5/cutr)*rr)/rr
          if(i.ne.1) rms=rms+(true-testv)**2
          if(i.eq.1) then
!             write(ounit,'(''1/r'',f8.4,4f12.6,30x,2f12.6)') rr,true,ewa,short_ewald,testv,testv-true,vs,vl
            write(ounit,'(''1/r'',f8.4,5f12.6,30x,2f12.6)') rr,true,ewa,testv,testv-true,vs,vl
          else
!            write(ounit,'(''1/r'',f8.4,6f12.6,2f9.3,2f12.6)') rr,true,ewa,short_ewald,testv,testv-true
!     &      ,1/rr,(true-true_s)/dx,(testv-test_s)/dx,vs,vl
!     write(ounit,'(''1/r'',f8.4,5f12.6,2f9.3,2f12.6)') rr,true,ewa,short_ewald,testv,testv-true
             write(ounit,'(''1/r'',f8.4,5f12.6,2f9.3,2f12.6)') rr,true,ewa,testv,testv-true &
            ,1/rr,(true-true_s)/dx,(testv-test_s)/dx,vs,vl
          endif
          true_s=true
        test_s=testv
        enddo
        rms=sqrt(rms/(npts-1))
        write(ounit,'(''Rms error of 1/r (prim cell) fit ='',d12.4)') rms
        if(rms.gt.1.d-3) write(ounit,'(''Warning: rms error of 1/r (prim cell) fit too large'',d12.4)') rms
        write(ounit,*) "End n-n Coulumb potential separation"
      endif




! e-e interactions (simulation cell) (we can reuse vbare_coul)
! put in uniform background by setting k=0 term to zero
      vbare_coul(1)=0.d0
      vbare_jas(1)=0.d0
      do k=2,ngnorm_sim_big
! Fourier transfom of 1/r
        vbare_coul(k)=2*twopi/(vcell_sim*gnorm_sim(k)**2)
! Fourier transform of -1/r*(1-exp(-r/f)) for Jastrow
        vbare_jas(k)=-vbare_coul(k)/(1+(f*gnorm_sim(k))**2)
      enddo

      if(ipr.ge.4) write(ounit,'(''start vbare_coul separate sim cell'')')
      call separate(vbare_coul,b0,lowest_pow,ngnorm_sim_big,igmult_sim,gnorm_sim,ngnorm_sim &
      ,cutr_sim,vcell_sim,ncoef_per,np,b_coul_sim,y_coul_sim,chisq,ifcon,isrange)
      if(ipr.ge.4) write(ounit,'(''finish vbare_coul separate sim cell'')')




      if(chisq.gt.0) then
        write(ounit,'(''Rms error in 1/r separation in simulation cell'',d12.5)') dsqrt(chisq)
       else
        write(ounit,'(''Warning: Rms error missing, chisq negative in 1/r simulation separate'',d12.4)') chisq
        if(chisq.lt.0.d0) call fatal_error ('chisq<0 in separate')
      endif

      if(ipr.ge.4) then
         write(ounit,'(/,''Separation of Coulomb interaction in simulation cell'')')
         write(ounit,'(''vbare_coul = '',20d12.4)') (vbare_coul(k),k=1,ngnorm_sim)
         write(ounit,'(''vbare_coul = '',20d12.4)') (vbare_coul(k),k=1,ngnorm_sim_big)
         write(ounit,'(''y_coul_sim = '',20d12.4)') (y_coul_sim(k),k=1,ngnorm_sim)
         write(ounit,'(''b_coul_sim = '',20d12.4)') (b_coul_sim(k),k=1,ncoef_per)
      endif


! debug e-e interaction (simulation cell)
! Note vbare_coul is used both for primitive and simulation cells
! Since Jastrow has singlularity at 0, cannot match there, so evaluate
! rms error only for latter 3/4 of interval
      if(ipr.ge.0) then
         write(ounit,*) "Test e-e coulumb separation"
         write(ounit,*) " "
         write(ounit, &
       '(''      r       "true"      ewald       test        test-true    1/r     d_true  d_test   vsrange   vlrange'')')
        lowest_pow=-1
        npts=101
        dx=cutr_sim/(npts-1)

        sum=0
        do i=1,npts
          rr=(i-1)*dx+1.d-20
          if(i.eq.1.or.i.eq.npts) then
            wt=0.5d0
           else
            wt=1
          endif
          sum=sum+wt*rr**2*vsrange(rr,cutr_sim,lowest_pow,ncoef_per,np,b_coul_sim)
        enddo
        constv=2*twopi*sum*dx/vcell_sim
        write(ounit,'(''const='',9f12.8)') constv

        rms=0.d0
        do i=1,npts
          r_tmp(1)=(i-1)*dx+1.d-20
          r_tmp(2)=0.d0
          r_tmp(3)=0.d0
          rr=sqrt(r_tmp(1)**2+r_tmp(2)**2+r_tmp(3)**2)
          vs=vsrange(rr,cutr_sim,lowest_pow,ncoef_per,np,b_coul_sim)
          vl=vlrange_old(r_tmp,gvec_sim,ngnorm_sim,igmult_sim,y_coul_sim)
          testv=vs+vl
!         true=vlrange_old(r_tmp,gvec_sim,ngnorm_sim_big,igmult_sim,vbare_coul)
          true=ewald_pot(r_tmp,rr,gvec_sim,gnorm_sim,ngnorm_sim_big,igmult_sim,vbare_coul,cutr_sim,vcell_sim)
          ewa=ewald_pot(r_tmp,rr,gvec_sim,gnorm_sim,ngnorm_sim,igmult_sim,vbare_coul,cutr_sim,vcell_sim)
          if(i.ne.1) rms=rms+(true-testv)**2
          if(i.eq.1) then
            write(ounit,'(''1/r'',f8.4,4f12.6,30x,2f12.6)') rr,true,ewa,testv,testv-true,vs,vl
           else
            write(ounit,'(''1/r'',f8.4,5f12.6,2f9.3,2f12.6)') rr,true,ewa,testv,testv-true &
            ,1/rr,(true-true_s)/dx,(testv-test_s)/dx,vs,vl
          endif
          true_s=true
        test_s=testv
        enddo
        rms=sqrt(rms/(npts-1))
        write(ounit,'(''Rms error of 1/r fit (sim cell) ='',d12.4)') rms
        if(rms.gt.1.d-3) write(ounit,'(''Warning: rms error of 1/r (sim cell) fit too large'',d12.4)') rms
        write(ounit,*) "End e-e coulumb separation"
      endif





!     e-e Jastrow
!      write(ounit,*) "nloc before computing e-ion local part", nloc
      if(iperiodic.eq.1) goto 88
      lowest_pow=0
      b0=0.5d0/f**2
      ifcon=1
      isrange=0

!      write(ounit,*) "end separate 00005"


      if(ipr.ge.4.) write(ounit,'(''start vbare_jas separate'')')
      call separate(vbare_jas,b0,lowest_pow,ngnorm_sim_big,igmult_sim, &
           gnorm_sim, ngnorm_sim, cutr_sim,vcell_sim,ncoef_per,np, &
           b_jas,y_jas,chisq,ifcon,isrange)
      if(ipr.ge.4.) write(ounit,'(''finnish vbare_jas separate'')')

      if(chisq.gt.0) then
        write(ounit,'(''Rms error in Jastrow separation'',d12.5)') dsqrt(chisq)
       else
        write(ounit,'(''Warning: Rms error missing, chisq negative in Jastrow separate'',d12.4)') chisq
        if(chisq.lt.0.d0) call fatal_error ('chisq<0 in separate')
      endif

      if(ipr.ge.4) then
         write(ounit,'(/,''Separation of Jastrow in simulation cell'')')
         write(ounit,'(''vbare_jas = '',20d12.4)') (vbare_jas(k),k=1,ngnorm_sim)
         write(ounit,'(''vbare_jas = '',20d12.4)') (vbare_jas(k),k=1,ngnorm_sim_big)
         write(ounit,'(''y_jas = '',20d12.4)') (y_jas(k),k=1,ngnorm_sim)
         write(ounit,'(''b_jas = '',20d12.4)') (b_jas(k),k=1,ncoef_per)
      endif

! debug e-e Jastrow
! Since Jastrow has singlularity at 0, cannot match there, so evaluate
! rms error only for latter 3/4 of interval
      if(ipr.ge.0) then
         write(ounit,'(/,''Debugging e-e Jastrow separation'')')
         write(ounit,'(''      r       "true"       test      test-true -1/r*(1-exp(-r/f) d_true d_test'')')
         lowest_pow=0
         npts=101
         dx=cutr_sim/(npts-1)
         rms=0.d0
         do i=1,npts
            r_tmp(1)=(i-1)*dx+1.d-20
            r_tmp(2)=0.d0
            r_tmp(3)=0.d0
            rr=sqrt(r_tmp(1)**2+r_tmp(2)**2+r_tmp(3)**2)
            testv=vsrange(rr,cutr_sim,lowest_pow,ncoef_per,np,b_jas)
            testv=testv+vlrange_old(r_tmp,gvec_sim,ngnorm_sim,igmult_sim,y_jas)
            true=vlrange_old(r_tmp,gvec_sim,ngnorm_sim_big,igmult_sim,vbare_jas)
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


! e-ion local pseudopotential

! One cannot numerically fourier transform a function with a long 1/r
! tail.  So, add in potential due to gaussian distribution, that cancels
! 1/r tail and subtract out the fourier components of this gaussian
! distribution after doing the numerical fourier transform.


      if(nloc.eq.0) goto 197


 88   if(ipr.ge.4) write(ounit,*) "compute vps_gauss and its fourier transform for range separation"
      do ict=1,nctype


! If this is done just to find fourier components of potential, the cut-off
!     radius can be rmax(ict), but if we are to use it also for one of the basis
!     functions for the optimal separation, then the cutoff has to be cutr.
!     alpha=5/rmax(ict)
         alpha=5/cutr
         r(1)=1.d-10
         vpot=0.d0
         dvpot=0.d0
         write(ounit,*) "lpot",lpot(ict)
         call gauss_pot(r(1),lpot(ict),ict,vpot,dvpot)
         vps_short(1)=vpot-znuc(ict)/r(1)+znuc(ict)*2*alpha/sqrt(pi)
         if(ipr.ge.4) then
            write(ounit,'(''alpha'',20d12.4)') alpha
            write(ounit,'(''arg'',20d12.4)') arg(ict)
         endif
!     write(ounit,'(''lalala CHECK rg, vpot, vps_short'',3d12.4)') r(ir), vpot, vps_short(1)

         do ir=2,nr_ps(ict)
!     write(ounit,*) 'ir',ir
            r(ir)=r0(ict)*(arg(ict)**(ir-1)-1.d0)
            vpot=0.d0
            dvpot=0.d0
            rg=max(1.0d-10,r(ir))
            call gauss_pot(rg,lpot(ict),ict,vpot,dvpot)
            vps_short(ir)=vpot-(znuc(ict)/rg)+(znuc(ict)*derf(alpha*rg)/rg)
!     write(ounit,'(''here CHECK ir rg, vpot, vps_short'',i20, 3d12.4)') ir, r(ir), vpot, vps_short(ir)
         enddo

         if(ipr.ge.4) write(ounit,'(/,''Grid parameters, r0,rmax,arg  '', 3d13.5, i2)') &
        r0(ict),r(nr_ps(ict)),arg(ict),nr_ps(ict)

      if(gnorm(ngnorm_big).gt.twopi/(100*(r(nr_ps(ict)-1)+r0(ict))*(arg(ict)-1))) then
        write(ounit,'(''**Warning: not enough grid pts in psp to do accurate numerical FT; &
        using a non-uniform grid only makes it worse'')')
        write(ounit,'(''gnorm(ngnorm_big)>twopi/(100*(r(nr_ps(ict)-1)+r0(ict))*(arg(ict)-1))'',9f9.4)') &
        gnorm(ngnorm_big),twopi/(100*(r(nr_ps(ict)-1)+r0(ict))*(arg(ict)-1)) &
        ,twopi/(100*r(nr_ps(ict))*(arg(ict)-1))
      endif

!      write(ounit,'(''CHECK vps_short'',d12.4)') vps_short(1:nr_ps(ict))

      call fourier_transform(r,arg(ict),r0(ict),nr_ps(ict),vps_short,vcell,gnorm,ngnorm_big, &
       vps_basis_fourier)

      if(ipr.ge.4) then
         write(ounit,'(''Start check vps_basis_fourier from vps_short'')')
         do ig=1,ngnorm_big,10
            write(ounit,'(''CHECK FT'',i5,g12.5,d12.4)') ig,gnorm(ig),vps_basis_fourier(ig)
         enddo
         write(ounit,'(''Finnish check vps_basis_fourier from vps_short'')')
      endif

      if(ipr.ge.4) then
         write(ounit,'(''Compare vps_short(ir) vs vps_shortf reconstructed'')')
         write(ounit,'(''vps_short   vps_shortf  diff'')')
      endif
      r_fk=0.d0
      do ir=2,nr_ps(ict)
         r_fk(1)=r0(ict)*(arg(ict)**(ir-1)-1.d0)
         vps_test=0.d0
         r_fk(1)=max(1.0d-10,r_fk(1))
         call fourier_it(r_fk,gvec,ngnorm_big,igmult,vps_basis_fourier, vps_test)
         if(ipr.ge.4) write(ounit,*) r_fk(1), vps_short(ir), vps_test, vps_short(ir)-vps_test
      enddo
      if(ipr.ge.4)       write(ounit,'(''finnish check vps_short reconstruted'')')




      if(ipr.ge.4) write(ounit,*) "Add in (4*pi*Z/V) to get correct background contribution"
! Add in (4*pi*Z/V) Integ d^3r (1-erf(alpha*r))/r to get correct background contribution
      vbare_psp(1)=vps_basis_fourier(1)+pi*znuc(ict)/(vcell*alpha**2)
      do ig=2,ngnorm_big
         g2=gnorm(ig)**2
         g2a=0.25d0*g2/alpha**2
         vbare_psp(ig)=vps_basis_fourier(ig)-2*twopi*znuc(ict)*exp(-g2a)/(vcell*g2)
      enddo

      if(ipr.ge.4) then
         do ig=1,ngnorm_big,10
            write(ounit,'(''CHECK FT2'',i5,g12.5,9d12.4)') ig,gnorm(ig),vbare_psp(ig),vps_basis_fourier(ig)
         enddo
      endif

! One cannot numerically fourier transform a function with a long 1/r
! tail.  So, add Z/r potential to cancel -Z/r tail and subtract out the
! fourier components of Z/r after doing the numerical fourier transform.
! It is OK that Z/r diverges at r, since it gets multiplied by r^2 when
! doing the FT.

!     do 140 ir=1,nr_ps(ict)
!     r(ir)=r0(ict)*(arg(ict)**(ir-1)-1.d0)
!     if(ir.eq.1) r(ir)=1.d-10
!     140   vps_short(ir)=vpseudo(ir,ict,lpot(ict))+znuc(ict)/r(ir)

!     do 150 ir=1,nr_ps(ict),50
! 150   write(ounit,'(''CHECK SR'',i5,g12.5,9d12.4)') ir,r(ir),vpseudo(ir,ict,lpot(ict)),vps_short(ir)

!     write(ounit,'(''ict='',9i5)') ict
!     write(ounit,'(''ict,nr_ps(ict)='',9i5)') ict,nr_ps(ict)
!     write(ounit,'(/,''Grid parameters, r0,rmax,arg  '',d10.4,2f10.5)')
!    & r0(ict),r(nr_ps(ict)),arg(ict)

!c    if(gnorm(ngnorm_big).gt.twopi/(100*(r(nr_ps(ict)-1)+r0(ict))*(arg(ict)-1))) then
!c      write(ounit,'(''**Warning: not enough grid pts in psp to do accurate numerical FT;
!c   &  using a non-uniform grid only makes it worse'')')
!c      write(ounit,'(''gnorm(ngnorm_big)>twopi/(100*(r(nr_ps(ict)-1)+r0(ict))*(arg(ict)-1))'',9f9.4)')
!c   &  gnorm(ngnorm_big),twopi/(100*(r(nr_ps(ict)-1)+r0(ict))*(arg(ict)-1))
!c   &  ,twopi/(100*r(nr_ps(ict))*(arg(ict)-1))
!c    endif

!     call fourier_transform(r,arg(ict),r0(ict),nr_ps(ict),vps_short,vcell,gnorm,ngnorm_big,
!    & vbare_psp)

!     do 160 ig=1,ngnorm_big,10
! 160   write(ounit,'(''CHECK FT'',i5,g12.5,d12.4)') ig,gnorm(ig),vbare_psp(ig)

!cSubtract out fourier components of Z/r
!cFourier coefs of Z/r or psp+Z/r go as 1/g^2, those of psp go as 1/g^5 because of deriv. discontinuity in TM psp
!     do 170 ig=2,ngnorm_big
! 170   vbare_psp(ig)=vbare_psp(ig)-2*twopi*znuc(ict)/(vcell*gnorm(ig)**2)

!     do 180 ig=1,ngnorm_big,10
! 180   write(ounit,'(''CHECK FT2'',i5,g12.5,d12.4)') ig,gnorm(ig),vbare_psp(ig)

! If isrange=0,1 set ifcon=1, if isrange=2,3 set ifcon=0
      lowest_pow=0
      b0=0
      ifcon=1
      isrange=1


      if(ipr.ge.4) then
         write(ounit,'(''start vbare_psp separate short range potential'')')
         write(ounit,'(''vbare_psp'', 20d12.4)') vbare_psp(1:ngnorm)
      endif
      call separate(vbare_psp,b0,lowest_pow,ngnorm_big,igmult,gnorm, &
           ngnorm,cutr,vcell,ncoef_per,np,b_psp(1,ict),y_psp(1,ict),chisq, &
           ifcon,isrange)


      if(ipr.ge.4) write(ounit,'(''finnish vbare_psp separate short range potential'')')



      if(chisq.gt.0) then
        write(ounit,'(''Rms error in pseudopotential separation'',d12.5)') dsqrt(chisq)
      else
         write(ounit,'(''Warning: Rms error missing, chisq negative in pseudopotential separate'',d12.4)') chisq
         if(chisq.lt.0.d0) call fatal_error ('chisq<0 in separate')
      endif

      if(ipr.ge.4) then
         write(ounit,'(/,''Separation of pseudopotential in primitive cell'')')
         write(ounit,'(''vbare_psp = '',20d12.4)') (vbare_psp(k),k=1,ngnorm)
         write(ounit,'(''vbare_psp = '',20d12.4)') (vbare_psp(k),k=1,ngnorm_big)
         write(ounit,'(''y_psp = '',20d12.4)') (y_psp(k,ict),k=1,ngnorm)
         write(ounit,'(''b_psp = '',20d12.4)') (b_psp(k,ict),k=1,ncoef_per)
      endif

! If sim cell is not primitive cell, vbare_coul has been overwritten, so restore it
! n-n, e-n interactions (primitive cell)
! put in uniform background by setting k=0 term to zero
      if(vcell_sim.ne.vcell) then
         vbare_coul(1)=0.d0
         do k=2,ngnorm_big
!     Fourier transfom of 1/r
            vbare_coul(k)=2*twopi/(vcell*gnorm(k)**2)
         enddo
      endif

! Note that 1/r and Jastrow have singlularities at r=0 and so at short r we cannot test
! the separation into short-range and long-range parts by just summing the
! bare long-range parts.  The psp does not have a singularity so we can test it everywhere.
      if(ipr.ge.0) then
         write(ounit,*)'testing vsrange ',isrange
        write(ounit,'(''      r       "true"      ewald       test        test-true   d_true  d_test   vsrange  vlrange'')')
!         write(ounit,'(''      r       "true"      ewald       test        test-true   d_true  d_test   vsrange  vlrange  w(r)'')')
        npts=101
        dx=cutr/(npts-1)
        rms=0.d0
        write(ounit,*)'when computing short range separate range is ',isrange

!        print*, 'cutr',cutr,'dx',dx
!        write(ounit,'(''start loop vsrange1'')')
!        write(ounit,'(''npts'',i20)') npts
        do i=1,npts
          r_tmp(1)=(i-1)*dx+1.d-20
          r_tmp(2)=0.d0
          r_tmp(3)=0.d0
          rr=sqrt(r_tmp(1)**2+r_tmp(2)**2+r_tmp(3)**2)
!          print*,'r_tmp',r_tmp
!          print*,'rr',rr

!     write(ounit,'(''i'',i20)') i
          if(isrange.eq.0) vs=vsrange(rr,cutr,lowest_pow,ncoef_per,np,b_psp(1,ict))
          if(isrange.eq.1) vs=vsrange1(rr,cutr,lowest_pow,ncoef_per,np,b_psp(1,ict),ict,lpot(ict))
          if(isrange.eq.2) vs=vsrange2(rr,cutr,lowest_pow,ncoef_per,np,b_psp(1,ict),ict,lpot(ict))
          if(isrange.eq.3) vs=vsrange3(rr,cutr,lowest_pow,ncoef_per,np,b_psp(1,ict),ict,lpot(ict))
!     write(ounit,'(''getting close'')')
          vl=vlrange_old(r_tmp,gvec,ngnorm,igmult,y_psp(1,ict))
          testv=vs+vl
!     true=vlrange_old(r_tmp,gvec,ngnorm_big,igmult,vbare_psp)
!     to compute w_approach in this case not sure if I should add corrections
!     call wf(rr,w_gauss(i),b_psp(1:ncoef_per,ict))
          if(rr.gt.cutr) write(ounit,*) "rr",rr,"pass the cutr limit", cutr
!     if(rr.lt.cutr) then
          true=ewald_pot_psp(r_tmp,rr,gvec,gnorm,ngnorm_big,igmult,vbare_coul,cutr,vcell,ict,lpot(ict),znuc(ict))
          ewa=ewald_pot_psp(r_tmp,rr,gvec,gnorm,ngnorm,igmult,vbare_coul,cutr,vcell,ict,lpot(ict),znuc(ict))
!          else
!             true=ewald_pot(r_tmp,rr,gvec,gnorm,ngnorm_big,igmult,vbare_coul,cutr,vcell)
!             ewa=ewald_pot(r_tmp,rr,gvec,gnorm,ngnorm,igmult,vbare_coul,cutr,vcell)
!          endif

          rms=rms+(true-testv)**2


          if(ipr.ge.4) then
             if(i.eq.1) then
                write(ounit,'(''rr'',f8.4,4f12.6,16x,9f12.6)') rr,true,ewa,testv,testv-true,vs,vl
!     write(ounit,'(''vps'',f8.4,4f12.6,16x,9f12.6)') rr,true,ewa,testv,testv-true,vs,vl, w_gauss(i)
             else
                write(ounit,'(''rr'',f8.4,4f12.6,2f8.3,9f12.6)') rr,true,ewa,testv,testv-true &
                     ,(true-true_s)/dx,(testv-test_s)/dx,vs,vl
!     write(ounit,'(''rr'',f8.4,4f12.6,2f8.3,9f12.6)') rr,true,ewa,testv,testv-true
!     &      ,(true-true_s)/dx,(testv-test_s)/dx,vs,vl, w_gauss(i)
             endif
          endif
          true_s=true
          test_s=testv
          if(ipr.ge.4) write(ounit,*) "rms current", dsqrt(rms/i)
       enddo
       if(ipr.ge.4) write(ounit,'(''Finish vsrange separatation'')')
       rms=sqrt(rms/npts)
       write(ounit,'(''Rms error of psp fit='',d12.4)') rms
       if(rms.gt.1.d-3) write(ounit,'(''Warning: rms error of psp fit too large'',d12.4)') rms
      endif
      enddo

      197 znuc_sum=0
      znuc2_sum=0
      do i=1,ncent
        znuc_sum=znuc_sum+znuc(iwctype(i))
        znuc2_sum=znuc2_sum+znuc(iwctype(i))**2
      enddo

!     call pot_nn_ewald_old
!     c_madelung=pecent*dist_nn/(znuc(1)*znuc(2)*ncent/2)
!     write(ounit,'(''pecent (old)='',f10.6)') pecent
!     write(ounit,'(''c_madelung_o='',f10.6)') c_madelung

      call pot_nn_ewald
!     c_madelung=pecent*dist_nn/(znuc(1)*znuc(2)*ncent/2)
      write(ounit,'(''pecent='',f12.6)') pecent
!     write(ounit,'(''c_madelung='',f10.6)') c_madelung

! images for periodic basis functions
      if(n_images.ge.1)then
!     n_images=26

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
                     if(ix.ne.0) ell(1,imcount)=1.d0*ix*rlatt_sim(1,1)
                     if(iy.ne.0) ell(2,imcount)=1.d0*iy*rlatt_sim(2,2)
                     if(iz.ne.0) ell(3,imcount)=1.d0*iz*rlatt_sim(3,3)
                  endif
               enddo
            enddo
         enddo
         if (ipr.ge.4) write(ounit,*) "Total number of peridic images counted: ",imcount

      endif

      return
      end
!-----------------------------------------------------------------------

      subroutine short_distance(vector,volume,dist_min,distcell)
! Written by Cyrus Umrigar
! distcell(i) is the perpendicular distance between cell faces parallel
! to the other 2 directions from i.
! dist_min is the shortest of these three.
! By choosing the range of the short-range part of the Ewald sums to be
! <= half the shortest perpendicular distance we ensure that the short-range
! part has zero or one terms.
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

      subroutine cross(v1,v2,v3)
! evaluates the cross-product of v1 and v2 and puts it in v3

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

      subroutine shells(cutg,glatt,gdist,igvec,gvec,gnorm,igmult,ngvec_big,ngnorm_big,ng1d,icell)
      use ewald_mod, only: NGVEC_BIGX
      use ewald_mod, only: NGVEC_SIM_BIGX
! Written by Cyrus Umrigar

! icell = 0  primitive cell
!         1  simulation cell

      use precision_kinds, only: dp
      use error, only: fatal_error
      implicit none

      integer :: i1, i2, i2min, i3, i3min
      integer :: icell, k, ngnorm_big, ngvec_big
      integer, dimension(3,*) :: igvec
      integer, dimension(*) :: igmult
      integer, dimension(*) :: ng1d
      real(dp) :: cutg, cutg2, glen2, gx, gy
      real(dp) :: gz
      real(dp), dimension(3,*) :: glatt
      real(dp), dimension(3) :: gdist
      real(dp), dimension(3,*) :: gvec
      real(dp), dimension(*) :: gnorm
      real(dp), dimension(NGVEC_SIM_BIGX) :: gnorm_tmp



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
              if(icell.eq.0 .and. ngvec_big.gt.NGVEC_BIGX) then
                call fatal_error ('ngvec_big > NGVEC_BIGX in shells')
               elseif(icell.eq.1 .and. ngvec_big.gt.NGVEC_SIM_BIGX) then
                call fatal_error ('ngvec_big > NGVEC_SIM_BIGX in shells')
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

      call sort(igvec,gvec,gnorm_tmp,gnorm,igmult,ngvec_big,ngnorm_big,icell)

      return
      end
!-----------------------------------------------------------------------

      subroutine sort(igvec,gvec,gnorm_tmp,gnorm,igmult,ngvec_big,ngnorm_big,icell)
      use ewald_mod, only: NGNORM_BIGX
      use ewald_mod, only: NGNORM_SIM_BIGX
! Written by Cyrus Umrigar
      use contrl_file,    only: ounit
      use control, only: ipr
      use error, only: fatal_error
      use precision_kinds, only: dp
      implicit none

      integer :: i, icell, icheck, icount, it
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
          if(icell.eq.0 .and. ngnorm_big.gt.NGNORM_BIGX) then
            if(ipr.ge.4) write(ounit,'(''ngnorm_big='',i8)') ngnorm_big
            call fatal_error ('ngnorm_big > NGNORM_BIGX in sort')
           elseif(icell.eq.1 .and. ngnorm_big.gt.NGNORM_SIM_BIGX) then
            if(ipr.ge.4) write(ounit,'(''ngnorm_sim_big='',i8)') ngnorm_big
            call fatal_error ('ngnorm_big > NGNORM_SIM_BIGX in sort')
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

!     j=0
!     do 100 i=1,ngnorm_big
!       do 100 im=1,igmult(i)
!         j=j+1
! 100     write(ounit,'(''CHECK '',2i4,9f10.4)')
!    &    i,igmult(i),gnorm(i),(gvec(k,j),k=1,3)
!     write(ounit,*)

      return
      end
!-----------------------------------------------------------------------

      subroutine k_vectors
! Written by Cyrus Umrigar
! Generate the unique k-vectors, i.e., those that are not related by a
! primitive cell reciprocal lattice vector and are not the inverses of
! existing vectors.  Note that for the moment we keep vectors that are
! related by primitive cell reciprocal lattice vectors to inverses of
! other vectors.  We should come back to the issue of whether that is
! a symmetry one could use later on.

      use ewald_mod, only: NSYM
      use periodic, only: cutg, cutg_sim
      use periodic, only: glatt_inv, gvec, gvec_sim, igvec_sim
      use periodic, only: k_inv, kvec
      use periodic, only: ngvec, ngvec_sim, nkvec
      use periodic, only: rknorm, rkvec, rkvec_shift, vcell
      use periodic, only: vcell_sim
      use precision_kinds, only: dp
      use contrl_file,    only: ounit
      use control, only: ipr
      use error, only: fatal_error
      implicit none

      integer :: i, ikv, j, k, l
      integer :: nkvec_tot
      real(dp) :: rnorm
      real(dp), dimension(3) :: rkvec_try
      real(dp), dimension(3) :: rkvec_latt
      real(dp), parameter :: eps = 1.d-6





      k_inv(1)=1
      do k=1,3
        kvec(k,1)=0
        rkvec(k,1)=rkvec_shift(k)
      enddo
!     write(ounit,'(''k-vec( 1)='',3i5,3f9.4)') (kvec(k,1),k=1,3),(rkvec(k,1),k=1,3)

      nkvec=0
! Warning: Need to think more about do loop limit
!     do 120 i=1,min(ngvec_sim,vcell_sim*NSYM/vcell)
      do i=1,min(ngvec_sim,nint(8*vcell_sim/vcell))
        do k=1,3
      20     rkvec_try(k)=rkvec_shift(k)+gvec_sim(k,i)
        enddo
! Check if after translation by primitive cell reciprocal lattice vector it is
! the same as an existing k-vector
        do j=2,ngvec
          do l=1,nkvec
            rnorm=0
            do k=1,3
              rnorm=rnorm+(rkvec_try(k)-gvec(k,j)-rkvec(k,l))**2
            enddo
            if(rnorm.lt.eps) goto 120
            rnorm=0
            do k=1,3
              rnorm=rnorm+(rkvec_try(k)+gvec(k,j)-rkvec(k,l))**2
            enddo
            if(rnorm.lt.eps) goto 120
      50   continue
          enddo
        enddo
! Check if after translation by primitive cell reciprocal lattice vector it is
! the inverse of an existing k-vector
        do j=2,ngvec
          do l=1,nkvec
            rnorm=0
            do k=1,3
              rnorm=rnorm+(rkvec_try(k)-gvec(k,j)+rkvec(k,l))**2
            enddo
            if(rnorm.lt.eps) then
              k_inv(l)=2
              goto 120
            endif
            rnorm=0
            do k=1,3
              rnorm=rnorm+(rkvec_try(k)+gvec(k,j)+rkvec(k,l))**2
            enddo
            if(rnorm.lt.eps) then
              k_inv(l)=2
              goto 120
            endif
          enddo
        enddo
! Voila, found a new one
        nkvec=nkvec+1
        k_inv(nkvec)=1
        rknorm(nkvec)=0
        do k=1,3
          kvec(k,nkvec)=igvec_sim(k,i)
          rkvec(k,nkvec)=rkvec_try(k)
          rknorm(nkvec)=rknorm(nkvec)+rkvec_try(k)**2
        enddo
        rknorm(nkvec)=sqrt(rknorm(nkvec))
!       write(ounit,'(''k-vec('',i2,'')='',3i5,3f9.4,f11.6)') nkvec,(kvec(k,nkvec),k=1,3),(rkvec(k,nkvec),k=1,3),rknorm(nkvec)
      120 continue
      enddo

! I could just get out of the above loop after finding vcell_sim/vcell
! but instead do check after loop to be safe.
      write(ounit,'(/,''k-vector k-inv      kvec               rkvec'')')
      nkvec_tot=0
      do i=1,nkvec
         nkvec_tot=nkvec_tot+k_inv(i)
         if(ipr.ge.4) then
            write(ounit,'(''k-vec('',i2,'')='',i2,2x,3i4,2x,3f14.10,f11.6)') i,k_inv(i),(kvec(k,i),k=1,3),(rkvec(k,i),k=1,3) &
       ,rknorm(i)
         endif
      enddo
      if(ipr.ge.4) write(ounit,'(''nkvec,nkvec_tot='',2i5)') nkvec,nkvec_tot

! Write out k-pts in reciprocal lattice units for input to pw program
      write(ounit,'(/,i2,'' k-vectors (shifted) in recip. latt. units for input to pw program'')') nkvec
      do ikv=1,nkvec
        do k=1,3
          rkvec_latt(k)=0
          do i=1,3
            rkvec_latt(k)=rkvec_latt(k)+glatt_inv(k,i)*rkvec(i,ikv)
          enddo
        enddo
        if(ipr.ge.4) write(ounit,'(''k-vec('',i2,'')='',i2,2x,3f14.10)') ikv,k_inv(ikv),(rkvec_latt(k),k=1,3)
      enddo

      if(nkvec_tot.ne.nint(vcell_sim/vcell)) then
        write(ounit,'(''Warning: nkvec != vcell_sim/vcell'',9i5)') nkvec_tot,nint(vcell_sim/vcell)
        write(ounit,'(''Possibly the primitive and simulation cells are not commensurate'')')
        if(nkvec_tot.lt.nint(vcell_sim/vcell)) &
        call fatal_error ('You probably need to increase limit of 120 loop if nkvec_tot < vcell_sim/vcell')
        if(nkvec_tot.gt.nint(vcell_sim/vcell)) &
        call fatal_error ('You probably need to increase cutg rel. to cutg_sim if nkvec_tot > vcell_sim/vcell')
      endif



      return
      end
!-----------------------------------------------------------------------

      subroutine fourier_transform(r,arg,r0,nr,vps_short,vcell,gnorm,ngnorm_big,vbare_psp)
! Written by Cyrus Umrigar and Claudia Filippi

! Note: vps_short overwritten
! g > 0 (4pi/vcell)*(int r*vps_short*sin(g*r)*dr)/g
! g = 0 (4pi/vcell)*(int r*2*vps_short*dr)

      use pseudo_mod, only: MPS_GRID
      use ewald_mod, only: NGNORM_BIGX
!      use constant, only: twopi
      use precision_kinds, only: dp
      implicit none

      integer :: ig, ir, ngnorm_big, nr
      real(dp) :: anorm, arg, dx, r0, rlogarg
      real(dp) :: vcell, pi,twopi
      real(dp), dimension(*) :: r
      real(dp), dimension(*) :: vps_short
      real(dp), dimension(*) :: gnorm
!      real(dp), dimension(MPS_GRID) :: y
      real(dp), dimension(nr) :: y
      real(dp), dimension(NGNORM_BIGX) :: vbare_psp


      pi=4.d0*datan(1.d0)
      twopi=2*pi



      anorm=2*twopi/vcell
!      anorm=4*twopi/vcell

! shifted exponential grid
      rlogarg=dlog(arg)
      do ir=1,nr
        vps_short(ir)=(r(ir)+r0)*vps_short(ir)*rlogarg
      enddo

      dx=1.d0
! g != 0 components
      do ig=2,ngnorm_big
         do  ir=1,nr
!     20         y(ir)=r(ir)*vps_short(ir)*sin(gnorm(ig)*r(ir))
            y(ir)=r(ir)*vps_short(ir)*sin(gnorm(ig))
         enddo
         call simson(y,vbare_psp(ig),dx,nr)
!     nrr=((nr-1)/4)*4+1
!     vbare_psp(ig)=bode(y,dx,nrr)
         vbare_psp(ig)=anorm*vbare_psp(ig)/gnorm(ig)
      enddo

! g=0 component
      do  ir=1,nr
        y(ir)=r(ir)*r(ir)*vps_short(ir)
      enddo
      call simson(y,vbare_psp(1),dx,nr)
      vbare_psp(1)=anorm*vbare_psp(1)

      return
      end
!-----------------------------------------------------------------------

      subroutine fourier_it(rvec,gvec,ngnorm_big,igmult,vps_k,vps_r)

      use ewald_mod, only: NGNORM_BIGX
      use contrl_file,    only: ounit
      use control,    only: ipr
      use precision_kinds, only: dp
      implicit none

      integer :: im, ivec, k, ngnorm_big
      integer, dimension(*) :: igmult
      real(dp) :: cos, product, vps_r
      real(dp), dimension(3) :: rvec
      real(dp), dimension(3,*) :: gvec
      real(dp), dimension(NGNORM_BIG) :: vps_k

!      wirte(ounit,*)  rvec

! Note: vpsp_k are:
!     g > 0 (4pi/vcell)*(int r*vps_short*sin(g*r)*dr)/g
!     g = 0 (4pi/vcell)*(int r*2*vps_short*dr)

      ivec=1
!     add other terms
!     since it has been assumed v_k = v_{-k}
      do k=2,ngnorm_big
         do im=1,igmult(k)
            ivec=ivec+1
            product=rvec(1)*gvec(1,ivec)+ &
                 rvec(2)*gvec(2,ivec)+ &
                 rvec(3)*gvec(3,ivec)
            vps_r=vps_r+cos(product)*vps_k(k)
         enddo
      enddo

!     a0 term k=0
      vps_r=2.d0*vps_r+vps_k(1)


      return
      end
!-----------------------------------------------------------------------
      subroutine separate(v,b0,lowest_pow,ngnorm_big,igmult,gnorm,ngnorm, &
       cutr,vcell,ncoef_per,np,b,y,chisq,ifcon,isrange)
! Written by Cyrus Umrigar and Claudia Filippi

      use ewald_mod, only: NCOEFX, NPX
!      use constant, only: twopi
      use precision_kinds, only: dp
      use contrl_file,    only: ounit
      use control,    only: ipr
      use error, only: fatal_error
      implicit none

      integer :: i, i0, ifcon, ig, info
      integer :: isrange, j, k, lowest_pow
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
         write(ounit,'(''ngnorm, ngnorm_big'',2i3)') ngnorm, ngnorm_big
      endif

      do k=ngnorm+1,ngnorm_big
        gr=gnorm(k)*cutr
        ig=k
        if(isrange.eq.0) then
           call integral_sin_poly(gr,lowest_pow,ncoef_per,np,anorm,c)
        elseif(isrange.eq.1) then
           call integral_sin_poly1(gr,ig,lowest_pow,ncoef_per,np,anorm,c)
        elseif(isrange.eq.2) then
           call integral_sin_poly2(gr,ig,lowest_pow,ncoef_per,np,anorm,c)
        elseif(isrange.eq.3) then
           call integral_sin_poly3(gr,ig,lowest_pow,ncoef_per,np,anorm,c)
        endif

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
!     if(chisq.gt.0) then
!       write(ounit,'(''Rms error '',d12.5)') dsqrt(chisq)
!      else
!       write(ounit,'(''Warning: Rms error missing, chisq negative in separate'',d12.4)') chisq
!       if(chisq.lt.0.d0) call fatal_error ('chisq<0 in separate')
!     endif

! this is cusp constraint
      if(ifcon.eq.1) b(1)=beta1+beta2*b(2)

! subtract effect of short range potential on fourier components

      do k=1,ngnorm
        gr=gnorm(k)*cutr
        ig=k
        if(isrange.eq.0) then
          call integral_sin_poly(gr,lowest_pow,ncoef_per,np,anorm,c)
         elseif(isrange.eq.1) then
          call integral_sin_poly1(gr,ig,lowest_pow,ncoef_per,np,anorm,c)
         elseif(isrange.eq.2) then
          call integral_sin_poly2(gr,ig,lowest_pow,ncoef_per,np,anorm,c)
         elseif(isrange.eq.3) then
          call integral_sin_poly3(gr,ig,lowest_pow,ncoef_per,np,anorm,c)
        endif
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
      20     c(i)=1.d0/(i+2+lowest_pow)
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

      subroutine integral_sin_poly1(g,ig,lowest_pow,n,np,anorm,c)
! Written by Cyrus Umrigar and Claudia Filippi
! anorm = 4*pi*cutr^3/volume
! g = g*cutr
! x = r/cutr
! output coefficients c

      use ewald_basis, only: vps_basis_fourier
      use precision_kinds, only: dp
      implicit none

      integer :: i, ig, j, k, lowest_pow
      integer :: n, np
      real(dp) :: anorm, dcmplx, dx, g
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
      20     c(i)=1.d0/(i+2+lowest_pow)
        enddo
      endif

! take care that expansion functions are h_i(x) = x**i*(1-x)**np
! Warning check if we need to go one more.
      do k=1,np
!       do 30 i=1,n+np-k
        do i=1,n+np-k-1
          c(i)=c(i)-c(i+1)
        enddo
      enddo

      if(n.gt.0) c(n)=vps_basis_fourier(ig)

!     write(ounit,'(''g,c1='',f5.1,9f9.5)') g,(c(i),i=1,n)

! Calculate c from numerical integral rather than recursion for small non-zero g's
       if(g.ne.0.d0 .and. g.lt.10.d0) then
       dx=1.d0/(NPTS-1)
!      do 36 i=1,n
       do i=1,n-1
         do j=1,NPTS
           x=(j-1)*dx
           if(g.gt.1.d-6) then
             y(j)=x**(i+lowest_pow)*(1-x)**np*sin(g*x)
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
!     do 40 i=1,n
      do i=1,n-1
        c(i)=anorm*c(i)
      enddo

      return
      end
!-----------------------------------------------------------------------

      subroutine integral_sin_poly2(g,ig,lowest_pow,n,np,anorm,c)
! Written by Cyrus Umrigar and Claudia Filippi
! (anorm/g) * integral_0^1 x*sin(g*x)*h(x) where
! h(x)= \sum_{i=1}^ncoef_per b_i (1-x^np)^{i+1}, x=r/cutr
! anorm = 4*pi*cutr^3/volume
! g = g*cutr
! x = r/cutr
! output coefficients c

      use ewald_mod, only: NCOEFX, NPX
      use ewald_basis, only: vps_basis_fourier
      use precision_kinds, only: dp
      implicit none

      integer :: i, ig, j, lowest_pow, n
      integer :: np
      real(dp) :: anorm, dcmplx, dx
!     real(dp) :: g, gi, sin, ti,et,em
      real(dp) :: g, gi, sin
      real(dp) :: x
      real(dp), dimension(*) :: c
      real(dp), dimension(NPX*(NCOEFX+1)) :: d
      integer, parameter :: NPTS = 1001
      real(dp), dimension(NPTS) :: y


      complex*16 ti,et,em



! integral of sin(g*x)*x**i for i=1 to np*(n+1)+1 and x from 0 to 1
      if(dabs(g).gt.1.d-10) then
        gi=1.d0/g
        ti=dcmplx(0.d0,-gi)
        et=dcmplx(dsin(g)*gi,-dcos(g)*gi)
        em=ti*(et-ti)
        do i=1,np*(n+1)+2
          if(i.gt.lowest_pow+1) d(i-lowest_pow-1)=dreal(em)
!         write(ounit,'(''g,et,i*em'',f7.3,9d12.4)') g,et,i*em
          em=ti*(et-i*em)
        enddo
       else
        do i=1,np*(n+1)+2
      20     d(i)=1.d0/(i+2+lowest_pow)
        enddo
      endif

!     write(ounit,'(''g,d='',f7.3,9f9.5)') g,(d(i),i=1,np*(n+1)+1)

! integral of sin(g*x)*x*(1-x^np)^{i+1} for i=1 to n and x from 0 to 1
!     do 30 i=1,n
      do i=1,n-1
        c(i)=0
        do j=0,i+1
          c(i)=c(i)+choose(i+1,j)*d(j*np+1)*(-1)**j
        enddo
      enddo

      if(n.gt.0) c(n)=vps_basis_fourier(ig)

!     write(ounit,'(''g1,c='',f5.1,9f9.5)') g,(c(i),i=1,n)

! Calculate c from numerical integral rather than recursion for small non-zero g's
       if(g.ne.0.d0 .and. g.lt.10.d0) then
       dx=1.d0/(NPTS-1)
!      do 36 i=1,n
       do i=1,n-1
         do j=1,NPTS
           x=(j-1)*dx
           if(g.gt.1.d-6) then
             y(j)=(1-x**np)**(i+1)*x*sin(g*x)
            else
             y(j)=(1-x**np)**(i+1)*x*x
           endif
         enddo
         c(i)=bode(y,dx,NPTS)
         if(g.gt.1.d-6) then
           c(i)=c(i)/g
         endif
       enddo

!     write(ounit,'(''g2,c='',f5.1,9f9.5)') g,(c(i),i=1,n)
      endif

! multiply by anorm
      do i=1,n-1
        c(i)=anorm*c(i)
      enddo

      return
      end
!-----------------------------------------------------------------------

      subroutine integral_sin_poly3(g,ig,lowest_pow,n,np,anorm,c)
! Written by Cyrus Umrigar and Claudia Filippi
! (anorm/g) * integral_0^1 x*sin(g*x)*h(x) where
! h(x)= \sum_{i=1}^ncoef_per b_i (1-x^np)^{i+1}, x=r/cutr
! anorm = 4*pi*cutr^3/volume
! g = g*cutr
! x = r/cutr
! output coefficients c

      use ewald_mod, only: NCOEFX, NPX
      use ewald_basis, only: vps_basis_fourier
      use precision_kinds, only: dp
      implicit none

      integer :: i, ig, j, lowest_pow, n
      integer :: np
      real(dp) :: anorm, dcmplx, dx
!     real(dp) :: g, gi, sin, ti,et,em
      real(dp) :: g, gi, sin
      real(dp) :: x
      real(dp), dimension(*) :: c
      real(dp), dimension(NPX*(NCOEFX+1)) :: d
      integer, parameter :: NPTS = 1001

      real(dp), dimension(NPTS) :: y


      complex*16 ti,et,em



! integral of sin(g*x)*x**i for i=1 to np*(n+1)+1 and x from 0 to 1
      if(dabs(g).gt.1.d-10) then
        gi=1.d0/g
        ti=dcmplx(0.d0,-gi)
        et=dcmplx(dsin(g)*gi,-dcos(g)*gi)
        em=ti*(et-ti)
        do i=1,np*(n+1)+2
          if(i.gt.lowest_pow+1) d(i-lowest_pow-1)=dreal(em)
!         write(ounit,'(''g,et,i*em'',f7.3,9d12.4)') g,et,i*em
          em=ti*(et-i*em)
        enddo
       else
        do i=1,np*(n+1)+2
      20     d(i)=1.d0/(i+2+lowest_pow)
        enddo
      endif

!     write(ounit,'(''g,d='',f7.3,9f9.5)') g,(d(i),i=1,np*(n+1)+1)

! integral of sin(g*x)*x*(1-x^np)^{i+1} for i=1 to n and x from 0 to 1
!     do 30 i=1,n
      do i=1,n-1
        c(i)=0
        do j=0,np
          c(i)=c(i)+choose(np,j)*d(j*(i+1)+1)*(-1)**j
        enddo
      enddo

      if(n.gt.0) c(n)=vps_basis_fourier(ig)

!     write(ounit,'(''g,c1='',f5.1,9f9.5)') g,(c(i),i=1,n)

! Calculate c from numerical integral rather than recursion for small non-zero g's
       if(g.ne.0.d0 .and. g.lt.10.d0) then
       dx=1.d0/(NPTS-1)
!      do 36 i=1,n
       do i=1,n-1
         do j=1,NPTS
           x=(j-1)*dx
           if(g.gt.1.d-6) then
             y(j)=(1-x**(i+1))**np*x*sin(g*x)
            else
             y(j)=(1-x**(i+1))**np*x*x
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
!     do 40 i=1,n
      do i=1,n-1
        c(i)=anorm*c(i)
      enddo

      return
      end
!-----------------------------------------------------------------------

      real*8 function choose(n,m)
! Written by Cyrus Umrigar
! Binomial coefficients ^nC_m
      use precision_kinds, only: dp
      implicit none

      integer :: i, m, n


      choose=1
      do i=1,m
        choose=choose*(n-i+1)/dfloat(i)
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

      real*8 function vsrange1(r,cutr,lowest_pow,ncoef_per,np,b,ict,l)
! Written by Cyrus Umrigar
! h(x)= \sum_{i=1}^ncoef_per b_i x^{i-1} (1-x)^np, x=r/cutr
!     use readps_tm_mod, only: splfit_tm
      use readps_gauss, only: gauss_pot
      use contrl_file,    only: ounit
      use precision_kinds, only: dp
      use system, only: znuc
      implicit none

      integer :: i, ict, l, ncoef_per, np, lowest_pow
      real(dp) :: cutr, r, vpot, dvpot, x, vsrange
      real(dp), dimension(*) :: b

      vpot=0.d0
      dvpot=0.d0

      r=max(1.0d-10,r)
      x=r/cutr

!      print*,'inside vsrange1'
!      print*, 'ncoef_per', ncoef_per
!      print*,'x',x,'r',r,'cutr',cutr
!     print*,"b",b(1:ncoef_per)
!      write(ounit,'(''inside vsrange1'')')
!      write(ounit,'(''ncoef_per'',i3)') ncoef_per
!      write(ounit,'(''b'',9f9.5)') b(1:ncoef_per)


!     This is added to guarantee it goes correctly to the limit matching condition at cutr
!     however it should not be used at larger range than rcut since is zero

      vsrange1=0.d0


!     to compute the expansion we got from the fit for w(r) potential no first coefficient
      vsrange=0.d0

!     do 10 i=1,ncoef_per
!     10   vsrange1=b(ncoef_per-i+1)+x*vsrange1

      do i=1,ncoef_per-1
!     print*,'b(i)', ncoef_per-i, b(ncoef_per-i)
!        write(ounit,'(''i,b(i)'',i3,9f9.5)') i,b(ncoef_per-i)
        vsrange=b(ncoef_per-i)+x*vsrange
      enddo

!      print*,"vsrange1  before scaling", vsrange
      vsrange=vsrange*(1-x)**np
!      write(ounit,*) " "
!      write(ounit,*) "vsrange1 incomplete",vsrange

!     print*,'splfit_tm input',cutr,r,l,ict,vpot




!      write(ounit,'(''ict,l,r,vsrange1,vpot'',2i3,9f9.5)') ict,l,r,vsrange1,vpot
!      call splfit_tm(r,l,ict,vpot)
!      vsrange1=vsrange+b(ncoef_per)*vpot
      call gauss_pot(r,l,ict,vpot,dvpot)
      vpot=vpot-znuc(ict)/r
      vsrange1=vsrange1+b(ncoef_per)*vpot

!     Adding term with correction (last term of expansion)
      vsrange1=vsrange1+vsrange



!     write(ounit,*) "b(n), vpot, b(n)*vpot",b(ncoef_per), vpot, b(ncoef_per)*vpot
!     write(ounit,*) "vsrange1 complete",vsrange1


      if(x.gt.1.d0) then
!        write(ounit,*) "Warning reached cutr threshold", vsrange1
         vsrange1=0.d0
         return
      endif


!      write(ounit,*) " "

!     if(lowest_pow.eq.-1) vsrange1=vsrange1/x

      return
      end
!-----------------------------------------------------------------------

      real*8 function vsrange2(r,cutr,lowest_pow,ncoef_per,np,b,ict,l)
! Written by Cyrus Umrigar
! h(x)= \sum_{i=1}^ncoef_per b_i (1-x^np)^{i+1}, x=r/cutr
!     use readps_tm_mod, only: splfit_tm
      use system, only: znuc
      use readps_gauss, only: gauss_pot
      use precision_kinds, only: dp
      implicit none

      integer :: i, ict, l, ncoef_per, np, lowest_pow
      real(dp) :: cutr, r, term, vpot,dvpot
      real(dp) :: x
      real(dp), dimension(*) :: b


      r=max(1.0d-10,r)
      x=r/cutr

      vsrange2=0
      if(x.gt.1.d0) return

      term=1-x**np
!     do 10 i=1,ncoef_per
!  10   vsrange2=b(ncoef_per-i+1)+term*vsrange2
      do i=1,ncoef_per-1
        vsrange2=b(ncoef_per-i)+term*vsrange2
      enddo

      vsrange2=vsrange2*term*term

!     call splfit_tm(r,l,ict,vpot)
      call gauss_pot(r,l,ict,vpot,dvpot)
      vpot=vpot-znuc(ict)/r
!     write(ounit,'(''ict,l,r,vsrange2,vpot'',2i3,9f9.5)') ict,l,r,vsrange2,vpot
      vsrange2=vsrange2+b(ncoef_per)*vpot

!     if(lowest_pow.eq.-1) vsrange2=vsrange2/x

      return
      end
!-----------------------------------------------------------------------

      real*8 function vsrange3(r,cutr,lowest_pow,ncoef_per,np,b,ict,l)
! Written by Cyrus Umrigar
! h(x)= \sum_{i=1}^ncoef_per b_i (1-x^{i+1})^np, x=r/cutr
!     use readps_tm_mod, only: splfit_tm
      use readps_gauss, only: gauss_pot
      use system, only: znuc
      use precision_kinds, only: dp
      implicit none

      integer :: i, ict, l, ncoef_per, np, lowest_pow
      real(dp) :: cutr, r, vpot, x, dvpot
      real(dp), dimension(*) :: b


      r=max(1.0d-10,r)
      x=r/cutr

      vsrange3=0
      if(x.gt.1.d0) return

!     do 10 i=1,ncoef_per
      do i=1,ncoef_per-1
        vsrange3=vsrange3+b(i)*(1-x**(i+1))**np
      enddo

!     call splfit_tm(r,l,ict,vpot)
      call gauss_pot(r,l,ict,vpot,dvpot)
      vpot=vpot-znuc(ict)/r
!     write(ounit,'(''ict,l,r,vsrange3,vpot'',2i3,9f9.5)') ict,l,r,vsrange3,vpot
      vsrange3=vsrange3+b(ncoef_per)*vpot

!     if(lowest_pow.eq.-1) vsrange3=vsrange3/x

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
!      gauss_exp=7.d0/cutr

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
      ewald_pot=2*ewald_pot+y(1)+derfc(gauss_exp*rr)/rr

      return
      end



!-----------------------------------------------------------------------


      real*8 function ewald_pot_psp(rvec,rr,gvec,gnorm,ngnorm,igmult,y,cutr,vcell,ict,l,z)
! Written by Cyrus Umrigar
!     use readps_tm_mod, only: splfit_tm
      use readps_gauss, only: gauss_pot
!      use constants, only: pi
      use system, only: znuc
      use precision_kinds, only: dp
      implicit none

      integer :: ict, im, ivec, k, l
      integer :: ngnorm
      integer, dimension(*) :: igmult
      real(dp) :: cos, cutr, expon, gaus_exp
      real(dp) :: product, rr, vcell, vpot, dvpot
      real(dp) :: z, pi
      real(dp), dimension(3) :: rvec
      real(dp), dimension(3,*) :: gvec
      real(dp), dimension(*) :: gnorm
      real(dp), dimension(*) :: y


      pi=4.d0*datan(1.d0)


      gaus_exp=5/cutr
      ivec=1
! The factor of 2 in the next line is just to compensate for the 2 in the
! last line, which is there because we keep only half the vectors in the star.
      ewald_pot_psp=-pi/(2*vcell*gaus_exp**2)
      do k=2,ngnorm
        expon=exp(-(gnorm(k)/(2*gaus_exp))**2)
        do im=1,igmult(k)
          ivec=ivec+1
          product=rvec(1)*gvec(1,ivec)+ &
                  rvec(2)*gvec(2,ivec)+ &
                  rvec(3)*gvec(3,ivec)
          ewald_pot_psp=ewald_pot_psp+cos(product)*y(k)*expon
        enddo
      enddo
!     call splfit_tm(rr,l,ict,vpot)
      call gauss_pot(rr,l,ict,vpot,dvpot)
      vpot=vpot-znuc(ict)/rr
!      wrtie(ounit,*) "inside ewald_psp y(1)",Y(1)
!      write(ounit,'(''rr,ewald_pot_psp'',f8.4,9f9.5)') rr,-z*(2*ewald_pot_psp+y(1)),vpot,-z*(2*ewald_pot_psp+y(1))+vpot
      ewald_pot_psp=-z*(2*ewald_pot_psp+y(1))+vpot

      return
      end
!-----------------------------------------------------------------------
      real*8 function vlrange_old(rvec,gvec,ngnorm,igmult,y)
      use ewald_mod, only: NGNORM_SIM_BIGX, NGVEC_SIM_BIGX
! Written by Cyrus Umrigar

      use precision_kinds, only: dp
      implicit none

      integer :: im, ivec, k, ngnorm
      integer, dimension(*) :: igmult
      real(dp) :: cos, product, vlrange
      real(dp), dimension(3) :: rvec
      real(dp), dimension(3,*) :: gvec
      real(dp), dimension(*) :: y


!     dimension rvec(3),gvec(3,NGVEC_SIM_BIGX),igmult(NGNORM_SIM_BIGX),y(NGNORM_SIM_BIGX)

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
      vlrange_old=2*vlrange+y(1)

      return
      end
!-----------------------------------------------------------------------

      real*8 function vlrange_nn_old2(ncent,znuc,iwctype,ngnorm,igmult,cos_g,sin_g,y)
      use system, only: nelec
! Written by Cyrus Umrigar

      use precision_kinds, only: dp
      implicit none

      integer :: i, im, ivec, k, ncent
      integer :: ngnorm
      integer, dimension(*) :: iwctype
      integer, dimension(*) :: igmult
      real(dp) :: cos_sum, sin_sum, vl, znuci
      real(dp), dimension(*) :: znuc
      real(dp), dimension(nelec,*) :: cos_g
      real(dp), dimension(nelec,*) :: sin_g
      real(dp), dimension(*) :: y



      ivec=1
      vl=0.d0
      do k=2,ngnorm
        do im=1,igmult(k)
          ivec=ivec+1
          cos_sum=0.d0
          sin_sum=0.d0
          do i=1,ncent
            znuci=znuc(iwctype(i))
            cos_sum=cos_sum+znuci*cos_g(i,ivec)
            sin_sum=sin_sum+znuci*sin_g(i,ivec)
          enddo
          vl=vl+y(k)*(cos_sum**2+sin_sum**2)
        enddo
      enddo
      vlrange_nn_old2=vl

      return
      end
!-----------------------------------------------------------------------

      real*8 function vlrange_ee_old2(nelec,ngnorm,igmult,cos_g,sin_g,y)
!      use system, only: nelec
! Written by Cyrus Umrigar

      use precision_kinds, only: dp
      implicit none
      integer nelec
      integer :: i, im, ivec, k, ngnorm
      integer, dimension(*) :: igmult
      real(dp) :: cos_sum, sin_sum, vl
      real(dp), dimension(nelec,*) :: cos_g
      real(dp), dimension(nelec,*) :: sin_g
      real(dp), dimension(*) :: y



      ivec=1
      vl=0.d0
      do k=2,ngnorm
        do im=1,igmult(k)
          ivec=ivec+1
          cos_sum=0.d0
          sin_sum=0.d0
          do i=1,nelec
            cos_sum=cos_sum+cos_g(i,ivec)
            sin_sum=sin_sum+sin_g(i,ivec)
          enddo
          vl=vl+y(k)*(cos_sum**2+sin_sum**2)
        enddo
      enddo
      vlrange_ee_old2=vl

      return
      end
!-----------------------------------------------------------------------

      real*8 function vlrange(ngnorm,igmult,cos1_sum,cos2_sum,sin1_sum,sin2_sum,y)
! Written by Cyrus Umrigar

      use precision_kinds, only: dp
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

      real*8 function vlrange_p(ngnorm,igmult,cos1_sum,cos2_sum,sin1_sum,sin2_sum)
! Written by Cyrus Umrigar

      use precision_kinds, only: dp
      implicit none

      integer :: im, ivec, k, ngnorm
      integer, dimension(*) :: igmult
      real(dp) :: vl
      real(dp), dimension(*) :: cos1_sum
      real(dp), dimension(*) :: cos2_sum
      real(dp), dimension(*) :: sin1_sum
      real(dp), dimension(*) :: sin2_sum



      ivec=1
      vl=0.5d0*(cos1_sum(1)*cos2_sum(1)+sin1_sum(1)*sin2_sum(1))
      do k=2,ngnorm
        do im=1,igmult(k)
          ivec=ivec+1
          vl=vl+(cos1_sum(ivec)*cos2_sum(ivec)+sin1_sum(ivec)*sin2_sum(ivec))
        enddo
      enddo
      vlrange_p=vl

      return
      end

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



      subroutine simson(y,s,h,n)


      use precision_kinds, only: dp
      implicit none

      real(dp), parameter :: c3=0.125d0
      real(dp), parameter :: c2=-0.625d0
      real(dp), parameter :: c1=1.375d0
      real(dp), parameter :: c0=-2.875d0
      real(dp), dimension (*) :: y
      real(dp) :: s, h, h3, odd, eve
      integer :: n, n1, i

      h3=0.333333333333333d0*h
      n1=n-1
      odd=0.d0
      eve=0.d0

      do i=1,n1,2
        odd=odd+y(i)
        eve=eve+y(i+1)
      enddo

      s=2.d0*(odd+2.d0*eve) - y(1)

      if(mod(n,2).eq.1) then
         s=(s+y(n)) * h3
      else
         s=(s+c3*y(n-3)+c2*y(n-2)+c1*y(n-1)+c0*y(n)) * h3
      endif



      return
      end



!-----------------------------------------------------------------------

      subroutine pot_nn_ewald_old
! Written by Cyrus Umrigar

      use system, only: znuc, cent, iwctype, ncent
      use multiple_geo, only: pecent
      use contrl_file,    only: ounit
      use ewald, only: b_coul, y_coul

      use periodic, only: cutr
      use periodic, only: gvec, igmult
      use periodic, only: ncoef_per, ngnorm
      use periodic, only: np, vcell
      use periodic, only: vcell_sim, znuc2_sum
      use find_pimage, only: find_image3,find_image_pbc
      use precision_kinds, only: dp
      implicit none

      integer :: i, j, k, lowest_pow
      real(dp) :: c0, rnorm, vl, vlr
!     real(dp) :: vs, vsrange, zprod
      real(dp) :: vs, zprod
      real(dp), dimension(3) :: r






      lowest_pow=-1
      c0=(b_coul(2)-np*b_coul(1))/2
      vs=c0*znuc2_sum
      vl=0
      do i=1,ncent
        do j=1,i
          zprod=znuc(iwctype(i))*znuc(iwctype(j))
          do k=1,3
            r(k)=cent(k,j)-cent(k,i)
          enddo
!     call find_image3(r,rnorm)
          call find_image_pbc(r,rnorm)
          if(i.ne.j) then
            vs=vs+zprod*vsrange(rnorm,cutr,lowest_pow,ncoef_per,np,b_coul)
          endif
          vlr=vlrange_old(r,gvec,ngnorm,igmult,y_coul)
          if(i.eq.j) vlr=0.5d0*vlr
      20     vl=vl+zprod*vlr
        enddo
      enddo
      pecent=vs+vl
      vs=vs*2/ncent
      vl=vl*2/ncent
      write(ounit,'(''v_nn,vs,vl,vs1,vl1='',9f12.8)') pecent*2/ncent,vs,vl,znuc2_sum*c0*2/ncent
      pecent=pecent*vcell_sim/vcell

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
      use periodic, only: np, vcell
      use periodic, only: vcell_sim, znuc2_sum, znuc_sum
      use precision_kinds, only: dp
      use ewald_mod, only: NGVECX
      use find_pimage, only: find_image3, find_image_pbc
      implicit none

      integer :: i, j, k, lowest_pow
      real(dp) :: c0, rnorm, vl
      real(dp) :: vs, zprod
      real(dp), dimension(3) :: r





! short-range sum
      lowest_pow=-1
      c0=(b_coul(2)-np*b_coul(1))/2
      vs=c0*znuc2_sum
      do i=1,ncent
        do j=1,i-1
          zprod=znuc(iwctype(i))*znuc(iwctype(j))
          do k=1,3
            r(k)=cent(k,j)-cent(k,i)
          enddo
!          print *, 'CIAO',i,j
!          print *, 'CIAO',(r(k),k=1,3)
!     call find_image3(r,rnorm)
          call find_image_pbc(r,rnorm)
!          print *, 'CIAO',(r(k),k=1,3)
!          print *, 'CIAO',sqrt(r(1)**2+r(2)**2+r(3)**2),rnorm
          vs=vs+zprod*vsrange(rnorm,cutr,lowest_pow,ncoef_per,np,b_coul)
!          print *, 'CIAO',vsrange(rnorm,cutr,lowest_pow,ncoef_per,np,b_coul)
        enddo
      enddo

! long-range sum
!     call cossin_old2(glatt,igvec,ngvec,cent,ncent,ng1d,cos_g,sin_g)
      call cossin_n(znuc,iwctype,glatt,igvec,ngvec,cent,ncent,ng1d,cos_n_sum,sin_n_sum)

!     vl=vlrange_nn_old2(ncent,znuc,iwctype,ngnorm,igmult,cos_g,sin_g,y_coul)
      vl=vlrange(ngnorm,igmult,cos_n_sum,cos_n_sum,sin_n_sum,sin_n_sum,y_coul)
!     vl=vl+0.5d0*y_coul(1)*znuc_sum**2

      pecent=vs+vl
      vs=vs*2/ncent
      vl=vl*2/ncent
      if(ipr.ge.4) then
         write(ounit,'(''v_nn,vs,vl,vs1,vl1='',9f12.8)') pecent*2/ncent,vs,vl,znuc2_sum*c0*2/ncent,y_coul(1)*znuc_sum**2/ncent
      endif
      pecent=pecent*vcell_sim/vcell

      return
      end


!-----------------------------------------------------------------------

      subroutine pot_en_coul_ewald(x,pe_en)
! Written by Cyrus Umrigar

      use contrl_file,    only: ounit
      use vmc_mod, only: nmat_dim2
      use system, only: znuc, cent, iwctype, ncent, ncent_tot
      use control, only: ipr
      use system, only: nelec
      use ewald, only: b_coul, y_coul, y_psp, b_psp, y_jas, b_jas
      use ewald, only: cos_n_sum, sin_n_sum, cos_e_sum, sin_e_sum, cos_e_sum_sim, sin_e_sum_sim, cos_p_sum, sin_p_sum
      use periodic, only: cutr, glatt
      use periodic, only: igmult, igvec
      use periodic, only: isrange, ncoef_per, ng1d, ngnorm
      use periodic, only: ngvec
      use periodic, only: np
      use periodic, only: znuc_sum
      use pseudo, only: lpot, nloc
      use distance_mod, only: r_en, rvec_en, r_ee, rvec_ee

      use ewald_mod, only: NGNORM_SIMX, NGVEC_SIMX, NCOEFX,NGNORMX, NGVECX
      use find_pimage, only: find_image3, find_image_pbc

      use precision_kinds, only: dp
      implicit none

      integer :: i, ict, j, k, lowest_pow
      real(dp) :: pe_en, vl
      real(dp) :: vs, vs_aux
!      real(dp) :: vsrange, vsrange1, vsrange2, vsrange3
      real(dp), dimension(3,*) :: x


!      write(ounit,*) "inside pot_en_ewald isrange", isrange
!      write(ounit,*) "inside pot_en_ewald nloc", nloc




! short-range sum
! Warning: I need to call the appropriate vsrange
      vs=0.d0
      lowest_pow=-1
      do i=1,ncent
!        ict=iwctype(i)
        do j=1,nelec
           do k=1,3
              rvec_en(k,j,i)=x(k,j)-cent(k,i)
!     write(ounit,'(''x, cent '',3d12.4)') x(k,j),cent(k,i)
           enddo
           call find_image_pbc(rvec_en(1,j,i),r_en(j,i))
!           vs=vs-znuc(ict)*vsrange(r_en(j,i),cutr,lowest_pow,ncoef_per,np,b_coul)
           vs=vs-znuc(iwctype(i))*vsrange(r_en(j,i),cutr,lowest_pow,ncoef_per,np,b_coul)
        enddo
      enddo


! long-range sum

      call cossin_e(glatt,igvec,ngvec,x,nelec,ng1d,cos_e_sum,sin_e_sum)

      vl=-2*vlrange(ngnorm,igmult,cos_n_sum,cos_e_sum,sin_n_sum,sin_e_sum,y_coul)



!      write(ounit,'(''before print vl'')')

      pe_en=vs+vl
!      write(ounit,'(''pe_en='',20d12.4)') pe_en

!      write(ounit,'(''vs1='',20f12.4)') vs
!      write(ounit,'(''vl='',20f12.4)') vl
!      write(ounit,'(''vl1='',20f12.4)') -y_coul(1)*znuc_sum


      vs=vs/nelec
      vl=vl/nelec
      if(ipr.ge.2) then
        if(nloc.eq.0) write(ounit,'(''v_en,vs,vl,vl1='',9f12.8)') pe_en/nelec,vs,vl,-y_coul(1)*znuc_sum
        if(nloc.ne.0) write(ounit,'(''v_en,vs,vl,vl1='',9f12.8)') pe_en/nelec,vs,vl,y_psp(1,1)*znuc_sum
      endif


!      write(ounit,'(''after print vl'')')


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
      use ewald, only: b_coul, y_coul, y_psp, b_psp, y_jas, b_jas
      use ewald, only: cos_n_sum, sin_n_sum, cos_e_sum, sin_e_sum, cos_e_sum_sim, sin_e_sum_sim, cos_p_sum, sin_p_sum
      use periodic, only: cutr, glatt
      use periodic, only: igmult, igvec
      use periodic, only: isrange, ncoef_per, ng1d, ngnorm
      use periodic, only: ngvec
      use periodic, only: np
      use periodic, only: znuc_sum
      use pseudo, only: lpot, nloc
      use distance_mod, only: r_en, rvec_en, r_ee, rvec_ee

      use ewald_mod, only: NGNORM_SIMX, NGVEC_SIMX, NCOEFX,NGNORMX, NGVECX
      use find_pimage, only: find_image3, find_image_pbc

      use precision_kinds, only: dp
      implicit none

      integer :: i, ict, j, k, lowest_pow
      real(dp) :: pe_en, vl
      real(dp) :: vs, vs_aux
!      real(dp) :: vsrange, vsrange1, vsrange2, vsrange3
      real(dp), dimension(3,*) :: x


!      write(ounit,*) "inside pot_en_ewald isrange", isrange
!      write(ounit,*) "inside pot_en_ewald nloc", nloc




! short-range sum
! Warning: I need to call the appropriate vsrange
      vs=0.d0
      do i=1,ncent
        ict=iwctype(i)
        do j=1,nelec
          do k=1,3
             rvec_en(k,j,i)=x(k,j)-cent(k,i)
!             write(ounit,'(''x, cent '',3d12.4)') x(k,j),cent(k,i)
         enddo
!         write(ounit,'(''rvec_en '',3d12.4)') rvec_en(:,j,i)
!     call find_image3(rvec_en(1,j,i),r_en(j,i))
         call find_image_pbc(rvec_en(1,j,i),r_en(j,i))
!          if(nloc.eq.0) then
            lowest_pow=-1
            vs=vs-znuc(ict)*vsrange(r_en(j,i),cutr,lowest_pow,ncoef_per,np,b_coul)
!            write(ounit,'(''vs if nloc==0'',1d12.4)') vs
!           else
!            lowest_pow=0
!     vs=vs+vsrange(r_en(j,i),cutr,lowest_pow,ncoef_per,np,b_psp(1,ict))
!            write(ounit,'(''before compute tutut  pot_en_ewald '')')
!            write(ounit,'(''ict'', i20)') ict
!            write(ounit,'(''j'', i20)') j
!            write(ounit,'(''r_en '',20d12.4)') r_en(j,i)
!            if(isrange.eq.0) vs=vs+vsrange(r_en(j,i),cutr,lowest_pow,ncoef_per,np,b_psp(1,ict))
!            if(isrange.eq.1) vs=vs+vsrange1(r_en(j,i),cutr,lowest_pow,ncoef_per,np,b_psp(1,ict),ict,lpot(ict))
!            if(isrange.eq.2) vs=vs+vsrange2(r_en(j,i),cutr,lowest_pow,ncoef_per,np,b_psp(1,ict),ict,lpot(ict))
!            if(isrange.eq.3) vs=vs+vsrange3(r_en(j,i),cutr,lowest_pow,ncoef_per,np,b_psp(1,ict),ict,lpot(ict))
!     here this function is
!            vs=vs+vsrange(r_en(j,i),cutr,lowest_pow,ncoef_per,np,b_coul)
!            write(ounit,'(''isrange, local_vs'',i20, d12.4)') isrange, vs
!            write(ounit,'(''after compute tutut pot_en_ewald '')')
!          endif
        enddo
      enddo



! long-range sum
!     call cossin_e(glatt,igvec,ngvec,xold,nelec,ng1d,cos_e_sum,sin_e_sum)
      call cossin_e(glatt,igvec,ngvec,x,nelec,ng1d,cos_e_sum,sin_e_sum)



      if(nloc.eq.0) then
        vl=-2*vlrange(ngnorm,igmult,cos_n_sum,cos_e_sum,sin_n_sum,sin_e_sum,y_coul)
!       vl=vl-y_coul(1)*znuc_sum*nelec
       else
!        call cossin_p(y_psp,iwctype,glatt,igvec,ngnorm,igmult,cent,ncent,ng1d,cos_p_sum,sin_p_sum)
!        vl=+2*vlrange_p(ngnorm,igmult,cos_p_sum,cos_e_sum,sin_p_sum,sin_e_sum)
!       vl=vl+y_psp(1,iwctype(i))*znuc_sum*nelec
        vl=-2*vlrange(ngnorm,igmult,cos_n_sum,cos_e_sum,sin_n_sum,sin_e_sum,y_coul)
!        write(ounit,*) "tralalata here we are vl vl vl"
      endif


!      write(ounit,'(''before print vl'')')

      pe_en=vs+vl
!      write(ounit,'(''pe_en='',20d12.4)') pe_en

!      write(ounit,'(''vs1='',20f12.4)') vs
!      write(ounit,'(''vl='',20f12.4)') vl
!      write(ounit,'(''vl1='',20f12.4)') -y_coul(1)*znuc_sum


      vs=vs/nelec
      vl=vl/nelec
      if(ipr.ge.2) then
        if(nloc.eq.0) write(ounit,'(''v_en,vs,vl,vl1='',9f12.8)') pe_en/nelec,vs,vl,-y_coul(1)*znuc_sum
        if(nloc.ne.0) write(ounit,'(''v_en,vs,vl,vl1='',9f12.8)') pe_en/nelec,vs,vl,y_psp(1,1)*znuc_sum
      endif


!      write(ounit,'(''after print vl'')')


      return
      end


!-----------------------------------------------------------------------

      subroutine pot_ee_ewald(x,pe_ee)
! Written by Cyrus Umrigar

      use contrl_file,    only: ounit
      use vmc_mod, only: nmat_dim2
      use control, only: ipr
      use system, only: nelec
      use ewald, only: b_coul_sim, y_coul_sim
      use ewald, only: cos_n_sum, sin_n_sum, cos_e_sum, sin_e_sum, cos_e_sum_sim, sin_e_sum_sim, cos_p_sum, sin_p_sum

      use periodic, only: cutr_sim, isrange
      use periodic, only: glatt_sim, igmult_sim, igvec_sim
      use periodic, only: ncoef_per, ng1d_sim
      use periodic, only: ngnorm_sim, ngvec_sim
      use periodic, only: np
      use distance_mod, only: r_en, rvec_en, r_ee, rvec_ee
      use ewald_mod, only: NGNORM_SIMX, NGVEC_SIMX, NCOEFX,NGNORMX, NGVECX
      use find_pimage, only: find_image3, find_image_pbc
      use pseudo, only: lpot, nloc

      use precision_kinds, only: dp
      implicit none

      integer :: i, ij, j, k, lowest_pow
      real(dp) :: c0, pe_ee, vl
      real(dp) :: vs
      real(dp), dimension(3,*) :: x


!      write(ounit,*) "inside pot_ee_ewald isrange", isrange
!      write(ounit,*) "inside pot_ee_ewald nloc", nloc


! short-range sum
      lowest_pow=-1
      c0=(b_coul_sim(2)-np*b_coul_sim(1))/2
      vs=c0*nelec
      ij=0
      do i=1,nelec
        do j=1,i-1
          ij=ij+1
          do k=1,3
            rvec_ee(k,ij)=x(k,i)-x(k,j)
          enddo
!     call find_image3(rvec_ee(1,ij),r_ee(ij))
          call find_image_pbc(rvec_ee(1,ij),r_ee(ij))
          vs=vs+vsrange(r_ee(ij),cutr_sim,lowest_pow,ncoef_per,np,b_coul_sim)
        enddo
      enddo




! long-range sum
!     call cossin_old2(glatt_sim,igvec_sim,ngvec_sim,xold,nelec,ng1d_sim,cos_g,sin_g)
!     call cossin_e(glatt_sim,igvec_sim,ngvec_sim,xold,nelec,ng1d_sim,cos_e_sum_sim,sin_e_sum_sim)
      call cossin_e(glatt_sim,igvec_sim,ngvec_sim,x,nelec,ng1d_sim,cos_e_sum_sim,sin_e_sum_sim)


!     vl=vlrange_ee_old2(nelec,ngnorm_sim,igmult_sim,cos_g,sin_g,y_coul_sim)
      vl=vlrange(ngnorm_sim,igmult_sim,cos_e_sum_sim,cos_e_sum_sim,sin_e_sum_sim,sin_e_sum_sim,y_coul_sim)
!     vl=vl+0.5d0*y_coul_sim(1)*nelec**2

!      write(ounit,'(''before call v_ee'')')
      pe_ee=vs+vl

      vs=vs*2/nelec
      vl=vl*2/nelec
      if(ipr.ge.2) write(ounit,'(''v_ee,vs,vl,vs1,vl1='',9f12.8)') pe_ee*2/nelec,vs,vl,c0*2,y_coul_sim(1)*nelec
!      write(ounit,'(''after call v_ee_ewald'')')

      return
      end
!-----------------------------------------------------------------------

      subroutine cossin_old2(glatt,igvec,ngvec,r,nr,ng1d,cos_g,sin_g)
      use ewald_mod, only: NG1DX
      use system, only: nelec
! Written by Cyrus Umrigar

      use precision_kinds, only: dp
      implicit none

      integer :: i, ir, k, n, ngvec
      integer :: nr
      integer, dimension(3,*) :: igvec
      integer, dimension(3) :: ng1d
      real(dp) :: cos, cos_tmp, dot, sin, sin_tmp
      real(dp), dimension(3,3) :: glatt
      real(dp), dimension(3,*) :: r
      real(dp), dimension(nelec,*) :: cos_g
      real(dp), dimension(nelec,*) :: sin_g
      real(dp), dimension(-NG1DX:NG1DX,3) :: cos_gr
      real(dp), dimension(-NG1DX:NG1DX,3) :: sin_gr



! Calculate cosines and sines for all positions and reciprocal lattice vectors
      do ir=1,nr
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
      20     sin_gr(-n,i)=-sin_gr(n,i)
        enddo
      enddo

      cos_g(ir,1)=1.d0
      sin_g(ir,1)=0.d0
      do i=2,ngvec
        cos_tmp=cos_gr(igvec(1,i),1)*cos_gr(igvec(2,i),2) &
               -sin_gr(igvec(1,i),1)*sin_gr(igvec(2,i),2)
        sin_tmp=sin_gr(igvec(1,i),1)*cos_gr(igvec(2,i),2) &
               +cos_gr(igvec(1,i),1)*sin_gr(igvec(2,i),2)
        cos_g(ir,i)=cos_tmp*cos_gr(igvec(3,i),3) &
                   -sin_tmp*sin_gr(igvec(3,i),3)
        sin_g(ir,i)=sin_tmp*cos_gr(igvec(3,i),3) &
                   +cos_tmp*sin_gr(igvec(3,i),3)
      enddo
      enddo

      return
      end
!-----------------------------------------------------------------------

      subroutine cossin_psi(glatt,gnorm,gvec,igvec,ngvec,r,nr,ng1d,cos_g,sin_g &
      ,dcos_g,dsin_g,ddcos_g,ddsin_g,g_shift,iflag)
      use ewald_mod, only: NG1DX
      use system, only: nelec
! Written by Cyrus Umrigar
! iflag = 0 Calculate cos(gr) and sin(gr) and first 2 derivs at electron positions.
!       = 1 Calculate cos(kr) and sin(kr) and first 2 derivs at electron positions.
! Needed for orbitals and their Laplacian.
! Presently using cossin_psi_g and cossin_psi_k instead.

      use error, only: fatal_error
      use precision_kinds, only: dp
      implicit none

      integer :: i, iflag, ir, k, n
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
      real(dp), dimension(*) :: g_shift
      real(dp), dimension(-NG1DX:NG1DX,3) :: cos_gr
      real(dp), dimension(-NG1DX:NG1DX,3) :: sin_gr



! Calculate cosines and sines for recip. lattice vectors along axes first.
      do ir=1,nr
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
      20     sin_gr(-n,i)=-sin_gr(n,i)
        enddo
      enddo

! If the calculation is for g-vectors then no shift; if for k-vectors there could be one.
      if(iflag.eq.0) then
        cos_tmp0=1.d0
        sin_tmp0=0.d0
       elseif(iflag.eq.1) then
        dot=0
        do k=1,3
          dot=dot+g_shift(k)*r(k,ir)
        enddo
        cos_tmp0=cos(dot)
        sin_tmp0=sin(dot)
       else
        call fatal_error ('iflag must be 0 or 1 in cossin_psi')
      endif

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
        ddcos_g(ir,i)=-gnorm(i)*gnorm(i)*cos_g(ir,i)
        ddsin_g(ir,i)=-gnorm(i)*gnorm(i)*sin_g(ir,i)
      enddo
      enddo

      return
      end
!-----------------------------------------------------------------------

!     subroutine cossin_psi_g(glatt,gnorm,igmult,ngnorm,gvec,igvec,ngvec,r,nr,ng1d,cos_g,sin_g
      subroutine cossin_psi_g(glatt,gnorm,igmult,ngnorm,gvec,igvec,ngvec,r,ir,ng1d,cos_g,sin_g &
      ,dcos_g,dsin_g,ddcos_g,ddsin_g,g_shift)
      use ewald_mod, only: NG1DX
      use system, only: nelec
! Written by Cyrus Umrigar
! Calculate cos(gr) and sin(gr) and first 2 derivs at electron positions.
! Needed for orbitals and their Laplacian.

      use precision_kinds, only: dp
      implicit none

      integer :: i, im, in, ir, k
      integer :: n, ngnorm, ngvec
      integer, dimension(*) :: igmult
      integer, dimension(3,*) :: igvec
      integer, dimension(3) :: ng1d
      real(dp) :: cos, cos_tmp, dot, sin, sin_tmp
      real(dp), dimension(3,3) :: glatt
      real(dp), dimension(*) :: gnorm
      real(dp), dimension(3,*) :: gvec
      real(dp), dimension(3) :: r
      real(dp), dimension(*) :: cos_g
      real(dp), dimension(*) :: sin_g
      real(dp), dimension(3,*) :: dcos_g
      real(dp), dimension(3,*) :: dsin_g
      real(dp), dimension(*) :: ddcos_g
      real(dp), dimension(*) :: ddsin_g
      real(dp), dimension(*) :: g_shift
      real(dp), dimension(-NG1DX:NG1DX,3) :: cos_gr
      real(dp), dimension(-NG1DX:NG1DX,3) :: sin_gr


!     dimension glatt(3,3),gnorm(*),igmult(*),gvec(3,*),igvec(3,*),r(3,*),ng1d(3)
!    &,cos_g(nelec,*),sin_g(nelec,*)
!    &,dcos_g(3,nelec,*),dsin_g(3,nelec,*)
!    &,ddcos_g(nelec,*),ddsin_g(nelec,*),g_shift(*)

! Calculate cosines and sines for recip. lattice vectors along axes first.
!     do 30 ir=1,nr
      do i=1,3
        dot=0
        do k=1,3
          dot=dot+glatt(k,i)*r(k)
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
      20     sin_gr(-n,i)=-sin_gr(n,i)
        enddo
      enddo

! If the calculation is for g-vectors then no shift; if for k-vectors there could be one.
!     if(iflag.eq.0) then
!       cos_tmp0=1.d0
!       sin_tmp0=0.d0
!      elseif(iflag.eq.1) then
!       dot=0
!       do 25 k=1,3
!  25     dot=dot+g_shift(k)*r(k,ir)
!       cos_tmp0=cos(dot)
!       sin_tmp0=sin(dot)
!      else
!       call fatal_error ('iflag must be 0 or 1 in cossin_psi')
!     endif

!     cos_g(1)=1.d0
!     sin_g(1)=0.d0
      i=0
      do in=1,ngnorm
        do im=1,igmult(in)
        i=i+1
        cos_tmp=cos_gr(igvec(1,i),1)*cos_gr(igvec(2,i),2) &
               -sin_gr(igvec(1,i),1)*sin_gr(igvec(2,i),2)
        sin_tmp=sin_gr(igvec(1,i),1)*cos_gr(igvec(2,i),2) &
               +cos_gr(igvec(1,i),1)*sin_gr(igvec(2,i),2)
        cos_g(i)=cos_tmp*cos_gr(igvec(3,i),3) &
                   -sin_tmp*sin_gr(igvec(3,i),3)
        sin_g(i)=sin_tmp*cos_gr(igvec(3,i),3) &
                   +cos_tmp*sin_gr(igvec(3,i),3)

!     do 30 i=1,ngvec
!       cos_tmp1=cos_tmp0*cos_gr(igvec(1,i),1)
!    &          -sin_tmp0*sin_gr(igvec(1,i),1)
!       sin_tmp1=sin_tmp0*cos_gr(igvec(1,i),1)
!    &          +cos_tmp0*sin_gr(igvec(1,i),1)
!       cos_tmp2=cos_tmp1*cos_gr(igvec(2,i),2)
!    &          -sin_tmp1*sin_gr(igvec(2,i),2)
!       sin_tmp2=sin_tmp1*cos_gr(igvec(2,i),2)
!    &          +cos_tmp1*sin_gr(igvec(2,i),2)
!       cos_g(i)=cos_tmp2*cos_gr(igvec(3,i),3)
!    &             -sin_tmp2*sin_gr(igvec(3,i),3)
!       sin_g(i)=sin_tmp2*cos_gr(igvec(3,i),3)
!    &             +cos_tmp2*sin_gr(igvec(3,i),3)
        do k=1,3
          dcos_g(k,i)=-gvec(k,i)*sin_g(i)
          dsin_g(k,i)= gvec(k,i)*cos_g(i)
        enddo
!       if(i.lt.5) write(ounit,'(''i,gnorm(in),cos_g(i),sin_g(i),dcos_g(k,i),dsin_g(k,i)'',2i5,9d12.4)')
!    & i,gnorm(in),cos_g(i),sin_g(i),(dcos_g(k,i),dsin_g(k,i),k=1,3)
        ddcos_g(i)=-gnorm(in)*gnorm(in)*cos_g(i)
        ddsin_g(i)=-gnorm(in)*gnorm(in)*sin_g(i)
        enddo
      enddo

      return
      end
!-----------------------------------------------------------------------

!     subroutine cossin_psi_k(glatt,gnorm,gvec,igvec,ngvec,r,nr,ng1d,cos_g,sin_g
      subroutine cossin_psi_k(glatt,gnorm,gvec,igvec,ngvec,r,ir,ng1d,cos_g,sin_g &
      ,dcos_g,dsin_g,ddcos_g,ddsin_g,g_shift)
      use system, only: nelec
! Written by Cyrus Umrigar
! Needed for orbitals and their Laplacian.
! For the k-vectors do it straightforwardly since there are few of them

      use precision_kinds, only: dp
      implicit none

      integer :: i, ir, k, ngvec
      integer, dimension(3,*) :: igvec
      integer, dimension(3) :: ng1d
      real(dp) :: cos, dot, sin
      real(dp), dimension(3,3) :: glatt
      real(dp), dimension(*) :: gnorm
      real(dp), dimension(3,*) :: gvec
      real(dp), dimension(3) :: r
      real(dp), dimension(*) :: cos_g
      real(dp), dimension(*) :: sin_g
      real(dp), dimension(3,*) :: dcos_g
      real(dp), dimension(3,*) :: dsin_g
      real(dp), dimension(*) :: ddcos_g
      real(dp), dimension(*) :: ddsin_g
      real(dp), dimension(*) :: g_shift


!     dimension glatt(3,3),gnorm(*),gvec(3,*),igvec(3,*),r(3,*),ng1d(3)
!    &,cos_g(nelec,*),sin_g(nelec,*)
!    &,dcos_g(3,nelec,*),dsin_g(3,nelec,*)
!    &,ddcos_g(nelec,*),ddsin_g(nelec,*),g_shift(*)

!     do 30 ir=1,nr
      do i=1,ngvec
        dot=0
        do k=1,3
          dot=dot+gvec(k,i)*r(k)
        enddo
        cos_g(i)=cos(dot)
        sin_g(i)=sin(dot)
        do k=1,3
          dcos_g(k,i)=-gvec(k,i)*sin_g(i)
          dsin_g(k,i)= gvec(k,i)*cos_g(i)
        enddo
!       if(i.lt.5) write(ounit,'(''i,gnorm(i),cos_g(i),sin_g(i),dcos_g(k,i),dsin_g(k,i)'',2i5,9d12.4)')
!    & i,gnorm(i),cos_g(i),sin_g(i),(dcos_g(k,i),dsin_g(k,i),k=1,3)
        ddcos_g(i)=-gnorm(i)*gnorm(i)*cos_g(i)
        ddsin_g(i)=-gnorm(i)*gnorm(i)*sin_g(i)
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
      20       sin_gr(-n,i,ir)=-sin_gr(n,i,ir)
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

      subroutine cossin_cent(znuc,glatt,igvec,ngvec,r,ng1d,cos_sum,sin_sum)
      use ewald_mod, only: NG1DX
      use system, only: ncent_tot
! Written by Cyrus Umrigar
! Calculate cos_sum and sin_sum for nuclei

      use precision_kinds, only: dp
      implicit none
      integer :: i, k, n, ngvec
      integer, dimension(3,*) :: igvec
      integer, dimension(3) :: ng1d
      real(dp) :: cos, cos_tmp, dot, sin, sin_tmp
      real(dp) :: znuc
      real(dp), dimension(3,3) :: glatt
      real(dp), dimension(3) :: r
      real(dp), dimension(*) :: cos_sum
      real(dp), dimension(*) :: sin_sum
      real(dp), dimension(-NG1DX:NG1DX,3) :: cos_gr
      real(dp), dimension(-NG1DX:NG1DX,3) :: sin_gr



! Calculate cosines and sines for all positions and reciprocal lattice vectors
        do i=1,3
          dot=0.d0
          do k=1,3
            dot=dot+glatt(k,i)*r(k)
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
   20       sin_gr(-n,i)=-sin_gr(n,i)
          enddo
        enddo


      do i=1,ngvec
        cos_sum(i)=0.d0
        sin_sum(i)=0.d0

          cos_tmp=cos_gr(igvec(1,i),1)*cos_gr(igvec(2,i),2) &
                 -sin_gr(igvec(1,i),1)*sin_gr(igvec(2,i),2)
          sin_tmp=sin_gr(igvec(1,i),1)*cos_gr(igvec(2,i),2) &
                 +cos_gr(igvec(1,i),1)*sin_gr(igvec(2,i),2)
          cos_sum(i)=cos_sum(i)+znuc* &
                     (cos_tmp*cos_gr(igvec(3,i),3) &
                     -sin_tmp*sin_gr(igvec(3,i),3))
          sin_sum(i)=sin_sum(i)+znuc* &
                     (sin_tmp*cos_gr(igvec(3,i),3) &
                     +cos_tmp*sin_gr(igvec(3,i),3))

      enddo

      return
      end

!-----------------------------------------------------------------------

      subroutine cossin_p(y_psp,iwctype,glatt,igvec,ngnorm,igmult,r,nr,ng1d,cos_sum,sin_sum)
      use ewald_mod, only: NGNORMX, NG1DX
      use system, only: ncent_tot, nctype_tot
! Written by Cyrus Umrigar
! Calculate cos_sum and sin_sum for pseudopotentials

      use precision_kinds, only: dp
      implicit none

      integer :: i, im, ir, k, n
      integer :: ngnorm, nr
      integer, dimension(*) :: iwctype
      integer, dimension(3,*) :: igvec
      integer, dimension(*) :: igmult
      integer, dimension(3) :: ng1d
      real(dp) :: cos_tmp, dot, sin_tmp
      real(dp), dimension(NGNORMX,nctype_tot) :: y_psp
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
      20       sin_gr(-n,i,ir)=-sin_gr(n,i,ir)
          enddo
        enddo
      enddo

      i=0
      do k=1,ngnorm
        do im=1,igmult(k)
        i=i+1
        cos_sum(i)=0
        sin_sum(i)=0
        do ir=1,nr
          cos_tmp=cos_gr(igvec(1,i),1,ir)*cos_gr(igvec(2,i),2,ir) &
                 -sin_gr(igvec(1,i),1,ir)*sin_gr(igvec(2,i),2,ir)
          sin_tmp=sin_gr(igvec(1,i),1,ir)*cos_gr(igvec(2,i),2,ir) &
                 +cos_gr(igvec(1,i),1,ir)*sin_gr(igvec(2,i),2,ir)
          cos_sum(i)=cos_sum(i)+y_psp(k,iwctype(ir))* &
                     (cos_tmp*cos_gr(igvec(3,i),3,ir) &
                     -sin_tmp*sin_gr(igvec(3,i),3,ir))
          sin_sum(i)=sin_sum(i)+y_psp(k,iwctype(ir))* &
                     (sin_tmp*cos_gr(igvec(3,i),3,ir) &
                     +cos_tmp*sin_gr(igvec(3,i),3,ir))
        enddo
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
      20       sin_gr(-n,i,ir)=-sin_gr(n,i,ir)
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


!---------------------------single electron--------------
      subroutine cossin_1e(glatt,igvec,ngvec,r,ng1d,cos_sum,sin_sum)
      use ewald_mod, only: NG1DX
      use system, only: nelec
! Written by Cyrus Umrigar
! Calculate cos_sum and sin_sum for electrons

      use precision_kinds, only: dp
      implicit none

      integer :: i, k, n, ngvec
      integer :: nr
      integer, dimension(3,*) :: igvec
      integer, dimension(3) :: ng1d
      real(dp) :: cos, cos_tmp, dot, sin, sin_tmp
      real(dp), dimension(3,3) :: glatt
      real(dp), dimension(3) :: r
      real(dp), dimension(*) :: cos_sum
      real(dp), dimension(*) :: sin_sum
      real(dp), dimension(-NG1DX:NG1DX,3) :: cos_gr
      real(dp), dimension(-NG1DX:NG1DX,3) :: sin_gr



! Calculate cosines and sines for all positions and reciprocal lattice vectors
      do i=1,3
         dot=0
         do k=1,3
            dot=dot+glatt(k,i)*r(k)
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
   20       sin_gr(-n,i)=-sin_gr(n,i)
         enddo
      enddo


      do i=1,ngvec
        cos_sum(i)=0
        sin_sum(i)=0

          cos_tmp=cos_gr(igvec(1,i),1)*cos_gr(igvec(2,i),2) &
                 -sin_gr(igvec(1,i),1)*sin_gr(igvec(2,i),2)
          sin_tmp=sin_gr(igvec(1,i),1)*cos_gr(igvec(2,i),2) &
                 +cos_gr(igvec(1,i),1)*sin_gr(igvec(2,i),2)
          cos_sum(i)=cos_sum(i)+ &
                     (cos_tmp*cos_gr(igvec(3,i),3) &
                     -sin_tmp*sin_gr(igvec(3,i),3))
          sin_sum(i)=sin_sum(i)+ &
                     (sin_tmp*cos_gr(igvec(3,i),3) &
                     +cos_tmp*sin_gr(igvec(3,i),3))

      enddo

      return
      end




!----------------------------------add for testing -------------------------------

      subroutine wf(x,wx,b)

      use periodic, only: cutr, np, npoly

      use precision_kinds, only: dp
      implicit none

      integer :: i,ncoef_per
      real(dp)  :: hf,xs, wx, x
      real(dp), dimension(*) :: b

      ncoef_per=npoly+1
      xs=x/cutr
      wx=0.d0

      do i=1,ncoef_per
         wx=wx+b(i)*((xs**i)*((1-xs)**np))
      enddo



      return
      end





!----------------------------------add for testing -------------------------------


! to keep adding
!--------------------



!-----------------------------------------------------------------------
      end module



