      module grid3d_orbitals
c Orbitals on a 3d grid with spline fit
c Written by A. Scemama, adapted from C. Umrigar's 2D routines
      use error, only: fatal_error
      use basis_fns_mod, only: basis_fns
      use grid3d, only: int_from_cart
      interface ! pspline interface
      subroutine fvtricub(ict,ivec,ivecd,
     >   fval,ii,jj,kk,xparam,yparam,zparam,
     >   hx,hxi,hy,hyi,hz,hzi,
     >   fin,inf2,inf3,nz)
        implicit none
        integer ict(10)
        integer ivec  
        integer ivecd  
        integer ii,jj,kk
        real*8 xparam,yparam,zparam
        real*8 hx,hy,hz 
        real*8 hxi,hyi,hzi 
        integer inf2,inf3,nz
        real fin(0:7,inf2,inf3,nz)  
        real*8 fval(10)     
      end subroutine
      subroutine mktricubw(x,nx,y,ny,z,nz,f,nf2,nf3,
     >                    ibcxmin,bcxmin,ibcxmax,bcxmax,inb1x,
     >                    ibcymin,bcymin,ibcymax,bcymax,inb1y,
     >                    ibczmin,bczmin,ibczmax,bczmax,inb1z,
     >                    wk,nwk,ilinx,iliny,ilinz,ier)
      integer nx, ny, nz
      real x(nx), y(ny), z(nz) 
      integer nf2,nf3
      real f(8,nf2,nf3,nz)            
      integer inb1x, inb1y, inb1z
      integer ibcxmin,ibcxmax
      integer ibcymin,ibcymax
      integer ibczmin,ibczmax
      real bcxmin(inb1x,nz),bcxmax(inb1x,nz)
      real bcymin(inb1y,nz),bcymax(inb1y,nz)
      real bczmin(inb1z,ny),bczmax(inb1z,ny)
      integer nwk   
      real wk(nwk)  
      integer ilinx 
      integer iliny 
      integer ilinz 
      integer ier
      end subroutine
      end interface
      
      contains

      subroutine setup_3dsplorb

      use grid_spline_mod, only: orb_num_spl
      use grid_mod, only: MXNSTEP, MXNSTEP3
      use grid_mod, only: cart_from_int
      use atom, only: cent, ncent
      use ghostatom, only: nghostcent
      use phifun, only: d2phin, dphin
      use phifun, only: phin
      use wfsec, only: iwf
      use coefs, only: coef, nbasis, norb
!      use contrl, only: idump, irstar, isite, nconf, nblk, nblkeq, nconf_new, nstep
      use control_vmc, only: vmc_idump, vmc_irstar, vmc_isite, vmc_nconf
      use control_vmc, only: vmc_nblk, vmc_nblkeq, vmc_nconf_new, vmc_nstep
      use phifun, only: d2phin, dphin, phin
      use grid3d_param, only: endpt, nstep3d, origin
      use distance_mod, only: r_en, rvec_en
      use precision_kinds, only: dp
      use contrl_file, only: ounit
      use basis_fns_mod, only: basis_fns
      use system, only: nelec
      implicit none

      integer :: i, ibcxmax, ibcxmin, ibcymax, ibcymin
      integer :: ibczmax, ibczmin, ic, ier
      integer :: ilinx, iliny, ilinz, iok
      integer :: iorb, ix, iy, iz
      integer :: m, memory, nwk
      integer, dimension(4) :: iaxis
      integer, dimension(3) :: ixyz
      real(dp) :: ddf, f, value
      real(dp), dimension(3) :: r
      real(dp), dimension(3) :: df

      real*4  bc(MXNSTEP,MXNSTEP,3:8,nelec/2+1), wk(80*MXNSTEP3)

c     Note:
c     The boundary condition array ranges from 3 to 8. This way, if we code
c     x=1, y=2 and z=3, the 3rd dimension is the sum of the values
c     corresponding to the axes defining the plane.
c     For the maximum boundary, add 3 :
c     xy_min = 1+2 = 3
c     xz_min = 1+3 = 4
c     yz_min = 2+3 = 5
c     xy_max = 1+2+3 = 6
c     xz_max = 1+3+3 = 7
c     yz_max = 2+3+3 = 8

      iwf=1
      iok=1

c     Check the sizes
c      if (norb.gt.MORB_OCC)
c     >  call fatal_error ('MORB_OCC too small. Recompile.')

c     We have no info on the derivatives, so use "not a knot" in the creation
c     of the fit
      ibcxmin=0
      ibcxmax=0
      ibcymin=0
      ibcymax=0
      ibczmin=0
      ibczmax=0

c     Evaluate the energy needed for the calculation
      memory=dfloat((norb+10)*4)
      memory=memory*dfloat(nstep3d(1)*nstep3d(2)*nstep3d(3))
      memory=memory*8.d-6

      write (45,*) 'Allocated memory for the 3D spline fits of the LCAO:',
     & memory, 'Mb'

      if ( vmc_irstar.ne.1 ) then

c      ----------------------------------------------------------------
c      Compute the orbitals values and gradients on the boundary points
c      ----------------------------------------------------------------

       do i=1,3
        iaxis(i)=i
       enddo

c      Permute the values in iaxis 3 times to get xyz, yzx, zxy
c      the 2 first indexes define the plane we compute,
c      the 3rd index defines the orthogonal direction,
c      the 4th index of is a temporary variable.
c      ixyz is used in the loops so we can loop over x, y or z
c      depending on the value of iaxis, in the purpose to use the values
c      ixyz(...) in the access to orb_num_spl.

c      write (45,*) 'Computation of the boundary...'
c      do icount=1,3

c       prepare for the computation of the minimum boundary
c       r(iaxis(3)) = origin(iaxis(3))
c       do i=1,3
c        ixyz(i)=1
c       enddo

c       do ii=0,3,3
c        compute the min boundary at ii=0 and the max boundary at ii=3

c        ixyz(iaxis(1)) = 1
c        do while (ixyz(iaxis(1)).le.nstep3d(iaxis(1)))
c         r(iaxis(1)) = cart_from_int (ixyz(iaxis(1)),iaxis(1))
c
c         ixyz(iaxis(2)) = 1
c         do while (ixyz(iaxis(2)).le.nstep3d(iaxis(2)))
c          r(iaxis(2)) = cart_from_int (ixyz(iaxis(2)),iaxis(2))

c          Calculate e-N inter-particle distances
c          do ic=1,ncent+nghostcent
c           do m=1,3
c            rvec_en(m,1,ic)=r(m)-cent(m,ic)
c           enddo
c           r_en(1,ic)=0.d0
c           do m=1,3
c            r_en(1,ic)=r_en(1,ic)+rvec_en(m,1,ic)**2
c           enddo
c           r_en(1,ic)=dsqrt(r_en(1,ic))
c          enddo
c
c          Calculate the value and the gradient of the orbital
c          call basis_fns(1,1,rvec_en,r_en,1)
c
c          do iorb=1,norb
c           orb_num_spl(1,ixyz(1),ixyz(2),ixyz(3),iorb)=0.d0
c           bc(ixyz(iaxis(1)),ixyz(iaxis(2)),iaxis(1)+iaxis(2)+ii,iorb)=0.d0

c           do m=1,nbasis
c            Value:
c            orb_num_spl(1,ixyz(1),ixyz(2),ixyz(3),iorb)=
c    >           orb_num_spl(1,ixyz(1),ixyz(2),ixyz(3),iorb)+
c    >           +coef(m,iorb,iwf)*phin(m,1)

c            Gradient:
c            bc(ixyz(iaxis(1)),ixyz(iaxis(2)),iaxis(1)+iaxis(2)+ii,iorb)=
c    >        bc(ixyz(iaxis(1)),ixyz(iaxis(2)),iaxis(1)+iaxis(2)+ii,iorb)
c    >          +coef(m,iorb,iwf)*dphin(m,1,iaxis(3))
c           enddo

c          enddo !iorb

c          ixyz(iaxis(2)) = ixyz(iaxis(2))+1
c         enddo ! ixyz

c         ixyz(iaxis(1)) = ixyz(iaxis(1))+1
c        enddo ! ixyz

c        r(iaxis(3)) = endpt(iaxis(3))
c        do i=1,3
c         ixyz(i)=nstep3d(i)
c        enddo
c       enddo ! ii

c       Permutation of the coordinates
c       iaxis(4) = iaxis(1)
c       do i=1,3
c        iaxis(i)=iaxis(i+1)
c       enddo

c      enddo ! icount

c      The derivatives are now computed on the sides of the box,
c      so we can use this info for the creation of the spline fit.
c      ibcxmin=1
c      ibcxmax=1
c      ibcymin=1
c      ibcymax=1
c      ibczmin=1
c      ibczmax=1


c      ----------------------------------------------------------------
c      Compute the orbitals on the remaining points
c      ----------------------------------------------------------------
       write (45,*) 'Computation of the grid points...'
       do ix=1, nstep3d(1)
        r(1) = cart_from_int (ix,1)

        do iy=1, nstep3d(2)
         r(2) = cart_from_int (iy,2)

         do iz=1, nstep3d(3)
          r(3) = cart_from_int (iz,3)


c         Calculate e-N inter-particle distances
          iok=1
          do ic=1,ncent+nghostcent
           do m=1,3
            rvec_en(m,1,ic)=r(m)-cent(m,ic)
           enddo
            r_en(1,ic)=0.d0
           do m=1,3
            r_en(1,ic)=r_en(1,ic)+rvec_en(m,1,ic)**2
           enddo
           r_en(1,ic)=dsqrt(r_en(1,ic))
           if (r_en(1,ic).lt.1.d-6) iok=0
          enddo

c         Check that no atom is exactly on a grid point
          if (iok.eq.0) then
            write(ounit,*) ''
            write(ounit,*) 'There is an atom exactly on one point of the grid.'
            write(ounit,*) 'Resubmit the job with the following parameters:'
 10         format (a7,3(2x,a2,1x,f8.2))
            write(ounit,10) '&3dgrid','x0', origin(1)+.01,
     >                            'y0', origin(2)+.01,
     >                            'z0', origin(3)+.01
            write(ounit,10) '&3dgrid','xn', endpt(1)+.01,
     >                            'yn', endpt(2)+.01,
     >                            'zn', endpt(3)+.01
            write(ounit,*) ''
            call fatal_error('aborted')
          endif

c         Calculate the value of the orbital
          call basis_fns(1,1,rvec_en,r_en,0)

          do iorb=1,norb
           orb_num_spl(1,ix,iy,iz,iorb)=0.d0
           do m=1,nbasis
            orb_num_spl(1,ix,iy,iz,iorb)=
     &        orb_num_spl(1,ix,iy,iz,iorb)+coef(m,iorb,iwf)*phin(m,1)
           enddo

          enddo

         enddo
        enddo
       enddo

       nwk=80*nstep3d(1)*nstep3d(2)*nstep3d(3)
       ier=0
       do iorb=1,norb

c       call r8mktricubw(cart_from_int(1,1),nstep3d(1),
        call mktricubw(cart_from_int(1,1),nstep3d(1),
     >                  cart_from_int(1,2),nstep3d(2),
     >                  cart_from_int(1,3),nstep3d(3),
     >                  orb_num_spl(1,1,1,1,iorb),
     >                  MXNSTEP,MXNSTEP,
     >                  ibcxmin,bc(1,1,2+3,iorb),
     >                  ibcxmax,bc(1,1,2+3+3,iorb),MXNSTEP,
     >                  ibcymin,bc(1,1,1+3,iorb),
     >                  ibcymax,bc(1,1,1+3+3,iorb),MXNSTEP,
     >                  ibczmin,bc(1,1,1+2,iorb),
     >                  ibczmax,bc(1,1,1+2+3,iorb),MXNSTEP,
     >                  wk,nwk,ilinx,iliny,ilinz,ier)
        if (ier.eq.1)
     >   call fatal_error ('Error in r8mktricubw')
        write (45,*) 'orbital ', iorb, 'splined'
       enddo
      endif ! ( irstar .ne. 0 )

c DEBUG
      goto 900
      r(1) = origin(1)
      do i=1,200
        r(2) = 1.
        r(3) = 1.
          do ic=1,ncent+nghostcent
           do m=1,3
            rvec_en(m,1,ic)=r(m)-cent(m,ic)
           enddo
            r_en(1,ic)=0.d0
           do m=1,3
            r_en(1,ic)=r_en(1,ic)+rvec_en(m,1,ic)**2
           enddo
           r_en(1,ic)=dsqrt(r_en(1,ic))
          enddo

          call basis_fns(1,1,rvec_en,r_en,0)

           value=0.
           do m=1,nbasis
            value=value+coef(m,norb,iwf)*phin(m,1)
           enddo

        ddf=0.
        f = 0.
        call spline_mo(r,norb,f,df,ddf,ier)
        print *, r(1), value, f
        r(1) = r(1) + 1./30.
      enddo
      stop 'end'
 900  continue
c DEBUG
      end ! subroutine setup_3dsplorb

c----------------------------------------------------------------------


      subroutine spline_mo(r,iorb,f,df,ddf,ier)
      use grid_spline_mod, only: orb_num_spl
      use grid_mod, only: MXNSTEP
      use grid_mod, only: IUNDEFINED
      use grid_mod, only: cart_from_int
      use insout, only: inout, inside
      use grid3d_param, only: nstep3d, step3d
      implicit none

      integer :: i

c     Input:
      integer   iorb    ! Index of the MO to spline
      real*8    r(3)    ! Cartesian coordinates

c     Output:
      integer   ier     ! error status
      real*8    f, ddf  ! Value, Laplacian
      real*8    df(3)   ! Gradients

c     Work:
      integer   ict(10)   ! Control of the spline subroutine
      integer   ix(3)     ! Integer coordinates
      real*8    fval(10)  ! values returned by the spline routine
      real*8    rscaled(3)! normalized displacement
      real*8    inv_step3d(3) ! Inverse of step sizes


      inout   = inout  +1.d0
      if (ier.eq.1) then
        return
      endif

      ict(1)=1

      do i=1,3
       ix(i) = int_from_cart(r(i),i)
      enddo

      if (  ( ix(1).eq.IUNDEFINED )
     >  .or.( ix(2).eq.IUNDEFINED )
     >  .or.( ix(3).eq.IUNDEFINED ) ) then
        ier = 1
      else

        inside = inside+1.d0
        do i=2,10
         ict(i)=0
        enddo

c       If ddf is equal to 1., compute the laplacian
        do i=1,3
         if (df(i).ne.0.d0) then
          ict(i+1)=1
         endif
        enddo

        if (ddf.ne.0.d0) then
         ict(5)=1
         ict(6)=1
         ict(7)=1
        endif

        do i=1,3
         inv_step3d(i) = 1.d0/step3d(i)
        enddo

        do i=1,3
         rscaled(i) = (r(i)-cart_from_int(ix(i),i))*inv_step3d(i)
        enddo

c       call r8fvtricub(ict,1,1,fval,
        call fvtricub(ict,1,1,fval,
     >                  ix(1),ix(2),ix(3),
     >                  rscaled(1),rscaled(2),rscaled(3),
     >                  step3d(1),inv_step3d(1),
     >                  step3d(2),inv_step3d(2),
     >                  step3d(3),inv_step3d(3),
     >                  orb_num_spl(1,1,1,1,iorb),
     >                  MXNSTEP,MXNSTEP,nstep3d(3))

        f     = fval(1)
        df(1) = fval(2)
        df(2) = fval(3)
        df(3) = fval(4)
        ddf   = fval(5) + fval(6) + fval(7)

      endif

      end ! subroutine spline_mo


c----------------------------------------------------------------------
c
c Lagrange interpolation routines

      subroutine setup_3dlagorb

      use grid_lagrange_mod, only: LAGSTART, LAGEND
      use grid_lagrange_mod, only: orb_num_lag
      use grid_mod, only: cart_from_int
      use vmc_mod, only: norb_tot
      use vmc_mod, only: norb_tot
      use atom, only: cent, ncent
      use wfsec, only: iwf
      use grid3d_param, only: nstep3d, endpt, origin
      use orbital_num_lag, only: denom
      use coefs, only: coef, nbasis, norb
      use ghostatom, only: nghostcent
      use control_vmc, only: vmc_irstar
      use phifun, only: phin, dphin, d2phin
      use distance_mod, only: r_en, rvec_en
      use precision_kinds, only: dp
      use contrl_file, only: ounit
      implicit none

      integer :: i, ic, idenom, ier, iok
      integer :: iorb, ix, iy, iz
      integer :: j, m, memory
      real(dp) :: value, value2
      real(dp), dimension(3) :: r
      real(dp), dimension(norb_tot) :: f




      character*(32) filename
      integer a,b,c

      iwf=1

c     Check the sizes
c      if (norb.gt.MORB_OCC)
c     >  call fatal_error ('MORB_OCC too small. Recompile.')


c     Evaluate the memory needed for the calculation
      memory=dfloat(norb)*5.d0
      memory=memory*dfloat(nstep3d(1)*nstep3d(2)*nstep3d(3))
      memory=memory*4.d-6

      write (45,*) 'Allocated memory for the 3D Lagrange fits of the LCAO:',
     & memory, 'Mb'

      if ( vmc_irstar.ne.1 ) then

c      ----------------------------------------------------------------
c      Compute the orbitals values, gradients and laplacians
c      ----------------------------------------------------------------

       write (45,*) 'Computation of the grid points...'
       do ix=1, nstep3d(1)
        r(1) = cart_from_int (ix,1)

        do iy=1, nstep3d(2)
         r(2) = cart_from_int (iy,2)

         do iz=1, nstep3d(3)
          r(3) = cart_from_int (iz,3)

c         Calculate e-N inter-particle distances
          iok=1
          do ic=1,ncent+nghostcent
           do m=1,3
            rvec_en(m,1,ic)=r(m)-cent(m,ic)
           enddo
            r_en(1,ic)=0.d0
           do m=1,3
            r_en(1,ic)=r_en(1,ic)+rvec_en(m,1,ic)**2
           enddo
           r_en(1,ic)=dsqrt(r_en(1,ic))
           if (r_en(1,ic).lt.1.d-6) iok=0
          enddo

c         Check that no atom is exactly on a grid point
          if (iok.eq.0) then
            write(ounit,*) ''
            write(ounit,*) 'There is an atom exactly on one point of the grid.'
            write(ounit,*) 'Resubmit the job with the following parameters:'
 10         format (a7,3(2x,a2,1x,f8.2))
            write(ounit,10) '&3dgrid','x0', origin(1)+.01,
     >                            'y0', origin(2)+.01,
     >                            'z0', origin(3)+.01
            write(ounit,10) '&3dgrid','xn', endpt(1)+.01,
     >                            'yn', endpt(2)+.01,
     >                            'zn', endpt(3)+.01
            write(ounit,*) ''
            call fatal_error('aborted')
          endif

c         Calculate the grids
          call basis_fns(1,1,rvec_en,r_en,2)

          do iorb=1,norb
           orb_num_lag(1,ix,iy,iz,iorb)=0.d0
           orb_num_lag(2,ix,iy,iz,iorb)=0.d0
           orb_num_lag(3,ix,iy,iz,iorb)=0.d0
           orb_num_lag(4,ix,iy,iz,iorb)=0.d0
           orb_num_lag(5,ix,iy,iz,iorb)=0.d0
           do m=1,nbasis
            orb_num_lag(1,ix,iy,iz,iorb)=
     >        orb_num_lag(1,ix,iy,iz,iorb)+coef(m,iorb,iwf)*phin(m,1)
            orb_num_lag(2,ix,iy,iz,iorb)=
     >        orb_num_lag(2,ix,iy,iz,iorb)+coef(m,iorb,iwf)*dphin(m,1,1)
            orb_num_lag(3,ix,iy,iz,iorb)=
     >        orb_num_lag(3,ix,iy,iz,iorb)+coef(m,iorb,iwf)*dphin(m,1,2)
            orb_num_lag(4,ix,iy,iz,iorb)=
     >        orb_num_lag(4,ix,iy,iz,iorb)+coef(m,iorb,iwf)*dphin(m,1,3)
            orb_num_lag(5,ix,iy,iz,iorb)=
     >        orb_num_lag(5,ix,iy,iz,iorb)+coef(m,iorb,iwf)*d2phin(m,1)
           enddo

          enddo

         enddo
        enddo
       enddo
      endif
c DEBUG
c      do igrid=1,5
c       do iorb=1,norb
c        do ix=1,nstep3d(1)
c         do iy=1,nstep3d(2)
c          do iz=1,nstep3d(3)
c           grid3d(ix,iy,iz) = orb_num_lag(igrid,ix,iy,iz,iorb)
c          enddo
c         enddo
c        enddo
c        write(filename,'(i1,1a,i2.2,5a)') igrid,'.',iorb,'.cube'
c        print *, filename
c        call write_cube(filename)
c       enddo
c      enddo
c      stop
c DEBUG

       do ix=1,3
        do i=LAGSTART,LAGEND
         idenom = 1
         do j=LAGSTART,LAGEND
          if (i.ne.j)
     >     idenom = idenom*(i-j)
         enddo
         denom(i,ix) = 1.d0/dfloat(idenom)
        enddo
       enddo

c DEBUG
      goto 900
      a=1
      b=2
      c=3
      r(a) = origin(a)
      do i=1,200
        r(b) = 1.
        r(c) = 1.
          do ic=1,ncent+nghostcent
           do m=1,3
            rvec_en(m,1,ic)=r(m)-cent(m,ic)
           enddo
            r_en(1,ic)=0.d0
           do m=1,3
            r_en(1,ic)=r_en(1,ic)+rvec_en(m,1,ic)**2
           enddo
           r_en(1,ic)=dsqrt(r_en(1,ic))
          enddo

          call basis_fns(1,1,rvec_en,r_en,0)

          value=0.
          do j=1,norb
           value2=0.
           do m=1,nbasis
            value2=value2 + coef(m,j,iwf)*phin(m,1)
           enddo
           value = value + value2
          enddo

        ier=0
        call lagrange_mose(1,r,f,ier)
        value2=0.
        do j=1,norb
         value2 = value2+f(j)
        enddo
        print *, r(a), value, value2
        r(a) = r(a) + 1./30.
      enddo
      stop 'end'
 900  continue
c DEBUG

      end

c-----------------------------------------------------------------------

      subroutine lagrange_mos(igrid,r,orb,iel,ier)
c Written by A Scemama
c Evaluate orbitals by a Lagrange interpolation (LAGMAX mesh pts)
c on an equally-spaced 3D grid.
c The mesh pts. on which the function values, f, are given, are assumed
c to be at 1,2,3,...nstep3d(1), and similarly for y and z.
c

      use grid_lagrange_mod, only: LAGMAX, LAGSTART, LAGEND
      use grid_lagrange_mod, only: orb_num_lag
      use grid_mod, only: cart_from_int
      use vmc_mod, only: norb_tot
      use insout, only: inout, inside
      use coefs, only: norb
      use grid3d_param, only: nstep3d, step3d
      use orbital_num_lag, only: denom

      use precision_kinds, only: dp
      use system, only: nelec
      implicit none

      integer :: i, i1, i2, i3, iel
      integer :: ier, igrid, iorb, j
      integer, dimension(3) :: ix
      real(dp) :: orb1, orb2
      real(dp), dimension(3) :: r
      real(dp), dimension(3) :: dr
      real(dp), dimension(nelec,norb_tot) :: orb
      real(dp), dimension(LAGSTART:LAGEND) :: xi



      real*8    num(LAGSTART:LAGEND,3)

      inout   = inout  +1.d0
      if (ier.eq.1) then
        return
      endif

      do i=1,3
       ix(i) = int_from_cart(r(i),i)
      enddo

      if (  ( ix(1).le.LAGMAX/2 )
     >  .or.( ix(2).le.LAGMAX/2 )
     >  .or.( ix(3).le.LAGMAX/2 )
     >  .or.( ix(1).ge.nstep3d(1)-LAGMAX/2 )
     >  .or.( ix(2).ge.nstep3d(2)-LAGMAX/2 )
     >  .or.( ix(3).ge.nstep3d(3)-LAGMAX/2 ) ) then
        ier = 1
      else
        inside = inside+1.d0
c Compute displacements
       do i=1,3
        dr(i) = (r(i)-cart_from_int(ix(i),i))/step3d(i)
       enddo

       do i1=1,3
         do i2=LAGSTART, LAGEND
           num(i2,i1)=denom(i2,i1)
           do j=LAGSTART, LAGEND
             if (j.ne.i2)
     >         num(i2,i1) = num(i2,i1)*(dr(i1) - dfloat(j))
           enddo
         enddo
       enddo

       do iorb=1,norb
         orb(iel,iorb) = 0.d0
         do i3=LAGSTART, LAGEND !z interpolation
           orb2 = 0.d0
           do i2=LAGSTART, LAGEND !y interpolation
             orb1 = 0.d0
             do i1=LAGSTART, LAGEND !x interpolation
               orb1 = orb1 + num(i1,1) *
     >           orb_num_lag(igrid,ix(1)+i1,ix(2)+i2,ix(3)+i3,iorb)
             enddo !x interpolation
             orb2 = orb2 + num(i2,2) * orb1
           enddo !y interpolation
           orb(iel,iorb) = orb(iel,iorb) + num(i3,3) * orb2
         enddo !z interpolation
       enddo ! iorb

      endif

      end

c----------------------------------------------------------------------

      subroutine lagrange_mos_2(igrid,r,orb,iel,ier)
c     Written by A Scemama
c     Evaluate orbitals by a Lagrange interpolation (LAGMAX mesh pts)
c     on an equally-spaced 3D grid.
c     The mesh pts. on which the function values, f, are given, are assumed
c     to be at 1,2,3,...nstep3d(1), and similarly for y and z.
c     
      
      use grid_lagrange_mod, only: LAGMAX, LAGSTART, LAGEND
      use grid_lagrange_mod, only: orb_num_lag
      use grid_mod, only: cart_from_int
      use vmc_mod, only: norb_tot
      use insout, only: inout, inside
      use coefs, only: norb
      use grid3d_param, only: nstep3d, step3d
      use orbital_num_lag, only: denom
      
      use precision_kinds, only: dp
      use system, only: nelec
      implicit none
      
      integer :: i, i1, i2, i3, iel
      integer :: ier, igrid, iorb, j
      integer, dimension(3) :: ix
      real(dp) :: orb1, orb2
      real(dp), dimension(3) :: r
      real(dp), dimension(3) :: dr
      real(dp), dimension(norb_tot,nelec) :: orb
      real(dp), dimension(LAGSTART:LAGEND) :: xi
      
      

      real*8    num(LAGSTART:LAGEND,3)

      inout   = inout  +1.d0
      if (ier.eq.1) then
         return
      endif
      
      do i=1,3
       ix(i) = int_from_cart(r(i),i)
      enddo

      if (  ( ix(1).le.LAGMAX/2 )
     >  .or.( ix(2).le.LAGMAX/2 )
     >  .or.( ix(3).le.LAGMAX/2 )
     >  .or.( ix(1).ge.nstep3d(1)-LAGMAX/2 )
     >  .or.( ix(2).ge.nstep3d(2)-LAGMAX/2 )
     >  .or.( ix(3).ge.nstep3d(3)-LAGMAX/2 ) ) then
        ier = 1
      else
        inside = inside+1.d0
c Compute displacements
       do i=1,3
        dr(i) = (r(i)-cart_from_int(ix(i),i))/step3d(i)
       enddo

       do i1=1,3
         do i2=LAGSTART, LAGEND
           num(i2,i1)=denom(i2,i1)
           do j=LAGSTART, LAGEND
             if (j.ne.i2)
     >         num(i2,i1) = num(i2,i1)*(dr(i1) - dfloat(j))
           enddo
         enddo
       enddo

       do iorb=1,norb
         orb(iorb,iel) = 0.d0
         do i3=LAGSTART, LAGEND !z interpolation
            orb2 = 0.d0
           do i2=LAGSTART, LAGEND !y interpolation
              orb1 = 0.d0
              do i1=LAGSTART, LAGEND !x interpolation
                 orb1 = orb1 + num(i1,1) *
     >                orb_num_lag(igrid,ix(1)+i1,ix(2)+i2,ix(3)+i3,iorb)
              enddo             !x interpolation
              orb2 = orb2 + num(i2,2) * orb1
           enddo                !y interpolation
           orb(iorb,iel) = orb(iorb,iel) + num(i3,3) * orb2
        enddo                   !z interpolation
      enddo                     ! iorb
      
      endif
      
      end

c----------------------------------------------------------------------

      subroutine lagrange_mos_grad(igrid,r,orb,iel,ier)
c Written by A Scemama
c Evaluate orbitals by a Lagrange interpolation (LAGMAX mesh pts)
c on an equally-spaced 3D grid.
c The mesh pts. on which the function values, f, are given, are assumed
c to be at 1,2,3,...nstep3d(1), and similarly for y and z.
c

      use grid_lagrange_mod, only: LAGMAX, LAGSTART, LAGEND
      use grid_lagrange_mod, only: orb_num_lag
      use grid_mod, only: cart_from_int
      use vmc_mod, only: norb_tot
      use insout, only: inout, inside
      use coefs, only: norb
      use grid3d_param, only: nstep3d, step3d
      use orbital_num_lag, only: denom

      use precision_kinds, only: dp
      use system, only: nelec
      implicit none

      integer :: i, i1, i2, i3, iaxis
      integer :: iel, ier, igrid, iorb
      integer :: j
      integer, dimension(3) :: ix
      real(dp) :: orb1, orb2
      real(dp), dimension(3) :: r
      real(dp), dimension(3) :: dr
      real(dp), dimension(norb_tot,nelec,3) :: orb
      real(dp), dimension(LAGSTART:LAGEND) :: xi


      real*8    num(LAGSTART:LAGEND,3)

      inout   = inout  +1.d0

      if (ier.eq.1) then
        return
      endif

      do i=1,3
       ix(i) = int_from_cart(r(i),i)
      enddo

      if (  ( ix(1).le.LAGMAX/2 )
     >  .or.( ix(2).le.LAGMAX/2 )
     >  .or.( ix(3).le.LAGMAX/2 )
     >  .or.( ix(1).ge.nstep3d(1)-LAGMAX/2 )
     >  .or.( ix(2).ge.nstep3d(2)-LAGMAX/2 )
     >  .or.( ix(3).ge.nstep3d(3)-LAGMAX/2 ) ) then
        ier = 1
      else
        inside = inside+1.d0
c Compute displacements
       do i=1,3
        dr(i) = (r(i)-cart_from_int(ix(i),i))/step3d(i)
       enddo

       do i1=1,3
         do i2=LAGSTART, LAGEND
           num(i2,i1)=denom(i2,i1)
           do j=LAGSTART, LAGEND
             if (j.ne.i2)
     >         num(i2,i1) = num(i2,i1)*(dr(i1) - dfloat(j))
           enddo
         enddo
       enddo

       iaxis = igrid-1
       do iorb=1,norb
         orb(iorb,iel,iaxis) = 0.d0
         do i3=LAGSTART, LAGEND !z interpolation
           orb2 = 0.d0
           do i2=LAGSTART, LAGEND !y interpolation
             orb1 = 0.d0
             do i1=LAGSTART, LAGEND !x interpolation
               orb1 = orb1 + num(i1,1) *
     >           orb_num_lag(igrid,ix(1)+i1,ix(2)+i2,ix(3)+i3,iorb)
             enddo !x interpolation
             orb2 = orb2 + num(i2,2) * orb1
           enddo !y interpolation
           orb(iorb,iel,iaxis) = orb(iorb,iel,iaxis) + num(i3,3) * orb2
         enddo !z interpolation
       enddo ! iorb

      endif

      end

c-----------------------------------------------------------------------

      subroutine lagrange_mose(igrid,r,orb,ier)
c Written by A Scemama
c Evaluate orbitals by a Lagrange interpolation (LAGMAX mesh pts)
c on an equally-spaced 3D grid.
c The mesh pts. on which the function values, f, are given, are assumed
c to be at 1,2,3,...nstep3d(1), and similarly for y and z.
c

      use grid_lagrange_mod, only: LAGMAX, LAGSTART, LAGEND
      use grid_lagrange_mod, only: orb_num_lag
      use grid_mod, only: cart_from_int
      use vmc_mod, only: norb_tot
      use insout, only: inout, inside
      use coefs, only: norb
      use grid3d_param, only: nstep3d, step3d
      use orbital_num_lag, only: denom

      use precision_kinds, only: dp
      implicit none

      integer :: i, i1, i2, i3, ier
      integer :: igrid, iorb, j
      integer, dimension(3) :: ix
      real(dp) :: orb1, orb2
      real(dp), dimension(3) :: r
      real(dp), dimension(3) :: dr
      real(dp), dimension(norb_tot) :: orb
      real(dp), dimension(LAGSTART:LAGEND) :: xi



      real*8    num(LAGSTART:LAGEND,3)

      inout   = inout  +1.d0
      if (ier.eq.1) then
        return
      endif

      do i=1,3
       ix(i) = int_from_cart(r(i),i)
      enddo

      if (  ( ix(1).le.LAGMAX/2 )
     >  .or.( ix(2).le.LAGMAX/2 )
     >  .or.( ix(3).le.LAGMAX/2 )
     >  .or.( ix(1).ge.nstep3d(1)-LAGMAX/2 )
     >  .or.( ix(2).ge.nstep3d(2)-LAGMAX/2 )
     >  .or.( ix(3).ge.nstep3d(3)-LAGMAX/2 ) ) then
        ier = 1
      else
        inside = inside+1.d0
c Compute displacements
       do i=1,3
        dr(i) = (r(i)-cart_from_int(ix(i),i))/step3d(i)
       enddo

       do i1=1,3
         do i2=LAGSTART, LAGEND
           num(i2,i1)=denom(i2,i1)
           do j=LAGSTART, LAGEND
             if (j.ne.i2)
     >         num(i2,i1) = num(i2,i1)*(dr(i1) - dfloat(j))
           enddo
         enddo
       enddo

       do iorb=1,norb
         orb(iorb) = 0.d0
         do i3=LAGSTART, LAGEND !z interpolation
           orb2 = 0.d0
           do i2=LAGSTART, LAGEND !y interpolation
             orb1 = 0.d0
             do i1=LAGSTART, LAGEND !x interpolation
               orb1 = orb1 + num(i1,1) *
     >           orb_num_lag(igrid,ix(1)+i1,ix(2)+i2,ix(3)+i3,iorb)
             enddo !x interpolation
             orb2 = orb2 + num(i2,2) * orb1
           enddo !y interpolation
           orb(iorb) = orb(iorb) + num(i3,3) * orb2
         enddo !z interpolation
       enddo ! iorb

      endif

      end

c----------------------------------------------------------------------

      subroutine lagrange_mos_grade(igrid,r,orb,ier)
c Written by A Scemama
c Evaluate orbitals by a Lagrange interpolation (LAGMAX mesh pts)
c on an equally-spaced 3D grid.
c The mesh pts. on which the function values, f, are given, are assumed
c to be at 1,2,3,...nstep3d(1), and similarly for y and z.
c

      use grid_lagrange_mod, only: LAGMAX, LAGSTART, LAGEND
      use grid_lagrange_mod, only: orb_num_lag
      use grid_mod, only: cart_from_int
      use vmc_mod, only: norb_tot
      use insout, only: inout, inside
      use coefs, only: norb
      use grid3d_param, only: nstep3d, step3d
      use orbital_num_lag, only: denom

      use precision_kinds, only: dp
      implicit none

      integer :: i, i1, i2, i3, iaxis
      integer :: ier, igrid, iorb, j
      integer, dimension(3) :: ix
      real(dp) :: orb1, orb2
      real(dp), dimension(3) :: r
      real(dp), dimension(3) :: dr
      real(dp), dimension(norb_tot, 3) :: orb
      real(dp), dimension(LAGSTART:LAGEND) :: xi


      real*8    num(LAGSTART:LAGEND,3)

      inout   = inout  +1.d0
      if (ier.eq.1) then
        return
      endif

      do i=1,3
       ix(i) = int_from_cart(r(i),i)
      enddo

      if (  ( ix(1).le.LAGMAX/2 )
     >  .or.( ix(2).le.LAGMAX/2 )
     >  .or.( ix(3).le.LAGMAX/2 )
     >  .or.( ix(1).ge.nstep3d(1)-LAGMAX/2 )
     >  .or.( ix(2).ge.nstep3d(2)-LAGMAX/2 )
     >  .or.( ix(3).ge.nstep3d(3)-LAGMAX/2 ) ) then
        ier = 1
      else
        inside = inside+1.d0
c Compute displacements
       do i=1,3
        dr(i) = (r(i)-cart_from_int(ix(i),i))/step3d(i)
       enddo

       iaxis = igrid-1
       do i1=1,3
         do i2=LAGSTART, LAGEND
           num(i2,i1)=denom(i2,i1)
           do j=LAGSTART, LAGEND
             if (j.ne.i2)
     >         num(i2,i1) = num(i2,i1)*(dr(i1) - dfloat(j))
           enddo
         enddo
       enddo

       do iorb=1,norb
         orb(iorb,iaxis) = 0.d0
         do i3=LAGSTART, LAGEND !z interpolation
           orb2 = 0.d0
           do i2=LAGSTART, LAGEND !y interpolation
             orb1 = 0.d0
             do i1=LAGSTART, LAGEND !x interpolation
               orb1 = orb1 + num(i1,1) *
     >           orb_num_lag(igrid,ix(1)+i1,ix(2)+i2,ix(3)+i3,iorb)
             enddo !x interpolation
             orb2 = orb2 + num(i2,2) * orb1
           enddo !y interpolation
           orb(iorb,iaxis) = orb(iorb,iaxis) + num(i3,3) * orb2
         enddo !z interpolation
       enddo ! iorb

      endif

      end
c-----------------------------------------------------------------------
      subroutine orb3d_dump(iu)
      use grid_mod, only: cart_from_int
      use coefs, only: norb
      use grid3d_param, only: endpt, nstep3d, origin, step3d
      use grid3dflag, only: i3dgrid, i3dlagorb, i3dsplorb

      implicit none

      integer :: i, iu, j





      if (i3dgrid.eq.0) return

      write (iu) norb
      write (iu) i3dsplorb, i3dlagorb
      write (iu) (origin(i), i=1,3)
      write (iu) (endpt(i), i=1,3)
      write (iu) (nstep3d(i), i=1,3)
      write (iu) (step3d(i), i=1,3)
      write (iu) ((cart_from_int(i,j), i=1,nstep3d(j)),j=1,3)

      if (i3dsplorb.ge.1) call splorb_dump(iu)
      if (i3dlagorb.ge.1) call lagorb_dump(iu)

      end
c-----------------------------------------------------------------------
      subroutine orb3d_rstrt(iu)

      use grid_mod, only: cart_from_int
      use coefs, only: norb

      use grid3d_param, only: endpt, nstep3d, origin, step3d
      use grid3dflag, only: i3dgrid, i3dlagorb, i3dsplorb

      implicit none

      integer :: i, iu, j





      if (i3dgrid.eq.0) return

      read (iu) norb
      read (iu) i3dsplorb, i3dlagorb
      read (iu) (origin(i), i=1,3)
      read (iu) (endpt(i), i=1,3)
      read (iu) (nstep3d(i), i=1,3)
      read (iu) (step3d(i), i=1,3)
      read (iu) ((cart_from_int(i,j), i=1,nstep3d(j)),j=1,3)

      if (i3dsplorb.ge.1) call splorb_rstrt(iu)
      if (i3dlagorb.ge.1) call lagorb_rstrt(iu)
      end
c-----------------------------------------------------------------------
      subroutine splorb_dump(iu)
      use grid_spline_mod, only: orb_num_spl
      use coefs, only: norb
      use grid3d_param, only: nstep3d
      implicit none

      integer :: i, iu, j, k, l
      integer :: m




      do i=1,8
       do m=1,norb
        write (iu) (((orb_num_spl(i,j,k,l,m), j=1,nstep3d(1)),
     >                k=1,nstep3d(2)), l=1,nstep3d(3))
       enddo
      enddo

      end
c-----------------------------------------------------------------------
      subroutine splorb_rstrt(iu)
      use grid_spline_mod, only: orb_num_spl
      use coefs, only: norb
      use grid3d_param, only: nstep3d
      implicit none

      integer :: i, iu, j, k, l
      integer :: m




      do i=1,8
       do m=1,norb
        read  (iu) (((orb_num_spl(i,j,k,l,m), j=1,nstep3d(1)),
     >                k=1,nstep3d(2)), l=1,nstep3d(3))
       enddo
      enddo
      end
c-----------------------------------------------------------------------
      subroutine lagorb_dump(iu)

      use grid_lagrange_mod, only: orb_num_lag
      use coefs, only: norb
      use grid3d_param, only: nstep3d

      implicit none

      integer :: i, iu, j, k, l
      integer :: m




      do i=1,5
       do m=1,norb
        write (iu) (((orb_num_lag(i,j,k,l,m), j=1,nstep3d(1)),
     >                k=1,nstep3d(2)), l=1,nstep3d(3))
       enddo
      enddo

      end
c-----------------------------------------------------------------------
      subroutine lagorb_rstrt(iu)

      use grid_lagrange_mod, only: orb_num_lag
      use coefs, only: norb
      use grid3d_param, only: nstep3d

      implicit none

      integer :: i, iu, j, k, l
      integer :: m




      do i=1,5
       do m=1,norb
        read  (iu) (((orb_num_lag(i,j,k,l,m), j=1,nstep3d(1)),
     >                k=1,nstep3d(2)), l=1,nstep3d(3))
       enddo
      enddo

      end
c-----------------------------------------------------------------------
      end module 
