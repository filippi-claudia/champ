c Orbitals on a 3d grid with spline fit
c Written by A. Scemama, adapted from C. Umrigar's 2D routines

      subroutine setup_3dsplorb

      use atom, only: cent, ncent

      use ghostatom, only: newghostype, nghostcent
      use ghostatom, only: newghostype, nghostcent
      use phifun, only: d2phin, d2phin_all, d3phin, dphin, n0_ibasis, n0_ic, n0_nbasis,
     &phin
      use wfsec, only: iwf, iwftype, nwftype
      use coefs, only: coef, nbasis, norb
      use contrl, only: idump, irstar, isite, n_conf, nblk, nblkeq, nconf_new, nstep

      use coefs, only: coef, nbasis, norb
      use phifun, only: d2phin, dphin, phin

      use contrl, only: idump, irstar, isite, n_conf, nblk, nblkeq, nconf_new, nstep
      implicit real*8(a-h,o-z)




      include 'force.h'
      include 'vmc.h'
      include '3dgrid.h'
      include '3dgrid_flags.h'
      include '3dgrid_spline.h'

      common /distance/ rshift(3,MELEC,MCENT),rvec_en(3,MELEC,MCENT)
     &,r_en(MELEC,MCENT),rvec_ee(3,MMAT_DIM2),r_ee(MMAT_DIM2)

      real*4  bc(MXNSTEP,MXNSTEP,3:8,MELEC/2+1), wk(80*MXNSTEP3)
      common /orbital_num_spl2/ bc, wk

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


      dimension r(3), iaxis(4), ixyz(3)

      dimension df(3)
      iwf=1
      iok=1

c     Check the sizes
      if (norb.gt.MORB_OCC) 
     >  call fatal_error ('MORB_OCC too small. Recompile.')

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

      if ( irstar.ne.1 ) then

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
c          call basis_fnse(1,rvec_en,r_en)
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
c    >          +coef(m,iorb,iwf)*dphin(iaxis(3),m,1)
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
            write(6,*) ''
            write(6,*) 'There is an atom exactly on one point of the grid.'
            write(6,*) 'Resubmit the job with the following parameters:'
 10         format (a7,3(2x,a2,1x,f8.2))
            write(6,10) '&3dgrid','x0', origin(1)+.01,
     >                            'y0', origin(2)+.01,
     >                            'z0', origin(3)+.01
            write(6,10) '&3dgrid','xn', endpt(1)+.01,
     >                            'yn', endpt(2)+.01,
     >                            'zn', endpt(3)+.01
            write(6,*) ''
            call fatal_error('aborted')
          endif

c         Calculate the value of the orbital 
          call basis_fnse_v(1,rvec_en,r_en)

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
       end do
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

          call basis_fnse_v(1,rvec_en,r_en)

           value=0.
           do m=1,nbasis
            value=value+
     &        coef(m,norb,iwf)*phin(m,1)
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
      use insout, only: inout, inside
      implicit real*8(a-h,o-z)
      include 'vmc.h'
      include '3dgrid.h'
      include '3dgrid_spline.h'

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

      use atom, only: cent, ncent
      use wfsec, only: iwf, iwftype, nwftype

      implicit real*8(a-h,o-z)

      include 'force.h'
      include 'vmc.h'
      include '3dgrid.h'
      include '3dgrid_flags.h'
      include '3dgrid_lagrange.h'

      common /distance/ rshift(3,MELEC,MCENT),rvec_en(3,MELEC,MCENT)
     &,r_en(MELEC,MCENT),rvec_ee(3,MMAT_DIM2),r_ee(MMAT_DIM2)

      character*(32) filename
      integer a,b,c
      dimension r(3)

      dimension f(MORB)
      iwf=1

c     Check the sizes
      if (norb.gt.MORB_OCC) 
     >  call fatal_error ('MORB_OCC too small. Recompile.')


c     Evaluate the memory needed for the calculation
      memory=dfloat(norb)*5.d0
      memory=memory*dfloat(nstep3d(1)*nstep3d(2)*nstep3d(3))
      memory=memory*4.d-6
     
      write (45,*) 'Allocated memory for the 3D Lagrange fits of the LCAO:',
     & memory, 'Mb'

      if ( irstar.ne.1 ) then

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
            write(6,*) ''
            write(6,*) 'There is an atom exactly on one point of the grid.'
            write(6,*) 'Resubmit the job with the following parameters:'
 10         format (a7,3(2x,a2,1x,f8.2))
            write(6,10) '&3dgrid','x0', origin(1)+.01,
     >                            'y0', origin(2)+.01,
     >                            'z0', origin(3)+.01
            write(6,10) '&3dgrid','xn', endpt(1)+.01,
     >                            'yn', endpt(2)+.01,
     >                            'zn', endpt(3)+.01
            write(6,*) ''
            call fatal_error('aborted')
          endif
         
c         Calculate the grids
          call basis_fnse_vgl(1,rvec_en,r_en)

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
     >        orb_num_lag(2,ix,iy,iz,iorb)+coef(m,iorb,iwf)*dphin(1,m,1)
            orb_num_lag(3,ix,iy,iz,iorb)=
     >        orb_num_lag(3,ix,iy,iz,iorb)+coef(m,iorb,iwf)*dphin(2,m,1)
            orb_num_lag(4,ix,iy,iz,iorb)=
     >        orb_num_lag(4,ix,iy,iz,iorb)+coef(m,iorb,iwf)*dphin(3,m,1)
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

          call basis_fnse_vgl(1,rvec_en,r_en)

          value=0.
          do j=1,norb
           value2=0.
           do m=1,nbasis
            value2=value2 + 
     &        coef(m,j,iwf)*phin(m,1)
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

      use insout, only: inout, inside
      use coefs, only: coef, nbasis, norb
      implicit real*8(a-h,o-z)


      include 'vmc.h'
      include 'force.h'
      include '3dgrid.h'
      include '3dgrid_lagrange.h'


      dimension r(3),dr(3),orb(MELEC,MORB),ix(3)
      dimension xi(LAGSTART:LAGEND)
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

      subroutine lagrange_mos_grad(igrid,r,orb,iel,ier)
c Written by A Scemama
c Evaluate orbitals by a Lagrange interpolation (LAGMAX mesh pts) 
c on an equally-spaced 3D grid.
c The mesh pts. on which the function values, f, are given, are assumed
c to be at 1,2,3,...nstep3d(1), and similarly for y and z.
c

      use insout, only: inout, inside
      use coefs, only: coef, nbasis, norb
      implicit real*8(a-h,o-z)


      include 'vmc.h'
      include 'force.h'
      include '3dgrid.h'
      include '3dgrid_lagrange.h'


      dimension r(3),dr(3),orb(3,MELEC,MORB),ix(3)
      dimension xi(LAGSTART:LAGEND)
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
         orb(iaxis,iel,iorb) = 0.d0
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
           orb(iaxis,iel,iorb) = orb(iaxis,iel,iorb) + num(i3,3) * orb2
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

      use insout, only: inout, inside
      use coefs, only: coef, nbasis, norb
      implicit real*8(a-h,o-z)


      include 'vmc.h'
      include 'force.h'
      include '3dgrid.h'
      include '3dgrid_lagrange.h'


      dimension r(3),dr(3),orb(MORB),ix(3)
      dimension xi(LAGSTART:LAGEND)
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

      use insout, only: inout, inside
      use coefs, only: coef, nbasis, norb
      implicit real*8(a-h,o-z)

      include 'vmc.h'
      include 'force.h'
      include '3dgrid.h'
      include '3dgrid_lagrange.h'

      dimension r(3),dr(3),orb(3,MORB),ix(3)
      dimension xi(LAGSTART:LAGEND)
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
         orb(iaxis,iorb) = 0.d0
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
           orb(iaxis,iorb) = orb(iaxis,iorb) + num(i3,3) * orb2
         enddo !z interpolation
       enddo ! iorb

      endif

      end
c-----------------------------------------------------------------------
      subroutine orb3d_dump(iu)
      use coefs, only: coef, nbasis, norb
      implicit real*8(a-h,o-z)

      include 'vmc.h'
      include 'force.h'
      include '3dgrid.h'
      include '3dgrid_flags.h'

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

      use coefs, only: coef, nbasis, norb

      implicit real*8(a-h,o-z)

      include 'vmc.h'
      include 'force.h'
      include '3dgrid.h'
      include '3dgrid_flags.h'

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
      use coefs, only: coef, nbasis, norb
      implicit real*8(a-h,o-z)

      include 'vmc.h'
      include '3dgrid.h'
      include '3dgrid_spline.h'
      include 'force.h'
 
      do i=1,8
       do m=1,norb
        write (iu) (((orb_num_spl(i,j,k,l,m), j=1,nstep3d(1)),
     >                k=1,nstep3d(2)), l=1,nstep3d(3))
       enddo
      enddo

      end
c-----------------------------------------------------------------------
      subroutine splorb_rstrt(iu)
      use coefs, only: coef, nbasis, norb
      implicit real*8(a-h,o-z)

      include 'vmc.h'
      include '3dgrid.h'
      include '3dgrid_spline.h'
      include 'force.h'

      do i=1,8
       do m=1,norb
        read  (iu) (((orb_num_spl(i,j,k,l,m), j=1,nstep3d(1)),
     >                k=1,nstep3d(2)), l=1,nstep3d(3))
       enddo
      enddo
      end
c-----------------------------------------------------------------------
      subroutine lagorb_dump(iu)
      use coefs, only: coef, nbasis, norb
      implicit real*8(a-h,o-z)

      include 'vmc.h'
      include '3dgrid.h'
      include '3dgrid_lagrange.h'
      include 'force.h'
 
      do i=1,5
       do m=1,norb
        write (iu) (((orb_num_lag(i,j,k,l,m), j=1,nstep3d(1)),
     >                k=1,nstep3d(2)), l=1,nstep3d(3))
       enddo
      enddo

      end
c-----------------------------------------------------------------------
      subroutine lagorb_rstrt(iu)
      use coefs, only: coef, nbasis, norb
      implicit real*8(a-h,o-z)

      include 'vmc.h'
      include '3dgrid.h'
      include '3dgrid_lagrange.h'
      include 'force.h'

      do i=1,5
       do m=1,norb
        read  (iu) (((orb_num_lag(i,j,k,l,m), j=1,nstep3d(1)),
     >                k=1,nstep3d(2)), l=1,nstep3d(3))
       enddo
      enddo

      end
c-----------------------------------------------------------------------
