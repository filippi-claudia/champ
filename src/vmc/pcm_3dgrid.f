c 3d grid module
c Created by Amovilli/Floris following subroutines by Scemama
c Example of input:
c
c &pcm nx_pcm 30 ny_pcm 30 nz_pcm 30
c &pcm dx_pcm .1 dy_pcm .1 dz_pcm .1
c
c The following are computed if they are not present in the input.
c &pcm x0_pcm 0. y0_pcm 0. z0_pcm 0.
c &pcm xn_pcm 1. yn_pcm 1. zn_pcm 1.
c----------------------------------------------------------------------
      subroutine pcm_setup_grid
      implicit real*8(a-h,o-z)

      include 'vmc.h'
      include 'pcm_3dgrid.h'

      common /atom/ znuc(MCTYPE),cent(3,MCENT),pecent
     &,iwctype(MCENT),nctype,ncent
      common /contrl/ nstep,nblk,nblkeq,nconf,nconf_new,isite,idump,irstar

CACTIVATE
c     if (irstar.ne.1) then
c Read the variables from the input

       call p2gtid('pcm:nx_pcm',ipcm_nstep3d(1),IUNDEFINED,1)
       call p2gtid('pcm:ny_pcm',ipcm_nstep3d(2),IUNDEFINED,1)
       call p2gtid('pcm:nz_pcm',ipcm_nstep3d(3),IUNDEFINED,1)

       call p2gtfd('pcm:dx_pcm',pcm_step3d(1),UNDEFINED,1)
       call p2gtfd('pcm:dy_pcm',pcm_step3d(2),UNDEFINED,1)
       call p2gtfd('pcm:dz_pcm',pcm_step3d(3),UNDEFINED,1)
 
       call p2gtfd('pcm:x0_pcm',pcm_origin(1),UNDEFINED,1)
       call p2gtfd('pcm:y0_pcm',pcm_origin(2),UNDEFINED,1)
       call p2gtfd('pcm:z0_pcm',pcm_origin(3),UNDEFINED,1)
 
       call p2gtfd('pcm:xn_pcm',pcm_endpt(1),UNDEFINED,1)
       call p2gtfd('pcm:yn_pcm',pcm_endpt(2),UNDEFINED,1)
       call p2gtfd('pcm:zn_pcm',pcm_endpt(3),UNDEFINED,1)

       call p2gtfd('pcm:shift',PCM_SHIFT,4.d0,1)

c Test if the input is consistent
       input_ok=1
       do i=1,3
        if(ipcm_nstep3d(i).eq.IUNDEFINED.and.pcm_step3d(i).eq.UNDEFINED) then
          write(6,*) 'ipcm_nstep3d(',i,') and pcm_step3d(',i,') are undefined'
          input_ok=0
        endif
        if(pcm_origin(i).gt.pcm_endpt(i)) then
         write (6,*) 'The pcm_origin coordinates have to be smaller than the end points'
         input_ok=0
        endif
       enddo

       if(input_ok.ne.1) call fatal_error('PCM_3DGRID: 3D Grid input inconsistent')

c Origin and end of the grid. If not in input, compute it from the atomic coordinates.
       do i=1,3
        if(pcm_origin(i).eq.UNDEFINED) then
          pcm_origin(i)=cent(i,1)
          do j=2,ncent
            if(cent(i,j).lt.pcm_origin(i)) pcm_origin(i) = cent(i,j)
          enddo
          pcm_origin(i)=pcm_origin(i)-PCM_SHIFT
        endif
 
        if(pcm_endpt(i).eq.UNDEFINED) then
          pcm_endpt(i)=cent(i,1)
          do j=2,ncent
            if(cent(i,j).gt.pcm_endpt(i)) pcm_endpt(i)=cent(i,j)
          enddo
          pcm_endpt(i)=pcm_endpt(i)+PCM_SHIFT
        endif
       enddo
 
c If the step is undefined, use the value in ipcm_nstep3d to compute it.
c Else, compute the value of ipcm_nstep3d
       do i=1, 3
        if(pcm_step3d(i).eq.UNDEFINED) then
         pcm_step3d(i)=(pcm_endpt(i)-pcm_origin(i))/(ipcm_nstep3d(i)-1)
        else
         ipcm_nstep3d(i)=int((pcm_endpt(i)-pcm_origin(i))/pcm_step3d(i))+1
        endif

        if (ipcm_nstep3d(i).gt.MGRID_PCM) then
         write(6,*) 'Warning: using ipcm_nstep3d(',i,') = ',MGRID_PCM
         ipcm_nstep3d(i)=MGRID_PCM
         pcm_endpt(i)=MGRID_PCM*pcm_step3d(i)+pcm_origin(i)
        endif
       enddo

c Prepare the integer->cartesian array
       do i=1,3
        do j=1,ipcm_nstep3d(i)
         pcm_cart_from_int(j,i)=pcm_origin(i)+(j-1)*pcm_step3d(i)
        enddo
       enddo     

c Update the end point 
       do i=1,3
         pcm_endpt(i)=pcm_cart_from_int(ipcm_nstep3d(i),i)
       enddo

CACTIVATE
c     endif

c     Print the parameters to the output file

      write(6,*) 
      write(6,'(''pcm 3D grid parameters'')')
      write(6,*) 
      write(6,'(''pcm origin and end points'')')
      write(6,'(3(F10.6, 3X))') ( pcm_origin(i), i=1,3 )
      write(6,'(3(F10.6, 3X))') ( pcm_endpt (i), i=1,3 )
      write(6,'(''pcm number of steps'')')
      write(6,'(3(I5, 3X))') ( ipcm_nstep3d(i), i=1,3 )
      write(6,'(''pcm step sizes'')')
      write(6,'(3(F10.6, 3X))') ( pcm_step3d (i), i=1,3 )
      write(6,*) 
      
      end ! subroutine pcm_setup_grid
c----------------------------------------------------------------------
      function ipcm_int_from_cart(value,iaxis)
      implicit real*8(a-h,o-z)

      include 'pcm_3dgrid.h'
      
      if (value.lt.pcm_origin(iaxis).or.value.ge.pcm_endpt(iaxis)) then
        ipcm_int_from_cart = IUNDEFINED
       else 
        ipcm_int_from_cart=int((value-pcm_origin(iaxis))/pcm_step3d(iaxis)+1.0)
      endif

      end ! subroutine ipcm_int_from_cart
c----------------------------------------------------------------------
c PCM on a 3d grid with spline fit
      subroutine pcm_setup_3dspl
      implicit real*8(a-h,o-z)

      include 'force.h'
      include 'vmc.h'
      include 'pcm_3dgrid.h'

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
      real*8  bc(MGRID_PCM,MGRID_PCM,3:8), wk(80*MGRID_PCM3)
      common /pcm_num_spl2/ bc, wk
      common /contrl/ nstep,nblk,nblkeq,nconf,nconf_new,isite,idump,irstar

      dimension r(3)

      iok=1
c We have no info on the derivatives, so use "not a knot" in the creation of the fit
      ibcxmin=0
      ibcxmax=0
      ibcymin=0
      ibcymax=0
      ibczmin=0
      ibczmax=0

c Evaluate the energy needed for the calculation
      memory=dfloat(8)
      memory=memory*dfloat(ipcm_nstep3d(1)*ipcm_nstep3d(2)*ipcm_nstep3d(3))
      memory=memory*8.d-6
      write(45,*) 'Allocated memory for the 3D spline fit of PCM:', memory, 'Mb'

      if(irstar.ne.1) then
       write(45,*) 'Computation of the grid points...'
       do 10 ix=1,ipcm_nstep3d(1)
          r(1)=pcm_cart_from_int(ix,1)

          do 10 iy=1,ipcm_nstep3d(2)
            r(2)=pcm_cart_from_int(iy,2)

            do 10 iz=1,ipcm_nstep3d(3)
              r(3) =pcm_cart_from_int(iz,3)
         
c Calculate the value of the pcm potential on the position [r(1),r(2),r(3)]
              call pcm_extpot_ene_elec(r,pepol_s,pepol_v)
  10          pcm_num_spl(1,ix,iy,iz)=pepol_s+pepol_v

       nwk=80*ipcm_nstep3d(1)*ipcm_nstep3d(2)*ipcm_nstep3d(3)
       ier=0
       call r8mktricubw(pcm_cart_from_int(1,1),ipcm_nstep3d(1),
     &                  pcm_cart_from_int(1,2),ipcm_nstep3d(2),
     &                  pcm_cart_from_int(1,3),ipcm_nstep3d(3),
     &                  pcm_num_spl,MGRID_PCM,MGRID_PCM,
     &                  ibcxmin,bc(1,1,2+3),
     &                  ibcxmax,bc(1,1,2+3+3),MGRID_PCM,
     &                  ibcymin,bc(1,1,1+3),
     &                  ibcymax,bc(1,1,1+3+3),MGRID_PCM,
     &                  ibczmin,bc(1,1,1+2),
     &                  ibczmax,bc(1,1,1+2+3),MGRID_PCM,
     &                  wk,nwk,ilinx,iliny,ilinz,ier)
        if(ier.eq.1) call fatal_error ('Error in r8mktricubw')
      endif ! (irstar.ne.0)

c DEBUG
c      do 30 ix=1,ipcm_nstep3d(1)
c         r(1)=pcm_cart_from_int(ix,1)
c         do 30 iy=1,ipcm_nstep3d(2)
c           r(2)=pcm_cart_from_int(iy,2)
c           do 30 iz=1,ipcm_nstep3d(3)
c             r(3) =pcm_cart_from_int(iz,3)
c         
c             call pcm_extpot_ene_elec(r,pepol_s,pepol_v)
c             call spline_pcm(r,value,ier)
c  30         write(6,*) 'DEBUG',pepol_s+pepol_v, value
c      stop

      end ! subroutine pcm_setup_3dspl

c----------------------------------------------------------------------
      subroutine spline_pcm(r,f,ier)
      implicit real*8(a-h,o-z)
      include 'vmc.h'
      include 'pcm_3dgrid.h'

      real*8 inside, inout  
      common /insout/ inside,inout  

c     Input:
      real*8    r(3)    ! Cartesian coordinates
c     Output:
      integer   ier     ! error status
      real*8    f       ! Value

c     Work:
      integer   ict(10)   ! Control of the spline subroutine
      integer   ix(3)     ! Integer coordinates
      real*8    fval(10)  ! values returned by the spline routine
      real*8    rscaled(3)! normalized displacement
      real*8    inv_step3d(3) ! Inverse of step sizes


      inout=inout+1.d0
      if (ier.eq.1) then
        return
      endif

      ict(1)=1
      do i=1,3
       ix(i) = ipcm_int_from_cart(r(i),i)
      enddo

      if ((ix(1).eq.IUNDEFINED ).or.(ix(2).eq.IUNDEFINED).or.(ix(3).eq.IUNDEFINED)) then
        ier=1
      else
        inside = inside+1.d0
        do i=2,10
         ict(i)=0
        enddo

        do i=1,3
         inv_step3d(i) = 1.d0/pcm_step3d(i)
        enddo
        do i=1,3
         rscaled(i)=(r(i)-pcm_cart_from_int(ix(i),i))*inv_step3d(i)
        enddo

       call r8fvtricub(ict,1,1,fval,
     &                  ix(1),ix(2),ix(3),
     &                  rscaled(1),rscaled(2),rscaled(3),
     &                  pcm_step3d(1),inv_step3d(1),
     &                  pcm_step3d(2),inv_step3d(2),
     &                  pcm_step3d(3),inv_step3d(3),
     &                  pcm_num_spl,
     &                  MGRID_PCM,MGRID_PCM,ipcm_nstep3d(3))
      
        f=fval(1)
      endif 
      
      end ! subroutine spline_pcm

c-----------------------------------------------------------------------
      subroutine pcm_3dgrid_dump(iu)
      implicit real*8(a-h,o-z)
      include 'vmc.h'
      include 'pcm_3dgrid.h'

      if (ipcm.eq.0.or.ipcm_grid.eq.0) return

      write (iu) (pcm_origin(i), i=1,3)
      write (iu) (pcm_endpt(i), i=1,3)
      write (iu) (ipcm_nstep3d(i), i=1,3)
      write (iu) (pcm_step3d(i), i=1,3)
      write (iu) ((pcm_cart_from_int(i,j), i=1,ipcm_nstep3d(j)),j=1,3)

      call splpcm_dump(iu)

      end
c-----------------------------------------------------------------------
      subroutine pcm_3dgrid_rstrt(iu)
      implicit real*8(a-h,o-z)
      include 'vmc.h'
      include 'pcm_3dgrid.h'

      if (ipcm.eq.0.or.ipcm_grid.eq.0) return

      read (iu) (pcm_origin(i), i=1,3)
      read (iu) (pcm_endpt(i), i=1,3)
      read (iu) (ipcm_nstep3d(i), i=1,3)
      read (iu) (pcm_step3d(i), i=1,3)
      read (iu) ((pcm_cart_from_int(i,j), i=1,ipcm_nstep3d(j)),j=1,3)
      
      call splpcm_rstrt(iu)
      end
c-----------------------------------------------------------------------
      subroutine splpcm_dump(iu)
      implicit real*8(a-h,o-z)
      include 'vmc.h'
      include 'pcm_3dgrid.h'
 
      do i=1,8
        write(iu)(((pcm_num_spl(i,j,k,l),j=1,ipcm_nstep3d(1)),k=1,ipcm_nstep3d(2)), l=1,ipcm_nstep3d(3))
      enddo

      end
c-----------------------------------------------------------------------
      subroutine splpcm_rstrt(iu)
      implicit real*8(a-h,o-z)
      include 'vmc.h'
      include 'pcm_3dgrid.h'

      do i=1,8
        read(iu)(((pcm_num_spl(i,j,k,l),j=1,ipcm_nstep3d(1)),k=1,ipcm_nstep3d(2)),l=1,ipcm_nstep3d(3))
      enddo
      end
c-----------------------------------------------------------------------
