      module misc_grdnts
c -----------------------------------------------------------------------
c   Various subroutines used when energy gradients for atoms are
c   calculated.
c
c   Written by Omar Valsson
c -----------------------------------------------------------------------
      use m_zmat_tools, only: zmat2cart_rc
      contains


c   Subroutine which at the start up prints out information about the
c   energy gradients (cartesian).
      subroutine inpwrt_grdnts_cart()
      use grdntsmv, only: igrdaidx, igrdcidx

      use grdntspar, only: delgrdxyz, ngradnts
      use contrl_file,    only: ounit
      implicit none

      integer :: ig





      write(ounit,*)
      write(ounit,'(''Correlated sampling used to calculate energy gradients (Cartesian)'')')
      write(ounit,'(''- Number of gradients: '',i0)') ngradnts
      write(ounit,'(''- Gradients calculated for:'')')
      do ig=1,ngradnts
        if(igrdcidx(ig).eq.1) then
          write(ounit,'(''   x coordinate of atom'',i0)') igrdaidx(ig)
        elseif(igrdcidx(ig).eq.2) then
          write(ounit,'(''   y coordinate of atom'',i0)') igrdaidx(ig)
        else
          write(ounit,'(''   z coordinate of atom'',i0)') igrdaidx(ig)
        endif
      enddo
      write(ounit,'(''- Displacement of x,y,z coordinates - delgrdxyz = '',f7.5)') delgrdxyz
      write(ounit,*)

      return
      end
c -----------------------------------------------------------------------

c   Subroutine which at the start up prints out information about the
c   energy gradients (z matrix/internal).
      subroutine inpwrt_grdnts_zmat()

      use grdntsmv, only: igrdaidx, igrdcidx
      use grdntspar, only: delgrdba, delgrdbl, delgrdda, ngradnts
      use zmatrix, only: izcmat
      use contrl_file,    only: ounit
      implicit none

      integer :: ig, na, nb, nc, nd




      write(ounit,*)
      write(ounit,'(''Correlated sampling used to calculate energy gradients (Z matrix)'')')
      write(ounit,'(''- Number of gradients: '',i0)') ngradnts
      write(ounit,'(''- Gradients calculated for:'')')
      do ig=1,ngradnts
        na=igrdaidx(ig)
        nb=izcmat(1,na)
        nc=izcmat(2,na)
        nd=izcmat(3,na)
        if(igrdcidx(ig).eq.1) then
          write(ounit,'(''   bond length between atoms '',i0,'' and '',i0)') na,nb
        elseif(igrdcidx(ig).eq.2) then
          write(ounit,'(''   bond angle between atoms '',i0,'','',i0,'' and '',i0)') na,nb,nc
        else
          write(ounit,'(''   dihedral angle between atoms '',i0,'','',i0,'','',i0,'' and '',i0)') na,nb,nc,nd
        endif
      enddo
      write(ounit,'(''- Displacement of bond lengths      -   delgrdbl = '',f7.5)') delgrdbl
      write(ounit,'(''- Displacement of bond angles       -   delgrdba = '',f7.5)') delgrdba
      write(ounit,'(''- Displacement of dihedral angles   -   delgrdda = '',f7.5)') delgrdda
      write(ounit,*)

      return
      end
c -----------------------------------------------------------------------

c   Subroutine which calculates and printouts energy gradients
c   for cartesian coordinates of atoms from energy differences
c   calculated using correlated smapling.
      subroutine finwrt_grdnts_cart(forces_ave,forces_err)
      use force_mod, only: MFORCE
      use system, only: iwctype, ncent
      use forcepar, only: nforce
      use grdntsmv, only: igrdaidx, igrdcidx, igrdmv

      use grdntspar, only: delgrdxyz, ngradnts
      use contrl_file,    only: ounit
      use precision_kinds, only: dp
      implicit none

      integer :: ic, ig, k, if
      real(dp) :: advance
      real(dp), dimension(MFORCE) :: forces_ave
      real(dp), dimension(MFORCE) :: forces_err
      real(dp), dimension(MFORCE) :: grdnts_ave
      real(dp), dimension(MFORCE) :: grdnts_err





      do if=1,nforce-2,2
        grdnts_ave((if+1)/2)=((-forces_ave(if)+forces_ave(if+1))/delgrdxyz)*0.5d0
        grdnts_err((if+1)/2)=(dsqrt(forces_err(if)**2+forces_err(if+1)**2)/delgrdxyz)*0.5d0
      enddo

      write(ounit,*)
      write(ounit,'(''Energy gradients (dE/d{x,y,z}):'')')
      do ig=1,ngradnts
        if(igrdcidx(ig).eq.1) then
          write(ounit,'(''  x coordinate of atom'',i3,'': '',f10.7,'' +- '',f10.7)')
     &      igrdaidx(ig),grdnts_ave(ig),grdnts_err(ig)
        elseif(igrdcidx(ig).eq.2) then
          write(ounit,'(''  y coordinate of atom'',i3,'': '',f10.7,'' +- '',f10.7)')
     &      igrdaidx(ig),grdnts_ave(ig),grdnts_err(ig)
        else
          write(ounit,'(''  z coordinate of atom'',i3,'': '',f10.7,'' +- '',f10.7)')
     &      igrdaidx(ig),grdnts_ave(ig),grdnts_err(ig)
        endif
      enddo

      write(ounit,*)

      ig=1
      write(ounit,'(''Energy gradients (Cartesian) - atom table:'')')
      write(ounit,'(1x,''Atom'',7x,''dE/dx'',21x,''dE/dy'',21x,''dE/dz'',15x,''iwctype'')')
      do ic=1,ncent
        write(ounit,'(1x,i3,3x)',advance='no') ic
        do k=1,3
          if(igrdmv(k,ic).eq.1) then
            write(ounit,'(f10.7,'' +- '',f9.7,3x)',advance='no') grdnts_ave(ig),grdnts_err(ig)
            ig=ig+1
          else
            write(ounit,'(a10,'' +- '',a9,3x)',advance='no') 'x.xxxxxxx','x.xxxxxxx'
          endif
        enddo
        write(ounit,'(i3)') iwctype(ic)
      enddo
      write(ounit,'('' ---------------------------------------------------------------------------------------'')')
      write(ounit,*)

      return
      end
c -----------------------------------------------------------------------

c   Subroutine which calculates and printouts energy gradients
c   for Z matrix coordinates of atoms from energy differences
c   calculated using correlated smapling.
      subroutine finwrt_grdnts_zmat(forces_ave,forces_err)
      use force_mod, only: MFORCE
      use system, only: iwctype, ncent
      use forcepar, only: nforce
      use grdntsmv, only: igrdaidx, igrdcidx, igrdmv

      use grdntspar, only: delgrdba, delgrdbl, delgrdda, ngradnts
      use zmatrix, only: izcmat
      use contrl_file,    only: ounit
      use precision_kinds, only: dp
      implicit none

      integer :: ic, ifr, ig, k, na
      integer :: nb, nc, nd
      real(dp) :: advance, delgrd_tmp
      real(dp), dimension(MFORCE) :: forces_ave
      real(dp), dimension(MFORCE) :: forces_err
      real(dp), dimension(MFORCE) :: grdnts_ave
      real(dp), dimension(MFORCE) :: grdnts_err






      do ifr=1,nforce-2,2
        if(igrdcidx((ifr+1)/2).eq.1) then
          delgrd_tmp=delgrdbl
        elseif(igrdcidx((ifr+1)/2).eq.2) then
          delgrd_tmp=delgrdba
        else
          delgrd_tmp=delgrdda
        endif
        grdnts_ave((ifr+1)/2)=((-forces_ave(ifr)+forces_ave(ifr+1))/delgrd_tmp)*0.5d0
        grdnts_err((ifr+1)/2)=(dsqrt(forces_err(ifr)**2+forces_err(ifr+1)**2)/delgrd_tmp)*0.5d0
      enddo

      write(ounit,*)
      write(ounit,'(''Energy gradients (dE/d{bl,ba,da}):'')')
      do ig=1,ngradnts
        na=igrdaidx(ig)
        nb=izcmat(1,na)
        nc=izcmat(2,na)
        nd=izcmat(3,na)
        if(igrdcidx(ig).eq.1) then
          write(ounit,'('' bond length between atoms '',i0,'' and '',i0,'': '',t50,f10.7,'' +- '',f10.7)')
     &            na,nb,grdnts_ave(ig),grdnts_err(ig)
        elseif(igrdcidx(ig).eq.2) then
          write(ounit,'('' bond angle between atoms '',i0,'','',i0,'' and '',i0,'': '',t50,f10.7,'' +- '',f10.7)')
     &            na,nb,nc,grdnts_ave(ig),grdnts_err(ig)
        else
          write(ounit,'('' dihedral angle between atoms '',i0,'','',i0,'','',i0,'' and '',i0,'': '',t50,f10.7,'' +- '',f10.7)')
     &            na,nb,nc,nd,grdnts_ave(ig),grdnts_err(ig)
        endif
      enddo

      write(ounit,*)

      ig=1
      write(ounit,'(''# Energy gradients (Z matrix) - atom table:'')')
      write(ounit,'(''#'',1x,''Atom'',3x,''izcmat'',11x,''dE/dbl'',20x,''dE/dba'',20x,''dE/dda'',15x,''iwctype'')')
      ic=1
      write(ounit,'(1x,i3,3x,3(2x,2x),2x)',advance='no') ic
      write(ounit,'(3(10x,4x,9x,3x))',advance='no')
      write(ounit,'(i3)') iwctype(ic)
      ic=2
      write(ounit,'(1x,i3,3x,1(i2,2x),2(2x,2x),2x)',advance='no') ic,izcmat(1,ic)
      k=1
      if(igrdmv(k,ic).eq.1) then
        write(ounit,'(f10.7,'' +- '',f9.7,3x)',advance='no') grdnts_ave(ig),grdnts_err(ig)
        ig=ig+1
      else
        write(ounit,'(a10,'' +- '',a9,3x)',advance='no') 'x.xxxxxxx','x.xxxxxxx'
      endif
      write(ounit,'(2(10x,4x,9x,3x))',advance='no')
      write(ounit,'(i3)') iwctype(ic)
      ic=3
      write(ounit,'(1x,i3,3x,2(i2,2x),1(2x,2x),2x)',advance='no') ic,(izcmat(k,ic),k=1,2)
      do k=1,2
        if(igrdmv(k,ic).eq.1) then
          write(ounit,'(f10.7,'' +- '',f9.7,3x)',advance='no') grdnts_ave(ig),grdnts_err(ig)
          ig=ig+1
        else
          write(ounit,'(a10,'' +- '',a9,3x)',advance='no') 'x.xxxxxxx','x.xxxxxxx'
        endif
      enddo
      write(ounit,'(1(10x,4x,9x,3x))',advance='no')
      write(ounit,'(i3)') iwctype(ic)
      do ic=4,ncent
        write(ounit,'(1x,i3,3x,3(i2,2x),2x)',advance='no') ic,(izcmat(k,ic),k=1,3)
        do k=1,3
          if(igrdmv(k,ic).eq.1) then
            write(ounit,'(f10.7,'' +- '',f9.7,3x)',advance='no') grdnts_ave(ig),grdnts_err(ig)
            ig=ig+1
          else
            write(ounit,'(a10,'' +- '',a9,3x)',advance='no') 'x.xxxxxxx','x.xxxxxxx'
          endif
        enddo
        write(ounit,'(i3)') iwctype(ic)
      enddo
      write(ounit,'(''# ---------------------------------------------------------------------------------------'')')
      write(ounit,*)


      open(2,file='GradientsZmat.data',status='unknown')

      ig=1
      write(2,'(''# Energy gradients (Z matrix) - atom table:'')')
      write(2,'(''#'',1x,''Atom'',3x,''izcmat'',11x,''dE/dbl'',20x,''dE/dba'',20x,''dE/dda'',15x,''iwctype'')')
      ic=1
      write(2,'(1x,i3,3x,3(2x,2x),2x)',advance='no') ic
      write(2,'(3(10x,4x,9x,3x))',advance='no')
      write(2,'(i3)') iwctype(ic)
      ic=2
      write(2,'(1x,i3,3x,1(i2,2x),2(2x,2x),2x)',advance='no') ic,izcmat(1,ic)
      k=1
      if(igrdmv(k,ic).eq.1) then
        write(2,'(f10.7,'' +- '',f9.7,3x)',advance='no') grdnts_ave(ig),grdnts_err(ig)
        ig=ig+1
      else
        write(2,'(a10,'' +- '',a9,3x)',advance='no') '0.0000000','0.0000000'
      endif
      write(2,'(2(10x,4x,9x,3x))',advance='no')
      write(2,'(i3)') iwctype(ic)
      ic=3
      write(2,'(1x,i3,3x,2(i2,2x),1(2x,2x),2x)',advance='no') ic,(izcmat(k,ic),k=1,2)
      do k=1,2
        if(igrdmv(k,ic).eq.1) then
          write(2,'(f10.7,'' +- '',f9.7,3x)',advance='no') grdnts_ave(ig),grdnts_err(ig)
          ig=ig+1
        else
          write(2,'(a10,'' +- '',a9,3x)',advance='no') '0.0000000','0.0000000'
        endif
      enddo
      write(2,'(1(10x,4x,9x,3x))',advance='no')
      write(2,'(i3)') iwctype(ic)
      do ic=4,ncent
        write(2,'(1x,i3,3x,3(i2,2x),2x)',advance='no') ic,(izcmat(k,ic),k=1,3)
        do k=1,3
          if(igrdmv(k,ic).eq.1) then
            write(2,'(f10.7,'' +- '',f9.7,3x)',advance='no') grdnts_ave(ig),grdnts_err(ig)
            ig=ig+1
          else
            write(2,'(a10,'' +- '',a9,3x)',advance='no') '0.0000000','0.0000000'
          endif
        enddo
        write(2,'(i3)') iwctype(ic)
      enddo
      write(2,'(''# ---------------------------------------------------------------------------------------'')')
      write(2,'(''# '')')

      close(2)


      return
      end
c -----------------------------------------------------------------------

c   Subroutine which calculates the displacement for energy gradients
c   using Z matrix (internal) coordinates
      subroutine grdzmat_displ(k_in,ic_in,ia_in,delfactor)
      use system, only: ncent, ncent_tot
      use forcestr, only: delc

      use grdntspar, only: delgrdba, delgrdbl, delgrdda
      use zmatrix, only: czcart, czint, czcart_ref, izcmat
      use contrl_file,    only: ounit
      use precision_kinds, only: dp
      implicit none

      integer :: ia_in, ic, ic_in, k, k_in
      real(dp) :: delfactor, delgrd_tmp
      real(dp), dimension(3,ncent_tot) :: czint_t1
      real(dp), dimension(3,ncent_tot) :: czcart_t1






      if(k_in.eq.1) then
        delgrd_tmp=delfactor*delgrdbl
      elseif(k_in.eq.2) then
        delgrd_tmp=delfactor*delgrdba
      else
        delgrd_tmp=delfactor*delgrdda
      endif

      do ic=1,ncent
        do k=1,3
          czint_t1(k,ic)=czint(k,ic)
          czcart_t1(k,ic)=0.0d0
        enddo
      enddo


      czint_t1(k_in,ic_in)=czint_t1(k_in,ic_in)+delgrd_tmp
      call zmat2cart_rc(ncent_tot,izcmat,czint_t1,czcart_t1,czcart_ref)

      do ic=1,ncent
        do k=1,3
          delc(k,ic,ia_in)=czcart_t1(k,ic)-czcart(k,ic)
        enddo
      enddo

      return
      end
c -----------------------------------------------------------------------


c   Subroutine which prints out at the start of a run
c   information regarding the Z matrix.
      subroutine inpwrt_zmatrix()
      use system, only: iwctype, ncent
      use zmatrix, only: czcart, czint, izcmat
      use contrl_file,    only: ounit
      implicit none

      integer :: ic, k


      write(ounit,'(''---------- Z matrix information ----------'')')
      write(ounit,'(''Geometry in internal coordinates:'')')
      write(ounit,'(1x,''Atom'',t9,''izcmat'',t24,''bond length'',t37,''bond angle'',t51,''dihed. angle '',t65,''iwctype'')')
      ic=1
      write(ounit,'(1x,i3,59x,i3)') ic,iwctype(ic)
      ic=2
      write(ounit,'(1x,i3,3x,1(i2,2x),10x,1(f11.7,3x),28x,i3)') ic,izcmat(1,ic),czint(1,ic),iwctype(ic)
      ic=3
      write(ounit,'(1x,i3,3x,2(i2,2x),6x,1(f11.7,3x),1(f11.6,3x),14x,i3)') ic,(izcmat(k,ic),k=1,2),(czint(k,ic),k=1,2),iwctype(ic)
      do ic=4,ncent
        write(ounit,'(1x,i3,3x,3(i2,2x),2x,1(f11.7,3x),2(f11.6,3x),i3)') ic,(izcmat(k,ic),k=1,3),(czint(k,ic),k=1,3),iwctype(ic)
      enddo
      write(ounit,'('' ----------------------------------------------------------------------'')')
      write(ounit,'(''Internal coordiantes converted back into cartesian coordinates:'')')
      write(ounit,'(1x,''Atom'',t12,''x'',t26,''y'',t40,''z'',t51,''iwctype'')')
      do ic=1,ncent
        write(ounit,'(1x,i3,3x,3(f11.7,3x),i3)') ic,(czcart(k,ic),k=1,3),iwctype(ic)
      enddo
      write(ounit,'('' --------------------------------------------------------'')')
      write(ounit,*)

      return
      end
c -----------------------------------------------------------------------



c   Subroutine which calculates and print outs the diagonal
c   part of the Hessian for Z matrix coordinates of atoms
c   from energy differences  calculated using correlated smapling.
      subroutine finwrt_diaghess_zmat(forces_ave,forces_err)
      use force_mod, only: MFORCE
      use system, only: iwctype, ncent
      use forcepar, only: nforce
      use grdntsmv, only: igrdaidx, igrdcidx, igrdmv

      use grdntspar, only: delgrdba, delgrdbl, delgrdda, ngradnts
      use zmatrix, only: izcmat
      use contrl_file,    only: ounit
      use precision_kinds, only: dp
      implicit none

      integer :: ic, ifr, ig, k, na
      integer :: nb, nc, nd
      real(dp) :: advance, delgrd_tmp
      real(dp), dimension(MFORCE) :: forces_ave
      real(dp), dimension(MFORCE) :: forces_err
      real(dp), dimension(MFORCE) :: diaghess_ave
      real(dp), dimension(MFORCE) :: diaghess_err









      do ifr=1,nforce-2,2
        if(igrdcidx((ifr+1)/2).eq.1) then
          delgrd_tmp=delgrdbl
        elseif(igrdcidx((ifr+1)/2).eq.2) then
          delgrd_tmp=delgrdba
        else
          delgrd_tmp=delgrdda
        endif
        diaghess_ave((ifr+1)/2)=((-forces_ave(ifr)-forces_ave(ifr+1))/delgrd_tmp**2)
        diaghess_err((ifr+1)/2)=(dsqrt(forces_err(ifr)**2+forces_err(ifr+1)**2)/delgrd_tmp**2)
      enddo

      write(ounit,*)
      write(ounit,'(''Diagonal of the Hessian (d^2E/d{bl,ba,da}^2):'')')
      do ig=1,ngradnts
        na=igrdaidx(ig)
        nb=izcmat(1,na)
        nc=izcmat(2,na)
        nd=izcmat(3,na)
        if(igrdcidx(ig).eq.1) then
          write(ounit,'('' bond length between atoms '',i0,'' and '',i0,'': '',t50,f10.7,'' +- '',f10.7)')
     &            na,nb,diaghess_ave(ig),diaghess_err(ig)
        elseif(igrdcidx(ig).eq.2) then
          write(ounit,'('' bond angle between atoms '',i0,'','',i0,'' and '',i0,'': '',t50,f10.7,'' +- '',f10.7)')
     &            na,nb,nc,diaghess_ave(ig),diaghess_err(ig)
        else
          write(ounit,'('' dihedral angle between atoms '',i0,'','',i0,'','',i0,'' and '',i0,'': '',t50,f10.7,'' +- '',f10.7)')
     &            na,nb,nc,nd,diaghess_ave(ig),diaghess_err(ig)
        endif
      enddo

      write(ounit,*)

      ig=1
      write(ounit,'(''# Diagonal of the Hessian (Z matrix) - atom table:'')')
      write(ounit,'(''#'',1x,''Atom'',3x,''izcmat'',11x,''d^2E/dbl^2'',20x,''d^2E/dba^2'',20x,''d^2E/dda^2'',15x,''iwctype'')')
      ic=1
      write(ounit,'(1x,i3,3x,3(2x,2x),2x)',advance='no') ic
      write(ounit,'(3(10x,4x,9x,3x))',advance='no')
      write(ounit,'(i3)') iwctype(ic)
      ic=2
      write(ounit,'(1x,i3,3x,1(i2,2x),2(2x,2x),2x)',advance='no') ic,izcmat(1,ic)
      k=1
      if(igrdmv(k,ic).eq.1) then
        write(ounit,'(f10.7,'' +- '',f9.7,3x)',advance='no') diaghess_ave(ig),diaghess_err(ig)
        ig=ig+1
      else
        write(ounit,'(a10,'' +- '',a9,3x)',advance='no') 'x.xxxxxxx','x.xxxxxxx'
      endif
      write(ounit,'(2(10x,4x,9x,3x))',advance='no')
      write(ounit,'(i3)') iwctype(ic)
      ic=3
      write(ounit,'(1x,i3,3x,2(i2,2x),1(2x,2x),2x)',advance='no') ic,(izcmat(k,ic),k=1,2)
      do k=1,2
        if(igrdmv(k,ic).eq.1) then
          write(ounit,'(f10.7,'' +- '',f9.7,3x)',advance='no') diaghess_ave(ig),diaghess_err(ig)
          ig=ig+1
        else
          write(ounit,'(a10,'' +- '',a9,3x)',advance='no') 'x.xxxxxxx','x.xxxxxxx'
        endif
      enddo
      write(ounit,'(1(10x,4x,9x,3x))',advance='no')
      write(ounit,'(i3)') iwctype(ic)
      do ic=4,ncent
        write(ounit,'(1x,i3,3x,3(i2,2x),2x)',advance='no') ic,(izcmat(k,ic),k=1,3)
        do k=1,3
          if(igrdmv(k,ic).eq.1) then
            write(ounit,'(f10.7,'' +- '',f9.7,3x)',advance='no') diaghess_ave(ig),diaghess_err(ig)
            ig=ig+1
          else
            write(ounit,'(a10,'' +- '',a9,3x)',advance='no') 'x.xxxxxxx','x.xxxxxxx'
          endif
        enddo
        write(ounit,'(i3)') iwctype(ic)
      enddo
      write(ounit,'(''# ---------------------------------------------------------------------------------------'')')
      write(ounit,*)

      open(2,file='HessianZmatDiag.data',status='unknown')

      ig=1
      write(2,'(''# Diagonal of the Hessian (Z matrix) - atom table:'')')
      write(2,'(''#'',1x,''Atom'',3x,''izcmat'',11x,''d^2E/dbl^2'',20x,''d^2E/dba^2'',20x,''d^2E/dda^2'',15x,''iwctype'')')
      ic=1
      write(2,'(1x,i3,3x,3(2x,2x),2x)',advance='no') ic
      write(2,'(3(10x,4x,9x,3x))',advance='no')
      write(2,'(i3)') iwctype(ic)
      ic=2
      write(2,'(1x,i3,3x,1(i2,2x),2(2x,2x),2x)',advance='no') ic,izcmat(1,ic)
      k=1
      if(igrdmv(k,ic).eq.1) then
        write(2,'(f10.7,'' +- '',f9.7,3x)',advance='no') diaghess_ave(ig),diaghess_err(ig)
        ig=ig+1
      else
        write(2,'(a10,'' +- '',a9,3x)',advance='no') '0.0000000','0.0000000'
      endif
      write(2,'(2(10x,4x,9x,3x))',advance='no')
      write(2,'(i3)') iwctype(ic)
      ic=3
      write(2,'(1x,i3,3x,2(i2,2x),1(2x,2x),2x)',advance='no') ic,(izcmat(k,ic),k=1,2)
      do k=1,2
        if(igrdmv(k,ic).eq.1) then
          write(2,'(f10.7,'' +- '',f9.7,3x)',advance='no') diaghess_ave(ig),diaghess_err(ig)
          ig=ig+1
        else
          write(2,'(a10,'' +- '',a9,3x)',advance='no') '0.0000000','0.0000000'
        endif
      enddo
      write(2,'(1(10x,4x,9x,3x))',advance='no')
      write(2,'(i3)') iwctype(ic)
      do ic=4,ncent
        write(2,'(1x,i3,3x,3(i2,2x),2x)',advance='no') ic,(izcmat(k,ic),k=1,3)
        do k=1,3
          if(igrdmv(k,ic).eq.1) then
            write(2,'(f10.7,'' +- '',f9.7,3x)',advance='no') diaghess_ave(ig),diaghess_err(ig)
            ig=ig+1
          else
            write(2,'(a10,'' +- '',a9,3x)',advance='no') '0.0000000','0.0000000'
          endif
        enddo
        write(2,'(i3)') iwctype(ic)
      enddo
      write(2,'(''# ---------------------------------------------------------------------------------------'')')
      write(2,'(''# '')')

      close(2)

      return
      end
c -----------------------------------------------------------------------
      subroutine transform_grad_zmat(force_cart)
      use vmc_mod, only: ncent3
      use system, only: cent, ncent, ncent_tot

      use grdntsmv, only: igrdmv
      use zmatrix, only: czint, czcart_ref, izcmat
      use force_analy, only: iuse_zmat
      use contrl_file,    only: ounit
      use precision_kinds, only: dp
      implicit none

      integer :: ic, ii, jc, jj, k
      integer :: l, lc, ncent_ind
      real(dp) :: czint_sav
      real(dp), dimension(3,ncent_tot) :: czcartp
      real(dp), dimension(3,ncent_tot) :: czcartm
      real(dp), dimension(ncent3,ncent3) :: bwilson
      real(dp), dimension(ncent3*ncent3) :: gmat
      real(dp), dimension(3,*) :: force_cart
      real(dp), dimension(ncent3) :: force_int
      real(dp), parameter :: eps = 1.d-5
      real(dp), parameter :: epsi = 0.5d0/eps






      do ic=1,3
        do k=1,3
          czcart_ref(k,ic)=cent(k,ic)
        enddo
      enddo

c     write(ounit,*) 'HELLO'
c     do 20 ic=1,ncent
c 20    write(ounit,'(3i5)') (izcmat(k,ic),k=1,3)

c     call cart2zmat(MCENT,czcart,izcmat,czint)
c     write(ounit,*) 'HELLO'
c     do 20 ic=1,ncent
c 20    write(ounit,'(1p3e12.5)') (czint(k,ic),k=1,3)
c     call zmat2cart(MCENT,izcmat,czint,czcart)

c     write(ounit,*) 'HELLO'
c     do 21 ic=1,ncent
c 21    write(ounit,'(1p3e12.5)') (czcart(k,ic),k=1,3)

c     return

      ii=0
      do ic=2,ncent
        do k=1,3
          if(ic.eq.2.and.k.ne.1) goto 100
          if(ic.eq.3.and.k.eq.3) goto 100

          ii=ii+1
          czint_sav=czint(k,ic)

          czint(k,ic)=czint_sav+eps
          call zmat2cart_rc(ncent,izcmat,czint,czcartp,czcart_ref)

          do lc=1,3
            do l=1,3
              czcart_ref(l,lc)=cent(l,lc)
            enddo
          enddo

          czint(k,ic)=czint_sav-eps
          call zmat2cart_rc(ncent,izcmat,czint,czcartm,czcart_ref)

          jj=0
          do jc=1,ncent
            do l=1,3
              jj=jj+1
              bwilson(ii,jj)=(czcartp(l,jc)-czcartm(l,jc))*epsi
            enddo
          enddo

          czint(k,ic)=czint_sav
 100  continue
        enddo
      enddo

      ncent_ind=3*ncent-6
      if(ncent.eq.2) ncent_ind=1
      ncent3=3*ncent
c     do ic=1,ncent_ind
c       write(ounit,*) 'Bmat',(bwilson(ic,jc),jc=1,ncent3)
c     enddo

      do ic=1,ncent_ind
        force_int(ic)=0
        jj=0
        do jc=1,ncent
          do l=1,3
          jj=jj+1
          force_int(ic)=force_int(ic)+bwilson(ic,jj)*force_cart(l,jc)
          enddo
        enddo
      enddo

      open(81,file='force_analytic_zmat',form='formatted',status='unknown')
      write(81,'(i5,1p3e14.5)') 1,0.d0,0.d0,0.d0
      write(81,'(i5,1p3e14.5)') 2,force_int(1)*igrdmv(1,2),0.d0,0.d0
      if(ncent.gt.2) write(81,'(i5,1p3e14.5)') 3,force_int(2)*igrdmv(1,3),force_int(3)*igrdmv(2,3),0.d0
      ii=3
      do ic=4,ncent_ind,3
        ii=ii+1
        write(81,'(i5,1p3e14.5)') ii,(force_int(ic+k)*igrdmv(k+1,ii),k=0,2)
      enddo
      close(81)

      if(iuse_zmat.eq.1) then
        force_cart(1,1)=0
        force_cart(2,1)=0
        force_cart(3,1)=0
        force_cart(1,2)=force_int(1)*igrdmv(1,2)
        force_cart(2,2)=0
        force_cart(3,2)=0
        if(ncent.gt.2) then
          force_cart(1,3)=force_int(2)*igrdmv(1,3)
          force_cart(2,3)=force_int(3)*igrdmv(2,3)
          force_cart(3,3)=0
        endif

        ii=3
        do ic=4,ncent_ind,3
          ii=ii+1
          do k=0,2
            force_cart(k+1,ii)=force_int(ic+k)*igrdmv(k+1,ii)
          enddo
        enddo
      endif

      return
      end
c -----------------------------------------------------------------------
      end module
