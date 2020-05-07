c -----------------------------------------------------------------------
c   Various subroutines used when energy gradients for atoms are 
c   calculated.
c
c   Written by Omar Valsson
c -----------------------------------------------------------------------


c   Subroutine which at the start up prints out information about the 
c   energy gradients (cartesian).
      subroutine inpwrt_grdnts_cart()
      use grdntsmv, only: igrdaidx, igrdcidx, igrdmv

      use grdntspar, only: delgrdba, delgrdbl, delgrdda, delgrdxyz, igrdtype, ngradnts

      implicit real*8(a-h,o-z)


      include 'force.h'
      include 'vmc.h'


      write(6,*)
      write(6,'(''Correlated sampling used to calculate energy gradients (Cartesian)'')')
      write(6,'(''- Number of gradients: '',i0)') ngradnts
      write(6,'(''- Gradients calculated for:'')')
      do 20 ig=1,ngradnts
        if(igrdcidx(ig).eq.1) then
          write(6,'(''   x coordinate of atom'',i0)') igrdaidx(ig)
        elseif(igrdcidx(ig).eq.2) then
          write(6,'(''   y coordinate of atom'',i0)') igrdaidx(ig)
        else
          write(6,'(''   z coordinate of atom'',i0)') igrdaidx(ig)
        endif
   20 continue
      write(6,'(''- Displacement of x,y,z coordinates - delgrdxyz = '',f7.5)') delgrdxyz
      write(6,*)
           
      return
      end
c -----------------------------------------------------------------------

c   Subroutine which at the start up prints out information about the 
c   energy gradients (z matrix/internal).
      subroutine inpwrt_grdnts_zmat()
      use grdntsmv, only: igrdaidx, igrdcidx, igrdmv

      use grdntspar, only: delgrdba, delgrdbl, delgrdda, delgrdxyz, igrdtype, ngradnts
      use zmatrix, only: czcart, czint, czcart_ref, izcmat, izmatrix

      implicit real*8(a-h,o-z)


      include 'force.h'
      include 'vmc.h'

      write(6,*)
      write(6,'(''Correlated sampling used to calculate energy gradients (Z matrix)'')')
      write(6,'(''- Number of gradients: '',i0)') ngradnts
      write(6,'(''- Gradients calculated for:'')')
      do 20 ig=1,ngradnts
        na=igrdaidx(ig)
        nb=izcmat(1,na)
        nc=izcmat(2,na)
        nd=izcmat(3,na)
        if(igrdcidx(ig).eq.1) then
          write(6,'(''   bond length between atoms '',i0,'' and '',i0)') na,nb
        elseif(igrdcidx(ig).eq.2) then
          write(6,'(''   bond angle between atoms '',i0,'','',i0,'' and '',i0)') na,nb,nc
        else
          write(6,'(''   dihedral angle between atoms '',i0,'','',i0,'','',i0,'' and '',i0)') na,nb,nc,nd
        endif
   20 continue
      write(6,'(''- Displacement of bond lengths      -   delgrdbl = '',f7.5)') delgrdbl
      write(6,'(''- Displacement of bond angles       -   delgrdba = '',f7.5)') delgrdba
      write(6,'(''- Displacement of dihedral angles   -   delgrdda = '',f7.5)') delgrdda
      write(6,*)
           
      return
      end
c -----------------------------------------------------------------------

c   Subroutine which calculates and printouts energy gradients 
c   for cartesian coordinates of atoms from energy differences 
c   calculated using correlated smapling.
      subroutine finwrt_grdnts_cart(forces_ave,forces_err)
      use atom, only: znuc, cent, pecent, iwctype, nctype, ncent
      use forcepar, only: deltot, istrech, nforce
      use grdntsmv, only: igrdaidx, igrdcidx, igrdmv

      use grdntspar, only: delgrdba, delgrdbl, delgrdda, delgrdxyz, igrdtype, ngradnts

      implicit real*8(a-h,o-z)



      include 'force.h'
      include 'vmc.h'

      dimension forces_ave(MFORCE),forces_err(MFORCE)
      dimension grdnts_ave(MFORCE),grdnts_err(MFORCE)     

      do 5 if=1,nforce-2,2
        grdnts_ave((if+1)/2)=((-forces_ave(if)+forces_ave(if+1))/delgrdxyz)*0.5d0
    5   grdnts_err((if+1)/2)=(dsqrt(forces_err(if)**2+forces_err(if+1)**2)/delgrdxyz)*0.5d0
 
      write(6,*)
      write(6,'(''Energy gradients (dE/d{x,y,z}):'')')
      do 10 ig=1,ngradnts
        if(igrdcidx(ig).eq.1) then
          write(6,'(''  x coordinate of atom'',i3,'': '',f10.7,'' +- '',f10.7)') 
     &      igrdaidx(ig),grdnts_ave(ig),grdnts_err(ig)
        elseif(igrdcidx(ig).eq.2) then
          write(6,'(''  y coordinate of atom'',i3,'': '',f10.7,'' +- '',f10.7)') 
     &      igrdaidx(ig),grdnts_ave(ig),grdnts_err(ig)
        else
          write(6,'(''  z coordinate of atom'',i3,'': '',f10.7,'' +- '',f10.7)') 
     &      igrdaidx(ig),grdnts_ave(ig),grdnts_err(ig)
        endif
   10 continue
      
      write(6,*)
      
      ig=1
      write(6,'(''Energy gradients (Cartesian) - atom table:'')')
      write(6,'(1x,''Atom'',7x,''dE/dx'',21x,''dE/dy'',21x,''dE/dz'',15x,''iwctype'')')
      do 20 ic=1,ncent
        write(6,'(1x,i3,3x)',advance='no') ic
        do 30 k=1,3          
          if(igrdmv(k,ic).eq.1) then
            write(6,'(f10.7,'' +- '',f9.7,3x)',advance='no') grdnts_ave(ig),grdnts_err(ig)
            ig=ig+1
          else
            write(6,'(a10,'' +- '',a9,3x)',advance='no') 'x.xxxxxxx','x.xxxxxxx'
          endif
   30   continue
   20   write(6,'(i3)') iwctype(ic)
      write(6,'('' ---------------------------------------------------------------------------------------'')')
      write(6,*)     

      return
      end     
c -----------------------------------------------------------------------

c   Subroutine which calculates and printouts energy gradients 
c   for Z matrix coordinates of atoms from energy differences 
c   calculated using correlated smapling.
      subroutine finwrt_grdnts_zmat(forces_ave,forces_err)
      use atom, only: znuc, cent, pecent, iwctype, nctype, ncent
      use forcepar, only: deltot, istrech, nforce
      use forcestr, only: delc
      use grdntsmv, only: igrdaidx, igrdcidx, igrdmv

      use grdntspar, only: delgrdba, delgrdbl, delgrdda, delgrdxyz, igrdtype, ngradnts
      use zmatrix, only: czcart, czint, czcart_ref, izcmat, izmatrix

      implicit real*8(a-h,o-z)

      include 'force.h'
      include 'vmc.h'

      dimension forces_ave(MFORCE),forces_err(MFORCE)
      dimension grdnts_ave(MFORCE),grdnts_err(MFORCE)     




      do 5 ifr=1,nforce-2,2
        if(igrdcidx((ifr+1)/2).eq.1) then
          delgrd_tmp=delgrdbl
        elseif(igrdcidx((ifr+1)/2).eq.2) then
          delgrd_tmp=delgrdba
        else        
          delgrd_tmp=delgrdda
        endif
        grdnts_ave((ifr+1)/2)=((-forces_ave(ifr)+forces_ave(ifr+1))/delgrd_tmp)*0.5d0
    5   grdnts_err((ifr+1)/2)=(dsqrt(forces_err(ifr)**2+forces_err(ifr+1)**2)/delgrd_tmp)*0.5d0
 
      write(6,*)
      write(6,'(''Energy gradients (dE/d{bl,ba,da}):'')')
      do 10 ig=1,ngradnts
        na=igrdaidx(ig)
        nb=izcmat(1,na)
        nc=izcmat(2,na)
        nd=izcmat(3,na)
        if(igrdcidx(ig).eq.1) then
          write(6,'('' bond length between atoms '',i0,'' and '',i0,'': '',t50,f10.7,'' +- '',f10.7)') 
     &            na,nb,grdnts_ave(ig),grdnts_err(ig)
        elseif(igrdcidx(ig).eq.2) then
          write(6,'('' bond angle between atoms '',i0,'','',i0,'' and '',i0,'': '',t50,f10.7,'' +- '',f10.7)') 
     &            na,nb,nc,grdnts_ave(ig),grdnts_err(ig)
        else
          write(6,'('' dihedral angle between atoms '',i0,'','',i0,'','',i0,'' and '',i0,'': '',t50,f10.7,'' +- '',f10.7)') 
     &            na,nb,nc,nd,grdnts_ave(ig),grdnts_err(ig)
        endif
   10 continue
      
      write(6,*)
      
      ig=1
      write(6,'(''# Energy gradients (Z matrix) - atom table:'')')
      write(6,'(''#'',1x,''Atom'',3x,''izcmat'',11x,''dE/dbl'',20x,''dE/dba'',20x,''dE/dda'',15x,''iwctype'')')
      ic=1
      write(6,'(1x,i3,3x,3(2x,2x),2x)',advance='no') ic
      write(6,'(3(10x,4x,9x,3x))',advance='no')
      write(6,'(i3)') iwctype(ic)
      ic=2
      write(6,'(1x,i3,3x,1(i2,2x),2(2x,2x),2x)',advance='no') ic,izcmat(1,ic)
      k=1      
      if(igrdmv(k,ic).eq.1) then
        write(6,'(f10.7,'' +- '',f9.7,3x)',advance='no') grdnts_ave(ig),grdnts_err(ig)
        ig=ig+1
      else
        write(6,'(a10,'' +- '',a9,3x)',advance='no') 'x.xxxxxxx','x.xxxxxxx'
      endif
      write(6,'(2(10x,4x,9x,3x))',advance='no')
      write(6,'(i3)') iwctype(ic)
      ic=3
      write(6,'(1x,i3,3x,2(i2,2x),1(2x,2x),2x)',advance='no') ic,(izcmat(k,ic),k=1,2)
      do 15 k=1,2
        if(igrdmv(k,ic).eq.1) then
          write(6,'(f10.7,'' +- '',f9.7,3x)',advance='no') grdnts_ave(ig),grdnts_err(ig)
          ig=ig+1
        else
          write(6,'(a10,'' +- '',a9,3x)',advance='no') 'x.xxxxxxx','x.xxxxxxx'
        endif
   15 continue
      write(6,'(1(10x,4x,9x,3x))',advance='no')
      write(6,'(i3)') iwctype(ic)
      do 20 ic=4,ncent
        write(6,'(1x,i3,3x,3(i2,2x),2x)',advance='no') ic,(izcmat(k,ic),k=1,3)
        do 30 k=1,3          
          if(igrdmv(k,ic).eq.1) then
            write(6,'(f10.7,'' +- '',f9.7,3x)',advance='no') grdnts_ave(ig),grdnts_err(ig)
            ig=ig+1
          else
            write(6,'(a10,'' +- '',a9,3x)',advance='no') 'x.xxxxxxx','x.xxxxxxx'
          endif
   30   continue
   20   write(6,'(i3)') iwctype(ic)
      write(6,'(''# ---------------------------------------------------------------------------------------'')')
      write(6,*)     


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
      do 55 k=1,2
        if(igrdmv(k,ic).eq.1) then
          write(2,'(f10.7,'' +- '',f9.7,3x)',advance='no') grdnts_ave(ig),grdnts_err(ig)
          ig=ig+1
        else
          write(2,'(a10,'' +- '',a9,3x)',advance='no') '0.0000000','0.0000000'
        endif
   55 continue
      write(2,'(1(10x,4x,9x,3x))',advance='no')
      write(2,'(i3)') iwctype(ic)
      do 60 ic=4,ncent
        write(2,'(1x,i3,3x,3(i2,2x),2x)',advance='no') ic,(izcmat(k,ic),k=1,3)
        do 70 k=1,3          
          if(igrdmv(k,ic).eq.1) then
            write(2,'(f10.7,'' +- '',f9.7,3x)',advance='no') grdnts_ave(ig),grdnts_err(ig)
            ig=ig+1
          else
            write(2,'(a10,'' +- '',a9,3x)',advance='no') '0.0000000','0.0000000'
          endif
   70   continue
   60   write(2,'(i3)') iwctype(ic)
      write(2,'(''# ---------------------------------------------------------------------------------------'')')
      write(2,'(''# '')')     

      close(2)


      return
      end     
c -----------------------------------------------------------------------

c   Subroutine which calculates the displacement for energy gradients
c   using Z matrix (internal) coordinates
      subroutine grdzmat_displ(k_in,ic_in,ia_in,delfactor)
      use atom, only: znuc, cent, pecent, iwctype, nctype, ncent
      use forcestr, only: delc
      use grdntsmv, only: igrdaidx, igrdcidx, igrdmv

      use grdntspar, only: delgrdba, delgrdbl, delgrdda, delgrdxyz, igrdtype, ngradnts
      use zmatrix, only: czcart, czint, czcart_ref, izcmat, izmatrix

      implicit real*8(a-h,o-z)



      include 'force.h'
      include 'vmc.h'

      dimension czint_t1(3,MCENT),czcart_t1(3,MCENT)


      if(k_in.eq.1) then
        delgrd_tmp=delfactor*delgrdbl
      elseif(k_in.eq.2) then
        delgrd_tmp=delfactor*delgrdba
      else 
        delgrd_tmp=delfactor*delgrdda
      endif
      
      do 20 ic=1,ncent
        do 20 k=1,3
          czint_t1(k,ic)=czint(k,ic)
   20     czcart_t1(k,ic)=0.0d0
      

      czint_t1(k_in,ic_in)=czint_t1(k_in,ic_in)+delgrd_tmp
      call zmat2cart_rc(MCENT,izcmat,czint_t1,czcart_t1,czcart_ref)

      do 30 ic=1,ncent
        do 30 k=1,3
   30     delc(k,ic,ia_in)=czcart_t1(k,ic)-czcart(k,ic)
      
      return
      end
c -----------------------------------------------------------------------

 
c   Subroutine which prints out at the start of a run 
c   information regarding the Z matrix.
      subroutine inpwrt_zmatrix()
      use atom, only: znuc, cent, pecent, iwctype, nctype, ncent
      use zmatrix, only: czcart, czint, czcart_ref, izcmat, izmatrix
      implicit real*8(a-h,o-z)
      include 'force.h'
      include 'vmc.h'

      write(6,'(''---------- Z matrix information ----------'')')
      write(6,'(''Geometry in internal coordinates:'')')
      write(6,'(1x,''Atom'',t9,''izcmat'',t24,''bond length'',t37,''bond angle'',t51,''dihed. angle '',t65,''iwctype'')')
      ic=1
      write(6,'(1x,i3,59x,i3)') ic,iwctype(ic)
      ic=2
      write(6,'(1x,i3,3x,1(i2,2x),10x,1(f11.7,3x),28x,i3)') ic,izcmat(1,ic),czint(1,ic),iwctype(ic)
      ic=3
      write(6,'(1x,i3,3x,2(i2,2x),6x,1(f11.7,3x),1(f11.6,3x),14x,i3)') ic,(izcmat(k,ic),k=1,2),(czint(k,ic),k=1,2),iwctype(ic)
      do 20 ic=4,ncent
   20   write(6,'(1x,i3,3x,3(i2,2x),2x,1(f11.7,3x),2(f11.6,3x),i3)') ic,(izcmat(k,ic),k=1,3),(czint(k,ic),k=1,3),iwctype(ic)
      write(6,'('' ----------------------------------------------------------------------'')')
      write(6,'(''Internal coordiantes converted back into cartesian coordinates:'')')
      write(6,'(1x,''Atom'',t12,''x'',t26,''y'',t40,''z'',t51,''iwctype'')')
      do 30 ic=1,ncent
   30   write(6,'(1x,i3,3x,3(f11.7,3x),i3)') ic,(czcart(k,ic),k=1,3),iwctype(ic)
      write(6,'('' --------------------------------------------------------'')')
      write(6,*)

      return
      end     
c -----------------------------------------------------------------------



c   Subroutine which calculates and print outs the diagonal
c   part of the Hessian for Z matrix coordinates of atoms 
c   from energy differences  calculated using correlated smapling.
      subroutine finwrt_diaghess_zmat(forces_ave,forces_err)
      use atom, only: znuc, cent, pecent, iwctype, nctype, ncent
      use forcepar, only: deltot, istrech, nforce
      use forcestr, only: delc
      use grdntsmv, only: igrdaidx, igrdcidx, igrdmv

      use grdntspar, only: delgrdba, delgrdbl, delgrdda, delgrdxyz, igrdtype, ngradnts
      use zmatrix, only: czcart, czint, czcart_ref, izcmat, izmatrix

      implicit real*8(a-h,o-z)




      include 'force.h'
      include 'vmc.h'

      dimension forces_ave(MFORCE),forces_err(MFORCE)
      dimension diaghess_ave(MFORCE),diaghess_err(MFORCE)     




      do 5 ifr=1,nforce-2,2
        if(igrdcidx((ifr+1)/2).eq.1) then
          delgrd_tmp=delgrdbl
        elseif(igrdcidx((ifr+1)/2).eq.2) then
          delgrd_tmp=delgrdba
        else        
          delgrd_tmp=delgrdda
        endif        
        diaghess_ave((ifr+1)/2)=((-forces_ave(ifr)-forces_ave(ifr+1))/delgrd_tmp**2)
    5   diaghess_err((ifr+1)/2)=(dsqrt(forces_err(ifr)**2+forces_err(ifr+1)**2)/delgrd_tmp**2)
 
      write(6,*)
      write(6,'(''Diagonal of the Hessian (d^2E/d{bl,ba,da}^2):'')')
      do 10 ig=1,ngradnts
        na=igrdaidx(ig)
        nb=izcmat(1,na)
        nc=izcmat(2,na)
        nd=izcmat(3,na)
        if(igrdcidx(ig).eq.1) then
          write(6,'('' bond length between atoms '',i0,'' and '',i0,'': '',t50,f10.7,'' +- '',f10.7)') 
     &            na,nb,diaghess_ave(ig),diaghess_err(ig)
        elseif(igrdcidx(ig).eq.2) then
          write(6,'('' bond angle between atoms '',i0,'','',i0,'' and '',i0,'': '',t50,f10.7,'' +- '',f10.7)') 
     &            na,nb,nc,diaghess_ave(ig),diaghess_err(ig)
        else
          write(6,'('' dihedral angle between atoms '',i0,'','',i0,'','',i0,'' and '',i0,'': '',t50,f10.7,'' +- '',f10.7)') 
     &            na,nb,nc,nd,diaghess_ave(ig),diaghess_err(ig)
        endif
   10 continue
      
      write(6,*)
      
      ig=1
      write(6,'(''# Diagonal of the Hessian (Z matrix) - atom table:'')')
      write(6,'(''#'',1x,''Atom'',3x,''izcmat'',11x,''d^2E/dbl^2'',20x,''d^2E/dba^2'',20x,''d^2E/dda^2'',15x,''iwctype'')')
      ic=1
      write(6,'(1x,i3,3x,3(2x,2x),2x)',advance='no') ic
      write(6,'(3(10x,4x,9x,3x))',advance='no')
      write(6,'(i3)') iwctype(ic)
      ic=2
      write(6,'(1x,i3,3x,1(i2,2x),2(2x,2x),2x)',advance='no') ic,izcmat(1,ic)
      k=1      
      if(igrdmv(k,ic).eq.1) then
        write(6,'(f10.7,'' +- '',f9.7,3x)',advance='no') diaghess_ave(ig),diaghess_err(ig)
        ig=ig+1
      else
        write(6,'(a10,'' +- '',a9,3x)',advance='no') 'x.xxxxxxx','x.xxxxxxx'
      endif
      write(6,'(2(10x,4x,9x,3x))',advance='no')
      write(6,'(i3)') iwctype(ic)
      ic=3
      write(6,'(1x,i3,3x,2(i2,2x),1(2x,2x),2x)',advance='no') ic,(izcmat(k,ic),k=1,2)
      do 15 k=1,2
        if(igrdmv(k,ic).eq.1) then
          write(6,'(f10.7,'' +- '',f9.7,3x)',advance='no') diaghess_ave(ig),diaghess_err(ig)
          ig=ig+1
        else
          write(6,'(a10,'' +- '',a9,3x)',advance='no') 'x.xxxxxxx','x.xxxxxxx'
        endif
   15 continue
      write(6,'(1(10x,4x,9x,3x))',advance='no')
      write(6,'(i3)') iwctype(ic)
      do 20 ic=4,ncent
        write(6,'(1x,i3,3x,3(i2,2x),2x)',advance='no') ic,(izcmat(k,ic),k=1,3)
        do 30 k=1,3          
          if(igrdmv(k,ic).eq.1) then
            write(6,'(f10.7,'' +- '',f9.7,3x)',advance='no') diaghess_ave(ig),diaghess_err(ig)
            ig=ig+1
          else
            write(6,'(a10,'' +- '',a9,3x)',advance='no') 'x.xxxxxxx','x.xxxxxxx'
          endif
   30   continue
   20   write(6,'(i3)') iwctype(ic)
      write(6,'(''# ---------------------------------------------------------------------------------------'')')
      write(6,*)    
  
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
      do 55 k=1,2
        if(igrdmv(k,ic).eq.1) then
          write(2,'(f10.7,'' +- '',f9.7,3x)',advance='no') diaghess_ave(ig),diaghess_err(ig)
          ig=ig+1
        else
          write(2,'(a10,'' +- '',a9,3x)',advance='no') '0.0000000','0.0000000'
        endif
   55 continue
      write(2,'(1(10x,4x,9x,3x))',advance='no')
      write(2,'(i3)') iwctype(ic)
      do 60 ic=4,ncent
        write(2,'(1x,i3,3x,3(i2,2x),2x)',advance='no') ic,(izcmat(k,ic),k=1,3)
        do 70 k=1,3          
          if(igrdmv(k,ic).eq.1) then
            write(2,'(f10.7,'' +- '',f9.7,3x)',advance='no') diaghess_ave(ig),diaghess_err(ig)
            ig=ig+1
          else
            write(2,'(a10,'' +- '',a9,3x)',advance='no') '0.0000000','0.0000000'
          endif
   70   continue
   60   write(2,'(i3)') iwctype(ic)
      write(2,'(''# ---------------------------------------------------------------------------------------'')')
      write(2,'(''# '')')    
  
      close(2)

      return
      end     
c -----------------------------------------------------------------------
      subroutine transform_grad_zmat(force_cart)
      use atom, only: znuc, cent, pecent, iwctype, nctype, ncent

      use grdntsmv, only: igrdaidx, igrdcidx, igrdmv
      use zmatrix, only: czcart, czint, czcart_ref, izcmat, izmatrix

      implicit real*8(a-h,o-z)

      include 'force.h'
      include 'vmc.h'

      parameter (eps=1.d-5,epsi=0.5d0/eps)
      parameter (MCENT3=3*MCENT)

      common /zmatrix_grad/ transform_grd(MCENT3,MCENT3)

      common /force_analy/ iforce_analy,iuse_zmat

      dimension czcartp(3,MCENT),czcartm(3,MCENT)
      dimension bwilson(MCENT3,MCENT3),gmat(MCENT3*MCENT3)

      dimension force_cart(3,*),force_int(MCENT3)

      do 10 ic=1,3
        do 10 k=1,3
   10     czcart_ref(k,ic)=cent(k,ic)

c     write(6,*) 'HELLO'
c     do 20 ic=1,ncent
c 20    write(6,'(3i5)') (izcmat(k,ic),k=1,3)

c     call cart2zmat(MCENT,czcart,izcmat,czint)
c     write(6,*) 'HELLO'
c     do 20 ic=1,ncent
c 20    write(6,'(1p3e12.5)') (czint(k,ic),k=1,3)
c     call zmat2cart(MCENT,izcmat,czint,czcart)
 
c     write(6,*) 'HELLO'
c     do 21 ic=1,ncent
c 21    write(6,'(1p3e12.5)') (czcart(k,ic),k=1,3)

c     return

      ii=0
      do 100 ic=2,ncent
        do 100 k=1,3
          if(ic.eq.2.and.k.ne.1) goto 100
          if(ic.eq.3.and.k.eq.3) goto 100

          ii=ii+1
          czint_sav=czint(k,ic)

          czint(k,ic)=czint_sav+eps
          call zmat2cart_rc(ncent,izcmat,czint,czcartp,czcart_ref)

          do 11 lc=1,3
            do 11 l=1,3
   11         czcart_ref(l,lc)=cent(l,lc)

          czint(k,ic)=czint_sav-eps
          call zmat2cart_rc(ncent,izcmat,czint,czcartm,czcart_ref)

          jj=0
          do 50 jc=1,ncent
            do 50 l=1,3
              jj=jj+1
              bwilson(ii,jj)=(czcartp(l,jc)-czcartm(l,jc))*epsi
  50      continue

          czint(k,ic)=czint_sav
 100  continue

      ncent_ind=3*ncent-6
      if(ncent.eq.2) ncent_ind=1
      ncent3=3*ncent
c     do ic=1,ncent_ind
c       write(6,*) 'Bmat',(bwilson(ic,jc),jc=1,ncent3)
c     enddo

      do 150 ic=1,ncent_ind
        force_int(ic)=0
        jj=0
        do 150 jc=1,ncent
          do 150 l=1,3
          jj=jj+1
 150      force_int(ic)=force_int(ic)+bwilson(ic,jj)*force_cart(l,jc)
          
      open(81,file='force_analytic_zmat',form='formatted',status='unknown')
      write(81,'(i5,1p3e14.5)') 1,0.d0,0.d0,0.d0
      write(81,'(i5,1p3e14.5)') 2,force_int(1)*igrdmv(1,2),0.d0,0.d0
      if(ncent.gt.2) write(81,'(i5,1p3e14.5)') 3,force_int(2)*igrdmv(1,3),force_int(3)*igrdmv(2,3),0.d0
      ii=3
      do 200 ic=4,ncent_ind,3
        ii=ii+1
 200    write(81,'(i5,1p3e14.5)') ii,(force_int(ic+k)*igrdmv(k+1,ii),k=0,2)
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
        do 300 ic=4,ncent_ind,3
          ii=ii+1
          do 300 k=0,2
 300        force_cart(k+1,ii)=force_int(ic+k)*igrdmv(k+1,ii)
      endif

      return
      end
c -----------------------------------------------------------------------
