!*********************************************************************
        subroutine qmmm_density_ini(znuc,cent,iwctype,mctype,mcent, &
      &  ncent)
!*********************************************************************

        use qmmm_density

        implicit none

        integer :: mctype,mcent,ncent
        integer :: iwctype(mcent)
        double precision :: znuc(mctype),cent(3,mcent),tmp
        integer :: i

        title_dens(1)='Variational Monte Carlo Density'
        title_dens(2)=''
        n_atomsd=ncent
        allocate(x_atomd(n_atomsd,3),id_atomd(n_atomsd))
        allocate(chrg_atomd(n_atomsd))

        cc_nuc(:)=0.d0
        cc_ele(:)=0.d0
        cc_ele2(:)=0.d0
        dipole(:)=0.d0
        dipole2(:)=0.d0
        cdipole=0
        tmp=0.d0
        do i=1,n_atomsd
          x_atomd(i,:)=cent(:,i)
!         write (*,*) i,iwctype(i),znuc(iwctype(i))
          id_atomd(i)=znuc(iwctype(i))
          if(znuc(iwctype(i)).gt.2) id_atomd(i)=znuc(iwctype(i))+2
          chrg_atomd(i)=0.d0
          cc_nuc(:)=cc_nuc(:)+znuc(iwctype(i))*cent(:,i)
          tmp=tmp+znuc(iwctype(i))
        enddo
        cc_nuc(:)=cc_nuc(:)/tmp
!       write (*,*) 'Center of nuclear charge:', cc_nuc(:),tmp

! define the box size
        n_xd=31;n_yd=31;n_zd=31
        x0d(1) = -7.d0
        x0d(2) = -7.d0
        x0d(3) = -7.d0
         deltad(1)= dabs(x0d(1)*2) / (n_xd-1)
         deltad(2)= dabs(x0d(2)*2) / (n_yd-1)
         deltad(3)= dabs(x0d(3)*2) / (n_zd-1)

!       n_xd=31;n_yd=31;n_zd=31
!       x0d(1) = 0.d0
!       x0d(2) = 0.d0
!       x0d(3) = 0.d0
!       deltad(1)=.516129
!       deltad(2)=.516129
!       deltad(3)=.516129

        allocate(dens(n_xd,n_yd,n_zd))
        allocate(sme(n_xd,n_yd,n_zd))
        dens(:,:,:)=0.d0
        sme(:,:,:)=0.d0
        outofbox=0.d0
        inofbox=0.d0
        outofboxs=0.d0
        inofboxs=0.d0

!        call qmmm_writecube("test.cube",title_dens,n_atomsd,x0d, &
!     &   n_xd,n_yd,n_zd,deltad,x_atomd,id_atomd,chrg_atomd,dens)

        return
        end

!*********************************************************************
        subroutine qmmm_density_accu(nelec,xold,weight)
!*********************************************************************

        use qmmm_density

        implicit none

        integer :: nelec
        double precision :: xold(3,nelec),weight,sumw,cc_dip(3)
        integer :: n, ix,iy,iz,jx,jy,jz,k,in
        double precision :: x_center(3),dist
        integer :: ind(3,125),maxind
        double precision :: value(125),ssme,win,wout

!       write (*,*) 'DENSITY_ACCU'
!       write (*,*) (xold(3,n), n=1,nelec)

        sumw=0.d0
        cc_dip(:)=0.d0
        do n=1,nelec
          ix = int ( (xold(1,n)-x0d(1)+deltad(1)/2)/(deltad(1)) ) + 1
          iy = int ( (xold(2,n)-x0d(2)+deltad(2)/2)/(deltad(2)) ) + 1
          iz = int ( (xold(3,n)-x0d(3)+deltad(3)/2)/(deltad(3)) ) + 1
          x_center(1) = x0d(1)+ deltad(1)*(ix-1)
!          write (*,*) xold(1,n)*Ang,ix,x_center
! calculate electronic center of charge
             cc_ele(:)=cc_ele(:)+xold(:,n)*weight
             cc_ele2(:)=cc_ele2(:)+(xold(:,n)*weight)**2
             cc_dip(:)=cc_dip(:)+xold(:,n)*weight
             sumw=sumw+weight
          if(ix.lt.1 .or. iy.lt.1 .or. iz.lt.1 .or. &
&            ix.gt.n_xd .or. iy.gt.n_yd .or. iz.gt.n_zd ) then
!            write (45,*) 'Walker out of the density box! ',xold(:,n)
             outofbox=outofbox+weight
             outofboxs=outofboxs+weight
          else
             inofbox=inofbox+weight
             dens(ix,iy,iz)=dens(ix,iy,iz)+weight
! smearing
! La densita' dell'elettrone viene sparsa sui primi vicini
! moltiplicata per un fattore gaussiano e poi rinormalizzata

             ssme=deltad(1)*1.0
             in=0
             if(1.ne.0) then
!                  write (*,*) ix,iy,iz
!                   write (*,*) "deltad",deltad(:)
               win=0.d0
               wout=0.d0
               value(:)=0.d0
               do jx=ix-2,ix+2
                 do jy=iy-2,iy+2
                   do jz=iz-2,iz+2
                     if(jx.lt.1 .or. jy.lt.1 .or. jz.lt.1 .or. &
&                    jx.gt.n_xd .or. jy.gt.n_yd .or. jz.gt.n_zd ) then
                       x_center(1) = x0d(1)+ deltad(1)*(jx-1)
                       x_center(2) = x0d(2)+ deltad(2)*(jy-1)
                       x_center(3) = x0d(3)+ deltad(3)*(jz-1)
                       dist=dsqrt(dot_product(xold(:,n)-x_center(:), &
&                                           xold(:,n)-x_center(:)))
                       wout=wout+exp(-((dist/ssme)**2))
                     else
                       in=in+1
                       ind(1,in)=jx
                       ind(2,in)=jy
                       ind(3,in)=jz
                       x_center(1) = x0d(1)+ deltad(1)*(jx-1)
                       x_center(2) = x0d(2)+ deltad(2)*(jy-1)
                       x_center(3) = x0d(3)+ deltad(3)*(jz-1)
                       dist=dsqrt(dot_product(xold(:,n)-x_center(:), &
&                                           xold(:,n)-x_center(:)))
                       value(in)=exp(-((dist/ssme)**2))
                       win=win+value(in)
!                   write (*,'(a5,4i4,3f12.8)') 'XXX', in, ind(1,in),ind(2,in),ind(3,in),dist,dist/ssme,value(in)
!                  write (*,'(7f8.3)') xold(:,n),x_center(:),dist
!                    write (*,'(7f12.8)') dist,value(in)
                     endif
                   enddo
                 enddo
               enddo
!              write (*,*) SUM(value(:),1),win,wout
               value(:)=value(:)/(win+wout)*weight
               inofboxs=inofboxs+weight*win/(win+wout)
               outofboxs=outofboxs+weight*wout/(win+wout)
              do k=1,in
                sme(ind(1,k),ind(2,k),ind(3,k)) = &
&                sme(ind(1,k),ind(2,k),ind(3,k)) + value(k)
              enddo
             endif
          endif
        enddo
        dipole(:)=dipole(:)+(cc_nuc(:)-cc_dip(:)/sumw)*nelec*2.5417
        dipole2(:)=dipole2(:)+ &
&          ( (cc_nuc(:)-cc_dip(:)/sumw)*nelec*2.5417 )**2
        cdipole=cdipole+1
!       write (*,*) 'DIPOLE ==',dipole(:)
!      write (999,*) (cc_nuc(:)-cc_dip(:)/sumw)*nelec*2.5417

        return
        end


!*********************************************************************
        subroutine qmmm_density_write(nelec,id)
!*********************************************************************

        use qmmm_density

        implicit none
        integer :: nelec,id
        double precision :: norm,deltav,totnorm,totnorms,tmp(3)
        character (len=80) :: file,files

        if(id.eq.0) then
          file = "density_vmc.cube"
          files = "densitys_vmc.cube"
        else
          file = "density_dmc.cube"
          files = "densitys_dmc.cube"
        endif

        deltav=deltad(1)*deltad(2)*deltad(3)
        write (*,*) 'deltav =', deltav

! Density
        norm =sum(sum(sum(dens(:,:,:),1),1),1)
        write(*,*) 'number of points inside the box:',inofbox
        write(*,*) 'number of points out of the box:',outofbox
        write (*,*) 'Normalizing with the total number of points'
        totnorm=inofbox+outofbox
        dens(:,:,:) = dens(:,:,:)/totnorm*nelec/deltav

        norm =sum(sum(sum(dens(:,:,:),1),1),1)
        write (*,*) 'Integrated density dens: ', norm*deltav

        call qmmm_writecube(file,title_dens,n_atomsd,x0d, &
     &    n_xd,n_yd,n_zd,deltad,x_atomd,id_atomd,chrg_atomd,dens)

! Smeared Density
!       norm =sum(sum(sum(sme(:,:,:),1),1),1)
        write (*,*) 'Smeared density:'
        write(*,*) 'number of points inside the box:',inofboxs
        write(*,*) 'number of points out of the box:',outofboxs
        write (*,*) 'Normalizing with the total number of points'
        totnorms=inofboxs+outofboxs
        sme(:,:,:) = sme(:,:,:)/totnorms*nelec/deltav

        norm =sum(sum(sum(sme(:,:,:),1),1),1)
        write (*,*) 'Integrated density sme: ', norm*deltav

        call qmmm_writecube(files,title_dens,n_atomsd,x0d, &
     &    n_xd,n_yd,n_zd,deltad,x_atomd,id_atomd,chrg_atomd,sme)

!...............Correct Dipole
        dipole(:)=dipole(:)/cdipole
        dipole2(:)=dipole2(:)/cdipole
        tmp(:)=dsqrt((dipole2(:)-dipole(:)**2)/cdipole)
        write (*,*) 'Mixed Average: Dipole moment (Debye) (x,y,z,|.|)'
        write(*,'(a12,4f12.6)') 'Dipole: ', dipole(:), &
     &               dsqrt(dot_product(dipole(:),dipole(:)))
        write(*,'(a12,4f12.6)') 'Dip. err.: ', tmp(:)
        write(*,'(a32,i16)') 'Total number of dipoles sampled:  ',cdipole
!.........................

!       cc_ele(:)=cc_ele(:)/totnorm
!       cc_ele2(:)=cc_ele2(:)/totnorm
!       tmp(:)=dsqrt((cc_ele2(:)-cc_ele(:)**2)/totnorm)
!       write (*,*) 'Center of electronic charge:', cc_ele(:)
!       dipole(:)=(cc_nuc(:)-cc_ele(:))*nelec*2.5417
!       tmp(:)=tmp(:)*nelec*2.5417
!       write (*,*) 'Mixed Average: Dipole moment (Debye) (x,y,z,|.|)'
!       write(*,'(a12,4f12.6)') 'Dipole: ', dipole(:), &
!&                    dsqrt(dot_product(dipole(:),dipole(:)))
!       write(*,'(a12,4f12.6)') 'Dip. err.: ', tmp(:)

        return

        end

!*********************************************************************
!       subroutine qmmm_prop_cc_nuc(znuc,cent,iwctype,mctype,mcent, &
!     &  ncent,cc_nuc)
!*********************************************************************

!       implicit none

!       integer :: mctype,mcent,ncent
!       integer :: iwctype(mcent)
!       double precision :: znuc(mctype),cent(3,mcent)
!       double precision :: cc_nuc(3), tmp
!       integer :: i,id
!
!       cc_nuc(:)=0.d0
!       tmp=0.d0
!       do i=1,ncent
!         id=znuc(iwctype(i))
!         if(znuc(iwctype(i)).gt.2) id=znuc(iwctype(i))+2
!         cc_nuc(:)=cc_nuc(:)+znuc(iwctype(i))*cent(:,i)
!         tmp=tmp+znuc(iwctype(i))
!       enddo
!       cc_nuc(:)=cc_nuc(:)/tmp
!       write (*,*) 'Center of nuclear charge:', cc_nuc(:),tmp

!       return
!       end

