c......................................................
      subroutine pcm_extpot_read(fcol,npmax)
c Written by Amovilli-Floris 
c...........................................................
c     read data for pcm calculations
c     comput nuclei-qpol interactions (penups,penupv)
c...........................................................

      use atom, only: znuc, cent, iwctype, ncent
      use pcm_cntrl, only: icall, ichpol, ipcm, isurf
      use pcm_unit, only: pcmfile_cavity, pcmfile_chs, pcmfile_chv
      use pcm_parms, only: ch, eps_solv, nch, nchs, nchs1, nchs2
      use pcm_parms, only: nchv, nesph, re
      use pcm_parms, only: surk, xe, xpol, ye, ze
      use pcm_ameta, only: amdlg, eta
      use pcm_pot, only: penupol, penups, penupv
      use pcm_fdc, only: feps, fs, rcol, rcolt

      use pcm_inda, only: inda

      implicit real*8(a-h,o-z)



      DATA PI/3.1415927D0/
      if(ipcm.eq.0)return

      if(isurf.eq.0)then
c Read geometrical properties of cavity
        open (50,file=pcmfile_cavity,form='formatted',status='unknown')
        rewind 50
        read(50,*)
        read(50,*)surk,nchs1,nchs2,nchs
        read(50,*)
        do i=1,nchs
          read(50,*)inda(i),(xpol(k,i),k=1,3),amdlg(i),(eta(k,i),k=1,3)
        enddo
        close(50)

c Read surface charges
        open (50,file=pcmfile_chs,form='formatted',status='unknown')
        rewind 50
        do i=1,nchs
          read(50,*)(xpol(k,i),k=1,3),ch(i)
        enddo
        close(50)

c Read volume charges
        open (50,file=pcmfile_chv,form='formatted',status='unknown')
        rewind 50
        read(50,*) nchv
        nch=nchs+nchv
        do i=nchs+1,nch
          read(50,*)(xpol(k,i),k=1,3),ch(i)
        enddo
        close(50)

       else
        call qpcm_surface(npmax)
        write (6,'(''pcm: a new solute cavity has been created'')')

c write geometrical properties of cavity
        open (50,file=pcmfile_cavity,form='formatted',status='unknown')
        rewind 50
        write (50,*) '  area tess.       nchs1    nchs2   nchs'
        write (50,*) surk,nchs1,nchs2,nchs
        write (50,*)'positions, madelung, normal versors on the surface'
        do j=1,nchs
          write (50,100) inda(j),(xpol(i,j),i=1,3),amdlg(j),(eta(i,j),i=1,3)
        enddo
 100    format (I4,8f10.5)
        close(50)

c write surface charges
        open (50,file=pcmfile_chs,form='formatted',status='unknown')
        rewind 50
        do i=1,nchs
          write(50,101)(xpol(k,i),k=1,3),0.d0
        enddo
 101    format (4f10.5)
        close(50)

c write volume charges
        open (50,file=pcmfile_chv,form='formatted',status='unknown')
        rewind 50
        write(50,*) 0
        close(50)

        stop
      endif
c...............................................
      icall=0
      fs=1.0d0/eps_solv
      feps=(eps_solv-1.0d0)/(4.0d0*pi*eps_solv)
      rcol=fcol*dsqrt(surk/pi)
      rcolt=dsqrt(surk/pi)
c...............................................
      write(6,'(''pcm nchs ='',i10)') nchs
      write(6,'(''pcm nchv ='',i10)') nchv
      write(6,'(''pcm surk ='',f12.6)') surk
c...............................................
      penups=0.0d0
      do i=1,ncent
        do j=1,nchs
          xx=(xpol(1,j)-cent(1,i))**2.0d0
          yy=(xpol(2,j)-cent(2,i))**2.0d0
          zz=(xpol(3,j)-cent(3,i))**2.0d0
          rnp2=xx+yy+zz
          rnp=dsqrt(rnp2)
          penups=penups+0.5d0*znuc(iwctype(i))*ch(j)/rnp
        enddo
      enddo
      call pcm_compute_penupv
      penupol=penups+penupv
      if (ichpol.eq.1)call chnu
c........................................................................      
      write(6,*)
      write(6,'(''pcm number of spheres ='',i5)') nesph
      write(6,'(''pcm cavity geometry'')')
      do i=1,nesph
      write(6,1000)i,xe(i),ye(i),ze(i),re(i)
      enddo
      write(6,*)
      write(6,'(''pcm epot nuclei-surface polarization charges ='',f12.6)') penups
      write(6,'(''pcm epot nuclei-volume polarization charges ='',f12.6)') penupv
      write(6,'(''pcm epot nuclei-polarization charges ='',f12.6)') penupol
      write(6,*)
c........................................................................      
 1000 format(I4,2x,3F12.5,2x,F12.5,2x,F12.5)
c........................................................................      
      return
      end
      
c......................................................
      subroutine pcm_qvol(n)
c Written by Amovilli-Floris

      use pcm_cntrl, only: ipcm
      use pcm_parms, only: nscv

      use pcm_fdc, only: fs, qvol
      implicit real*8(a-h,o-z)



      if(ipcm.eq.0) return

      tmp=(1.0d0-fs)/nscv
      qvol=tmp/n

      return
      end

      subroutine chnu 
C     ***************************************************************
C     contribution from nuclei to polarization charghes
C     ***************************************************************

      use pcm, only: MCHS, MCHV, MSPHERE
      use vmc, only: MELEC, MORB, MBASIS, MDET, MCENT, MCTYPE, MCTYP3X
      use vmc, only: NSPLIN, nrad, MORDJ, MORDJ1, MMAT_DIM, MMAT_DIM2, MMAT_DIM20
      use vmc, only: radmax, delri
      use vmc, only: NEQSX, MTERMS
      use vmc, only: MCENT3, NCOEF, MEXCIT
      use atom, only: znuc, cent, pecent, iwctype, nctype, ncent
      use pcm_parms, only: ch, eps_solv, iscov, nch, nchs, nchs1, nchs2
      use pcm_parms, only: nchv, ncopcm, nesph, nscv, nvopcm, re, re2
      use pcm_parms, only: retk, surk, xe, xpol, ye, ze

      use pcm_ameta, only: amdlg, eta
      use pcm_ah, only: ahca, bh
      use pcm_fdc, only: feps, fs, qfree, qvol, rcol, rcolt, rcolv
      use pcm_inda, only: inda

      implicit real*8(a-h,o-z)





      real*8 ah_vec,det
C     
C    
      dimension ah_vec(MCHS*MCHS),bhn(MCHS),ah(MCHS,MCHS)

c
      DATA PI/3.1415927D0/
c............................................................
c     The matrix ah is computed and inverted 
c     correzione distinta in base ad inda (modifica rispetto 
c     alla versione 5)
c..............................................................
      do k=1,nchs
      cork=1.0d0/(2.0d0*amdlg(k))
c     sum(k)=0.0d0
      do l=1,nchs
      if (l.eq.k)then
      ah(k,l)=1.0d0+0.5d0*(1.d0-eps_solv)/eps_solv
      else
      xx=xpol(1,k)-xpol(1,l) 
      yy=xpol(2,k)-xpol(2,l) 
      zz=xpol(3,k)-xpol(3,l) 
      s1=xx*eta(1,k)
      s2=yy*eta(2,k)
      s3=zz*eta(3,k)
      ss=s1+s2+s3
      xx2=xx**2.0d0
      yy2=yy**2.0d0
      zz2=zz**2.0d0
      rkl2=xx2+yy2+zz2
      rkl3=rkl2**1.5d0
c     aa=sur(l)*ss/rkl3
c     sum(k)=sum(k)+aa
c      ah(k,l)=feps*surk*ss/rkl3
       corr=cork
       if(inda(l).ne.inda(k))corr=1.0d0
       ah(k,l)=feps*surk*ss*corr/rkl3
      endif
      enddo
c     ah(k,l)=ah(k,l)*2.0d0*pi/sum(k)
c     write(6,*)k,'sum =',sum(k)
      enddo
      do i=1,nchs
        do j=1,nchs
          ij=(j-1)*nchs+i
          ah_vec(ij)=ah(i,j)
        enddo
      enddo
c     MODIFIED (avoid the use of automatic array)............................
      call qpcm_matinv(ah_vec,nchs,det)
      do i=1,nchs
        do j=1,nchs
          ij=(j-1)*nchs+i
          ah(i,j)=ah_vec(ij)
        enddo
      enddo
c     call qpcm_matinv(ah,nchs,det)
c..............................................................
c     calcolo di bhn
c..............................................................
      do k=1,nchs
      enk=0.0d0
      do l=1,ncent
      xx=xpol(1,k)-cent(1,l) 
      yy=xpol(2,k)-cent(2,l) 
      zz=xpol(3,k)-cent(3,l) 
      rr2=xx**2+yy**2+zz**2
      rr3=rr2**1.5d0
      cc1=xx*eta(1,k)
      cc2=yy*eta(2,k)
      cc3=zz*eta(3,k)
      cc=cc1+cc2+cc3 
      enk=enk+znuc(iwctype(l))*cc/rr3
      enddo
      bhn(k)=-feps*surk*enk
      bh(k)=bhn(k)
      enddo
      do k=1,nchs
      do l=1,nchs
      ahca(k,l)=ah(k,l)
      enddo
      enddo
c........................................................
 1000 format(8F12.5)
c........................................................
      return
      END


      subroutine qpcm_efield(nelec,coord)
C     ***************************************************************
c     1) electrons are classified with respect to the cavity
c
c    fac =1   e- is in the cavity
c    fac=fs   e- is out of the cavity 
c
c    quopcm= escaped electron charge out of the cavity
c
c     2) volume charges out of the cavity are sampled
c       A volume point charge is associated to each escaped 
c       electron
c
c     3) for the accepted configuration, the normal component 
c        of the electron field  at the point on the surface is computed 
C     ***************************************************************

      use pcm, only: MCHS, MCHV, MSPHERE
      use vmc, only: MELEC, MORB, MBASIS, MDET, MCENT, MCTYPE, MCTYP3X
      use vmc, only: NSPLIN, nrad, MORDJ, MORDJ1, MMAT_DIM, MMAT_DIM2, MMAT_DIM20
      use vmc, only: radmax, delri
      use vmc, only: NEQSX, MTERMS
      use vmc, only: MCENT3, NCOEF, MEXCIT
      use atom, only: znuc, cent, iwctype, ncent
      use pcm_hpsi, only: enfpcm, qopcm
      use pcm_xv_new, only: xv_new
      use pcm_cntrl, only: icall, ichpol, ipcm, ipcmprt, isurf
      use pcm_parms, only: ch, eps_solv, iscov, nch, nchs, nchs1, nchs2
      use pcm_parms, only: nchv, ncopcm, nesph, nscv, nvopcm, re, re2
      use pcm_parms, only: retk, surk, xe, xpol, ye, ze

      use pcm_ameta, only: eta
      use pcm_fdc, only: fs, qvol, rcol, rcolt
      implicit real*8(a-h,o-z)




      dimension coord(3,*)
      dimension fac(100)
C     ***************************************************************
      DATA PI/3.1415927D0/,GC/1.9872159D0/,AV/0.60228D0/
      data hatokc/627.509541d0/  
c............................................................
     
      write (6,*) 'hello! hello! qpcm_efield called -------'

      ncopcm=ncopcm+1
      iscv=mod(ncopcm,iscov)
      qopcm=0.0d0
      do 1 i=1,nchs
   1    enfpcm(i)=0

      do i=1,nelec
        in=0
        fac(i)=1.0d0
        do j=1,nesph
          xx=(coord(1,i)-xe(j))**2.0d0 
          yy=(coord(2,i)-ye(j))**2.0d0 
          zz=(coord(3,i)-ze(j))**2.0d0 
          res2=xx+yy+zz
          if(res2.le.re2(j)) in=in+1
        enddo
        if (in.eq.0)then
          fac(i)=fs
          qopcm=qopcm+1.d0
        endif
      enddo
c ATTENZION modificato 30/3/2009
      if (ipcm.eq.3) return
c..............................................................
c    enfpcm(k) normal componenent of e- field on the point k of
c            cavity surface 
c..............................................................
c             (first set of points 1---->nchs1)
c..............................................................
      do k=1,nchs1
        eek=0.0d0
        do i=1,nelec
          xx=xpol(1,k)-coord(1,i) 
          yy=xpol(2,k)-coord(2,i) 
          zz=xpol(3,k)-coord(3,i) 
          rr2=xx**2+yy**2+zz**2
          rr1=dsqrt(rr2)
          cc1=xx*eta(1,k)
          cc2=yy*eta(2,k)
          cc3=zz*eta(3,k)
          cc=(cc1+cc2+cc3)/rr1
c..............................................................
c  correction for collision between e- and polarization charges 
c  n.b. the fields of e- which are out of the cavity 
c       are scaled 
c..............................................................
          if (rr1.lt.rcol) rr2=rcol**2.0d0
          eek=eek-cc*fac(i)/rr2
        enddo
        enfpcm(k)=eek
      enddo
c..............................................................
c             (second set of points nchs1+1---->nchs)
c..............................................................
c     rcol=fcol*dsqrt(surk/pi)
      do k=nchs1+1,nchs
        eek=0.0d0
        do i=1,nelec
          xx=xpol(1,k)-coord(1,i) 
          yy=xpol(2,k)-coord(2,i) 
          zz=xpol(3,k)-coord(3,i) 
          rr2=xx**2+yy**2+zz**2
          rr1=dsqrt(rr2)
          cc1=xx*eta(1,k)
          cc2=yy*eta(2,k)
          cc3=zz*eta(3,k)
          cc=(cc1+cc2+cc3)/rr1
          if(rr1.lt.rcolt.and.cc.gt.0.0d0) then
            if(fac(i).ne.1.0d0)then
              qopcm=qopcm-1.d0
              fac(i)=1.0d0
            endif
          endif
c..............................................................
c  correction for collision between e- and polarization charges 
c  n.b. the fields of e- which are out of the cavity 
c       are scaled 
c..............................................................
          if (rr1.lt.rcol) rr2=rcol**2.0d0
          eek=eek-cc*fac(i)/rr2
        enddo
        enfpcm(k)=eek
      enddo
c..............................................................
c     samples volume charges
c............................................................
      if(iscv.eq.0)then
        do i=1,nelec
	  if(fac(i).ne.1.0d0)then
	    nvopcm=nvopcm+1
	    xv_new(1,nvopcm)=coord(1,i)
	    xv_new(2,nvopcm)=coord(2,i)
	    xv_new(3,nvopcm)=coord(3,i)
	  endif
	enddo
      endif
      return
      end

      subroutine qpcm_update_vol(iupdate)
c............................................................
c      update of volume charges and penupol
c............................................................

      use pcm, only: MCHS, MCHV, MSPHERE
      use vmc, only: MELEC, MORB, MBASIS, MDET, MCENT, MCTYPE, MCTYP3X
      use vmc, only: NSPLIN, nrad, MORDJ, MORDJ1, MMAT_DIM, MMAT_DIM2, MMAT_DIM20
      use vmc, only: radmax, delri
      use vmc, only: NEQSX, MTERMS
      use vmc, only: MCENT3, NCOEF, MEXCIT
      use atom, only: znuc, cent, iwctype, ncent
      use pcm_xv_new, only: xv_new
      use pcm_cntrl, only: ipcm
      use pcm_parms, only: ch, eps_solv, iscov, nch, nchs, nchs1, nchs2
      use pcm_parms, only: nchv, ncopcm, nesph, nscv, nvopcm, re, re2
      use pcm_parms, only: retk, surk, xe, xpol, ye, ze

      use pcm_pot, only: penupol, penupv
      use pcm_fdc, only: qvol
      implicit real*8(a-h,o-z)




      iupdate=0
      if(ipcm.eq.0.or.ipcm.eq.3) return

      nsco=ncopcm/iscov
c     write(6,*) 'HELLO',nsco,ncopcm
      if(nsco.eq.nscv)then	
         iupdate=1
	 nchv=nvopcm
	 nch=nchs+nchv
	 do ko=1,nchv
	   kn=nchs+ko
  	   xpol(1,kn)=xv_new(1,ko)
	   xpol(2,kn)=xv_new(2,ko)
	   xpol(3,kn)=xv_new(3,ko)
           ch(kn)=qvol
	 enddo

c        write(6,*)
c        write(6,*)'update of volume charges: nchv =',nchv,'nscv =',nscv
c        write(6,*)
         ncopcm=0
         nvopcm=0
      endif
c..............................................................
      return
      end

      subroutine pcm_compute_penupv
c............................................................
c     compute penupv of volume charges 
c............................................................

      use pcm, only: MCHS, MCHV, MSPHERE
      use atom, only: znuc, cent, iwctype, ncent
      use pcm_cntrl, only: ipcm
      use pcm_parms, only: ch, eps_solv, iscov, nch, nchs, nchs1, nchs2
      use pcm_parms, only: nchv, ncopcm, nesph, nscv, nvopcm, re, re2
      use pcm_parms, only: retk, surk, xe, xpol, ye, ze

      use pcm_pot, only: penupv
      implicit real*8(a-h,o-z)


       

      if(ipcm.eq.0) return

      penupv=0.0d0
      do i=1,ncent
        do j=nchs+1,nch
          xx=(xpol(1,j)-cent(1,i))**2.0d0
          yy=(xpol(2,j)-cent(2,i))**2.0d0
          zz=(xpol(3,j)-cent(3,i))**2.0d0
          rnp2=xx+yy+zz
          rnp=dsqrt(rnp2)
          penupv=penupv+0.5d0*znuc(iwctype(i))*ch(j)/rnp
         enddo
       enddo
c..............................................................
      return
      end

      subroutine pcm_write_chvol

      use pcm, only: MCHS, MCHV, MSPHERE
      use pcm_parms, only: ch, eps_solv, iscov, nch, nchs, nchs1, nchs2
      use pcm_parms, only: nchv, ncopcm, nesph, nscv, nvopcm, re, re2
      use pcm_parms, only: retk, surk, xe, xpol, ye, ze

      implicit real*8(a-h,o-z)

      open (52,file='chvol_new',form='formatted',status='unknown')
      rewind 52

      nchtmp=nch-nchs
      write(52,*) nchtmp
      do ko=1,nchtmp
        kn=nchs+ko
        write(52 ,1000)(xpol(k,kn),k=1,3),ch(kn)
      enddo
c..............................................................
 1000 format(4F15.7)
      return
      end

      subroutine pcm_extpot_ene(coord,nelec,pepcms,pepcmv)
c Written by Amovilli-Floris 
c......................................................
c       Calculate e-qpol interactions (pcm)
c       and adds nuclei-qpol interactions
c......................................................
      use pcm, only: MCHS, MCHV, MSPHERE
      use pcm_cntrl, only: icall
      use pcm_parms, only: ch, eps_solv, iscov, nch, nchs, nchs1, nchs2
      use pcm_parms, only: nchv, ncopcm, nesph, nscv, nvopcm, re, re2
      use pcm_parms, only: retk, surk, xe, xpol, ye, ze

      use pcm_pot, only: penups, penupv
      use pcm_grid3d_contrl, only: ipcm_3dgrid
      implicit real*8(a-h,o-z)



      include 'pcm_3dgrid.h'

      dimension coord(3,*)
      DATA PI/3.1415927D0/,GC/1.9872159D0/,AV/0.60228D0/
      data hatokc/627.509541d0/
      icall=icall+1

      pepcms=penups
      pepcmv=penupv
      do 200 i=1,nelec
        if(ipcm_3dgrid.gt.0) then
          ier=0
          call spline_pcm(coord(1,i),pepol_sv,ier)
          if(ier.eq.1) then
            call pcm_extpot_ene_elec(coord(1,i),pepol_s,pepol_v)
            pepol_sv=pepol_s+pepol_v
          endif
          pepcms=pepcms+pepol_sv
         else
          call pcm_extpot_ene_elec(coord(1,i),pepol_s,pepol_v)
          pepcms=pepcms+pepol_s
          pepcmv=pepcmv+pepol_v
        endif
  200 continue

      if (icall.eq.1)then
        write(6,*)
        write(6,'(''nchs='',i6,'' nch='',i6,'' ch(nch)='',1p1d12.4)') nchs,nch,ch(nch)
        write(6,'(''pcm epot solute-polarization charges ='',f12.6)') pepcms+pepcmv
        write(6,*)
c       write(6,*)
c       write(6,*)'pepcms=',pepcms,'penups=',penups
c       write(6,*)'pepcmv=',pepcmv,'penupv=',penupv
c       if(ipcm_3dgrid.gt.0) then
c     write(6,'(''pcm epot solute-polarization charges ='',f12.6)') pepcms     
c       write(6,*)
c       else
c     write(6,'(''pcm epot solute-polarization charges ='',f12.6)') pepcms+pepcmv
c       write(6,*)
c       endif
      endif
 1000 format(6F15.8)
      return
      end

      subroutine pcm_extpot_ene_elec(x,pepol_s,pepol_v)
c Written by Amovilli-Floris
c......................................................
c       Calculate e-qpol interactions (pcm)
c       and adds nuclei-qpol interactions
c......................................................

      use pcm, only: MCHS, MCHV, MSPHERE
      use pcm_parms, only: ch, eps_solv, iscov, nch, nchs, nchs1, nchs2
      use pcm_parms, only: nchv, ncopcm, nesph, nscv, nvopcm, re, re2
      use pcm_parms, only: retk, surk, xe, xpol, ye, ze

      use pcm_fdc, only: rcol, rcolv
      implicit real*8(a-h,o-z)



      dimension x(3)
      DATA PI/3.1415927D0/,GC/1.9872159D0/,AV/0.60228D0/

      pepol_s=0.0d0
c......................................................
c     interaction with surface point charges
c......................................................
      do j=1,nchs
        xx=(x(1)-xpol(1,j))**2.0d0
        yy=(x(2)-xpol(2,j))**2.0d0
        zz=(x(3)-xpol(3,j))**2.0d0
        r2=xx+yy+zz
        repol=dsqrt(r2)
c......................................................
c    corrections for collisions electrons-qpol
c......................................................
        if (repol.lt.rcol) repol=rcol
        pepol_s=pepol_s-0.5d0*ch(j)/repol
      enddo

c......................................................
c     interaction with volume point charges
c......................................................
      pepol_v=0.0d0
      do j=nchs+1,nch
        xx=(x(1)-xpol(1,j))**2.0d0
        yy=(x(2)-xpol(2,j))**2.0d0
        zz=(x(3)-xpol(3,j))**2.0d0
        r2=xx+yy+zz
        repol=dsqrt(r2)
c......................................................
c    corrections for collisions electrons-qpol

c......................................................
        if (repol.lt.rcolv) repol=rcolv
        pepol_v=pepol_v-0.5d0*ch(j)/repol
      enddo

      return
      end
c......................................................
c    AVERAGES   subroutines
c......................................................

      subroutine pcm_init(iflag)

      use pcm, only: MCHS, MCHV, MSPHERE
      use pcm_cntrl, only: icall, ichpol, ipcm, ipcmprt, isurf
      use pcm_parms, only: ch, eps_solv, iscov, nch, nchs, nchs1, nchs2
      use pcm_parms, only: nchv, ncopcm, nesph, nscv, nvopcm, re, re2
      use pcm_parms, only: retk, surk, xe, xpol, ye, ze
      use pcm_averages, only: spcmsum, spcmcum, spcmcm2, vpcmsum, vpcmcum, vpcmcm2
      use pcm_averages, only: qopcm_sum, qopcm_cum, qopcm_cm2
      use pcm_averages, only: enfpcm_sum, enfpcm_cum, enfpcm_cm2

      implicit real*8(a-h,o-z)



      if(ipcm.eq.0) return
      spcmsum=0
      vpcmsum=0
      qopcm_sum=0

      do i=1,nchs
      enfpcm_sum(i)=0.0d0
      enddo

      if(iflag.gt.0) return
      spcmcum=0
      vpcmcum=0
      qopcm_cum=0

      spcmcm2=0
      vpcmcm2=0
      qopcm_cm2=0

      do i=1,nchs
      enfpcm_cum(i)=0.0d0
      enfpcm_cm2(i)=0.0d0
      enddo

      ncopcm=0

      return
      end

c-----------------------------------------------------------------------
      subroutine pcm_dump(iu)

      use pcm, only: MCHS, MCHV, MSPHERE
      use pcm_cntrl, only: ipcm
      use pcm_parms, only: ch, eps_solv, iscov, nch, nchs, nchs1, nchs2
      use pcm_parms, only: nchv, ncopcm, nesph, nscv, nvopcm, re, re2
      use pcm_parms, only: retk, surk, xe, xpol, ye, ze
      use pcm_averages, only: spcmcum, spcmcm2, vpcmcum, vpcmcm2
      use pcm_averages, only: qopcm_cum, qopcm_cm2
      use pcm_averages, only: enfpcm_cum, enfpcm_cm2

      implicit real*8(a-h,o-z)


      if(ipcm.eq.0) return
      write(iu) spcmcum,spcmcm2
      write(iu) vpcmcum,vpcmcm2
      write(iu) qopcm_cum,qopcm_cm2

      do i=1,nchs
      write(iu) enfpcm_cum(i),enfpcm_cm2(i)
      enddo

      return
      end
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
      subroutine pcm_rstrt(iu)

      use pcm, only: MCHS, MCHV, MSPHERE
      use pcm_cntrl, only: ipcm
      use pcm_parms, only: ch, eps_solv, iscov, nch, nchs, nchs1, nchs2
      use pcm_parms, only: nchv, ncopcm, nesph, nscv, nvopcm, re, re2
      use pcm_parms, only: retk, surk, xe, xpol, ye, ze
      use pcm_averages, only: spcmcum, spcmcm2, vpcmcum, vpcmcm2
      use pcm_averages, only: qopcm_cum, qopcm_cm2
      use pcm_averages, only: enfpcm_cum, enfpcm_cm2

      implicit real*8(a-h,o-z)


      if(ipcm.eq.0) return
      read(iu) spcmcum,spcmcm2
      read(iu) vpcmcum,vpcmcm2
      read(iu) qopcm_cum,qopcm_cm2

      do i=1,nchs
      read(iu) enfpcm_cum(i),enfpcm_cm2(i)
      enddo

      return
      end
c......................................................
c     end AVERAGES subroutines
c......................................................

      subroutine qpcm_charges(enfpcm_ave,enfpcm_err,qpol,sqpol2)

      use pcm, only: MCHS, MCHV, MSPHERE
      use pcm_parms, only: ch, eps_solv, iscov, nch, nchs, nchs1, nchs2
      use pcm_parms, only: nchv, ncopcm, nesph, nscv, nvopcm, re, re2
      use pcm_parms, only: retk, surk, xe, xpol, ye, ze
      use pcm_ah, only: ahca, bh
      use pcm_fdc, only: feps, qvol

      use pcm_inda, only: inda

      implicit real*8(a-h,o-z)



      dimension ch_new(MCHS),enfpcm_ave(*),enfpcm_err(*)
      dimension sch(MCHS),sch2(MCHS)
      dimension qpolsp(MSPHERE),sqpol2_sp(MSPHERE)

c...................................................................
c     only for ground states
c     We don't need to add qvol contribution to En because we have 
c     scaled contributions from e- out of the cavity
c...................................................................
c...................................................................
c     computes bh(k)
c...................................................................

      open (51,file='chsurf_new',form='formatted',status='unknown')
      open (54,file='field',form='formatted',status='unknown')
      rewind 51
      rewind 54

      do k=1,nchs
      bhe=-feps*surk*enfpcm_ave(k)
      bh(k)=bh(k)+bhe
      enddo

      qpol=0.0d0
      sqpol2=0.0d0
      do i=1,nchs
      ch_new(i)=0.0d0
      sch2(i)=0.0d0
	   do j=1,nchs
           ch_new(i)=ch_new(i)+ahca(i,j)*bh(j)
           sch2(i)=sch2(i)+(ahca(i,j)*feps*surk*enfpcm_err(i))**2.0d0
           enddo
      sch(i)=dsqrt(sch2(i)) 
      qpol=qpol+ch_new(i)
      sqpol2=sqpol2+sch2(i)
      enddo
c
c    computes qpol for each portion of sphere 
c
      do k=1,nesph
      qpolsp(k)=0.0D0                
      sqpol2_sp(k)=0.0D0
      enddo
      do i=1,nchs
      do k=1,nesph
      xx=(xpol(1,i)-xe(k))**2.0d0
      yy=(xpol(2,i)-ye(k))**2.0d0
      zz=(xpol(3,i)-ze(k))**2.0d0
      rr=xx+yy+zz
      rr=dsqrt(rr)
      write(54,1001)i,inda(i),k,re(k),rr
      if(inda(i).eq.k)then
      qpolsp(k)=qpolsp(k)+ch_new(i)
      sqpol2_sp(k)=sqpol2_sp(k)+sch2(i)
      endif
      enddo
      enddo
      write(6,*)
      write(6,*) 'pol. charges on spheres'
      do i=1,nesph
      sqpol_sp=dsqrt(sqpol2_sp(i))
      write(6,1200)i,xe(i),ye(i),ze(i),re(i),qpolsp(i),sqpol_sp
      enddo
      write(6,*)
c
c     write in 'chsurf_new' the new charges
c
      do i=1,nchs
        write(51,1000)(xpol(k,i),k=1,3),ch_new(i)
        write(54,1000)(xpol(k,i),k=1,3),enfpcm_ave(i),enfpcm_err(i),ch_new(i),sch(i)
      enddo
c..................................................

 1000 format(9F14.7)
 1001 format(3i7,5F14.7)
 1200 format(I4,9F14.7)

      return
      end

c................................................................
c     da modificare per stati eccitati
c................................................................
      subroutine qpcm_charges2(enfpcm_ave,enfpcm_err,qpol)

      use pcm, only: MCHS, MCHV, MSPHERE
      use vmc, only: MELEC, MORB, MBASIS, MDET, MCENT, MCTYPE, MCTYP3X
      use vmc, only: NSPLIN, nrad, MORDJ, MORDJ1, MMAT_DIM, MMAT_DIM2, MMAT_DIM20
      use vmc, only: radmax, delri
      use vmc, only: NEQSX, MTERMS
      use vmc, only: MCENT3, NCOEF, MEXCIT
      use pcm_parms, only: ch, nch, nchs, nchs1, nchs2
      use pcm_parms, only: nesph, re
      use pcm_parms, only: surk, xe, xpol, ye, ze

      use pcm_ameta, only: amdlg, eta
      use pcm_ah, only: ahca, bh
      use pcm_fdc, only: feps, qvol
      implicit real*8(a-h,o-z)





      dimension ch_new(MCHS)
      dimension env(MCHS)

c...................................................................
c     computes normal component electric field due to qvol
c...................................................................
      do k=1,nchs
      env(k)=0.0d0
      do l=nchs+1,nch
      xx=xpol(1,k)-xpol(1,l) 
      yy=xpol(2,k)-xpol(2,l) 
      zz=xpol(3,k)-xpol(3,l) 
      rr2=xx**2+yy**2+zz**2
      rr3=rr2**1.5d0
      cc1=xx*eta(1,k)
      cc2=yy*eta(2,k)
      cc3=zz*eta(3,k)
      cc=cc1+cc2+cc3 
      env(k)=env(k)+ch(l)*cc/rr3
      enddo
      enddo

c...................................................................
c     computes bh(k)
c...................................................................

      do k=1,nchs
      bhe=-feps*surk*(enfpcm_ave(k)+env(k))
      bh(k)=bh(k)+bhe
      enddo

      qpol=0.0d0
      do i=1,nchs
      ch_new(i)=0.0d0
	   do j=1,nchs
           ch_new(i)=ch_new(i)+ahca(i,j)*bh(j)
           enddo
      qpol=qpol+ch_new(i)
      enddo
      
c
c     write in 'chsurf_new' the new charges
c
      do i=1,nchs
      write(51,1000)(xpol(k,i),k=1,3),ch_new(i)
      write(54,1000)(xpol(k,i),k=1,3),enfpcm_ave(i)
      enddo
c..................................................

 1000 format(8F14.7)

      return
      end

c      
C     ***************************************************************
      subroutine qpcm_surface (npmax)
C     ***************************************************************
C     ***************************************************************
c     This subroutine computes the coordinates of point charges 
c     on the cavity surface
C     ***************************************************************

      use atom, only: cent, iwctype, ncent
      use spc, only: nsf, num
      use spc1, only: csf, qsf, rsf
      use spc2, only: nxyz, sfxyz, usf
      use pcm_parms, only: nch, nchs, nchs1, nchs2
      use pcm_parms, only: nesph, re
      use pcm_parms, only: surk, xe, xpol, ye, ze

      use pcm_ameta, only: amdlg, eta
      use pcm_inda, only: inda

      implicit real*8(a-h,o-z)



c
c....................................................................
c....................................................................
      dimension dist(2000000),ijpair(2000000,2)
      dimension inda1(5000),xpolt(3,5000)
      dimension indat(5000),itoro(5000)
      dimension amdlgt(5000),etat(3,5000)
c
c............................TRIAL VERSION-------------------
      character*80 row1,row2
      open (70,file='surf_info',form='formatted',status='unknown')
      rewind 70
      nxyz=0
c............................TRIAL VERSION-------------------

      open (53,file='surface',form='formatted',status='unknown')
      rewind 53

      pi=4.d0*datan(1.d0)
      nsf=nesph
      do i=1,nsf
      qsf(i,1)=xe(i)
      qsf(i,2)=ye(i)
      qsf(i,3)=ze(i)
      rsf(i)=re(i)
      enddo
      rmax=0.d0
      do i=1,nsf
      if (rsf(i).gt.rmax) rmax=rsf(i)
      enddo
      do i=1,nsf
      qi=dble(npmax)
      qi=qi*(rsf(i)/rmax)**2+0.5d0
      num(i)=qi
      enddo
      area=4.d0*pi*rmax**2/dble(npmax)
c------------------------TRIAL VERSION----------------------------
      read (70,*) icheck_sphere
      if (icheck_sphere.eq.0) then

c---------------prep and prep1 run only once----------------------
          call prep
          call prep1
c---------------storage of info from prep1------------------------
          rewind 70
          icheck_sphere=1
          write (70,*) icheck_sphere
          write (70,700) nsf,(num(i),i=1,nsf)
          do i=1,nsf
             do j=1,num(i)
                write (70,701) (csf(j,k,i),k=1,4)
             enddo
          enddo

      else

          if (ncent.ne.nsf) then
             write (6,*) 'FATAL: ncent not equal to nsf'
             write (6,*) 'not allowed in geometry optimization '
             stop
          endif
c---------------prep and prep1 bypassed---------------------------
          read  (70,700) nsf,(num(i),i=1,nsf)
          do i=1,nsf
             do j=1,num(i)
                read  (70,701) (csf(j,k,i),k=1,4)
             enddo
          enddo
c------------------uptade cavity_spheres (old file is lost!)------
       open (71,file='cavity_spheres',form='formatted',status='unknown')
          rewind 71
          read (71,702) row1
          do i=1,nsf+1
            read (71,702) row2
          enddo
          rewind 71
          write (71,702) row1
          do i=1,ncent
             write (71,703) (cent(j,i),j=1,3),rsf(i),iwctype(i)
c------------------uptade sphere center coordinates---------------
             do j=1,3
                qsf(i,j)=cent(j,i)
             enddo
          enddo
          write (71,702) row2
          close (71)

      endif
 700  format (20i5)
 701  format (8f15.8)
 702  format (a80)
 703  format (4f15.8,i3)
c------------------------TRIAL VERSION----------------------------
      do i=1,nsf
      do ii=1,num(i)
      kode=0
      x=qsf(i,1)+csf(ii,1,i)
      y=qsf(i,2)+csf(ii,2,i)
      z=qsf(i,3)+csf(ii,3,i)
      w=csf(ii,4,i)
      do j=1,nsf
      if (j.ne.i) then
      xx=x-qsf(j,1)
      yy=y-qsf(j,2)
      zz=z-qsf(j,3)
      rr=dsqrt(xx**2+yy**2+zz**2)
      if (rr.lt.rsf(j)) then
      kode=1
      goto 1
      endif
      endif
      enddo
 1    continue
      if (kode.eq.0) then
      nxyz=nxyz+1
      sfxyz(nxyz,1)=x
      sfxyz(nxyz,2)=y
      sfxyz(nxyz,3)=z
      sfxyz(nxyz,4)=w
      ux=csf(ii,1,i)
      uy=csf(ii,2,i)
      uz=csf(ii,3,i)
      uu=dsqrt(ux**2+uy**2+uz**2)
      usf(nxyz,1)=ux/uu
      usf(nxyz,2)=uy/uu
      usf(nxyz,3)=uz/uu
      inda1(nxyz)=i
      itoro(nxyz)=0
      endif
      enddo
      enddo
      deltaij=rmax/50.d0
      do i=1,1000
      dist(i)=0.d0
      enddo
      imax=0
      ipair=0
      cutoff=rmax/dsqrt(dble(npmax))*2.7d0
      do i=2,nxyz
      do j=1,i-1
      rij=dsqrt((sfxyz(i,1)-sfxyz(j,1))**2+(sfxyz(i,2)-sfxyz(j,2))**2+
     +(sfxyz(i,3)-sfxyz(j,3))**2)
      if (rij.lt.cutoff) then
      ipair=ipair+1
      ijpair(ipair,1)=i
      ijpair(ipair,2)=j
      endif
      rd=rij/deltaij
      ird=rd+1
      if (imax.lt.ird) imax=ird
      dist(ird)=dist(ird)+1.d0
      enddo
      enddo
      do ij=1,ipair
      i1=ijpair(ij,1)
      i2=ijpair(ij,2)
      x=(sfxyz(i1,1)+sfxyz(i2,1))/2.d0
      y=(sfxyz(i1,2)+sfxyz(i2,2))/2.d0
      z=(sfxyz(i1,3)+sfxyz(i2,3))/2.d0
      w=(sfxyz(i1,4)+sfxyz(i2,4))/2.d0
      sfxyz(i1,1)=x
      sfxyz(i1,2)=y
      sfxyz(i1,3)=z
      sfxyz(i1,4)=w
      sfxyz(i2,1)=x
      sfxyz(i2,2)=y
      sfxyz(i2,3)=z
      sfxyz(i2,4)=w
      ux=usf(i1,1)+usf(i2,1)
      uy=usf(i1,2)+usf(i2,2)
      uz=usf(i1,3)+usf(i2,3)
      uu=dsqrt(ux**2+uy**2+uz**2)
      prdm=rannyu(0)
      if (prdm.le.0.5d0) then
      usf(i1,1)=ux/uu
      usf(i1,2)=uy/uu
      usf(i1,3)=uz/uu
      usf(i2,1)=0.d0
      usf(i2,2)=0.d0
      usf(i2,3)=0.d0
      itoro(i1)=1
      else
      usf(i2,1)=ux/uu
      usf(i2,2)=uy/uu
      usf(i2,3)=uz/uu
      usf(i1,1)=0.d0
      usf(i1,2)=0.d0
      usf(i1,3)=0.d0
      itoro(i2)=1
      endif
      enddo
      icount=0
      icount2=0
      do i=1,nxyz
      ux=usf(i,1)
      uy=usf(i,2)
      uz=usf(i,3)
      uu=dsqrt(ux**2+uy**2+uz**2)
      ith=itoro(i)
      if (uu.gt.1.d-5.and.ith.eq.0) then
      icount=icount+1
      xpol(1,icount)=sfxyz(i,1)
      xpol(2,icount)=sfxyz(i,2)
      xpol(3,icount)=sfxyz(i,3)
      amdlg(icount)=sfxyz(i,4)
      eta(1,icount)=usf(i,1)
      eta(2,icount)=usf(i,2)
      eta(3,icount)=usf(i,3)
      inda(icount)=inda1(i)
      endif
      if (uu.gt.1.d-5.and.ith.eq.1) then
      icount2=icount2+1
      xpolt(1,icount2)=sfxyz(i,1)
      xpolt(2,icount2)=sfxyz(i,2)
      xpolt(3,icount2)=sfxyz(i,3)
      amdlgt(icount2)=sfxyz(i,4)
      etat(1,icount2)=usf(i,1)
      etat(2,icount2)=usf(i,2)
      etat(3,icount2)=usf(i,3)
      indat(icount2)=inda1(i)
      endif
      enddo
      nchs1=icount
      nchs2=icount2
      nchs=nchs1+nchs2
      do i=icount+1,nchs
      iold=i-icount
      xpol(1,i)= xpolt(1,iold)
      xpol(2,i)= xpolt(2,iold)
      xpol(3,i)= xpolt(3,iold)
      amdlg(i)= amdlgt(iold)
      eta(1,i)= etat(1,iold)
      eta(2,i)= etat(2,iold)
      eta(3,i)= etat(3,iold)
      inda(i)=indat(iold)
      enddo
      surk=area
c..................................................
c    print results in the file 'surface'
c..................................................
      write (53,*) '  area tess.       nchs1    nchs2   nchs'
      write (53,*) area,icount,icount2,nchs,' --> (nch) '
      write (53,*) 'positions, madelung, normal versors on the surface'
      do j=1,nchs
      write (53,100) inda(j),(xpol(i,j),i=1,3),amdlg(j),(eta(i,j),i=1,3)
      enddo
c..................................................
 100  format (I4,8f10.5)
c..................................................
      return
      end


      subroutine prep
      use spc, only: nsf, num
      use spc1, only: csf, rsf
      implicit real*8(a-h,o-z)


      dimension x(1000),y(1000),z(1000)
      dimension xm(1000),ym(1000),zm(1000)
      dimension el1(1000),el2(1000)
      dimension elm1(1000),elm2(1000)
      maxit=100
      epsilon=0.d0
      pi2=8.d0*datan(1.d0)
      do isf=1,nsf
      np=num(isf)
      r0=rsf(isf)
      np5=np
      nf=1
      cutoff=2.d0*r0/dble(np5+1)
      dd=0.d0
      if (cutoff.gt.dd) cutoff=dd
      z1=-r0
      z2=r0
      ww=r0/dble(np*np)
      imax=1000
      dz=(z2-z1)/dble(imax)
      zz=z1-dz
      ncyc=0
      energy=1.d10
      pigreco=pi2*0.5d0
      pi25=pi2/dble(nf)
      icount=0
      do ncyc=1,maxit
      do ncyc1=1,maxit
      icount=icount+1
      pol=0.d0
      do i=1,np5
300   rrand=rannyu(0)
      a1=rrand
c     theta=pigreco*a1
c     qm=r0*dcos(theta)
      cos_theta=-1.d0+2.d0*a1
      qm=r0*cos_theta
      qf=pi25*rannyu(0)
      el1(i)=qm
      el2(i)=qf
c     f1=r0*dsin(theta)
      sin_theta=dsqrt(dabs(1.d0-cos_theta**2))
      f1=r0*sin_theta
      xi=f1*dcos(qf)
      yi=f1*dsin(qf)
      zi=qm
      if (i.gt.1) then
      do j=1,i-1
      rij=dsqrt((x(j)-xi)**2+(y(j)-yi)**2+(z(j)-zi)**2)
      if (rij.lt.cutoff) goto 300
      rij=dsqrt((x(j+np5)-xi)**2+(y(j+np5)-yi)**2+(z(j+np5)-zi)**2)
      if (rij.lt.cutoff) goto 300
      rij=dsqrt((x(j+np5*(nf-1))-xi)**2+
     +(y(j+np5*(nf-1))-yi)**2+(z(j+np5*(nf-1))-zi)**2)
      if (rij.lt.cutoff) goto 300
      enddo
      endif
      x(i)=xi
      y(i)=yi
      z(i)=zi
      do k=1,nf-1
      qf=qf+pi25
      x(i+k*np5)=f1*dcos(qf)
      y(i+k*np5)=f1*dsin(qf)
      z(i+k*np5)=qm
      pol=pol+epsilon*qm
      enddo
      pol=pol+epsilon*z(i)
      enddo
      call rep (x,y,z,np,e0)
      e0=e0*ww+pol
      if (e0.lt.energy) then
      energy=e0
      do i=1,np
      elm1(i)=el1(i)
      elm2(i)=el2(i)
      xm(i)=x(i)
      ym(i)=y(i)
      zm(i)=z(i)
      enddo
      endif
      enddo
      enddo
      do i=1,np
      csf(i,1,isf)=xm(i)
      csf(i,2,isf)=ym(i)
      csf(i,3,isf)=zm(i)
      enddo
      enddo
      return
      end

      subroutine rep (x,y,z,n,e)
      implicit real*8 (a-h,o-z)
      dimension x(1),y(1),z(1)
      e=0.d0
      do i=1,n-1
      do j=i+1,n
      rij=(x(i)-x(j))**2+(y(i)-y(j))**2+(z(i)-z(j))**2
      rij=dsqrt(rij)
      if (rij.eq.0.d0) goto 10
      e=e+1.d0/rij
      enddo
      enddo
      return
 10   e=1.d10
      return
      end

      subroutine prep1
      use pcm, only: MCHS
      use spc, only: nsf, num
      use spc1, only: csf, rsf
      implicit real*8(a-h,o-z)


      dimension x(1000),y(1000),z(1000)
      dimension xm(1000),ym(1000),zm(1000)
      dimension el1(1000),el2(1000)
      dimension elm1(1000),elm2(1000)
      maxit=300
      epsilon=0.d0
      pi2=8.d0*datan(1.d0)
      do isf=1,nsf
      np=num(isf)
      r0=rsf(isf)
      np5=np
      nf=1
      energy=1.d8
      do i=1,np
      xm(i)=csf(i,1,isf)
      ym(i)=csf(i,2,isf)
      zm(i)=csf(i,3,isf)
      enddo
      r0=dsqrt(xm(1)**2+ym(1)**2+zm(1)**2)
      cutoff=2.d0*r0/dble(np5+1)
      dd=0.d0
      if (cutoff.gt.dd) cutoff=dd
      z1=-r0
      z2=r0
      ww=r0/dble(np*np)
c----------------------------------------------------------
      call rep  (xm,ym,zm,np,erep)
      do i=1,np
      sumz=sumz+zm(i)
      enddo
      energy=erep*ww+epsilon*sumz
      amadelung=erep*ww
c----------------------------------------------------------
      imax=1000
      pi2=8.d0*datan(1.d0)
      dz=(z2-z1)/dble(imax)
      zz=z1-dz
      ncyc=0
      pigreco=pi2*0.5d0
      pi25=pi2/dble(nf)
      ff1=pigreco/700.d0
      icount=0
      do ncyc=1,maxit
      do ncyc1=1,maxit
      icount=icount+1
      pol=0.d0
      do i=1,np5
300   continue
      f1=ff1*(rannyu(0)-0.5d0)
      f2=ff1*(rannyu(0)-0.5d0)
      c1=dcos(f1)
      s1=dsin(f1)
      c2=dcos(f2)
      s2=dsin(f2)
      xx=xm(i)*c1-ym(i)*s1
      yy=xm(i)*s1+ym(i)*c1
      zz=zm(i)
      xi=xx*c2+zz*s2
      yi=yy
      zi=-xx*s2+zz*c2
      x(i)=xi
      y(i)=yi
      z(i)=zi
      ro=dsqrt(xi**2+yi**2)
      qf=datan(yi/xi)
      do k=1,nf-1
      qf=qf+pi25
      x(i+k*np5)=ro*dcos(qf)
      y(i+k*np5)=ro*dsin(qf)
      z(i+k*np5)=zi
      pol=pol+epsilon*zi
      enddo
      pol=pol+epsilon*z(i)
      enddo
      call rep (x,y,z,np,erep)
      e0=erep*ww+pol
      if (e0.lt.energy) then
      energy=e0
      amadelung=erep*ww
      do i=1,np
      xm(i)=x(i)
      ym(i)=y(i)
      zm(i)=z(i)
      enddo
      endif
      enddo
c     write (3,*) isf,ncyc,amadelung
      enddo
      do i=1,np
      amdg=0.d0
      do j=1,np
      if (j.ne.i) then
      rij=dsqrt((xm(i)-xm(j))**2+(ym(i)-ym(j))**2+(zm(i)-zm(j))**2)
      amdg=amdg+1.d0/rij
      endif
      enddo
      amdg=amdg*r0/dble(2*np)
c     write (6,100) xm(i),ym(i),zm(i)
      csf(i,1,isf)=xm(i)
      csf(i,2,isf)=ym(i)
      csf(i,3,isf)=zm(i)
      csf(i,4,isf)=amdg
      enddo
c     write (6,*) 'energy',energy
      enddo
      return
      end

c........................................................
c.............IMPORTANT!!!!!!!!!!!!!!!!!!!!!!!!!!........
c...........REVISION DONE ON OCTOBER 2013 ...............
c.............IMPORTANT!!!!!!!!!!!!!!!!!!!!!!!!!!........
c........................................................
c
c     subroutine qpcm_matinv(a,nsub,det)
c     implicit real*8 (a-h,o-z)
c     include 'pcm.h'
c
c routine to calculate inverse and determinant of matrix a
c assumed to be dimensioned a(nsub,nsub).
c the matrix a is replaced by its inverse.
c
c     dimension a(nsub,nsub)
c     dimension ipvt(nsub),work(nsub),dt(2)
c     dimension ipvt(1000000),work(1000000),dt(2)
c
c     call dgefa(a,nsub,nsub,ipvt,info)
c     if(info.gt.0) then
c       write(6,'(''qpcm_MATINV: u(k,k)=0 with k= '',i5)') info
c       call fatal_error('MATINV: info ne 0 in dgefa')
c     endif
c     call dgedi(a,nsub,nsub,ipvt,dt,work,11)
c
c     det = dt(1)*10.0**dt(2)
c     return
c
c     end
c
c........................................................
c........................................................
c...................NEW!!!!!!!!!!!!!!!!!.................
c........................................................
c........................................................

      subroutine qpcm_matinv(a,nsub,determinant)
      use pcm, only: MCHS
      implicit real*8 (a-h,o-z)

c routine to calculate inverse and determinant of matrix a
c assumed to be dimensioned a(nsub,nsub).
c the matrix a is replaced by its inverse.

      dimension a(nsub,nsub)
      dimension ipvt(MCHS),work(MCHS),det(2)

      call dgetrf(nsub,nsub,a,nsub,ipvt,info)
      if(info.gt.0) then
        write(6,'(''MATINV: u(k,k)=0 with k= '',i5)') info
        call fatal_error('MATINV: info ne 0 in dgetrf')
      endif

      det(1) = 1.0d0
      det(2) = 0.0d0
      ten = 10.0d0
      do 50 i = 1, nsub
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
   50 continue
   60 continue

      determinant = det(1)*10.0**det(2)

      call dgetri(nsub,a,nsub,ipvt,work,MCHS,info)

      return

      end

      subroutine sigma_R (efield,q_strech)
C     ***************************************************************
C     march/2014: compute surface charges for a stretched cavity
C     ***************************************************************
      use pcm, only: MCHS
      use atom, only: znuc, cent, pecent, iwctype, nctype, ncent
      use pcm_parms, only: ch, eps_solv, iscov, nch, nchs, nchs1, nchs2
      use pcm_parms, only: nchv, ncopcm, nesph, nscv, nvopcm, re, re2
      use pcm_parms, only: retk, surk, xe, xpol, ye, ze

      use pcm_ameta, only: amdlg, eta
      use pcm_fdc, only: feps, fs, qfree, qvol, rcol, rcolt, rcolv
      use pcm_inda, only: inda

      implicit real*8(a-h,o-z)



C     
C    

      dimension ah_vec(MCHS*MCHS),ah(MCHS,MCHS)
      dimension q_strech(MCHS),efield(MCHS)

c
      DATA PI/3.1415927D0/

      do k=1,nchs
        cork=1.0d0/(2.0d0*amdlg(k))
         do l=1,nchs
          if (l.eq.k)then
           ah(k,l)=1.0d0+0.5d0*(1.d0-eps_solv)/eps_solv
          else
           xx=xpol(1,k)-xpol(1,l) 
           yy=xpol(2,k)-xpol(2,l) 
           zz=xpol(3,k)-xpol(3,l) 
           s1=xx*eta(1,k)
           s2=yy*eta(2,k)
           s3=zz*eta(3,k)
           ss=s1+s2+s3
           xx2=xx**2.0d0
           yy2=yy**2.0d0
           zz2=zz**2.0d0
           rkl2=xx2+yy2+zz2
           rkl3=rkl2**1.5d0
           corr=cork
           if(inda(l).ne.inda(k))corr=1.0d0
           ah(k,l)=feps*surk*ss*corr/rkl3
          endif
        enddo
      enddo
      do i=1,nchs
        do j=1,nchs
          ij=(j-1)*nchs+i
          ah_vec(ij)=ah(i,j)
        enddo
      enddo
      call qpcm_matinv(ah_vec,nchs,det)
      do i=1,nchs
        do j=1,nchs
          ij=(j-1)*nchs+i
          ah(i,j)=ah_vec(ij)
        enddo
      enddo

c     surface charges recomputed for the new stretched cavity

      do i=1,nchs
         q_strech(i)=0.d0
         do j=1,nchs 
           q_strech(i)=q_strech(i)-feps*surk*ah(i,j)*efield(j)
         enddo
      enddo

      return
      END

