      module basis_fnse_vg_mod
      contains
      subroutine basis_fnse_vg(k,rvec_en,r_en)
c Written by Claudia Filippi by modifying basis_fns
c routine to calculate basis functions and derivatives for electron k
c vg -> value,gradient

      use numbas_mod, only: MRWF
      use atom, only: iwctype, ncent, ncent_tot
      use ghostatom, only: nghostcent
      use numbas, only: iwrwf, nrbas, numr

      use phifun, only: dphin, n0_nbasis
      use phifun, only: phin
      use wfsec, only: iwf
      use basis, only: n1s, n2p, n3dzr, n3dx2, n3dxy, n3dxz, n3dyz
      use basis, only: n4fxxx, n4fyyy, n4fzzz, n4fxxy, n4fxxz, n4fyyx, n4fyyz
      use basis, only: n4fzzx, n4fzzy, n4fxyz
      use const, only: nelec

      use precision_kinds, only: dp

      use basis_fns_vgl_mod, only: n0_inc
      use splfit_mod, only: splfit
      implicit none

      integer :: i, iabs, ic, ider, irb
      integer :: j, ju, k, l
      integer :: ll
      real(dp) :: cd1, cd2, cf, cf2, cf3
      real(dp) :: cp, cs, r, r2
      real(dp) :: ri, ri2, ri3, rk
      real(dp) :: rt3, rt3b2, term, wf1
      real(dp) :: x2y, x2y2, x2z, xvec
      real(dp) :: xvec3, xy, xyz, xz
      real(dp) :: y2x, y2z, yvec, yvec3
      real(dp) :: yz, z2x, z2y, zr
      real(dp) :: zvec, zvec3
      real(dp), dimension(4, MRWF) :: wfv
      real(dp), dimension(3) :: xc
      real(dp), dimension(3, nelec, ncent_tot) :: rvec_en
      real(dp), dimension(nelec, ncent_tot) :: r_en
      real(dp), parameter :: one = 1.d0
      real(dp), parameter :: two = 2.d0
      real(dp), parameter :: three = 3.d0
      real(dp), parameter :: four = 4.d0
      real(dp), parameter :: half = 0.5d0






      data rt3,rt3b2/1.732050808d0,0.866025404d0/
c cs=1/sqrt(4*pi), cp=sqrt(3/(4*pi)), cd1=sqrt(5/(4*pi)), cd2=sqrt(15/(4*pi))
      data cs,cp,cd1,cd2/0.28209479d0,0.48860251d0,
     &0.63078313d0,1.0925484d0/
c cf=sqrt(7/(4*pi)),cf2=cf*sqrt(5),cf3=cf*sqrt(15)
      data cf,cf2,cf3/0.746352665180231d0,1.66889529453114d0,
     &2.89061144264055d0/

      l=0
      n0_nbasis(k)=0

c loop through centers

      do ic=1,ncent+nghostcent
      ll=0

      i=iwctype(ic)

c get distance to center

      xc(1)=rvec_en(1,k,ic)
      xc(2)=rvec_en(2,k,ic)
      xc(3)=rvec_en(3,k,ic)
      r=r_en(k,ic)
      r2=r*r
      ri=one/r
      ri2=ri**2
      ri3=ri2*ri


c analytical orbital

      if(numr.gt.0) then


c numerical orbitals

      ider=1
      rk=r
      do irb=1,nrbas(i)
      call splfit(rk,irb,i,iwf,wfv(1,irb),ider)
      enddo

c s states

      if (iabs(n1s(i)).lt.1) goto 620
      ju=iabs(n1s(i))
      do j=1,ju
      l=l+1
      ll=ll+1
      irb=iwrwf(ll,i)
      phin(l,k)=cs*wfv(1,irb)
      dphin(1,l,k)=cs*xc(1)*ri*wfv(2,irb)
      dphin(2,l,k)=cs*xc(2)*ri*wfv(2,irb)
      dphin(3,l,k)=cs*xc(3)*ri*wfv(2,irb)
      call n0_inc(l,k,ic)
      enddo
  620 continue

c p states

      if (iabs(n2p(1,i)).lt.1) goto 640
      ju=iabs(n2p(1,i))
      do j=1,ju
      l=l+1
      ll=ll+1
      irb=iwrwf(ll,i)
      xvec=xc(1)*ri
      yvec=xc(2)*ri
      zvec=xc(3)*ri
      wf1=wfv(2,irb)-ri*wfv(1,irb)
      phin(l,k)=cp*xvec*wfv(1,irb)
      dphin(1,l,k)=cp*xvec*xvec*wf1+cp*ri*wfv(1,irb)
      dphin(2,l,k)=cp*xvec*yvec*wf1
      dphin(3,l,k)=cp*xvec*zvec*wf1
      call n0_inc(l,k,ic)
      enddo
  640 continue

      if (iabs(n2p(2,i)).lt.1) goto 660
      ju=iabs(n2p(2,i))
      do j=1,ju
      l=l+1
      ll=ll+1
      irb=iwrwf(ll,i)
      xvec=xc(1)*ri
      yvec=xc(2)*ri
      zvec=xc(3)*ri
      wf1=wfv(2,irb)-ri*wfv(1,irb)
      phin(l,k)=cp*yvec*wfv(1,irb)
      dphin(1,l,k)=cp*yvec*xvec*wf1
      dphin(2,l,k)=cp*yvec*yvec*wf1+cp*ri*wfv(1,irb)
      dphin(3,l,k)=cp*yvec*zvec*wf1
      call n0_inc(l,k,ic)
      enddo
  660 continue

      if (iabs(n2p(3,i)).lt.1) goto 680
      ju=iabs(n2p(3,i))
      do j=1,ju
      l=l+1
      ll=ll+1
      irb=iwrwf(ll,i)
      xvec=xc(1)*ri
      yvec=xc(2)*ri
      zvec=xc(3)*ri
      wf1=wfv(2,irb)-ri*wfv(1,irb)
      phin(l,k)=cp*zvec*wfv(1,irb)
      dphin(1,l,k)=cp*zvec*xvec*wf1
      dphin(2,l,k)=cp*zvec*yvec*wf1
      dphin(3,l,k)=cp*zvec*zvec*wf1+cp*ri*wfv(1,irb)
      call n0_inc(l,k,ic)
      enddo
  680 continue

c d states

      if (iabs(n3dzr(i)).lt.1) goto 700
      ju=iabs(n3dzr(i))
      do j=1,ju
      l=l+1
      ll=ll+1
      irb=iwrwf(ll,i)
      zr=half*(three*xc(3)**2-r2)*ri2
      xvec=xc(1)*ri
      yvec=xc(2)*ri
      zvec=xc(3)*ri
      wf1=wfv(2,irb)-two*ri*wfv(1,irb)
      phin(l,k)=cd1*zr*wfv(1,irb)
      dphin(1,l,k)=cd1*xvec*(zr*wf1-wfv(1,irb)*ri)
      dphin(2,l,k)=cd1*yvec*(zr*wf1-wfv(1,irb)*ri)
      dphin(3,l,k)=cd1*zvec*(zr*wf1+two*wfv(1,irb)*ri)
      call n0_inc(l,k,ic)
      enddo

  700 continue
      if (iabs(n3dx2(i)).lt.1) goto 720
      ju=iabs(n3dx2(i))
      do j=1,ju
      l=l+1
      ll=ll+1
      irb=iwrwf(ll,i)
      x2y2=half*(xc(1)**2-xc(2)**2)*ri2
      xvec=xc(1)*ri
      yvec=xc(2)*ri
      zvec=xc(3)*ri
      wf1=wfv(2,irb)-two*ri*wfv(1,irb)
      phin(l,k)=cd2*x2y2*wfv(1,irb)
      dphin(1,l,k)=cd2*xvec*(x2y2*wf1+wfv(1,irb)*ri)
      dphin(2,l,k)=cd2*yvec*(x2y2*wf1-wfv(1,irb)*ri)
      dphin(3,l,k)=cd2*zvec* x2y2*wf1
      call n0_inc(l,k,ic)
      enddo

  720 continue
      if (iabs(n3dxy(i)).lt.1) goto 740
      ju=iabs(n3dxy(i))
      do j=1,ju
      l=l+1
      ll=ll+1
      irb=iwrwf(ll,i)
      xy=xc(1)*xc(2)*ri2
      xvec=xc(1)*ri
      yvec=xc(2)*ri
      zvec=xc(3)*ri
      wf1=wfv(2,irb)-two*ri*wfv(1,irb)
      phin(l,k)=cd2*xy*wfv(1,irb)
      dphin(1,l,k)=cd2*(xvec*xy*wf1+yvec*wfv(1,irb)*ri)
      dphin(2,l,k)=cd2*(yvec*xy*wf1+xvec*wfv(1,irb)*ri)
      dphin(3,l,k)=cd2* zvec*xy*wf1
      call n0_inc(l,k,ic)
      enddo

  740 continue
      if (iabs(n3dxz(i)).lt.1) goto 760
      ju=iabs(n3dxz(i))
      do j=1,ju
      l=l+1
      ll=ll+1
      irb=iwrwf(ll,i)
      xz=xc(1)*xc(3)*ri2
      xvec=xc(1)*ri
      yvec=xc(2)*ri
      zvec=xc(3)*ri
      wf1=wfv(2,irb)-two*ri*wfv(1,irb)
      phin(l,k)=cd2*xz*wfv(1,irb)
      dphin(1,l,k)=cd2*(xvec*xz*wf1+zvec*wfv(1,irb)*ri)
      dphin(2,l,k)=cd2* yvec*xz*wf1
      dphin(3,l,k)=cd2*(zvec*xz*wf1+xvec*wfv(1,irb)*ri)
      call n0_inc(l,k,ic)
      enddo

  760 continue
      if (iabs(n3dyz(i)).lt.1) goto 780
      ju=iabs(n3dyz(i))
      do j=1,ju
      l=l+1
      ll=ll+1
      irb=iwrwf(ll,i)
      yz=xc(2)*xc(3)*ri2
      xvec=xc(1)*ri
      yvec=xc(2)*ri
      zvec=xc(3)*ri
      wf1=wfv(2,irb)-two*ri*wfv(1,irb)
      phin(l,k)=cd2*yz*wfv(1,irb)
      dphin(1,l,k)=cd2* xvec*yz*wf1
      dphin(2,l,k)=cd2*(yvec*yz*wf1+zvec*wfv(1,irb)*ri)
      dphin(3,l,k)=cd2*(zvec*yz*wf1+yvec*wfv(1,irb)*ri)
      call n0_inc(l,k,ic)
      enddo
  780 continue

      if (iabs(n4fxxx(i)).lt.1) goto 790
      ju=iabs(n4fxxx(i))
      do j=1,ju
      l=l+1
      ll=ll+1
      irb=iwrwf(ll,i)
      xvec=xc(1)*ri
      yvec=xc(2)*ri
      zvec=xc(3)*ri
      xvec3=xvec**3
      phin(l,k)=cf*xvec3*wfv(1,irb)
      term=cf*xvec3*(wfv(2,irb)-3*wfv(1,irb)*ri)
      dphin(1,l,k)=cf*3*xvec**2*wfv(1,irb)*ri+term*xvec
      dphin(2,l,k)=term*yvec
      dphin(3,l,k)=term*zvec
      call n0_inc(l,k,ic)
      enddo
 790  continue
      if (iabs(n4fyyy(i)).lt.1) goto 805
      ju=iabs(n4fyyy(i))
      do j=1,ju
      l=l+1
      ll=ll+1
      irb=iwrwf(ll,i)
      xvec=xc(1)*ri
      yvec=xc(2)*ri
      zvec=xc(3)*ri
      yvec3=yvec**3
      phin(l,k)=cf*yvec3*wfv(1,irb)
      term=cf*yvec3*(wfv(2,irb)-3*wfv(1,irb)*ri)
      dphin(1,l,k)=term*xvec
      dphin(2,l,k)=cf*3*yvec**2*wfv(1,irb)*ri+term*yvec
      dphin(3,l,k)=term*zvec
      call n0_inc(l,k,ic)
      enddo
 805  continue
      if (iabs(n4fzzz(i)).lt.1) goto 810
      ju=iabs(n4fzzz(i))
      do j=1,ju
      l=l+1
      ll=ll+1
      irb=iwrwf(ll,i)
      xvec=xc(1)*ri
      yvec=xc(2)*ri
      zvec=xc(3)*ri
      zvec3=zvec**3
      phin(l,k)=cf*zvec3*wfv(1,irb)
      term=cf*zvec3*(wfv(2,irb)-3*wfv(1,irb)*ri)
      dphin(1,l,k)=term*xvec
      dphin(2,l,k)=term*yvec
      dphin(3,l,k)=cf*3*zvec**2*wfv(1,irb)*ri+term*zvec
      call n0_inc(l,k,ic)
      enddo
 810  continue
      if (iabs(n4fxxy(i)).lt.1) goto 820
      ju=iabs(n4fxxy(i))
      do j=1,ju
      l=l+1
      ll=ll+1
      irb=iwrwf(ll,i)
      xvec=xc(1)*ri
      yvec=xc(2)*ri
      zvec=xc(3)*ri
      x2y=xvec**2*yvec
      phin(l,k)=cf2*x2y*wfv(1,irb)
      term=cf2*x2y*(wfv(2,irb)-3*wfv(1,irb)*ri)
      dphin(1,l,k)=cf2*2*xvec*yvec*wfv(1,irb)*ri+term*xvec
      dphin(2,l,k)=cf2*xvec*xvec*wfv(1,irb)*ri+term*yvec
      dphin(3,l,k)=term*zvec
      call n0_inc(l,k,ic)
      enddo
 820  continue
      if (iabs(n4fxxz(i)).lt.1) goto 825
      ju=iabs(n4fxxz(i))
      do j=1,ju
      l=l+1
      ll=ll+1
      irb=iwrwf(ll,i)
      xvec=xc(1)*ri
      yvec=xc(2)*ri
      zvec=xc(3)*ri
      x2z=xvec**2*zvec
      phin(l,k)=cf2*x2z*wfv(1,irb)
      term=cf2*x2z*(wfv(2,irb)-3*wfv(1,irb)*ri)
      dphin(1,l,k)=cf2*2*xvec*zvec*wfv(1,irb)*ri+term*xvec
      dphin(2,l,k)=term*yvec
      dphin(3,l,k)=cf2*xvec*xvec*wfv(1,irb)*ri+term*zvec
      call n0_inc(l,k,ic)
      enddo
 825  continue
      if (iabs(n4fyyx(i)).lt.1) goto 835
      ju=iabs(n4fyyx(i))
      do j=1,ju
      l=l+1
      ll=ll+1
      irb=iwrwf(ll,i)
      xvec=xc(1)*ri
      yvec=xc(2)*ri
      zvec=xc(3)*ri
      y2x=yvec**2*xvec
      phin(l,k)=cf2*y2x*wfv(1,irb)
      term=cf2*y2x*(wfv(2,irb)-3*wfv(1,irb)*ri)
      dphin(1,l,k)=cf2*yvec*yvec*wfv(1,irb)*ri+term*xvec
      dphin(2,l,k)=cf2*2*yvec*xvec*wfv(1,irb)*ri+term*yvec
      dphin(3,l,k)=term*zvec
      call n0_inc(l,k,ic)
      enddo
 835  continue
      if (iabs(n4fyyz(i)).lt.1) goto 845
      ju=iabs(n4fyyz(i))
      do j=1,ju
      l=l+1
      ll=ll+1
      irb=iwrwf(ll,i)
      xvec=xc(1)*ri
      yvec=xc(2)*ri
      zvec=xc(3)*ri
      y2z=yvec**2*zvec
      phin(l,k)=cf2*y2z*wfv(1,irb)
      term=cf2*y2z*(wfv(2,irb)-3*wfv(1,irb)*ri)
      dphin(1,l,k)=term*xvec
      dphin(2,l,k)=cf2*2*yvec*zvec*wfv(1,irb)*ri+term*yvec
      dphin(3,l,k)=cf2*yvec*yvec*wfv(1,irb)*ri+term*zvec
      call n0_inc(l,k,ic)
      enddo
 845  continue
      if (iabs(n4fzzx(i)).lt.1) goto 855
      ju=iabs(n4fzzx(i))
      do j=1,ju
      l=l+1
      ll=ll+1
      irb=iwrwf(ll,i)
      xvec=xc(1)*ri
      yvec=xc(2)*ri
      zvec=xc(3)*ri
      z2x=zvec**2*xvec
      phin(l,k)=cf2*z2x*wfv(1,irb)
      term=cf2*z2x*(wfv(2,irb)-3*wfv(1,irb)*ri)
      dphin(1,l,k)=cf2*zvec*zvec*wfv(1,irb)*ri+term*xvec
      dphin(2,l,k)=term*yvec
      dphin(3,l,k)=cf2*2*zvec*xvec*wfv(1,irb)*ri+term*zvec
      call n0_inc(l,k,ic)
      enddo
 855  continue
      if (iabs(n4fzzy(i)).lt.1) goto 865
      ju=iabs(n4fzzy(i))
      do j=1,ju
      l=l+1
      ll=ll+1
      irb=iwrwf(ll,i)
      xvec=xc(1)*ri
      yvec=xc(2)*ri
      zvec=xc(3)*ri
      z2y=zvec**2*yvec
      phin(l,k)=cf2*z2y*wfv(1,irb)
      term=cf2*z2y*(wfv(2,irb)-3*wfv(1,irb)*ri)
      dphin(1,l,k)=term*xvec
      dphin(2,l,k)=cf2*zvec*zvec*wfv(1,irb)*ri+term*yvec
      dphin(3,l,k)=cf2*2*zvec*yvec*wfv(1,irb)*ri+term*zvec
      call n0_inc(l,k,ic)
      enddo
 865  continue
      if (iabs(n4fxyz(i)).lt.1) goto 875
      ju=iabs(n4fxyz(i))
      do j=1,ju
      l=l+1
      ll=ll+1
      irb=iwrwf(ll,i)
      xvec=xc(1)*ri
      yvec=xc(2)*ri
      zvec=xc(3)*ri
      xyz=xvec*yvec*zvec
      phin(l,k)=cf3*xyz*wfv(1,irb)
      term=cf3*xyz*(wfv(2,irb)-3*wfv(1,irb)*ri)
      dphin(1,l,k)=cf3*yvec*zvec*wfv(1,irb)*ri+term*xvec
      dphin(2,l,k)=cf3*xvec*zvec*wfv(1,irb)*ri+term*yvec
      dphin(3,l,k)=cf3*xvec*yvec*wfv(1,irb)*ri+term*zvec
      call n0_inc(l,k,ic)
      enddo
 875  continue

      else
       stop
      endif

      enddo

      return
      end
      end module
