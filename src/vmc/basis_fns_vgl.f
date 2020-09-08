      subroutine basis_fns_vgl(x,rvec_en,r_en)
c Written by Cyrus Umrigar and Claudia Filippi, starting from Kevin Schmidt routine
c routine to calculate the values of the basis functions and their derivatives
c vgl -> value, gradient, laplacian

      use numbas_mod, only: MRWF
      use vmc, only: MELEC, MCENT
      use atom, only: iwctype, ncent
      use ghostatom, only: nghostcent
      use const, only: pi, nelec
      use numbas, only: iwrwf, nrbas, numr
      use numbas1, only: iwlbas, nbastyp
      use phifun, only: d2phin, d2phin_all, d3phin, dphin, n0_nbasis
      use phifun, only: phin
      use wfsec, only: iwf
      use force_analy, only: iforce_analy
      use basis, only: zex, n1s, n2s, n2p, n3s, n3p, n3dzr, n3dx2, n3dxy, n3dxz, n3dyz
      use basis, only: n4s, n4p

      implicit real*8(a-h,o-z)

      include 'pseudo.h'

      parameter (one=1.d0,two=2.d0,three=3.d0,four=4.d0)
      parameter (five=5.d0,six=6.d0,seven=7.d0,eight=8.d0)
      parameter (ten=10.d0,half=.5d0)
      parameter (twelve=12.d0)





      dimension x(3,*),rvec_en(3,MELEC,MCENT),r_en(MELEC,MCENT)
      dimension dy(3),ddy(3,3),dlapy(3)

      dimension wfv(4,MELEC,MRWF)
      dimension xc(3,MELEC),r(MELEC),ri(MELEC),ri2(MELEC)
      dimension ri3(MELEC),r2(MELEC)
      data rt3,rt3b2/1.732050808d0,0.866025404d0/
c cs=1/sqrt(4*pi), cp=sqrt(3/(4*pi)), cd1=sqrt(5/(4*pi)), cd2=sqrt(15/(4*pi))
      data cs,cp,cd1,cd2/0.28209479d0,0.48860251d0,
     &0.63078313d0,1.0925484d0/
c cf=sqrt(7/(4*pi)),cf2=cf*sqrt(5),cf3=cf*sqrt(15)
      data cf,cf2,cf3/0.746352665180231d0,1.66889529453114d0,
     &2.89061144264055d0/

      l=0

c loop through centers
      do 1 j=1,nelec
   1    n0_nbasis(j)=0

      do 990 ic=1,ncent+nghostcent
      ll=0

      i=iwctype(ic)

c get distance to center

      do 20 j=1,nelec
      xc(1,j)=rvec_en(1,j,ic)
      xc(2,j)=rvec_en(2,j,ic)
      xc(3,j)=rvec_en(3,j,ic)
      r(j)=r_en(j,ic)
      r2(j)=r(j)*r(j)
      ri(j)=one/r(j)
      ri2(j)=ri(j)*ri(j)
   20 ri3(j)=ri2(j)*ri(j)

c analytical orbital (Slater basis)

      if(numr.eq.0) then

      if (iabs(n1s(i)).lt.1) goto 50

c 1s states

      ju=n1s(i)
      do 40 j=1,ju
      l=l+1
      do 40 k=1,nelec
      ex=dexp(-zex(l,iwf)*r(k))
      phin(l,k)=ex
      x1=-zex(l,iwf)*ri(k)
      dphin(1,l,k)=x1*xc(1,k)*ex
      dphin(2,l,k)=x1*xc(2,k)*ex
      dphin(3,l,k)=x1*xc(3,k)*ex
      d2phin(l,k)=(two*x1+zex(l,iwf)**2)*ex
   40 continue
   50 continue

c 2s states
      if (iabs(n2s(i)).lt.1) goto 70
      ju=n2s(i)
      do 60 j=1,ju
      l=l+1
      do 60 k=1,nelec
      ex=dexp(-zex(l,iwf)*r(k))
      phin(l,k)=r(k)*ex
      x1=ri(k)-zex(l,iwf)
      x2=(zex(l,iwf)**2-ri2(k)-zex(l,iwf)*ri(k))*ri(k)
      dphin(1,l,k)=x1*xc(1,k)*ex
      dphin(2,l,k)=x1*xc(2,k)*ex
      dphin(3,l,k)=x1*xc(3,k)*ex
      d2phin(l,k)=(two*ri(k)-four*zex(l,iwf)+zex(l,iwf)**2*r(k))*ex
   60 continue
   70 continue

c 2p states
      if (iabs(n2p(1,i)).lt.1) goto 90
      ju=n2p(1,i)
      do 80 j=1,ju
      l=l+1
      do 80 k=1,nelec
      ex=dexp(-zex(l,iwf)*r(k))
      phin(l,k)=xc(1,k)*ex
      x1=-zex(l,iwf)*ri(k)
      x2=zex(l,iwf)*(zex(l,iwf)*ri2(k)+ri3(k))
      dphin(1,l,k)=ex*(one-zex(l,iwf)*xc(1,k)**2*ri(k))
      dphin(2,l,k)=x1*xc(2,k)*phin(l,k)
      dphin(3,l,k)=x1*xc(3,k)*phin(l,k)
      d2phin(l,k)=(-four*zex(l,iwf)*ri(k)+zex(l,iwf)**2)*phin(l,k)
   80 continue
   90 continue

      if (iabs(n2p(2,i)).lt.1) goto 110
      ju=n2p(2,i)
      do 100 j=1,ju
      l=l+1
      do 100 k=1,nelec
      ex=dexp(-zex(l,iwf)*r(k))
      phin(l,k)=xc(2,k)*ex
      x1=-zex(l,iwf)*ri(k)
      x2=zex(l,iwf)*(zex(l,iwf)*ri2(k)+ri3(k))
      dphin(1,l,k)=x1*xc(1,k)*phin(l,k)
      dphin(2,l,k)=ex*(one-zex(l,iwf)*xc(2,k)**2*ri(k))
      dphin(3,l,k)=x1*xc(3,k)*phin(l,k)
      d2phin(l,k)=(-four*zex(l,iwf)*ri(k)+zex(l,iwf)**2)*phin(l,k)
  100 continue
  110 continue

      if (iabs(n2p(3,i)).lt.1) goto 130
      ju=n2p(3,i)
      do 120 j=1,ju
      l=l+1
      do 120 k=1,nelec
      ex=dexp(-zex(l,iwf)*r(k))
      phin(l,k)=xc(3,k)*ex
      x1=-zex(l,iwf)*ri(k)
      x2=zex(l,iwf)*(zex(l,iwf)*ri2(k)+ri3(k))
      dphin(1,l,k)=x1*xc(1,k)*phin(l,k)
      dphin(2,l,k)=x1*xc(2,k)*phin(l,k)
      dphin(3,l,k)=ex*(one-zex(l,iwf)*xc(3,k)**2*ri(k))
      d2phin(l,k)=(-four*zex(l,iwf)*ri(k)+zex(l,iwf)**2)*phin(l,k)
  120 continue
  130 continue

c 3s state

      if (iabs(n3s(i)).lt.1) goto 155
      ju=n3s(i)
      do 150 j=1,ju
      l=l+1
      do 150 k=1,nelec
      ex=dexp(-zex(l,iwf)*r(k))
      phin(l,k)=r(k)**2*ex
      x1=(two-zex(l,iwf)*r(k))*ex
      dphin(1,l,k)=x1*xc(1,k)
      dphin(2,l,k)=x1*xc(2,k)
      dphin(3,l,k)=x1*xc(3,k)
      d2phin(l,k)=(six-six*zex(l,iwf)*r(k)+zex(l,iwf)**2*r(k)**2)*ex
  150 continue
  155 continue

c 3p states

      if (iabs(n3p(1,i)).lt.1) goto 170
      ju=n3p(1,i)
      do 160 j=1,ju
      l=l+1
      do 160 k=1,nelec
      ex=dexp(-zex(l,iwf)*r(k))
      phin(l,k)=xc(1,k)*r(k)*ex
      x1=(ri(k)-zex(l,iwf))*xc(1,k)*ex
      dphin(1,l,k)=(xc(1,k)**2*ri(k)-zex(l,iwf)*xc(1,k)**2+r(k))*ex
      dphin(2,l,k)=x1*xc(2,k)
      dphin(3,l,k)=x1*xc(3,k)
      d2phin(l,k)=(four*ri(k)-six*zex(l,iwf)+zex(l,iwf)**2*r(k))
     &                 *xc(1,k)*ex
  160 continue
  170 continue

      if (iabs(n3p(2,i)).lt.1) goto 190
      ju=n3p(2,i)
      do 180 j=1,ju
      l=l+1
      do 180 k=1,nelec
      ex=dexp(-zex(l,iwf)*r(k))
      phin(l,k)=xc(2,k)*r(k)*ex
      x1=(ri(k)-zex(l,iwf))*xc(2,k)*ex
      dphin(1,l,k)=x1*xc(1,k)
      dphin(2,l,k)=(xc(2,k)**2*ri(k)-zex(l,iwf)*xc(2,k)**2+r(k))*ex
      dphin(3,l,k)=x1*xc(3,k)
      d2phin(l,k)=(four*ri(k)-six*zex(l,iwf)+zex(l,iwf)**2*r(k))
     &                 *xc(2,k)*ex
  180 continue
  190 continue

      if (iabs(n3p(3,i)).lt.1) goto 210
      ju=n3p(3,i)
      do 200 j=1,ju
      l=l+1
      do 200 k=1,nelec
      ex=dexp(-zex(l,iwf)*r(k))
      phin(l,k)=xc(3,k)*r(k)*ex
      x1=(ri(k)-zex(l,iwf))*xc(3,k)*ex
      dphin(1,l,k)=x1*xc(1,k)
      dphin(2,l,k)=x1*xc(2,k)
      dphin(3,l,k)=(xc(3,k)**2*ri(k)-zex(l,iwf)*xc(3,k)**2+r(k))*ex
      d2phin(l,k)=(four*ri(k)-six*zex(l,iwf)+zex(l,iwf)**2*r(k))
     &                 *xc(3,k)*ex
  200 continue
  210 continue

c 3d states

      if (iabs(n3dzr(i)).lt.1) goto 230
      do 220 j=1,n3dzr(i)
      l=l+1
      do 220 k=1,nelec
      ex=half*dexp(-zex(l,iwf)*r(k))
      x2=three*xc(3,k)**2-r(k)**2
      phin(l,k)=x2*ex
      x1=(-zex(l,iwf)*ri(k)*x2-two)*ex
      dphin(1,l,k)=x1*xc(1,k)
      dphin(2,l,k)=x1*xc(2,k)
      dphin(3,l,k)=(four-zex(l,iwf)*ri(k)*x2)*xc(3,k)*ex
      d2phin(l,k)=(zex(l,iwf)**2-six*zex(l,iwf)*ri(k))*phin(l,k)
  220 continue
  230 continue
      if (iabs(n3dx2(i)).lt.1) goto 250
      do 240 j=1,n3dx2(i)
      l=l+1
      do 240 k=1,nelec
      ex=rt3b2*dexp(-zex(l,iwf)*r(k))
      x2=xc(1,k)**2-xc(2,k)**2
      phin(l,k)=x2*ex
      dphin(1,l,k)=(two-zex(l,iwf)*x2*ri(k))*xc(1,k)*ex
      dphin(2,l,k)=-(two+zex(l,iwf)*x2*ri(k))*xc(2,k)*ex
      dphin(3,l,k)=-zex(l,iwf)*x2*ri(k)*xc(3,k)*ex
      d2phin(l,k)=(zex(l,iwf)**2-six*zex(l,iwf)*ri(k))*phin(l,k)
  240 continue

  250 continue
      if (iabs(n3dxy(i)).lt.1) goto 270
      do 260 j=1,n3dxy(i)
      l=l+1
      do 260 k=1,nelec
      ex=rt3*dexp(-zex(l,iwf)*r(k))
      phin(l,k)=xc(1,k)*xc(2,k)*ex
      dphin(1,l,k)=(one-zex(l,iwf)*xc(1,k)**2*ri(k))*xc(2,k)*ex
      dphin(2,l,k)=(one-zex(l,iwf)*xc(2,k)**2*ri(k))*xc(1,k)*ex
      dphin(3,l,k)=-zex(l,iwf)*xc(3,k)*ri(k)*phin(l,k)
      d2phin(l,k)=(zex(l,iwf)**2-six*zex(l,iwf)*ri(k))*phin(l,k)
  260 continue

  270 continue
      if (iabs(n3dxz(i)).lt.1) goto 290
      do 280 j=1,n3dxz(i)
      l=l+1
      do 280 k=1,nelec
      ex=rt3*dexp(-zex(l,iwf)*r(k))
      phin(l,k)=xc(1,k)*xc(3,k)*ex
      dphin(1,l,k)=(one-zex(l,iwf)*xc(1,k)**2*ri(k))*xc(3,k)*ex
      dphin(2,l,k)=-zex(l,iwf)*xc(2,k)*ri(k)*phin(l,k)
      dphin(3,l,k)=(one-zex(l,iwf)*xc(3,k)**2*ri(k))*xc(1,k)*ex
      d2phin(l,k)=(zex(l,iwf)**2-six*zex(l,iwf)*ri(k))*phin(l,k)
  280 continue

  290 continue
      if (iabs(n3dyz(i)).lt.1) goto 310
      do 300 j=1,n3dyz(i)
      l=l+1
      do 300 k=1,nelec
      ex=rt3*dexp(-zex(l,iwf)*r(k))
      phin(l,k)=xc(2,k)*xc(3,k)*ex
      dphin(1,l,k)=-zex(l,iwf)*xc(1,k)*ri(k)*phin(l,k)
      dphin(2,l,k)=(one-zex(l,iwf)*xc(2,k)**2*ri(k))*xc(3,k)*ex
      dphin(3,l,k)=(one-zex(l,iwf)*xc(3,k)**2*ri(k))*xc(2,k)*ex
      d2phin(l,k)=(zex(l,iwf)**2-six*zex(l,iwf)*ri(k))*phin(l,k)
  300 continue
  310 continue

c 4s state

      if (iabs(n4s(i)).lt.1) goto 330
      do 320 j=1,n4s(i)
      l=l+1
      do 320 k=1,nelec
      ex=dexp(-zex(l,iwf)*r(k))
      phin(l,k)=r(k)**3*ex
      x1=(three*r(k)-zex(l,iwf)*r(k)**2)*ex
      dphin(1,l,k)=x1*xc(1,k)
      dphin(2,l,k)=x1*xc(2,k)
      dphin(3,l,k)=x1*xc(3,k)
      d2phin(l,k)=(zex(l,iwf)**2*r(k)**3-eight*zex(l,iwf)*r(k)**2
     &             +twelve*r(k))*ex
  320 continue
  330 continue

c 4p states

      if (iabs(n4p(1,i)).lt.1) goto 350
      do 340 j=1,n4p(1,i)
      l=l+1
      do 340 k=1,nelec
      ex=dexp(-zex(l,iwf)*r(k))
      phin(l,k)=xc(1,k)*r(k)**2*ex
      x1=(two-zex(l,iwf)*r(k))*xc(1,k)*ex
      dphin(1,l,k)=(-zex(l,iwf)*r(k)*xc(1,k)**2+two*xc(1,k)**2+r(k)**2)
     &*ex
      dphin(2,l,k)=x1*xc(2,k)
      dphin(3,l,k)=x1*xc(3,k)
      d2phin(l,k)=(zex(l,iwf)**2*r(k)**2-eight*zex(l,iwf)*r(k)+ten)
     &*xc(1,k)*ex
  340 continue

  350 continue
      if (iabs(n4p(2,i)).lt.1) goto 370
      do 360 j=1,n4p(2,i)
      l=l+1
      do 360 k=1,nelec
      ex=dexp(-zex(l,iwf)*r(k))
      phin(l,k)=xc(2,k)*r(k)**2*ex
      x1=(two-zex(l,iwf)*r(k))*xc(2,k)*ex
      dphin(1,l,k)=x1*xc(1,k)
      dphin(2,l,k)=(-zex(l,iwf)*r(k)*xc(2,k)**2+two*xc(2,k)**2+r(k)**2)
     &*ex
      dphin(3,l,k)=x1*xc(3,k)
      d2phin(l,k)=(zex(l,iwf)**2*r(k)**2-eight*zex(l,iwf)*r(k)+ten)
     &*xc(2,k)*ex
  360 continue

  370 continue
      if (iabs(n4p(3,i)).lt.1) goto 390
      do 380 j=1,n4p(3,i)
      l=l+1
      do 380 k=1,nelec
      ex=dexp(-zex(l,iwf)*r(k))
      phin(l,k)=xc(3,k)*r(k)**2*ex
      x1=(two-zex(l,iwf)*r(k))*xc(3,k)*ex
      dphin(1,l,k)=x1*xc(1,k)
      dphin(2,l,k)=x1*xc(2,k)
      dphin(3,l,k)=(-zex(l,iwf)*r(k)*xc(3,k)**2+two*xc(3,k)**2+r(k)**2)
     &*ex
      d2phin(l,k)=(zex(l,iwf)**2*r(k)**2-eight*zex(l,iwf)*r(k)+ten)
     &*xc(3,k)*ex
  380 continue
  390 continue

      else

c numerical orbitals

      ider=1
      do 600 k=1,nelec
      rk=r(k)
      do 600 irb=1,nrbas(i)
      call splfit(rk,irb,i,iwf,wfv(1,k,irb),ider)
  600 continue

      k0=0
      l0=l
      do 950 k=1,nelec

      l=l0
      ll=0

      iwlbas0=0
      do 800 j=1,nbastyp(i)
      l=l+1
      ll=ll+1
      irb=iwrwf(ll,i)
      if(iwlbas(ll,i).ne.iwlbas0.or.k.ne.k0) then
        k0=k
        iwlbas0=iwlbas(ll,i)
        call slm(iwlbas0,xc(1,k),r2(k),y,dy,ddy,ddy_lap,dlapy,iforce_analy)
      endif
      call phi_combine(iwlbas0,xc(1,k),ri(k),ri2(k),wfv(1,k,irb),y,dy,ddy,ddy_lap,dlapy,
     &                 phin(l,k),dphin(1,l,k),d2phin(l,k),d2phin_all(1,1,l,k),d3phin(1,l,k),iforce_analy)
      call n0_inc(l,k,ic)
 800  continue

 950  continue

c end of numerical orbitals
      endif
  990 continue

      return
      end
c-------------------------------------------------------------------
      subroutine phi_combine(l,xc,ri,ri2,wfv,y,dy,ddy,ddy_lap,dlapy,phi,dphi,d2phi,d2phi_all,d3phi,iforce_analy)
      implicit real*8 (a-h,o-z)
      dimension xc(3),xcri(3),wfv(4),wfvn(4),dy(3),ddy(3,3),dphi(3),d2phi_all(3,3),d3phi(3),dlapy(3)

      parameter (two=2.d0,three=3.d0,four=4.d0,five=5.d0,six=6.d0)

      xcri(1)=xc(1)*ri
      xcri(2)=xc(2)*ri
      xcri(3)=xc(3)*ri
      if(l.eq.1) then
        wfvn(1)=wfv(1)
        wfvn(2)=wfv(2)
        wfvn(3)=wfv(3)
      elseif(l.lt.5) then
        wfvn(3)=ri*(wfv(3)+two*ri*(wfv(1)*ri-wfv(2)))
        wfvn(2)=-wfv(1)*ri2+wfv(2)*ri
        wfvn(1)=wfv(1)*ri
      elseif(l.lt.10) then
        wfvn(3)=ri2*(wfv(3)+two*ri*(three*wfv(1)*ri-two*wfv(2)))
        wfvn(2)=(-two*wfv(1)*ri+wfv(2))*ri2
        wfvn(1)=wfv(1)*ri2
      elseif(l.lt.20) then
        ri3=ri*ri2
        wfvn(3)=ri3*(wfv(3)+six*ri*(two*wfv(1)*ri-wfv(2)))
        wfvn(2)=(-three*wfv(1)*ri+wfv(2))*ri3
        wfvn(1)=wfv(1)*ri3
      else
       stop 'to fix for >f functions'
      endif

      phi=y*wfvn(1)
      d2phi=y*wfvn(3)+y*two*ri*wfvn(2)+ddy_lap*wfvn(1)
      dum=0.d0
      do 5 jj=1,3
        dphi(jj)=y*xcri(jj)*wfvn(2)+dy(jj)*wfvn(1)
        dum=dum+dy(jj)*xcri(jj)
   5  continue
      d2phi=d2phi+two*dum*wfvn(2)

      if(iforce_analy.eq.1) then

      if(l.eq.1) then
        wfvn(4)=wfv(4)
      elseif(l.lt.5) then
        wfvn(4)=-ri*wfvn(3)+ri*(wfv(4)-two*ri*(two*wfv(1)*ri2-two*wfv(2)*ri+wfv(3)))
      elseif(l.lt.10) then
        wfvn(4)=-two*ri*wfvn(3)+ri2*(wfv(4)-two*ri*(six*wfv(1)*ri2-five*wfv(2)*ri+two*wfv(3)))
      elseif(l.lt.20) then
        wfvn(4)=-three*ri*wfvn(3)+ri3*(wfv(4)-six*ri*(four*wfv(1)*ri2-three*wfv(2)*ri+wfv(3)))
      else
       stop 'to fix for >f functions'
      endif

      do 20 jj=1,3
        dum1=0
        do 10 ii=1,3
  10      dum1=dum1+ddy(jj,ii)*xcri(ii)
  20    d3phi(jj)=wfvn(4)*y*xcri(jj)
     &           +wfvn(3)*(dy(jj)+two*xcri(jj)*(y*ri+dum))
     &           +wfvn(2)*(xcri(jj)*(ddy_lap-two*ri*(dum+y*ri))+two*(dum1+two*dy(jj)*ri))
     &           +wfvn(1)*dlapy(jj)

      do 40 jj=1,3
        do 30 ii=jj,3
          prod=xcri(jj)*xcri(ii)
          d2phi_all(ii,jj)=ddy(ii,jj)*wfvn(1)+wfvn(2)*(dy(ii)*xcri(jj)+dy(jj)*xcri(ii)-y*ri*prod)+wfvn(3)*y*prod
  30      d2phi_all(jj,ii)=d2phi_all(ii,jj)
  40    d2phi_all(jj,jj)=d2phi_all(jj,jj)+y*ri*wfvn(2)
      endif

      return
      end
c-------------------------------------------------------------------
      subroutine phie_combine(l,ri,ri2,wfv,y,phi)
      implicit real*8 (a-h,o-z)

      parameter (two=2.d0,three=3.d0,five=5.d0,six=6.d0)

      ri3=ri*ri2
      if(l.eq.1) then
        wfvn=wfv
      elseif(l.lt.5) then
        wfvn=wfv*ri
      elseif(l.lt.10) then
        wfvn=wfv*ri2
      elseif(l.lt.20) then
        wfvn=wfv*ri3
      else
       stop 'to fix for >f functions'
      endif

      phi=y*wfvn

      return
      end
c-------------------------------------------------------------------
      subroutine n0_inc(l,k,ic)

      use phifun, only: dphin, n0_ibasis, n0_ic, n0_nbasis
      use phifun, only: phin
      implicit real*8(a-h,o-z)



      if(abs(phin(l,k))+abs(dphin(1,l,k))+abs(dphin(2,l,k))+abs(dphin(3,l,k)).gt.1.d-20)then
       n0_nbasis(k)=n0_nbasis(k)+1
       n0_ibasis(n0_nbasis(k),k)=l
       n0_ic(n0_nbasis(k),k)=ic
      endif

      return
      end
