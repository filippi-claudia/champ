      subroutine basis_fnse_vgl(k,rvec_en,r_en)
c Written by Cyrus Umrigar and Claudia Filippi, starting from Kevin Schmidt routine
c routine to calculate the values of the basis functions and their derivatives
      use atom, only: znuc, cent, pecent, iwctype, nctype, ncent

      use ghostatom, only: newghostype, nghostcent
      implicit real*8(a-h,o-z)

      include 'vmc.h'
      include 'force.h'
      include 'basis.h'
      include 'numbas.h'
      include 'pseudo.h'

      parameter (one=1.d0,two=2.d0,three=3.d0,four=4.d0)
      parameter (five=5.d0,six=6.d0,seven=7.d0,eight=8.d0)
      parameter (ten=10.d0,half=.5d0)
      parameter (twelve=12.d0)

c positive n1s,n2s,... Slater orbitals; positive nsa, npa, nda asymptotic orbitals
c n1s < 0   exp(-a*r^2)
c nsa < 0   r^2*exp(-a*r^2)
c npa < 0   (xyz)*exp(-a*r^2)
c nda < 0   (zr,x2y2...)*exp(-a*r^2)

      common /const/ pi,hb,etrial,delta,deltai,fbias,nelec,imetro,ipr
      common /phifun/ phin(MBASIS,MELEC),dphin(3,MBASIS,MELEC)
     &,d2phin(MBASIS,MELEC),d2phin_all(3,3,MBASIS,MELEC),d3phin(3,MBASIS,MELEC)
     &,n0_nbasis(MELEC),n0_ibasis(MBASIS,MELEC),n0_ic(MBASIS,MELEC)
      common /numbas/ arg(MCTYPE),r0(MCTYPE)
     &,rwf(MRWF_PTS,MRWF,MCTYPE,MWF),d2rwf(MRWF_PTS,MRWF,MCTYPE,MWF)
     &,numr,nrbas(MCTYPE),igrid(MCTYPE),nr(MCTYPE),iwrwf(MBASIS,MCTYPE)

      common /wfsec/ iwftype(MFORCE),iwf,nwftype

      dimension rvec_en(3,MELEC,MCENT),r_en(MELEC,MCENT)

      dimension wfv(4,MELEC,MRWF),xc(3)

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

      do 900 ic=1,ncent+nghostcent
      ll=0

      i=iwctype(ic)

c get distance to center

      xc(1)=rvec_en(1,k,ic)
      xc(2)=rvec_en(2,k,ic)
      xc(3)=rvec_en(3,k,ic)
      r=r_en(k,ic)
      r2=r*r
      ri=one/r
      ri2=ri*ri
      ri3=ri2*ri
      ri4=ri3*ri
      ri5=ri4*ri

c     write(6,*) (xc(k,1),k=1,3)

c analytical orbital

      if(numr.eq.0) then

      if (iabs(n1s(i)).lt.1) goto 50

c 1s states

      ju=n1s(i)
      do 40 j=1,ju
      l=l+1
      ex=dexp(-zex(l,iwf)*r)
      phin(l,k)=ex
      x1=-zex(l,iwf)*ri
      dphin(1,l,k)=x1*xc(1)*ex
      dphin(2,l,k)=x1*xc(2)*ex
      dphin(3,l,k)=x1*xc(3)*ex
      d2phin(l,k)=(two*x1+zex(l,iwf)**2)*ex
   40 continue
      ju=-ju
      do 45 j=1,ju
      l=l+1
      ex=dexp(-zex(l,iwf)*r2)
      phin(l,k)=ex
      x1=-two*zex(l,iwf)
      dphin(1,l,k)=x1*xc(1)*ex
      dphin(2,l,k)=x1*xc(2)*ex
      dphin(3,l,k)=x1*xc(3)*ex
      d2phin(l,k)=x1*(three+x1*r2)*ex
   45 continue
   50 continue
      if (iabs(n2s(i)).lt.1) goto 70

c 2s states

      ju=n2s(i)
      do 60 j=1,ju
      l=l+1
      ex=dexp(-zex(l,iwf)*r)
      phin(l,k)=r*ex
      x1=ri-zex(l,iwf)
      x2=(zex(l,iwf)**2-ri2-zex(l,iwf)*ri)*ri
      dphin(1,l,k)=x1*xc(1)*ex
      dphin(2,l,k)=x1*xc(2)*ex
      dphin(3,l,k)=x1*xc(3)*ex
      d2phin(l,k)=(two*ri-four*zex(l,iwf)+zex(l,iwf)**2*r)*ex
   60 continue
      if(ju.lt.0) call fatal_error('BAS: negative n2s not implemented')
      ju=-ju
      do 65 j=1,ju
      l=l+1
   65 continue
   70 continue

c 2p states

      if (iabs(n2p(1,i)).lt.1) goto 90
      ju=n2p(1,i)
      do 80 j=1,ju
      l=l+1
      ex=dexp(-zex(l,iwf)*r)
      phin(l,k)=xc(1)*ex
      x1=-zex(l,iwf)*ri
      x2=zex(l,iwf)*(zex(l,iwf)*ri2+ri3)
      dphin(1,l,k)=ex*(one-zex(l,iwf)*xc(1)**2*ri)
      dphin(2,l,k)=x1*xc(2)*phin(l,k)
      dphin(3,l,k)=x1*xc(3)*phin(l,k)
      d2phin(l,k)=(-four*zex(l,iwf)*ri+zex(l,iwf)**2)*phin(l,k)
   80 continue
      if(ju.lt.0) call fatal_error('BAS: negative n2p not implemented')
      ju=-ju
      do 85 j=1,ju
      l=l+1
   85 continue

   90 continue
      if (iabs(n2p(2,i)).lt.1) goto 110
      ju=n2p(2,i)
      do 100 j=1,ju
      l=l+1
      ex=dexp(-zex(l,iwf)*r)
      phin(l,k)=xc(2)*ex
      x1=-zex(l,iwf)*ri
      x2=zex(l,iwf)*(zex(l,iwf)*ri2+ri3)
      dphin(1,l,k)=x1*xc(1)*phin(l,k)
      dphin(2,l,k)=ex*(one-zex(l,iwf)*xc(2)**2*ri)
      dphin(3,l,k)=x1*xc(3)*phin(l,k)
      d2phin(l,k)=(-four*zex(l,iwf)*ri+zex(l,iwf)**2)*phin(l,k)
  100 continue
      if(ju.lt.0) call fatal_error('BAS: negative n2p not implemented')
      ju=-ju
      do 105 j=1,ju
      l=l+1
  105 continue

  110 continue
      if (iabs(n2p(3,i)).lt.1) goto 130
      ju=n2p(3,i)
      do 120 j=1,ju
      l=l+1
      ex=dexp(-zex(l,iwf)*r)
      phin(l,k)=xc(3)*ex
      x1=-zex(l,iwf)*ri
      x2=zex(l,iwf)*(zex(l,iwf)*ri2+ri3)
      dphin(1,l,k)=x1*xc(1)*phin(l,k)
      dphin(2,l,k)=x1*xc(2)*phin(l,k)
      dphin(3,l,k)=ex*(one-zex(l,iwf)*xc(3)**2*ri)
      d2phin(l,k)=(-four*zex(l,iwf)*ri+zex(l,iwf)**2)*phin(l,k)
  120 continue
      if(ju.lt.0) call fatal_error('BAS: negative n2p not implemented')
      ju=-ju
      do 125 j=1,ju
      l=l+1
  125 continue
  130 continue

c 3s state

      if (iabs(n3s(i)).lt.1) goto 155
      ju=n3s(i)
      do 150 j=1,ju
      l=l+1
      ex=dexp(-zex(l,iwf)*r)
      phin(l,k)=r**2*ex
      x1=(two-zex(l,iwf)*r)*ex
      dphin(1,l,k)=x1*xc(1)
      dphin(2,l,k)=x1*xc(2)
      dphin(3,l,k)=x1*xc(3)
      d2phin(l,k)=(six-six*zex(l,iwf)*r+zex(l,iwf)**2*r**2)*ex
  150 continue
      if(ju.lt.0) call fatal_error('BAS: negative n3s not implemented')
      ju=-ju
      do 152 j=1,ju
      l=l+1
  152 continue
  155 continue

c 3p states

      if (iabs(n3p(1,i)).lt.1) goto 170
      ju=n3p(1,i)
      do 160 j=1,ju
      l=l+1
      ex=dexp(-zex(l,iwf)*r)
      phin(l,k)=xc(1)*r*ex
      x1=(ri-zex(l,iwf))*xc(1)*ex
      dphin(1,l,k)=(xc(1)**2*ri-zex(l,iwf)*xc(1)**2+r)*ex
      dphin(2,l,k)=x1*xc(2)
      dphin(3,l,k)=x1*xc(3)
      d2phin(l,k)=(four*ri-six*zex(l,iwf)+zex(l,iwf)**2*r)
     &                 *xc(1)*ex
  160 continue
      if(ju.lt.0) call fatal_error('BAS: negative n3p not implemented')
      ju=-ju
      do 165 j=1,ju
      l=l+1
  165 continue

  170 continue
      if (iabs(n3p(2,i)).lt.1) goto 190
      ju=n3p(2,i)
      do 180 j=1,ju
      l=l+1
      ex=dexp(-zex(l,iwf)*r)
      phin(l,k)=xc(2)*r*ex
      x1=(ri-zex(l,iwf))*xc(2)*ex
      dphin(1,l,k)=x1*xc(1)
      dphin(2,l,k)=(xc(2)**2*ri-zex(l,iwf)*xc(2)**2+r)*ex
      dphin(3,l,k)=x1*xc(3)
      d2phin(l,k)=(four*ri-six*zex(l,iwf)+zex(l,iwf)**2*r)
     &                 *xc(2)*ex
  180 continue
      if(ju.lt.0) call fatal_error('BAS: negative n3p not implemented')
      ju=-ju
      do 185 j=1,ju
      l=l+1
  185 continue

  190 continue
      if (iabs(n3p(3,i)).lt.1) goto 210
      ju=n3p(3,i)
      do 200 j=1,ju
      l=l+1
      ex=dexp(-zex(l,iwf)*r)
      phin(l,k)=xc(3)*r*ex
      x1=(ri-zex(l,iwf))*xc(3)*ex
      dphin(1,l,k)=x1*xc(1)
      dphin(2,l,k)=x1*xc(2)
      dphin(3,l,k)=(xc(3)**2*ri-zex(l,iwf)*xc(3)**2+r)*ex
      d2phin(l,k)=(four*ri-six*zex(l,iwf)+zex(l,iwf)**2*r)
     &                 *xc(3)*ex
  200 continue
      if(ju.lt.0) call fatal_error('BAS: negative n3p not implemented')
      ju=-ju
      do 205 j=1,ju
      l=l+1
  205 continue
  210 continue

c 3d states

      if (iabs(n3dzr(i)).lt.1) goto 230
      do 220 j=1,n3dzr(i)
      l=l+1
      ex=half*dexp(-zex(l,iwf)*r)
      x2=three*xc(3)**2-r**2
      phin(l,k)=x2*ex
      x1=(-zex(l,iwf)*ri*x2-two)*ex
      dphin(1,l,k)=x1*xc(1)
      dphin(2,l,k)=x1*xc(2)
      dphin(3,l,k)=(four-zex(l,iwf)*ri*x2)*xc(3)*ex
      d2phin(l,k)=(zex(l,iwf)**2-six*zex(l,iwf)*ri)*phin(l,k)
  220 continue

  230 continue
      if (iabs(n3dx2(i)).lt.1) goto 250
      do 240 j=1,n3dx2(i)
      l=l+1
      ex=rt3b2*dexp(-zex(l,iwf)*r)
      x2=xc(1)**2-xc(2)**2
      phin(l,k)=x2*ex
      dphin(1,l,k)=(two-zex(l,iwf)*x2*ri)*xc(1)*ex
      dphin(2,l,k)=-(two+zex(l,iwf)*x2*ri)*xc(2)*ex
      dphin(3,l,k)=-zex(l,iwf)*x2*ri*xc(3)*ex
      d2phin(l,k)=(zex(l,iwf)**2-six*zex(l,iwf)*ri)*phin(l,k)
  240 continue

  250 continue
      if (iabs(n3dxy(i)).lt.1) goto 270
      do 260 j=1,n3dxy(i)
      l=l+1
      ex=rt3*dexp(-zex(l,iwf)*r)
      phin(l,k)=xc(1)*xc(2)*ex
      dphin(1,l,k)=(one-zex(l,iwf)*xc(1)**2*ri)*xc(2)*ex
      dphin(2,l,k)=(one-zex(l,iwf)*xc(2)**2*ri)*xc(1)*ex
      dphin(3,l,k)=-zex(l,iwf)*xc(3)*ri*phin(l,k)
      d2phin(l,k)=(zex(l,iwf)**2-six*zex(l,iwf)*ri)*phin(l,k)
  260 continue

  270 continue
      if (iabs(n3dxz(i)).lt.1) goto 290
      do 280 j=1,n3dxz(i)
      l=l+1
      ex=rt3*dexp(-zex(l,iwf)*r)
      phin(l,k)=xc(1)*xc(3)*ex
      dphin(1,l,k)=(one-zex(l,iwf)*xc(1)**2*ri)*xc(3)*ex
      dphin(2,l,k)=-zex(l,iwf)*xc(2)*ri*phin(l,k)
      dphin(3,l,k)=(one-zex(l,iwf)*xc(3)**2*ri)*xc(1)*ex
      d2phin(l,k)=(zex(l,iwf)**2-six*zex(l,iwf)*ri)*phin(l,k)
  280 continue

  290 continue
      if (iabs(n3dyz(i)).lt.1) goto 310
      do 300 j=1,n3dyz(i)
      l=l+1
      ex=rt3*dexp(-zex(l,iwf)*r)
      phin(l,k)=xc(2)*xc(3)*ex
      dphin(1,l,k)=-zex(l,iwf)*xc(1)*ri*phin(l,k)
      dphin(2,l,k)=(one-zex(l,iwf)*xc(2)**2*ri)*xc(3)*ex
      dphin(3,l,k)=(one-zex(l,iwf)*xc(3)**2*ri)*xc(2)*ex
      d2phin(l,k)=(zex(l,iwf)**2-six*zex(l,iwf)*ri)*phin(l,k)
  300 continue
  310 continue

c 4s state

      if (iabs(n4s(i)).lt.1) goto 330
      do 320 j=1,n4s(i)
      l=l+1
      ex=dexp(-zex(l,iwf)*r)
      phin(l,k)=r**3*ex
      x1=(three*r-zex(l,iwf)*r**2)*ex
      dphin(1,l,k)=x1*xc(1)
      dphin(2,l,k)=x1*xc(2)
      dphin(3,l,k)=x1*xc(3)
      d2phin(l,k)=(zex(l,iwf)**2*r**3-eight*zex(l,iwf)*r**2
     &             +twelve*r)*ex
  320 continue
  330 continue

c 4p states

      if (iabs(n4p(1,i)).lt.1) goto 350
      do 340 j=1,n4p(1,i)
      l=l+1
      ex=dexp(-zex(l,iwf)*r)
      phin(l,k)=xc(1)*r**2*ex
      x1=(two-zex(l,iwf)*r)*xc(1)*ex
      dphin(1,l,k)=(-zex(l,iwf)*r*xc(1)**2+two*xc(1)**2+r**2)
     &*ex
      dphin(2,l,k)=x1*xc(2)
      dphin(3,l,k)=x1*xc(3)
      d2phin(l,k)=(zex(l,iwf)**2*r**2-eight*zex(l,iwf)*r+ten)
     &*xc(1)*ex
  340 continue

  350 continue
      if (iabs(n4p(2,i)).lt.1) goto 370
      do 360 j=1,n4p(2,i)
      l=l+1
      ex=dexp(-zex(l,iwf)*r)
      phin(l,k)=xc(2)*r**2*ex
      x1=(two-zex(l,iwf)*r)*xc(2)*ex
      dphin(1,l,k)=x1*xc(1)
      dphin(2,l,k)=(-zex(l,iwf)*r*xc(2)**2+two*xc(2)**2+r**2)
     &*ex
      dphin(3,l,k)=x1*xc(3)
      d2phin(l,k)=(zex(l,iwf)**2*r**2-eight*zex(l,iwf)*r+ten)
     &*xc(2)*ex
  360 continue

  370 continue
      if (iabs(n4p(3,i)).lt.1) goto 390
      do 380 j=1,n4p(3,i)
      l=l+1
      ex=dexp(-zex(l,iwf)*r)
      phin(l,k)=xc(3)*r**2*ex
      x1=(two-zex(l,iwf)*r)*xc(3)*ex
      dphin(1,l,k)=x1*xc(1)
      dphin(2,l,k)=x1*xc(2)
      dphin(3,l,k)=(-zex(l,iwf)*r*xc(3)**2+two*xc(3)**2+r**2)
     &*ex
      d2phin(l,k)=(zex(l,iwf)**2*r**2-eight*zex(l,iwf)*r+ten)
     &*xc(3)*ex
  380 continue
  390 continue

c s asymptotic/gaussian states

      if (iabs(nsa(i)).lt.1) goto 410
      ju=nsa(i)
      do 400 j=1,ju
      l=l+1
      beta=betaq/zex(l,iwf)-one
      bb1=beta*(beta-one)
      tzb=two*zex(l,iwf)*beta
      z2=zex(l,iwf)**2
      r1=r+one
      r1i=one/r1
      ex=dexp(-zex(l,iwf)*r)
      phin(l,k)=r1**beta*ex
      x1=(beta*r1i-zex(l,iwf))*ri
      x2=bb1*r1i**2-tzb*r1i+z2
      dphin(1,l,k)=x1*xc(1)*phin(l,k)
      dphin(2,l,k)=x1*xc(2)*phin(l,k)
      dphin(3,l,k)=x1*xc(3)*phin(l,k)
      d2phin(l,k)=(two*x1+x2)*phin(l,k)
  400 continue
      ju=-ju
      do 405 j=1,ju
      l=l+1
      ex=dexp(-zex(l,iwf)*r2)
      phin(l,k)=r2*ex
      x1=-two*zex(l,iwf)
      dphin(1,l,k)=(two+x1*r2)*xc(1)*ex
      dphin(2,l,k)=(two+x1*r2)*xc(2)*ex
      dphin(3,l,k)=(two+x1*r2)*xc(3)*ex
      d2phin(l,k)=(six+x1*r2*(seven+x1*r2))*ex
  405 continue
  410 continue

c p asymptotic/gaussian states

      if (iabs(npa(1,i)).lt.1) goto 430
      ju=npa(1,i)
      do 420 j=1,ju
      l=l+1
      beta=betaq/zex(l,iwf)-one
      b1=beta-one
      b1b2=b1*(beta-two)
      tb1=two*b1
      tb1z=tb1*zex(l,iwf)
      tz=two*zex(l,iwf)
      z2=zex(l,iwf)**2
      r1=r+one
      r1i=one/r1
      ex=r1**b1*dexp(-zex(l,iwf)*r)
      phin(l,k)=xc(1)*ex
      x1=(b1*r1i-zex(l,iwf))*ri
      x2=b1b2*r1i**2+tb1*r1i*ri-tz*ri-tb1z*r1i+z2
      dphin(1,l,k)=ex*(one+x1*xc(1)**2)
      dphin(2,l,k)=x1*xc(2)*phin(l,k)
      dphin(3,l,k)=x1*xc(3)*phin(l,k)
      d2phin(l,k)=(two*x1+x2)*phin(l,k)
  420 continue
      ju=-ju
      do 425 j=1,ju
      l=l+1
      ex=dexp(-zex(l,iwf)*r2)
      phin(l,k)=xc(1)*ex
      x1=-two*zex(l,iwf)
      dphin(1,l,k)=(one+x1*xc(1)*xc(1))*ex
      dphin(2,l,k)=x1*xc(2)*phin(l,k)
      dphin(3,l,k)=x1*xc(3)*phin(l,k)
      d2phin(l,k)=x1*(five+x1*r2)*phin(l,k)
  425 continue

  430 continue
      if (iabs(npa(2,i)).lt.1) goto 450
      ju=npa(2,i)
      do 440 j=1,ju
      l=l+1
      beta=betaq/zex(l,iwf)-one
      b1=beta-one
      b1b2=b1*(beta-two)
      tb1=two*b1
      tb1z=tb1*zex(l,iwf)
      tz=two*zex(l,iwf)
      z2=zex(l,iwf)**2
      r1=r+one
      r1i=one/r1
      ex=r1**b1*dexp(-zex(l,iwf)*r)
      phin(l,k)=xc(2)*ex
      x1=(b1*r1i-zex(l,iwf))*ri
      x2=b1b2*r1i**2+tb1*r1i*ri-tz*ri-tb1z*r1i+z2
      dphin(1,l,k)=x1*xc(1)*phin(l,k)
      dphin(2,l,k)=ex*(one+x1*xc(2)**2)
      dphin(3,l,k)=x1*xc(3)*phin(l,k)
      d2phin(l,k)=(two*x1+x2)*phin(l,k)
  440 continue
      ju=-ju
      do 445 j=1,ju
      l=l+1
      ex=dexp(-zex(l,iwf)*r2)
      phin(l,k)=xc(2)*ex
      x1=-two*zex(l,iwf)
      dphin(1,l,k)=x1*xc(1)*phin(l,k)
      dphin(2,l,k)=(one+x1*xc(2)*xc(2))*ex
      dphin(3,l,k)=x1*xc(3)*phin(l,k)
      d2phin(l,k)=x1*(five+x1*r2)*phin(l,k)
  445 continue

  450 continue
      if (iabs(npa(3,i)).lt.1) goto 470
      ju=npa(3,i)
      do 460 j=1,ju
      l=l+1
      beta=betaq/zex(l,iwf)-one
      b1=beta-one
      b1b2=b1*(beta-two)
      tb1=two*b1
      tb1z=tb1*zex(l,iwf)
      tz=two*zex(l,iwf)
      z2=zex(l,iwf)**2
      r1=r+one
      r1i=one/r1
      ex=r1**b1*dexp(-zex(l,iwf)*r)
      phin(l,k)=xc(3)*ex
      x1=(b1*r1i-zex(l,iwf))*ri
      x2=b1b2*r1i**2+tb1*r1i*ri-tz*ri-tb1z*r1i+z2
      dphin(1,l,k)=x1*xc(1)*phin(l,k)
      dphin(2,l,k)=x1*xc(2)*phin(l,k)
      dphin(3,l,k)=ex*(one+x1*xc(3)**2)
      d2phin(l,k)=(two*x1+x2)*phin(l,k)
  460 continue
      ju=-ju
      do 465 j=1,ju
      l=l+1
      ex=dexp(-zex(l,iwf)*r2)
      phin(l,k)=xc(3)*ex
      x1=-two*zex(l,iwf)
      dphin(1,l,k)=x1*xc(1)*phin(l,k)
      dphin(2,l,k)=x1*xc(2)*phin(l,k)
      dphin(3,l,k)=(one+x1*xc(3)*xc(3))*ex
      d2phin(l,k)=x1*(five+x1*r2)*phin(l,k)
  465 continue
  470 continue

c d asymptotic/gaussian states

      if (ndzra(i).eq.0) goto 490
      ju=ndzra(i)
      do 480 j=1,ju
      l=l+1
      beta=betaq/zex(l,iwf)-one
      b2=beta-two
      b2b3=b2*(beta-three)
      fb2=four*b2
      fz=four*zex(l,iwf)
      tb2z=two*b2*zex(l,iwf)
      r1=r+one
      r1i=one/r1
      ex=half*r1**b2*dexp(-zex(l,iwf)*r)
      xa=three*xc(3)**2-r**2
      phin(l,k)=xa*ex
      x1=(b2*r1i-zex(l,iwf))*ri
      x2=b2b3*r1i**2+two*ri**2+fb2*ri*r1i+fz*ri-tb2z*r1i+z2
      dphin(1,l,k)=ex*(x1*xa-two)*xc(1)
      dphin(2,l,k)=ex*(x1*xa-two)*xc(2)
      dphin(3,l,k)=ex*(x1*xa+four)*xc(3)
      d2phin(l,k)=(two*x1+x2)*phin(l,k)
  480 continue
      ju=-ju
      do 485 j=1,ju
      l=l+1
      ex=dexp(-zex(l,iwf)*r2)
      phin(l,k)=(three*xc(3)**2-r**2)*ex
      x1=-two*zex(l,iwf)
      dphin(1,l,k)=x1*xc(1)*phin(l,k)-two*xc(1)*ex
      dphin(2,l,k)=x1*xc(2)*phin(l,k)-two*xc(2)*ex
      dphin(3,l,k)=x1*xc(3)*phin(l,k)+four*xc(3)*ex
      d2phin(l,k)=x1*(seven+x1*r2)*phin(l,k)
  485 continue

  490 continue
      if (ndx2a(i).eq.0) goto 510
      ju=ndx2a(i)
      do 500 j=1,ju
      l=l+1
      beta=betaq/zex(l,iwf)-one
      b2=beta-two
      b2b3=b2*(beta-three)
      fb2=four*b2
      fz=four*zex(l,iwf)
      tb2z=two*b2*zex(l,iwf)
      r1=r+one
      r1i=one/r1
      ex=rt3b2*r1**b2*dexp(-zex(l,iwf)*r)
      xa=xc(1)**2-xc(2)**2
      phin(l,k)=xa*ex
      x1=(b2*r1i-zex(l,iwf))*ri
      x2=b2b3*r1i**2+two*ri**2+fb2*ri*r1i+fz*ri-tb2z*r1i+z2
      dphin(1,l,k)=ex*(x1*xa+two)*xc(1)
      dphin(2,l,k)=ex*(x1*xa-two)*xc(2)
      dphin(3,l,k)=ex*x1*xa*xc(3)
      d2phin(l,k)=(two*x1+x2)*phin(l,k)
  500 continue
      ju=-ju
      do 505 j=1,ju
      l=l+1
      ex=dexp(-zex(l,iwf)*r2)
      phin(l,k)=(xc(1)**2-xc(2)**2)*ex
      x1=-two*zex(l,iwf)
      dphin(1,l,k)=x1*xc(1)*phin(l,k)+two*xc(1)*ex
      dphin(2,l,k)=x1*xc(2)*phin(l,k)-two*xc(2)*ex
      dphin(3,l,k)=x1*xc(3)*phin(l,k)
      d2phin(l,k)=x1*(seven+x1*r2)*phin(l,k)
  505 continue

  510 continue
      if (ndxya(i).eq.0) goto 530
      ju=ndxya(i)
      do 520 j=1,ju
      l=l+1
      beta=betaq/zex(l,iwf)-one
      b2=beta-two
      b2b3=b2*(beta-three)
      fb2=four*b2
      fz=four*zex(l,iwf)
      tb2z=two*b2*zex(l,iwf)
      r1=r+one
      r1i=one/r1
      ex=rt3*r1**b2*dexp(-zex(l,iwf)*r)
      xa=xc(1)*xc(2)
      phin(l,k)=xa*ex
      x1=(b2*r1i-zex(l,iwf))*ri
      x2=b2b3*r1i**2+two*ri**2+fb2*ri*r1i+fz*ri-tb2z*r1i+z2
      dphin(1,l,k)=ex*(x1*xa*xc(1)+xc(2))
      dphin(2,l,k)=ex*(x1*xa*xc(2)+xc(1))
      dphin(3,l,k)=ex*x1*xa*xc(3)
      d2phin(l,k)=(two*x1+x2)*phin(l,k)
  520 continue
      ju=-ju
      do 525 j=1,ju
      l=l+1
      ex=dexp(-zex(l,iwf)*r2)
      phin(l,k)=xc(1)*xc(2)*ex
      x1=-two*zex(l,iwf)
      dphin(1,l,k)=x1*xc(1)*phin(l,k)+xc(2)*ex
      dphin(2,l,k)=x1*xc(2)*phin(l,k)+xc(1)*ex
      dphin(3,l,k)=x1*xc(3)*phin(l,k)
      d2phin(l,k)=x1*(seven+x1*r2)*phin(l,k)
  525 continue

  530 continue
      if (ndxza(i).eq.0) goto 550
      ju=ndxza(i)
      do 540 j=1,ju
      l=l+1
      beta=betaq/zex(l,iwf)-one
      b2=beta-two
      b2b3=b2*(beta-three)
      fb2=four*b2
      fz=four*zex(l,iwf)
      tb2z=two*b2*zex(l,iwf)
      r1=r+one
      r1i=one/r1
      ex=rt3*r1**b2*dexp(-zex(l,iwf)*r)
      xa=xc(1)*xc(3)
      phin(l,k)=xa*ex
      x1=(b2*r1i-zex(l,iwf))*ri
      x2=b2b3*r1i**2+two*ri**2+fb2*ri*r1i+fz*ri-tb2z*r1i+z2
      dphin(1,l,k)=ex*(x1*xa*xc(1)+xc(3))
      dphin(2,l,k)=ex*x1*xa*xc(2)
      dphin(3,l,k)=ex*(x1*xa*xc(3)+xc(1))
      d2phin(l,k)=(two*x1+x2)*phin(l,k)
  540 continue
      ju=-ju
      do 545 j=1,ju
      l=l+1
      ex=dexp(-zex(l,iwf)*r2)
      phin(l,k)=xc(1)*xc(3)*ex
      x1=-two*zex(l,iwf)
      dphin(1,l,k)=x1*xc(1)*phin(l,k)+xc(3)*ex
      dphin(2,l,k)=x1*xc(2)*phin(l,k)
      dphin(3,l,k)=x1*xc(3)*phin(l,k)+xc(1)*ex
      d2phin(l,k)=x1*(seven+x1*r2)*phin(l,k)
  545 continue

  550 continue
      if (ndyza(i).eq.0) goto 570
      ju=ndyza(i)
      do 560 j=1,ju
      l=l+1
      beta=betaq/zex(l,iwf)-one
      b2=beta-two
      b2b3=b2*(beta-three)
      fb2=four*b2
      fz=four*zex(l,iwf)
      tb2z=two*b2*zex(l,iwf)
      r1=r+one
      r1i=one/r1
      ex=rt3*r1**b2*dexp(-zex(l,iwf)*r)
      xa=xc(2)*xc(3)
      phin(l,k)=xa*ex
      x1=(b2*r1i-zex(l,iwf))*ri
      x2=b2b3*r1i**2+two*ri**2+fb2*ri*r1i+fz*ri-tb2z*r1i+z2
      dphin(1,l,k)=ex*x1*xa*xc(1)
      dphin(2,l,k)=ex*(x1*xa*xc(2)+xc(3))
      dphin(3,l,k)=ex*(x1*xa*xc(3)+xc(2))
      d2phin(l,k)=(two*x1+x2)*phin(l,k)
  560 continue
      ju=-ju
      do 565 j=1,ju
      l=l+1
      ex=dexp(-zex(l,iwf)*r2)
      phin(l,k)=xc(2)*xc(3)*ex
      x1=-two*zex(l,iwf)
      dphin(1,l,k)=x1*xc(1)*phin(l,k)
      dphin(2,l,k)=x1*xc(2)*phin(l,k)+xc(3)*ex
      dphin(3,l,k)=x1*xc(3)*phin(l,k)+xc(2)*ex
      d2phin(l,k)=x1*(seven+x1*r2)*phin(l,k)
  565 continue
  570 continue

      else

c numerical orbitals

      ider=1
      rk=r
      do 600 irb=1,nrbas(i)
      call splfit(rk,irb,i,iwf,wfv(1,k,irb),ider)
  600 continue

c s states

      if (iabs(n1s(i)).lt.1) goto 620
      ju=iabs(n1s(i))
      do 610 j=1,ju
      l=l+1
      ll=ll+1
      irb=iwrwf(ll,i)
      phin(l,k)=cs*wfv(1,k,irb)
      dphin(1,l,k)=cs*xc(1)*ri*wfv(2,k,irb)
      dphin(2,l,k)=cs*xc(2)*ri*wfv(2,k,irb)
      dphin(3,l,k)=cs*xc(3)*ri*wfv(2,k,irb)
      d2phin(l,k)=cs*(wfv(3,k,irb)+two*ri*wfv(2,k,irb))
      call n0_inc(l,k,ic)
  610 continue
  620 continue

c p states

      if (iabs(n2p(1,i)).lt.1) goto 640
      ju=iabs(n2p(1,i))
      do 630 j=1,ju
      l=l+1
      ll=ll+1
      irb=iwrwf(ll,i)
      xvec=xc(1)*ri
      yvec=xc(2)*ri
      zvec=xc(3)*ri
      wf1=wfv(2,k,irb)-ri*wfv(1,k,irb)
      phin(l,k)=cp*xvec*wfv(1,k,irb)
      dphin(1,l,k)=cp*xvec*xvec*wf1+cp*ri*wfv(1,k,irb)
      dphin(2,l,k)=cp*xvec*yvec*wf1
      dphin(3,l,k)=cp*xvec*zvec*wf1
      d2phin(l,k)=cp*xvec*(wfv(3,k,irb)+two*ri*wf1)
      call n0_inc(l,k,ic)
  630 continue
  640 continue

      if (iabs(n2p(2,i)).lt.1) goto 660
      ju=iabs(n2p(2,i))
      do 650 j=1,ju
      l=l+1
      ll=ll+1
      irb=iwrwf(ll,i)
      xvec=xc(1)*ri
      yvec=xc(2)*ri
      zvec=xc(3)*ri
      wf1=wfv(2,k,irb)-ri*wfv(1,k,irb)
      phin(l,k)=cp*yvec*wfv(1,k,irb)
      dphin(1,l,k)=cp*yvec*xvec*wf1
      dphin(2,l,k)=cp*yvec*yvec*wf1+cp*ri*wfv(1,k,irb)
      dphin(3,l,k)=cp*yvec*zvec*wf1
      d2phin(l,k)=cp*yvec*(wfv(3,k,irb)+two*ri*wf1)
      call n0_inc(l,k,ic)
  650 continue
  660 continue

      if (iabs(n2p(3,i)).lt.1) goto 680
      ju=iabs(n2p(3,i))
      do 670 j=1,ju
      l=l+1
      ll=ll+1
      irb=iwrwf(ll,i)
      xvec=xc(1)*ri
      yvec=xc(2)*ri
      zvec=xc(3)*ri
      wf1=wfv(2,k,irb)-ri*wfv(1,k,irb)
      phin(l,k)=cp*zvec*wfv(1,k,irb)
      dphin(1,l,k)=cp*zvec*xvec*wf1
      dphin(2,l,k)=cp*zvec*yvec*wf1
      dphin(3,l,k)=cp*zvec*zvec*wf1+cp*ri*wfv(1,k,irb)
      d2phin(l,k)=cp*zvec*(wfv(3,k,irb)+two*ri*wf1)
      call n0_inc(l,k,ic)
  670 continue
  680 continue

c d states

      if (iabs(n3dzr(i)).lt.1) goto 700
      ju=iabs(n3dzr(i))
      do 690 j=1,ju
      l=l+1
      ll=ll+1
      irb=iwrwf(ll,i)
      zr=half*(three*xc(3)**2-r2)*ri2
      xvec=xc(1)*ri
      yvec=xc(2)*ri
      zvec=xc(3)*ri
      wf1=wfv(2,k,irb)-two*ri*wfv(1,k,irb)
      wf2=wfv(2,k,irb)-three*ri*wfv(1,k,irb)
      phin(l,k)=cd1*zr*wfv(1,k,irb)
      dphin(1,l,k)=cd1*xvec*(zr*wf1-wfv(1,k,irb)*ri)
      dphin(2,l,k)=cd1*yvec*(zr*wf1-wfv(1,k,irb)*ri)
      dphin(3,l,k)=cd1*zvec*(zr*wf1+two*wfv(1,k,irb)*ri)
      d2phin(l,k)=cd1*zr*(wfv(3,k,irb)+two*ri*wf2)
      call n0_inc(l,k,ic)
  690 continue

  700 continue
      if (iabs(n3dx2(i)).lt.1) goto 720
      ju=iabs(n3dx2(i))
      do 710 j=1,ju
      l=l+1
      ll=ll+1
      irb=iwrwf(ll,i)
      x2y2=half*(xc(1)**2-xc(2)**2)*ri2
      xvec=xc(1)*ri
      yvec=xc(2)*ri
      zvec=xc(3)*ri
      wf1=wfv(2,k,irb)-two*ri*wfv(1,k,irb)
      wf2=wfv(2,k,irb)-three*ri*wfv(1,k,irb)
      phin(l,k)=cd2*x2y2*wfv(1,k,irb)
      dphin(1,l,k)=cd2*xvec*(x2y2*wf1+wfv(1,k,irb)*ri)
      dphin(2,l,k)=cd2*yvec*(x2y2*wf1-wfv(1,k,irb)*ri)
      dphin(3,l,k)=cd2*zvec* x2y2*wf1
      d2phin(l,k)=cd2*x2y2*(wfv(3,k,irb)+two*ri*wf2)
      call n0_inc(l,k,ic)
  710 continue

  720 continue
      if (iabs(n3dxy(i)).lt.1) goto 740
      ju=iabs(n3dxy(i))
      do 730 j=1,ju
      l=l+1
      ll=ll+1
      irb=iwrwf(ll,i)
      xy=xc(1)*xc(2)*ri2
      xvec=xc(1)*ri
      yvec=xc(2)*ri
      zvec=xc(3)*ri
      wf1=wfv(2,k,irb)-two*ri*wfv(1,k,irb)
      wf2=wfv(2,k,irb)-three*ri*wfv(1,k,irb)
      phin(l,k)=cd2*xy*wfv(1,k,irb)
      dphin(1,l,k)=cd2*(xvec*xy*wf1+yvec*wfv(1,k,irb)*ri)
      dphin(2,l,k)=cd2*(yvec*xy*wf1+xvec*wfv(1,k,irb)*ri)
      dphin(3,l,k)=cd2* zvec*xy*wf1
      d2phin(l,k)=cd2*xy*(wfv(3,k,irb)+two*ri*wf2)
      call n0_inc(l,k,ic)
  730 continue

  740 continue
      if (iabs(n3dxz(i)).lt.1) goto 760
      ju=iabs(n3dxz(i))
      do 750 j=1,ju
      l=l+1
      ll=ll+1
      irb=iwrwf(ll,i)
      xz=xc(1)*xc(3)*ri2
      xvec=xc(1)*ri
      yvec=xc(2)*ri
      zvec=xc(3)*ri
      wf1=wfv(2,k,irb)-two*ri*wfv(1,k,irb)
      wf2=wfv(2,k,irb)-three*ri*wfv(1,k,irb)
      phin(l,k)=cd2*xz*wfv(1,k,irb)
      dphin(1,l,k)=cd2*(xvec*xz*wf1+zvec*wfv(1,k,irb)*ri)
      dphin(2,l,k)=cd2* yvec*xz*wf1
      dphin(3,l,k)=cd2*(zvec*xz*wf1+xvec*wfv(1,k,irb)*ri)
      d2phin(l,k)=cd2*xz*(wfv(3,k,irb)+two*ri*wf2)
      call n0_inc(l,k,ic)
  750 continue

  760 continue
      if (iabs(n3dyz(i)).lt.1) goto 780
      ju=iabs(n3dyz(i))
      do 770 j=1,ju
      l=l+1
      ll=ll+1
      irb=iwrwf(ll,i)
      yz=xc(2)*xc(3)*ri2
      xvec=xc(1)*ri
      yvec=xc(2)*ri
      zvec=xc(3)*ri
      wf1=wfv(2,k,irb)-two*ri*wfv(1,k,irb)
      wf2=wfv(2,k,irb)-three*ri*wfv(1,k,irb)
      phin(l,k)=cd2*yz*wfv(1,k,irb)
      dphin(1,l,k)=cd2* xvec*yz*wf1
      dphin(2,l,k)=cd2*(yvec*yz*wf1+zvec*wfv(1,k,irb)*ri)
      dphin(3,l,k)=cd2*(zvec*yz*wf1+yvec*wfv(1,k,irb)*ri)
      d2phin(l,k)=cd2*yz*(wfv(3,k,irb)+two*ri*wf2)
      call n0_inc(l,k,ic)
  770 continue

  780 continue
      if (iabs(n4fxxx(i)).lt.1) goto 790
      ju=iabs(n4fxxx(i))
      do 785 j=1,ju
      l=l+1
      ll=ll+1
      irb=iwrwf(ll,i)
      xvec=xc(1)*ri
      yvec=xc(2)*ri
      zvec=xc(3)*ri
      xvec3=xvec**3
      phin(l,k)=cf*xvec3*wfv(1,k,irb)
      term=cf*xvec3*(wfv(2,k,irb)-3*wfv(1,k,irb)*ri)
      dphin(1,l,k)=cf*3*xvec**2*wfv(1,k,irb)*ri+term*xvec
      dphin(2,l,k)=term*yvec
      dphin(3,l,k)=term*zvec
      d2phin(l,k)=cf*xvec3*(wfv(3,k,irb)+2*ri*wfv(2,k,irb))
     &            +6*(cf*xvec*wfv(1,k,irb)-2*phin(l,k))*ri2
      call n0_inc(l,k,ic)
 785  continue

 790  continue
      if (iabs(n4fyyy(i)).lt.1) goto 805
      ju=iabs(n4fyyy(i))
      do 795 j=1,ju
      l=l+1
      ll=ll+1
      irb=iwrwf(ll,i)
      xvec=xc(1)*ri
      yvec=xc(2)*ri
      zvec=xc(3)*ri
      yvec3=yvec**3
      phin(l,k)=cf*yvec3*wfv(1,k,irb)
      term=cf*yvec3*(wfv(2,k,irb)-3*wfv(1,k,irb)*ri)
      dphin(1,l,k)=term*xvec
      dphin(2,l,k)=cf*3*yvec**2*wfv(1,k,irb)*ri+term*yvec
      dphin(3,l,k)=term*zvec
      d2phin(l,k)=cf*yvec3*(wfv(3,k,irb)+2*ri*wfv(2,k,irb))
     &            +6*(cf*yvec*wfv(1,k,irb)-2*phin(l,k))*ri2
      call n0_inc(l,k,ic)
 795  continue

 805  continue
      if (iabs(n4fzzz(i)).lt.1) goto 810
      ju=iabs(n4fzzz(i))
      do 808 j=1,ju
      l=l+1
      ll=ll+1
      irb=iwrwf(ll,i)
      xvec=xc(1)*ri
      yvec=xc(2)*ri
      zvec=xc(3)*ri
      zvec3=zvec**3
      phin(l,k)=cf*zvec3*wfv(1,k,irb)
      term=cf*zvec3*(wfv(2,k,irb)-3*wfv(1,k,irb)*ri)
      dphin(1,l,k)=term*xvec
      dphin(2,l,k)=term*yvec
      dphin(3,l,k)=cf*3*zvec**2*wfv(1,k,irb)*ri+term*zvec
      d2phin(l,k)=cf*zvec3*(wfv(3,k,irb)+2*ri*wfv(2,k,irb))
     &            +6*(cf*zvec*wfv(1,k,irb)-2*phin(l,k))*ri2
      call n0_inc(l,k,ic)
 808  continue

 810  continue
      if (iabs(n4fxxy(i)).lt.1) goto 820
      ju=iabs(n4fxxy(i))
      do 815 j=1,ju
      l=l+1
      ll=ll+1
      irb=iwrwf(ll,i)
      xvec=xc(1)*ri
      yvec=xc(2)*ri
      zvec=xc(3)*ri
      x2y=xvec**2*yvec
      phin(l,k)=cf2*x2y*wfv(1,k,irb)
      term=cf2*x2y*(wfv(2,k,irb)-3*wfv(1,k,irb)*ri)
      dphin(1,l,k)=cf2*2*xvec*yvec*wfv(1,k,irb)*ri+term*xvec
      dphin(2,l,k)=cf2*xvec*xvec*wfv(1,k,irb)*ri+term*yvec
      dphin(3,l,k)=term*zvec
      d2phin(l,k)=cf2*x2y*(wfv(3,k,irb)+2*ri*wfv(2,k,irb))
     &            +2*(cf2*yvec*wfv(1,k,irb)-6*phin(l,k))*ri2
      call n0_inc(l,k,ic)
 815  continue

 820  continue
      if (iabs(n4fxxz(i)).lt.1) goto 825
      ju=iabs(n4fxxz(i))
      do 818 j=1,ju
      l=l+1
      ll=ll+1
      irb=iwrwf(ll,i)
      xvec=xc(1)*ri
      yvec=xc(2)*ri
      zvec=xc(3)*ri
      x2z=xvec**2*zvec
      phin(l,k)=cf2*x2z*wfv(1,k,irb)
      term=cf2*x2z*(wfv(2,k,irb)-3*wfv(1,k,irb)*ri)
      dphin(1,l,k)=cf2*2*xvec*zvec*wfv(1,k,irb)*ri+term*xvec
      dphin(2,l,k)=term*yvec
      dphin(3,l,k)=cf2*xvec*xvec*wfv(1,k,irb)*ri+term*zvec
      d2phin(l,k)=cf2*x2z*(wfv(3,k,irb)+2*ri*wfv(2,k,irb))
     &            +2*(cf2*zvec*wfv(1,k,irb)-6*phin(l,k))*ri2
      call n0_inc(l,k,ic)
 818  continue

 825  continue
      if (iabs(n4fyyx(i)).lt.1) goto 835
      ju=iabs(n4fyyx(i))
      do 830 j=1,ju
      l=l+1
      ll=ll+1
      irb=iwrwf(ll,i)
      xvec=xc(1)*ri
      yvec=xc(2)*ri
      zvec=xc(3)*ri
      y2x=yvec**2*xvec
      phin(l,k)=cf2*y2x*wfv(1,k,irb)
      term=cf2*y2x*(wfv(2,k,irb)-3*wfv(1,k,irb)*ri)
      dphin(1,l,k)=cf2*yvec*yvec*wfv(1,k,irb)*ri+term*xvec
      dphin(2,l,k)=cf2*2*yvec*xvec*wfv(1,k,irb)*ri+term*yvec
      dphin(3,l,k)=term*zvec
      d2phin(l,k)=cf2*y2x*(wfv(3,k,irb)+2*ri*wfv(2,k,irb))
     &            +2*(cf2*xvec*wfv(1,k,irb)-6*phin(l,k))*ri2
      call n0_inc(l,k,ic)
 830  continue

 835  continue
      if (iabs(n4fyyz(i)).lt.1) goto 845
      ju=iabs(n4fyyz(i))
      do 840 j=1,ju
      l=l+1
      ll=ll+1
      irb=iwrwf(ll,i)
      xvec=xc(1)*ri
      yvec=xc(2)*ri
      zvec=xc(3)*ri
      y2z=yvec**2*zvec
      phin(l,k)=cf2*y2z*wfv(1,k,irb)
      term=cf2*y2z*(wfv(2,k,irb)-3*wfv(1,k,irb)*ri)
      dphin(1,l,k)=term*xvec
      dphin(2,l,k)=cf2*2*yvec*zvec*wfv(1,k,irb)*ri+term*yvec
      dphin(3,l,k)=cf2*yvec*yvec*wfv(1,k,irb)*ri+term*zvec
      d2phin(l,k)=cf2*y2z*(wfv(3,k,irb)+2*ri*wfv(2,k,irb))
     &            +2*(cf2*zvec*wfv(1,k,irb)-6*phin(l,k))*ri2
      call n0_inc(l,k,ic)
 840  continue

 845  continue
      if (iabs(n4fzzx(i)).lt.1) goto 855
      ju=iabs(n4fzzx(i))
      do 850 j=1,ju
      l=l+1
      ll=ll+1
      irb=iwrwf(ll,i)
      xvec=xc(1)*ri
      yvec=xc(2)*ri
      zvec=xc(3)*ri
      z2x=zvec**2*xvec
      phin(l,k)=cf2*z2x*wfv(1,k,irb)
      term=cf2*z2x*(wfv(2,k,irb)-3*wfv(1,k,irb)*ri)
      dphin(1,l,k)=cf2*zvec*zvec*wfv(1,k,irb)*ri+term*xvec
      dphin(2,l,k)=term*yvec
      dphin(3,l,k)=cf2*2*zvec*xvec*wfv(1,k,irb)*ri+term*zvec
      d2phin(l,k)=cf2*z2x*(wfv(3,k,irb)+2*ri*wfv(2,k,irb))
     &            +2*(cf2*xvec*wfv(1,k,irb)-6*phin(l,k))*ri2
      call n0_inc(l,k,ic)
 850  continue

 855  continue
      if (iabs(n4fzzy(i)).lt.1) goto 865
      ju=iabs(n4fzzy(i))
      do 860 j=1,ju
      l=l+1
      ll=ll+1
      irb=iwrwf(ll,i)
      xvec=xc(1)*ri
      yvec=xc(2)*ri
      zvec=xc(3)*ri
      z2y=zvec**2*yvec
      phin(l,k)=cf2*z2y*wfv(1,k,irb)
      term=cf2*z2y*(wfv(2,k,irb)-3*wfv(1,k,irb)*ri)
      dphin(1,l,k)=term*xvec
      dphin(2,l,k)=cf2*zvec*zvec*wfv(1,k,irb)*ri+term*yvec
      dphin(3,l,k)=cf2*2*zvec*yvec*wfv(1,k,irb)*ri+term*zvec
      d2phin(l,k)=cf2*z2y*(wfv(3,k,irb)+2*ri*wfv(2,k,irb))
     &            +2*(cf2*yvec*wfv(1,k,irb)-6*phin(l,k))*ri2
      call n0_inc(l,k,ic)
 860  continue


 865  continue
      if (iabs(n4fxyz(i)).lt.1) goto 875
      ju=iabs(n4fxyz(i))
      do 870 j=1,ju
      l=l+1
      ll=ll+1
      irb=iwrwf(ll,i)
      xvec=xc(1)*ri
      yvec=xc(2)*ri
      zvec=xc(3)*ri
      xyz=xvec*yvec*zvec
      phin(l,k)=cf3*xyz*wfv(1,k,irb)
      term=cf3*xyz*(wfv(2,k,irb)-3*wfv(1,k,irb)*ri)
      dphin(1,l,k)=cf3*yvec*zvec*wfv(1,k,irb)*ri+term*xvec
      dphin(2,l,k)=cf3*xvec*zvec*wfv(1,k,irb)*ri+term*yvec
      dphin(3,l,k)=cf3*xvec*yvec*wfv(1,k,irb)*ri+term*zvec
      d2phin(l,k)=cf3*xyz*(wfv(3,k,irb)+2*ri*wfv(2,k,irb))
     &            -12*phin(l,k)*ri2
      call n0_inc(l,k,ic)
 870  continue

 875  continue
      endif

  900 continue

      return
      end
