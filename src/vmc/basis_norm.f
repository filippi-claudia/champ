      subroutine basis_norm(iwf,anorm,iflag)
c Written by Cyrus Umrigar and Claudia Filippi, starting from Kevin Schmidt routine
c Set normalization of basis fns.
      use force_mod, only: MFORCE, MFORCE_WT_PRD, MWF
      use vmc, only: MELEC, MORB, MBASIS, MDET, MCENT, MCTYPE, MCTYP3X
      use vmc, only: NSPLIN, nrad, MORDJ, MORDJ1, MMAT_DIM, MMAT_DIM2, MMAT_DIM20
      use vmc, only: radmax, delri
      use vmc, only: NEQSX, MTERMS
      use vmc, only: MCENT3, NCOEF, MEXCIT
      use atom, only: iwctype, ncent

      use ghostatom, only: nghostcent
      use const, only: pi
      use coefs, only: coef, nbasis, norb
      use basis, only: zex, n1s, n2s, n2p, n3s, n3p, n3dzr, n3dx2, n3dxy, n3dxz, n3dyz
      use basis, only: n4s, n4p

      use orbval, only: ddorb, dorb, nadorb, ndetorb, orb
      implicit real*8(a-h,o-z)



      parameter (one=1.d0,two=2.d0,three=3.d0)
      parameter (five=5.d0,six=6.d0,seven=7.d0)
      parameter (third=1.d0/3.d0)


      dimension anorm(*)

      do 5 j=1,nbasis
    5   anorm(j)=one

      l=0
      do 10 ic=1,ncent+nghostcent
      i=iwctype(ic)

      if (iabs(n1s(i)).lt.1) goto 20
      ju=n1s(i)
      do 25 j=1,ju
      l=l+1
      anorm(l)=dsqrt(zex(l,iwf)**3/pi)
   25 continue
      ju=-ju
      do 27 j=1,ju
      l=l+1
      anorm(l)=one
   27 continue
   20 continue

      if (iabs(n2s(i)).lt.1) goto 30
      ju=n2s(i)
      do 35 j=1,ju
      l=l+1
      anorm(l)=dsqrt(zex(l,iwf)**5/(three*pi))
   35 continue
      ju=-ju
      do 37 j=1,ju
      l=l+1
      anorm(l)=one
   37 continue
   30 continue

      do 40 k=1,3
      if (iabs(n2p(k,i)).lt.1) goto 40
      ju=n2p(k,i)
      do 45 j=1,ju
      l=l+1
      anorm(l)=dsqrt(zex(l,iwf)**5/pi)
   45 continue
      ju=-ju
      do 47 j=1,ju
      l=l+1
      anorm(l)=one
   47 continue
   40 continue

      if (iabs(n3s(i)).lt.1) goto 100
      ju=n3s(i)
      do 110 j=1,ju
      l=l+1
      anorm(l)=third*dsqrt(two*zex(l,iwf)**7/(five*pi))
  110 continue
      ju=-ju
      do 115 j=1,ju
      l=l+1
      anorm(l)=one
  115 continue
  100 continue

      do 120 k=1,3
      if (iabs(n3p(k,i)).lt.1) goto 130
      ju=n3p(k,i)
      do 140 j=1,ju
      l=l+1
      anorm(l)=third*dsqrt(six*zex(l,iwf)**7/(five*pi))
  140 continue
      ju=-ju
      do 145 j=1,ju
      l=l+1
      anorm(l)=one
  145 continue
  130 continue
  120 continue

      if (iabs(n3dzr(i)).lt.1) goto 150
      ju=n3dzr(i)
      do 160 j=1,ju
      l=l+1
      anorm(l)=third*dsqrt(two*zex(l,iwf)**7/pi)
  160 continue
      ju=-ju
      do 165 j=1,ju
      l=l+1
      anorm(l)=one
  165 continue
  150 continue
      if (iabs(n3dx2(i)).lt.1) goto 170
      ju=n3dx2(i)
      do 180 j=1,ju
      l=l+1
      anorm(l)=third*dsqrt(two*zex(l,iwf)**7/pi)
  180 continue
      ju=-ju
      do 185 j=1,ju
      l=l+1
      anorm(l)=one
  185 continue
  170 continue
      if (iabs(n3dxy(i)).lt.1) goto 190
      ju=n3dxy(i)
      do 200 j=1,ju
      l=l+1
      anorm(l)=third*dsqrt(two*zex(l,iwf)**7/pi)
  200 continue
      ju=-ju
      do 205 j=1,ju
      l=l+1
      anorm(l)=one
  205 continue
  190 continue
      if (iabs(n3dxz(i)).lt.1) goto 210
      ju=n3dxz(i)
      do 220 j=1,ju
      l=l+1
      anorm(l)=third*dsqrt(two*zex(l,iwf)**7/pi)
  220 continue
      ju=-ju
      do 225 j=1,ju
      l=l+1
      anorm(l)=one
  225 continue
  210 continue
      if (iabs(n3dyz(i)).lt.1) goto 230
      ju=n3dyz(i)
      do 240 j=1,ju
      l=l+1
      anorm(l)=third*dsqrt(two*zex(l,iwf)**7/pi)
  240 continue
      ju=-ju
      do 245 j=1,ju
      l=l+1
      anorm(l)=one
  245 continue
  230 continue

      if (n4s(i).lt.1) goto 250
      do 260 j=1,n4s(i)
      l=l+1
      anorm(l)=third*dsqrt(zex(l,iwf)**9/(seven*five*pi))
  260 continue
  250 continue

      do 270 k=1,3
      if (n4p(k,i).lt.1) goto 280
      do 290 j=1,n4p(k,i)
      l=l+1
      anorm(l)=third*dsqrt(three*zex(l,iwf)**9/(seven*five*pi))
  290 continue
  280 continue
  270 continue
   10 continue

      if(iflag.gt.0) return

      do 300 iorb=1,norb+nadorb
        do 300 j=1,nbasis
  300     coef(j,iorb,iwf)=coef(j,iorb,iwf)*anorm(j)

      return
      end
