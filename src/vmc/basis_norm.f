      subroutine basis_norm(iwf,anorm,iflag)
c Written by Cyrus Umrigar and Claudia Filippi, starting from Kevin Schmidt routine
c Set normalization of basis fns.
      use atom, only: iwctype, ncent

      use ghostatom, only: nghostcent
      use const, only: pi
      use coefs, only: coef, nbasis, norb
      use basis, only: zex, n1s, n2s, n2p, n3s, n3p, n3dzr, n3dx2, n3dxy, n3dxz, n3dyz
      use basis, only: n4s, n4p
      use orbval, only: nadorb

      use precision_kinds, only: dp
      implicit none

      integer :: i, iabs, ic, iflag, iorb
      integer :: iwf, j, ju, k
      integer :: l

      real(dp), dimension(*) :: anorm
      real(dp), parameter :: one = 1.d0
      real(dp), parameter :: two = 2.d0
      real(dp), parameter :: three = 3.d0
      real(dp), parameter :: five = 5.d0
      real(dp), parameter :: six = 6.d0
      real(dp), parameter :: seven = 7.d0
      real(dp), parameter :: third = 1.d0/3.d0






      do j=1,nbasis
        anorm(j)=one
      enddo

      l=0
      do ic=1,ncent+nghostcent
      i=iwctype(ic)

      if (iabs(n1s(i)).lt.1) goto 20
      ju=n1s(i)
      do j=1,ju
      l=l+1
      anorm(l)=dsqrt(zex(l,iwf)**3/pi)
      enddo
      ju=-ju
      do j=1,ju
      l=l+1
      anorm(l)=one
      enddo
   20 continue

      if (iabs(n2s(i)).lt.1) goto 30
      ju=n2s(i)
      do j=1,ju
      l=l+1
      anorm(l)=dsqrt(zex(l,iwf)**5/(three*pi))
      enddo
      ju=-ju
      do j=1,ju
      l=l+1
      anorm(l)=one
      enddo
   30 continue

      do k=1,3
      if (iabs(n2p(k,i)).lt.1) goto 40
      ju=n2p(k,i)
      do j=1,ju
      l=l+1
      anorm(l)=dsqrt(zex(l,iwf)**5/pi)
      enddo
      ju=-ju
      do j=1,ju
      l=l+1
      anorm(l)=one
      enddo
   40 continue
      enddo

      if (iabs(n3s(i)).lt.1) goto 100
      ju=n3s(i)
      do j=1,ju
      l=l+1
      anorm(l)=third*dsqrt(two*zex(l,iwf)**7/(five*pi))
      enddo
      ju=-ju
      do j=1,ju
      l=l+1
      anorm(l)=one
      enddo
  100 continue

      do k=1,3
      if (iabs(n3p(k,i)).lt.1) goto 130
      ju=n3p(k,i)
      do j=1,ju
      l=l+1
      anorm(l)=third*dsqrt(six*zex(l,iwf)**7/(five*pi))
      enddo
      ju=-ju
      do j=1,ju
      l=l+1
      anorm(l)=one
      enddo
  130 continue
      enddo

      if (iabs(n3dzr(i)).lt.1) goto 150
      ju=n3dzr(i)
      do j=1,ju
      l=l+1
      anorm(l)=third*dsqrt(two*zex(l,iwf)**7/pi)
      enddo
      ju=-ju
      do j=1,ju
      l=l+1
      anorm(l)=one
      enddo
  150 continue
      if (iabs(n3dx2(i)).lt.1) goto 170
      ju=n3dx2(i)
      do j=1,ju
      l=l+1
      anorm(l)=third*dsqrt(two*zex(l,iwf)**7/pi)
      enddo
      ju=-ju
      do j=1,ju
      l=l+1
      anorm(l)=one
      enddo
  170 continue
      if (iabs(n3dxy(i)).lt.1) goto 190
      ju=n3dxy(i)
      do j=1,ju
      l=l+1
      anorm(l)=third*dsqrt(two*zex(l,iwf)**7/pi)
      enddo
      ju=-ju
      do j=1,ju
      l=l+1
      anorm(l)=one
      enddo
  190 continue
      if (iabs(n3dxz(i)).lt.1) goto 210
      ju=n3dxz(i)
      do j=1,ju
      l=l+1
      anorm(l)=third*dsqrt(two*zex(l,iwf)**7/pi)
      enddo
      ju=-ju
      do j=1,ju
      l=l+1
      anorm(l)=one
      enddo
  210 continue
      if (iabs(n3dyz(i)).lt.1) goto 230
      ju=n3dyz(i)
      do j=1,ju
      l=l+1
      anorm(l)=third*dsqrt(two*zex(l,iwf)**7/pi)
      enddo
      ju=-ju
      do j=1,ju
      l=l+1
      anorm(l)=one
      enddo
  230 continue

      if (n4s(i).lt.1) goto 250
      do j=1,n4s(i)
      l=l+1
      anorm(l)=third*dsqrt(zex(l,iwf)**9/(seven*five*pi))
      enddo
  250 continue

      do k=1,3
      if (n4p(k,i).lt.1) goto 280
      do j=1,n4p(k,i)
      l=l+1
      anorm(l)=third*dsqrt(three*zex(l,iwf)**9/(seven*five*pi))
      enddo
  280 continue
      enddo
      enddo

      if(iflag.gt.0) return

      ! potentially resize the 2nd dim of
      ! coefs to norb+nadorb
    !   call resize_tensor(coef, norb+nadorb, 2)

      do iorb=1,norb+nadorb
        do j=1,nbasis
          coef(j,iorb,iwf)=coef(j,iorb,iwf)*anorm(j)
        enddo
      enddo

      return
      end
