      module pw_find_image
      contains
      subroutine check_lattice(rlatt,cutr,isim_cell)
c Written by Cyrus Umrigar
c Checks to see if the lattice vectors specified are the smallest
c ones possible.  This is necessary for the simple heuristic algorithm
c Mathew Foulkes suggested to find the image particle (if any) that lies in the
c inscribing sphere of the nearest Wigner Seitz cell.  Also set cutjas to
c 1/2 the shortest simulation cell lattice vector (inscribing sphere radius.
c If the input cutjas is smaller, it will be reset to the smaller value in read_input.
c However, cutjas is used not only for r_ee but also r_en, so for that purpose
c we should use shortest primitive cell lattice vector or sum over atoms in sim cell.
c Warning:  I need to fix the above:
c Also return rlenmin to set cutr to 1/2 the shortest lattice vector.  I think that is
c good enough -- no need to use 1/2 the shortest perpendicular distance.

      use jaspar6, only: cutjas
      use precision_kinds, only: dp
      use contrl_file,    only: ounit
      implicit none

      integer :: i, i1, i2, i3, imax
      integer :: imin, isim_cell, k
      real(dp) :: cutr, rlen, rlenmax, rlenmin
      real(dp), dimension(3,3) :: rlatt
      real(dp), parameter :: eps = 1.d-12





      rlenmax=0
      rlenmin=9.d99
      do i=1,3
        rlen=0
        do k=1,3
          rlen=rlen+rlatt(k,i)**2
        enddo
        if (rlen.gt.rlenmax) then
          rlenmax=max(rlen,rlenmax)
          imax=i
        endif
        if (rlen.lt.rlenmin) then
          rlenmin=min(rlen,rlenmin)
          imin=i
        endif
      enddo
      rlenmax=sqrt(rlenmax)
      rlenmin=sqrt(rlenmin)
      cutr=rlenmin/2

c Warning: setting cutjas=rlenmin/2 for sim cell is OK for B terms, but not for A and C.
      if(isim_cell.eq.0) then
        write(ounit,'(''primitive  cell lattice vector'',i3,'' is longest ; length='',f8.3)')
     &  imax,rlenmax
        write(ounit,'(''primitive  cell lattice vector'',i3,'' is shortest; length='',f8.3)')
     &  imin,rlenmin
       else
        write(ounit,'(''simulation cell lattice vector'',i3,'' is longest ; length='',f8.3)')
     &  imax,rlenmax
        write(ounit,'(''simulation cell lattice vector'',i3,'' is shortest; length='',f8.3)')
     &  imin,rlenmin
        cutjas=rlenmin/2
      endif

c     if(cutjas.gt.rlenmin/2) then
c       write(ounit,'(''Warning: input cutjas > half shortest lattice vector;
c    &  cutjas reset from'',f9.5,'' to'',f9.5)') cutjas,rlenmin/2
c       cutjas=rlenmin/2
c     endif
c     cutjas=rlenmin/2

      do i1=-1,1
        do i2=-1,1
          do i3=-1,1
            if((imax.eq.1.and.i1.ne.0).or.(imax.eq.2.and.i2.ne.0)
     &      .or.(imax.eq.3.and.i3.ne.0)) then
              rlen=0
              do k=1,3
                rlen=rlen+(i1*rlatt(k,1)+i2*rlatt(k,2)+i3*rlatt(k,3))**2
              enddo
              rlen=sqrt(rlen)
              if (rlen.lt.rlenmax-eps) then
                write(ounit,*) 'found shorter lattice vector'
                write(ounit,'(''i1,i2,i3,rlen='',3i3,f8.3)') i1,i2,i3,rlen
                write(ounit,'(''new rlatt='',3f8.3)')
     &          i1*rlatt(1,1)+i2*rlatt(1,2)+i3*rlatt(1,3),
     &          i1*rlatt(2,1)+i2*rlatt(2,2)+i3*rlatt(2,3),
     &          i1*rlatt(3,1)+i2*rlatt(3,2)+i3*rlatt(3,3)
                call fatal_error ('one can find shorter lattice vectors: see check_lattice')
              endif
            endif
          enddo
        enddo
      enddo

      return
      end
c-----------------------------------------------------------------------

      subroutine reduce_sim_cell(r,rlatt,rlatt_inv)
c Written by Cyrus Umrigar
c For any electron position, replace it by the equivalent
c position in the simulation cell centered at the origin.
c r       = position in cartesian coords
c r_basis = position in lattice coords
c rlatt   = lattice vectors
c r       = rlatt * r_basis
c r_basis = rlatt_inv * r

      use grid3d_param, only: origin
      use precision_kinds, only: dp
      implicit none

      integer :: i, k

      real(dp), dimension(3) :: r
      real(dp), dimension(3) :: r_basis
      real(dp), dimension(3,3) :: rlatt
      real(dp), dimension(3,3) :: rlatt_inv



c Find vector in basis coordinates
      do k=1,3
        r_basis(k)=0
        do i=1,3
          r_basis(k)=r_basis(k)+rlatt_inv(k,i)*r(i)
        enddo
        r_basis(k)=r_basis(k)-nint(r_basis(k))
      enddo

c     write(ounit,'(''r_basis'',9f9.4)') r_basis

c Convert back to cartesian coodinates
      do k=1,3
        r(k)=0
        do i=1,3
          r(k)=r(k)+rlatt(k,i)*r_basis(i)
        enddo
      enddo

      return
      end
c-----------------------------------------------------------------------

      subroutine find_sim_cell(r,rlatt_inv,r_basis,i_basis)
c Written by Cyrus Umrigar
c For any electron position, find its lattice coordinates
c r       = position in cartesian coords
c r_basis = position in lattice coords
c i_basis = which simulation cell it is in
c rlatt   = lattice vectors
c r       = rlatt * r_basis
c r_basis = rlatt_inv * r

      use precision_kinds, only: dp
      implicit none

      integer :: i, k
      integer, dimension(3) :: i_basis

      real(dp), dimension(3) :: r
      real(dp), dimension(3) :: r_basis
      real(dp), dimension(3,3) :: rlatt_inv


c Find vector in basis coordinates
      do k=1,3
        r_basis(k)=0
        do i=1,3
          r_basis(k)=r_basis(k)+rlatt_inv(k,i)*r(i)
        enddo
        i_basis(k)=nint(r_basis(k))
      enddo

      return
      end
c-----------------------------------------------------------------------

      subroutine find_image(r,rlatt,rlatt_inv)
c Written by Cyrus Umrigar
c For any vector (from one particle to another) it finds the
c image that is closest.
      use contrl_file,    only: ounit
      use precision_kinds, only: dp
      implicit none

      integer :: i, i1, i2, i3, k
      integer, dimension(3) :: i_sav
      integer, dimension(3) :: isign
      real(dp) :: r2, r_try2
      real(dp), dimension(3) :: r
      real(dp), dimension(3) :: r_basis
      real(dp), dimension(3,3) :: rlatt
      real(dp), dimension(3,3) :: rlatt_inv
      real(dp), dimension(3) :: r1_try
      real(dp), dimension(3) :: r2_try
      real(dp), dimension(3) :: r3_try


c Starting from a vector, which is a diff. of 2 vectors, each of which
c have been reduced to the central lattice cell, calculate
c a) its length
c b) sign along each of lattice directions

      r2=0
      do k=1,3
        r2=r2+r(k)**2
        r_basis(k)=0
        do i=1,3
          r_basis(k)=r_basis(k)+rlatt_inv(k,i)*r(i)
        enddo
        if(abs(r_basis(k)).gt.1.d0) write(ounit,'(''**Warning, abs(r_basis)>1'')')
        isign(k)=nint(sign(1.d0,r_basis(k)))
      enddo

      do k=1,3
        i_sav(k)=0
      enddo

c Check just 8, rather than 27, trapezoids
      do i1=0,isign(1),isign(1)
        do k=1,3
          r1_try(k)=r(k)-i1*rlatt(k,1)
        enddo
        do i2=0,isign(2),isign(2)
          do k=1,3
            r2_try(k)=r1_try(k)-i2*rlatt(k,2)
          enddo
          do i3=0,isign(3),isign(3)
            r_try2=0
            do k=1,3
              r3_try(k)=r2_try(k)-i3*rlatt(k,3)
              r_try2=r_try2+r3_try(k)**2
            enddo
          if(r_try2.lt.r2) then
            i_sav(1)=i1
            i_sav(2)=i2
            i_sav(3)=i3
            r2=r_try2
          endif
          enddo
        enddo
      enddo

c Replace r by its shortest image
      do i=1,3
        do k=1,3
          r(k)=r(k)-i_sav(i)*rlatt(k,i)
        enddo
      enddo

c     write(ounit,'(''rnew'',9f10.5)') (r(k),k=1,3),sqrt(r2)

c debug
c     r2_tmp=0
c     do 80 k=1,3
c  80  r2_tmp=r2_tmp+r(k)**2
c     if(r2_tmp.ne.r2) write(ounit,'(''r2,r2_tmp'',3d12.4)') r2,r2_tmp,r2-r2_tmp

      return
      end
c-----------------------------------------------------------------------

      subroutine find_image2(r,rlatt,r_basis1,r_basis2,i_basis1,i_basis2)
c Written by Cyrus Umrigar
c For any electron positions in lattice coordinates, it finds the
c image that is closest.
c Needs precomputed r_basis1,r_basis2,i_basis1,i_basis2.

      use precision_kinds, only: dp
      implicit none

      integer :: i, i1, i2, i3, k
      integer, dimension(3) :: i_basis1
      integer, dimension(3) :: i_basis2
      integer, dimension(3) :: i_sav
      integer, dimension(3) :: isign
      real(dp) :: r2, r_try2
      real(dp), dimension(3) :: r
      real(dp), dimension(3,3) :: rlatt
      real(dp), dimension(3) :: r_basis1
      real(dp), dimension(3) :: r_basis2
      real(dp), dimension(3) :: r1_try
      real(dp), dimension(3) :: r2_try
      real(dp), dimension(3) :: r3_try


c Find length of original vector and sign along each of lattice directions
      r2=0
      do k=1,3
      r2=r2+r(k)**2
        isign(k)=int(sign(1.d0,r_basis2(k)-r_basis1(k)-i_basis2(k)+i_basis1(k)))
      enddo

      do k=1,3
        i_sav(k)=0
      enddo

c Check just 8, rather than 27, trapezoids (not needed for orthorhombic lattice)
      do i1=0,isign(1),isign(1)
        do k=1,3
          r1_try(k)=r(k)-rlatt(k,1)*(i1+i_basis2(1)-i_basis1(1))
        enddo
        do i2=0,isign(2),isign(2)
          do k=1,3
            r2_try(k)=r1_try(k)-rlatt(k,2)*(i2+i_basis2(2)-i_basis1(2))
          enddo
          do i3=0,isign(3),isign(3)
            r_try2=0
            do k=1,3
              r3_try(k)=r2_try(k)-rlatt(k,3)*(i3+i_basis2(3)-i_basis1(3))
              r_try2=r_try2+r3_try(k)**2
            enddo
          if(r_try2.lt.r2) then
            i_sav(1)=i1+i_basis2(1)-i_basis1(1)
            i_sav(2)=i2+i_basis2(2)-i_basis1(2)
            i_sav(3)=i3+i_basis2(3)-i_basis1(3)
            r2=r_try2
          endif
          enddo
        enddo
      enddo

c Replace r by its shortest image
      do i=1,3
        do k=1,3
          r(k)=r(k)-rlatt(k,i)*i_sav(i)
        enddo
      enddo

      return
      end
c-----------------------------------------------------------------------

      subroutine find_image3(r,rnorm)
c Written by Cyrus Umrigar
c For any vector r (from one particle to another) it replaces the vector
c by its closest image and finds its norm

      use periodic, only: rlatt, rlatt_inv
      use precision_kinds, only: dp
      implicit none

      integer :: i, i1, i2, i3, k
      integer, dimension(3) :: i_sav
      integer, dimension(3) :: isign
      real(dp) :: r2, r_try2, rnorm
      real(dp), dimension(3) :: r
      real(dp), dimension(3) :: r_basis
      real(dp), dimension(3) :: r1_try
      real(dp), dimension(3) :: r2_try
      real(dp), dimension(3) :: r3_try
      real(dp), dimension(3) :: rsav





c Warning: tempor
      do k=1,3
        rsav(k)=r(k)
      enddo

c a) reduce vector to central cell by expressing vector in lattice coordinates and
c    removing nint of it in each direction
c b) sign along each of lattice directions of vector reduced to central cell
      do k=1,3
        r_basis(k)=0
        do i=1,3
          r_basis(k)=r_basis(k)+rlatt_inv(k,i)*r(i)
        enddo
        r_basis(k)=r_basis(k)-nint(r_basis(k))
        isign(k)=nint(sign(1.d0,r_basis(k)))
      enddo

c Convert back to cartesian coodinates and find squared length
      r2=0
      do k=1,3
        r(k)=0
        do i=1,3
          r(k)=r(k)+rlatt(k,i)*r_basis(i)
        enddo
        r2=r2+r(k)**2
      enddo

      do k=1,3
        i_sav(k)=0
      enddo

c Check just 8, rather than 27, trapezoids (not needed for orthorhombic lattice)
      do i1=0,isign(1),isign(1)
c     do 60 i1=-1,1,1
        do k=1,3
          r1_try(k)=r(k)-i1*rlatt(k,1)
        enddo
        do i2=0,isign(2),isign(2)
c       do 60 i2=-1,1,1
          do k=1,3
            r2_try(k)=r1_try(k)-i2*rlatt(k,2)
          enddo
          do i3=0,isign(3),isign(3)
c         do 60 i3=-1,1,1
            r_try2=0
            do k=1,3
              r3_try(k)=r2_try(k)-i3*rlatt(k,3)
              r_try2=r_try2+r3_try(k)**2
            enddo
          if(r_try2.lt.r2) then
            i_sav(1)=i1
            i_sav(2)=i2
            i_sav(3)=i3
            r2=r_try2
          endif
          enddo
        enddo
      enddo

c Replace r by its shortest image
      rnorm=0
      do k=1,3
        do i=1,3
          r(k)=r(k)-rlatt(k,i)*i_sav(i)
        enddo
        rnorm=rnorm+r(k)**2
      enddo
      rnorm=sqrt(rnorm)

c     if(rnorm.gt.5.d0) write(ounit,'(''long'',6i2,10f8.4)')
c    &(isign(k),k=1,3),(i_sav(k),k=1,3),rnorm,(r(k),k=1,3),(rsav(k),k=1,3),(r_basis(k),k=1,3)

      return
      end
c-----------------------------------------------------------------------

      subroutine find_image4(rshift,r,rnorm)
c Written by Cyrus Umrigar
c For any vector r (from one particle to another) it replaces the vector
c by its closest image and finds its norm and the shift needed.

      use periodic, only: rlatt, rlatt_inv
      use precision_kinds, only: dp
      implicit none

      integer :: i, i1, i2, i3, k
      integer, dimension(3) :: i_sav
      integer, dimension(3) :: isign
      real(dp) :: r2, r_try2, rnorm
      real(dp), dimension(3) :: r
      real(dp), dimension(3) :: r_basis
      real(dp), dimension(3) :: rshift
      real(dp), dimension(3) :: r1_try
      real(dp), dimension(3) :: r2_try
      real(dp), dimension(3) :: r3_try





c a) reduce vector to central cell by expressing vector in lattice coordinates and
c    removing nint of it in each direction
c b) sign along each of lattice directions of vector reduced to central cell
c Note: rhift is just a work array here; calculated for real only at end.
      do k=1,3
        r_basis(k)=0
        do i=1,3
          r_basis(k)=r_basis(k)+rlatt_inv(k,i)*r(i)
        enddo
        rshift(k)=r_basis(k)-nint(r_basis(k))
        isign(k)=nint(sign(1.d0,rshift(k)))
      enddo

c Convert back to cartesian coodinates and find squared length
      r2=0
      do k=1,3
        r(k)=0
        do i=1,3
          r(k)=r(k)+rlatt(k,i)*rshift(i)
        enddo
        r2=r2+r(k)**2
      enddo

      do k=1,3
        i_sav(k)=0
      enddo

c Check just 8, rather than 27, trapezoids (not needed for orthorhombic lattice)
      do i1=0,isign(1),isign(1)
        do k=1,3
          r1_try(k)=r(k)-i1*rlatt(k,1)
        enddo
        do i2=0,isign(2),isign(2)
          do k=1,3
            r2_try(k)=r1_try(k)-i2*rlatt(k,2)
          enddo
          do i3=0,isign(3),isign(3)
            r_try2=0
            do k=1,3
              r3_try(k)=r2_try(k)-i3*rlatt(k,3)
              r_try2=r_try2+r3_try(k)**2
            enddo
          if(r_try2.lt.r2) then
            i_sav(1)=i1
            i_sav(2)=i2
            i_sav(3)=i3
            r2=r_try2
          endif
          enddo
        enddo
      enddo

c Replace r by its shortest image and calculate rshift
      rnorm=0
      do k=1,3
        rshift(k)=0
        do i=1,3
          rshift(k)=rshift(k)+rlatt(k,i)*(nint(r_basis(i))+i_sav(i))
          r(k)=r(k)-rlatt(k,i)*i_sav(i)
        enddo
        rnorm=rnorm+r(k)**2
      enddo
      rnorm=sqrt(rnorm)

      return
      end
      end module
