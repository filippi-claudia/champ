! Modified by Edgar Josue Landinez Borda
      module find_pimage
      contains

      subroutine check_lattice(rlatt,cutr,isim_cell)
! Written by Cyrus Umrigar
! Checks to see if the lattice vectors specified are the smallest
! ones possible.  This is necessary for the simple heuristic algorithm
! Mathew Foulkes suggested to find the image particle (if any) that lies in the
! inscribing sphere of the nearest Wigner Seitz cell.  
! Warning:  I need to fix the above:
! Also return rlenmin to set cutr to 1/2 the shortest lattice vector.  I think that is
! good enough -- no need to use 1/2 the shortest perpendicular distance.

      use precision_kinds, only: dp
      use contrl_file,    only: ounit
      use error, only: fatal_error

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

      write(ounit,'(''primitive  cell lattice vector'',i3,'' is longest ; length='',f8.3)') &
      imax,rlenmax
      write(ounit,'(''primitive  cell lattice vector'',i3,'' is shortest; length='',f8.3)') &
      imin,rlenmin

      do i1=-1,1
        do i2=-1,1
          do i3=-1,1
            if((imax.eq.1.and.i1.ne.0).or.(imax.eq.2.and.i2.ne.0) &
            .or.(imax.eq.3.and.i3.ne.0)) then
              rlen=0
              do k=1,3
                rlen=rlen+(i1*rlatt(k,1)+i2*rlatt(k,2)+i3*rlatt(k,3))**2
              enddo
              rlen=sqrt(rlen)
              if (rlen.lt.rlenmax-eps) then
                write(ounit,*) 'found shorter lattice vector'
                write(ounit,'(''i1,i2,i3,rlen='',3i3,f8.3)') i1,i2,i3,rlen
                write(ounit,'(''new rlatt='',3f8.3)') &
                i1*rlatt(1,1)+i2*rlatt(1,2)+i3*rlatt(1,3), &
                i1*rlatt(2,1)+i2*rlatt(2,2)+i3*rlatt(2,3), &
                i1*rlatt(3,1)+i2*rlatt(3,2)+i3*rlatt(3,3)
                call fatal_error ('one can find shorter lattice vectors: see check_lattice')
              endif
            endif
          enddo
        enddo
      enddo

      return
      end
!-----------------------------------------------------------------------

      subroutine find_image_pbc(r,rnorm)
! Written by Edgar Landinez
! Simple algorithm for PBC minimum image convention
! to get the minimum distnace between two particles and it's norm


      use periodic, only: rlatt, rlatt_inv
      use precision_kinds, only: dp
      implicit none

      integer :: i, k
      real(dp) :: rnorm
      real(dp), dimension(3) :: s
      real(dp), dimension(3) :: r

!     minimum image in relative coordiantes space (a cube of length 1)
!     rlatt or rlatt is assumend to be rlatt=(a,b,c)
!     (a,b,c) the box vectors (the input should be always consistent with this )

      do k=1,3
        s(k)=0.d0
        do i=1,3
          s(k)=s(k)+rlatt_inv(k,i)*r(i)
        enddo
        s(k)=s(k)-nint(s(k))
      enddo

! resotring coordinates in real space
      do k=1,3
        r(k)=0.d0
        do i=1,3
          r(k)=r(k)+rlatt(k,i)*s(i)
       enddo
      enddo

! compute norm of the distance
      rnorm=0.d0
      do k=1,3
         rnorm=rnorm+(r(k)*r(k))
      enddo
      rnorm=dsqrt(rnorm)


      return
      end
!-----------------------------------------------------------------------
      end module
