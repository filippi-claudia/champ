      module matinv_mod
      use error,   only: fatal_error
      interface !LAPACK interface
        SUBROUTINE dgetrf( M, N, A, LDA, IPIV, INFO )
!*  -- LAPACK computational routine --
!*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
          INTEGER            INFO, LDA, M, N
          INTEGER            IPIV( * )
          DOUBLE PRECISION   A( LDA, * )
        END SUBROUTINE
        SUBROUTINE DGETRI( N, A, LDA, IPIV, WORK, LWORK, INFO )
!*  -- LAPACK routine (version 3.1) --
!*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!*     November 2006
          INTEGER            INFO, LDA, LWORK, N
          INTEGER            IPIV( * )
          DOUBLE PRECISION   A( LDA, * ), WORK( * )
        END SUBROUTINE
      end interface
      contains
      subroutine matinv(a,nsub,determinant)
      use contrl_file, only: ounit
      use precision_kinds, only: dp
      use system,  only: nelec
      implicit none

      integer :: i, info, nsub
      integer, dimension(nsub) :: ipvt
      real(dp) :: aux, determinant, deti
      real(dp) :: ten
      real(dp), dimension(nsub, nsub) :: a
      real(dp), dimension(nelec) :: work
      real(dp), dimension(nsub, nsub) :: a_save

      real(dp), dimension(2) :: det
      real(dp), parameter :: eps = 10.d0**(-40)
      integer :: ii,jj
! routine to calculate inverse and determinant of matrix a
! assumed to be dimensioned a(nsub,nsub).
! the matrix a is replaced by its inverse.

      if(nsub.eq.0) return

      if(nsub.eq.1) then
        determinant=a(1,1)
        a(1,1)=1.0d0/a(1,1)
       elseif(nsub.eq.2) then
        determinant=a(1,1)*a(2,2)-a(1,2)*a(2,1)
        deti=1.d0/determinant
        aux=a(1,1)
        a(1,1)= a(2,2)*deti
        a(2,2)= aux   *deti
        a(2,1)=-a(2,1)*deti
        a(1,2)=-a(1,2)*deti
       else

        ! Save a copy of the matrix before dgetrf overwrites it
        a_save = a

        call dgetrf(nsub,nsub,a,nsub,ipvt,info)

        if(info.gt.0) then
          open(101, file='matinv_err.log', position='append', action='write')
          write(101,'(a)')       '================================================='
          write(101,'(a)')       'MATINV: dgetrf returned singular matrix'
          write(101,'(a,i5)')    'MATINV: u(k,k)=0 at k = ', info
          write(101,'(a,i5)')    'MATINV: matrix size nsub = ', nsub
          write(101,'(a)')       'MATINV: Input matrix (before factorization):'
          do ii=1,nsub
            write(101,'(100(ES14.6,1x))') (a_save(ii,jj), jj=1,nsub)
          end do
          write(101,'(a)')       'MATINV: Checking for identical rows:'
          do ii=1,nsub
            do jj=ii+1,nsub
              if(maxval(abs(a_save(ii,:)-a_save(jj,:))) .lt. 1.d-12) then
                write(101,'(a,i4,a,i4,a,ES14.6)') '  Rows ', ii, ' and ', jj, &
                  ' are identical (maxdiff=', maxval(abs(a_save(ii,:)-a_save(jj,:))), ')'
              end if
            end do
          end do
          write(101,'(a)')       'MATINV: Checking for identical columns:'
          do ii=1,nsub
            do jj=ii+1,nsub
              if(maxval(abs(a_save(:,ii)-a_save(:,jj))) .lt. 1.d-12) then
                write(101,'(a,i4,a,i4,a,ES14.6)') '  Cols ', ii, ' and ', jj, &
                  ' are identical (maxdiff=', maxval(abs(a_save(:,ii)-a_save(:,jj))), ')'
              end if
            end do
          end do
          write(101,'(a)')       'MATINV: Checking for zero rows/columns:'
          do ii=1,nsub
            if(maxval(abs(a_save(ii,:))) .lt. 1.d-20) then
              write(101,'(a,i4,a,ES14.6)') '  Row ', ii, ' is zero (max=', maxval(abs(a_save(ii,:))), ')'
            end if
            if(maxval(abs(a_save(:,ii))) .lt. 1.d-20) then
              write(101,'(a,i4,a,ES14.6)') '  Col ', ii, ' is zero (max=', maxval(abs(a_save(:,ii))), ')'
            end if
          end do
          write(101,'(a)')       'MATINV: Checking for NaN entries:'
          do ii=1,nsub
            do jj=1,nsub
              if(a_save(ii,jj) .ne. a_save(ii,jj)) then
                write(101,'(a,i4,a,i4)') '  NaN at row ', ii, ' col ', jj
              end if
            end do
          end do
          write(101,'(a)')       '================================================='
          call flush(101)
          close(101)
          call fatal_error('MATINV: info ne 0 in dgetrf')

        endif

        det(1) = 1.0d0
        det(2) = 0.0d0
        ten = 10.0d0
        do i = 1, nsub
          if (ipvt(i) .ne. i) det(1) = -det(1)
          det(1) = a(i,i)*det(1)
          if (det(1) .eq. 0.0d0) exit
          do
            if (dabs(det(1)) .ge. 1.0d0) exit
            det(1) = ten*det(1)
            det(2) = det(2) - 1.0d0
          end do
          do
            if (dabs(det(1)) .lt. ten) exit
            det(1) = det(1)/ten
            det(2) = det(2) + 1.0d0
          end do
        enddo

        determinant = det(1)*10.0**det(2)

        call dgetri(nsub,a,nsub,ipvt,work,nelec,info)

      endif

      return

      end
      end module
