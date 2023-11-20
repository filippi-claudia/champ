module cuspexact4_mod
contains
      subroutine cuspexact4(iprin,iadiag)
! Written by Cyrus Umrigar
      use cuspmat4, only: d,iwc4,nterms
      use jastrow, only: c,nordc
      use precision_kinds, only: dp
      use system,  only: nctype
      implicit none

      integer :: i, iadiag, iprin, it, j
      real(dp) :: sum



! The last 2 columns are what we care about in the foll. table
!------------------------------------------------------------------------------
! ord  # of   cumul.    # terms  # terms   # 3-body  Cumul. #      Cumul # indep
!      terms  # terms   even t   odd t      terms    3-body terms  3-body terms
!  n  (n+1)* (n^3+5n)/6         int((n+1)/2            nterms
!    (n+2)/2  +n^2+n           *int((n+2)/2
!------------------------------------------------------------------------------
!  1     3       3        2         1          0         0              0
!  2     6       9        4         2          2         2              0
!  3    10      19        6         4          4         6              2
!  4    15      34        9         6          7        13              7
!  5    21      55       12         9         10        23             15
!  6    28      83       16        12         14        37             27
!  7    36     119       20        16         18        55             43
!------------------------------------------------------------------------------

! Dependent coefs. fixed by e-e and e-n cusp conditions resp. are;
! order:   2  3  4  5  6  7  2  3  4  5  6  7
! coefs:   1  4 10 19 32 49  2  6 12 22 35 53

! So the terms varied for a 5th, 6th order polynomial are:
!    3   5   7 8 9    11    13 14 15 16 17 18    20 21    23 (iwjasc(iparm),iparm=1,nparmc)
!    3   5   7 8 9    11    13 14 15 16 17 18    20 21    23 24 25 26 27 28 29 30 31    33 34    36 37
!                                                                 (iwjasc(iparm),iparm=1,nparmc)


! All the dependent variables, except one (the one from the 2nd order
! e-n cusp) depend only on independent variables.  On the other hand
! the one from the 2nd order e-n cusp depends only on other dependent
! variables.

      do it=1,nctype

! Set dep. variables from e-e cusp
        do i=1,nordc-1
          sum=0
          do j=1,nterms
            if(j.ne.iwc4(i)) sum=sum+d(i,j)*c(j,it,iadiag)
          enddo
          c(iwc4(i),it,iadiag)=-sum/d(i,iwc4(i))
        enddo

! Set dep. variables from 3rd and higher order e-n cusp
        do i=nordc+1,2*(nordc-1)
          sum=0
          do j=1,nterms
            if(j.ne.iwc4(i)) sum=sum+d(i,j)*c(j,it,iadiag)
          enddo
          c(iwc4(i),it,iadiag)=-sum/d(i,iwc4(i))
        enddo

! Set dep. variables from 2nd order e-n cusp
        if(nordc.gt.1) then
          i=nordc
          sum=0
          do j=1,nterms
            if(j.ne.iwc4(i)) sum=sum+d(i,j)*c(j,it,iadiag)
          enddo
          c(iwc4(i),it,iadiag)=-sum/d(i,iwc4(i))
        endif

      enddo

      return
      end
end module
