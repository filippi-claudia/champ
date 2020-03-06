      subroutine cuspexact4(iprin,iadiag)
c Written by Cyrus Umrigar
      use atom, only: znuc, cent, pecent, iwctype, nctype, ncent

      use jaspar3, only: a, b, c, fck, nord, scalek

      implicit real*8(a-h,o-z)


      include 'vmc.h'
      include 'force.h'

c For Jastrow4 NEQSX=2*(MORDJ-1) is sufficient.
c For Jastrow3 NEQSX=2*MORDJ should be sufficient.
c I am setting NEQSX=6*MORDJ simply because that is how it was for
c Jastrow3 for reasons I do not understand.
c     parameter(NEQSX=2*(MORDJ-1),MTERMS=55)
      parameter(NEQSX=6*MORDJ,MTERMS=55)

c The last 2 columns are what we care about in the foll. table
c------------------------------------------------------------------------------
c ord  # of   cumul.    # terms  # terms   # 3-body  Cumul. #      Cumul # indep
c      terms  # terms   even t   odd t      terms    3-body terms  3-body terms
c  n  (n+1)* (n^3+5n)/6         int((n+1)/2            nterms
c    (n+2)/2  +n^2+n           *int((n+2)/2
c------------------------------------------------------------------------------
c  1     3       3        2         1          0         0              0
c  2     6       9        4         2          2         2              0
c  3    10      19        6         4          4         6              2
c  4    15      34        9         6          7        13              7
c  5    21      55       12         9         10        23             15
c  6    28      83       16        12         14        37             27
c  7    36     119       20        16         18        55             43
c------------------------------------------------------------------------------

c Dependent coefs. fixed by e-e and e-n cusp conditions resp. are;
c order:   2  3  4  5  6  7  2  3  4  5  6  7
c coefs:   1  4 10 19 32 49  2  6 12 22 35 53

c So the terms varied for a 5th, 6th order polynomial are:
c    3   5   7 8 9    11    13 14 15 16 17 18    20 21    23 (iwjasc(iparm),iparm=1,nparmc)
c    3   5   7 8 9    11    13 14 15 16 17 18    20 21    23 24 25 26 27 28 29 30 31    33 34    36 37 (iwjasc(iparm),iparm=1,nparmc)


c All the dependent variables, except one (the one from the 2nd order
c e-n cusp) depend only on independent variables.  On the other hand
c the one from the 2nd order e-n cusp depends only on other dependent
c variables.

      common /jaspar4/ a4(MORDJ1,MCTYPE,MWF),norda,nordb,nordc
      common /cuspmat4/ d(NEQSX,MTERMS),icusp(NEQSX),nterms

      do 100 it=1,nctype

c Set dep. variables from e-e cusp
        do 20 i=1,nordc-1
          sum=0
          do 10 j=1,nterms
   10       if(j.ne.icusp(i)) sum=sum+d(i,j)*c(j,it,iadiag)
   20     c(icusp(i),it,iadiag)=-sum/d(i,icusp(i))

c Set dep. variables from 3rd and higher order e-n cusp
        do 40 i=nordc+1,2*(nordc-1)
          sum=0
          do 30 j=1,nterms
   30       if(j.ne.icusp(i)) sum=sum+d(i,j)*c(j,it,iadiag)
   40     c(icusp(i),it,iadiag)=-sum/d(i,icusp(i))

c Set dep. variables from 2nd order e-n cusp
        if(nordc.gt.1) then
          i=nordc
          sum=0
          do 50 j=1,nterms
   50       if(j.ne.icusp(i)) sum=sum+d(i,j)*c(j,it,iadiag)
          c(icusp(i),it,iadiag)=-sum/d(i,icusp(i))
        endif

 100  continue

      return
      end
