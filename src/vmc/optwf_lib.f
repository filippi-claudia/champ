c-----------------------------------------------------------------------
      subroutine chlsky(a,n,np,ierr)
c chlsky: purpose: cholesky decomposition of a and determinant
c in: matrix a of order n stored with physical dimension np
c out: lower triangular matrix stored in lower portion of a
c note: lower triangular portion of original a is overwritten
      use precision_kinds, only: dp
      use contrl_file,    only: ounit
      implicit none

      integer :: ierr, ip, j, jm1, k
      integer :: n, np
      real(dp) :: det, deta, diag_prod, s, sum
      real(dp), dimension(np,np) :: a
      real(dp), parameter :: ZERO = 0
      real(dp), parameter :: ONE = 1

      diag_prod=1
      do j=1,n
         diag_prod=diag_prod*a(j,j)
      enddo

      det=1
      ierr=0
      do j=1,n
c        diag_prod=diag_prod*a(j,j)
         if(j.gt.1) then
            jm1=j-1
            do k=j,n
               sum=0
               do ip=1,jm1
                  sum=sum+a(k,ip)*a(j,ip)
               enddo
               a(k,j)=a(k,j)-sum
            enddo
         endif

         det=det*a(j,j)
         if(a(j,j).le.ZERO) then
           write(ounit,'(''Warning: '',i2,'' element of a is <0'',d9.2)') j,a(j,j)
           ierr=j
           return
         endif

         s=ONE/sqrt(a(j,j))
         do k=j,n
            a(k,j)=a(k,j)*s
         enddo
      enddo

      deta=1
      do j=1,n
         deta=deta*a(j,j)**2
      enddo

c     write(ounit,'(''diag_prod,det,deta='',9d12.4)') diag_prod,det,deta

      return
      end
c-----------------------------------------------------------------------
      subroutine lxb(a,n,np,b)
c lxb: purpose: solve equation Lx=b, for lower triangular matrix L=a
c Golub and van Loan: algorithm 3.1.3
c in:  a = matrix of order n with physical dimension np
c      b = vector of order n
c out: b = solution of eqs.; overwites original b

      use precision_kinds, only: dp

      implicit none

      integer :: j, k, n, np
      real(dp), dimension(np,np) :: a
      real(dp), dimension(np) :: b

      do j=1,n-1
         b(j)=b(j)/a(j,j)
         do k=j+1,n
            b(k)=b(k)-b(j)*a(k,j)
         enddo
      enddo
      b(n)=b(n)/a(n,n)
      return
      end
c-----------------------------------------------------------------------
      subroutine uxb(a,n,np,b)
c uxb: purpose: solve equation Ux=b, for upper triangular matrix U=a
c Golub and van Loan: algorithm 3.1.4
c in:  a = matrix of order n with physical dimension np
c      b = vector of order n
c out: b = solution of eqs.; overwites original b

      use precision_kinds, only: dp

      implicit none

      integer :: j, k, n, np
      real(dp), dimension(np,np) :: a
      real(dp), dimension(np) :: b

      do j=n,2,-1
         b(j)=b(j)/a(j,j)
         do k=1,j-1
            b(k)=b(k)-b(j)*a(k,j)
         enddo
      enddo
      b(1)=b(1)/a(1,1)
      return
      end
c-----------------------------------------------------------------------
      subroutine sort(n,arr,indx)
C     index sort from numerical recepires
C     sorts in ascending order

      use precision_kinds, only: dp

      implicit none

      integer, parameter :: m = 7
      integer, parameter :: nstack = 50
      integer n, indx(n)
      real(dp) :: a, arr(n)
      integer i,indxt,ir,itemp,j,jstack,k,l,istack(nstack)

      do j=1,n
        indx(j)=j
      enddo
      jstack=0
      l=1
      ir=n
1     if(ir-l.lt.m)then
        do j=l+1,ir
          indxt=indx(j)
          a=arr(indxt)
          do i=j-1,l,-1
            if(arr(indx(i)).le.a)goto 2
            indx(i+1)=indx(i)
          enddo
          i=l-1
2         indx(i+1)=indxt
        enddo
        if(jstack.eq.0)return
        ir=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
      else
        k=(l+ir)/2
        itemp=indx(k)
        indx(k)=indx(l+1)
        indx(l+1)=itemp
        if(arr(indx(l)).gt.arr(indx(ir)))then
          itemp=indx(l)
          indx(l)=indx(ir)
          indx(ir)=itemp
        endif
        if(arr(indx(l+1)).gt.arr(indx(ir)))then
          itemp=indx(l+1)
          indx(l+1)=indx(ir)
          indx(ir)=itemp
        endif
        if(arr(indx(l)).gt.arr(indx(l+1)))then
          itemp=indx(l)
          indx(l)=indx(l+1)
          indx(l+1)=itemp
        endif
        i=l+1
        j=ir
        indxt=indx(l+1)
        a=arr(indxt)
3       continue
          i=i+1
        if(arr(indx(i)).lt.a)goto 3
4       continue
          j=j-1
        if(arr(indx(j)).gt.a)goto 4
        if(j.lt.i)goto 5
        itemp=indx(i)
        indx(i)=indx(j)
        indx(j)=itemp
        goto 3
5       indx(l+1)=indx(j)
        indx(j)=indxt
        jstack=jstack+2
        if(jstack.gt.nstack) stop 'nstack too small in indexx'
        if(ir-i+1.ge.j-l)then
          istack(jstack)=ir
          istack(jstack-1)=i
          ir=j-1
        else
          istack(jstack)=j-1
          istack(jstack-1)=l
          l=i
        endif
      endif
      goto 1
      end
c-----------------------------------------------------------------------
