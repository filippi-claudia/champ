c-----------------------------------------------------------------------
      subroutine chlsky(a,n,np,ierr)
c chlsky: purpose: cholesky decomposition of a and determinant
c in: matrix a of order n stored with physical dimension np
c out: lower triangular matrix stored in lower portion of a
c note: lower triangular portion of original a is overwritten
      implicit real*8(a-h,o-z)
      dimension a(np,np)
      parameter (ZERO=0,ONE=1)

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
           write(6,'(''Warning: '',i2,'' element of a is <0'',d9.2)') j,a(j,j)
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

c     write(6,'(''diag_prod,det,deta='',9d12.4)') diag_prod,det,deta

      return
      end
c-----------------------------------------------------------------------
      subroutine lxb(a,n,np,b)
c lxb: purpose: solve equation Lx=b, for lower triangular matrix L=a
c Golub and van Loan: algorithm 3.1.3
c in:  a = matrix of order n with physical dimension np
c      b = vector of order n
c out: b = solution of eqs.; overwites original b
      implicit real*8(a-h,o-z)
      dimension a(np,np),b(np)
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
      implicit real*8(a-h,o-z)
      dimension a(np,np),b(np)
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
      subroutine svdchk(a,np,n,b,x,nk,ftol,iprt)
c solve system of linear equations via singular value decomposition
c solve a*x=b using SVD
c Test result
      implicit double precision (a-h,o-z)
      include 'svd.h'
      dimension a(np,n)
      dimension b(n)
      dimension x(n)
c
      dimension a_save(MXDIM,MXDIM)
c
      if(n.gt.MXDIM)then
       print *,'n >  MXDIM ',n, MXDIM 
       stop 'dimension too small in svdchk!'
      endif
      do i=1,n
       do j=1,n
        a_save(j,i)=a(j,i)
       enddo
      enddo
c
      call linsvd(a,np,n,b,x,nk,ftol,iprt)
c
      rmax=0.d0
      cn1=0.d0
      cn2=0.d0
      do i=1,n
       s=0.d0
       do j=1,n
        s=s+a_save(j,i)*x(j)
       enddo
       r=s-b(i)
       if(abs(r).gt.rmax)then
        rmax=r
        imax=i
       endif
       cn1=cn1+abs(r)
       cn2=cn2+r**2
      enddo
      cn2=sqrt(cn2)
      do i=1,n
       do j=1,n
        a(j,i)=a_save(j,i)
       enddo
      enddo
      if((iprt.gt.0).or.(rmax.gt.ftol)) then
       write(6,'(/,''svdchk: max deviation  '',e20.8,i5)') rmax,imax
       write(6,'(''svdchk: residual norms '',2e20.8)') cn1,cn2
      endif

      end
c-----------------------------------------------------------------------
      subroutine linsvd(a,np,n,b,x,nk,ftol,iprt)
c solve a*x=b using SVD
c Input: matrix a  : replaced with a^-1 on output 
c        phys dim np, real dim n
c        rhs vector b 
c Output:
c        vector x
c        nk : number of eigenvalues kept
c Threshold:  ftol
      implicit double precision (a-h,o-z)
      include 'svd.h'
c     
      dimension a(np,n)
      dimension b(n)
      dimension x(n)
C
      dimension evect(MXDIM,MXDIM)
      dimension eval(MXDIM)
      dimension tmp(MXDIM**2+MBUF) 
C
      ierr=0

c     ---symmetrise a ---
      qmx=0.d0
      imx=0
      jmx=0
      do i=2,n
       do j=1,i-1
        d=dabs(a(j,i)-a(i,j))
        if(d.gt.qmx) then
         qmx=d
         imx=i
         jmx=j
        endif
        a(j,i)=(a(j,i)+a(i,j))/2.d0
       enddo
      enddo
      if(iprt.gt.1) then
       print 68,qmx,imx,jmx
 68    format('linsvd: matrix symmetrised, max. dev. ',e22.12,' at ',2i4
     $      ) 
      endif
      if(qmx.gt.ftol) then
       call warn('linsvd: non-symmetric matrix !')
      endif

      do i=1,n
       do j=1,n
        evect(j,i)=a(j,i)
       enddo
      enddo
      lwork=n**2+MBUF
      call  DSYEV('V', 'U', n, evect(1,1), MXDIM, eval(1), tmp(1),
     &     lwork,info)
      if(info.ne.0) then
       print *,'Error in DSYEV: info = ',info
       call fatal('linsvd: can not diagonalise matrix')
      endif

 69   format('linsvd: threshold for keeping eigenvalues is ',e22.12)
 70   format('linsvd: eigenvalues (',i4,')')
 71   format('_ev_  ',i4,e22.12)
 72   format('linsvd: eigenvectors')
 73   format('_evec_ ',e22.12,' : ',100f16.8)
C
C     --- inversion ----
      icount=0
      do i=1,n
       if(iprt.gt.2)then
        print *,'SVD EV ',i,eval(i)
       endif
       if(dabs(eval(i)).gt.ftol) then
        eval(i) = 1.d0 / eval(i)
        icount = icount + 1
       else
        eval(i) = 0.d0
       endif
      enddo
      if(iprt.gt.1) then
       print 74,icount,n
      endif
      nk=icount
 74   format('linsvd: keeping ',i5,' eigenvalues out of ',i5)

      do i=1,n
       do j=1,n
        s=0.d0
        do k=1,n
         s = s + evect(j,k) * eval(k) * evect(i,k)
        enddo
        a(j,i) = s
       enddo
      enddo

      do i=1,n
       s=0.d0
       do j=1,n
        s = s + a(j,i)*b(j)
       enddo
       x(i) = s
      enddo

      end 
c-----------------------------------------------------------------------
      subroutine ipagt2(p,iorb,ia,ib)
      implicit double precision (a-h,o-z)
c alpha/beta occupation from string
      character p*(*)
      character x
      x=p(iorb:iorb)
      if(x.eq.'0') then
       ia=0
       ib=0
      elseif(x.eq.'2') then
       ia=1
       ib=1
      elseif(x.eq.'+') then
       ia=1
       ib=0
      elseif(x.eq.'-') then
       ia=0
       ib=1
      else
       print *,'ipaget: illegal occupation pattern'
       stop 'ERROR in ipaget'
      endif
      end
c----------------------------------------------------------------------
      subroutine sort(n,arr,indx)
C     index sort from numerical recepires
C     sorts in ascending order
      implicit double precision(a-h,o-z)
      integer n,indx(n),m,nstack
      double precision arr(n)
      parameter (m=7,nstack=50)
      integer i,indxt,ir,itemp,j,jstack,k,l,istack(nstack)
      double precision a
      do 11 j=1,n
        indx(j)=j
11    continue
      jstack=0
      l=1
      ir=n
1     if(ir-l.lt.m)then
        do 13 j=l+1,ir
          indxt=indx(j)
          a=arr(indxt)
          do 12 i=j-1,l,-1
            if(arr(indx(i)).le.a)goto 2
            indx(i+1)=indx(i)
12        continue
          i=l-1
2         indx(i+1)=indxt
13      continue
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
        if(jstack.gt.nstack)pause 'nstack too small in indexx'
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
      subroutine idxsort(n,arr,indx)
C$    index sort from numerical recepires
C$    sorts in ascending order
      implicit double precision(a-h,o-z)
      integer n,indx(n),m,nstack
      double precision arr(n)
      parameter (m=7,nstack=50)
      integer i,indxt,ir,itemp,j,jstack,k,l,istack(nstack)
      double precision a
      do 11 j=1,n
        indx(j)=j
11    continue
      jstack=0
      l=1
      ir=n
1     if(ir-l.lt.m)then
        do 13 j=l+1,ir
          indxt=indx(j)
          a=arr(indxt)
          do 12 i=j-1,l,-1
            if(arr(indx(i)).le.a)goto 2
            indx(i+1)=indx(i)
12        continue
          i=l-1
2         indx(i+1)=indxt
13      continue
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
        if(jstack.gt.nstack)pause 'nstack too small in indexx'
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

