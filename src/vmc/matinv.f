      subroutine matinv(a,nsub,determinant)
      use vmc, only: MELEC, MORB, MBASIS, MDET, MCENT, MCTYPE, MCTYP3X
      use vmc, only: NSPLIN, nrad, MORDJ, MORDJ1, MMAT_DIM, MMAT_DIM2, MMAT_DIM20
      use vmc, only: radmax, delri
      use vmc, only: NEQSX, MTERMS
      use vmc, only: MCENT3, NCOEF, MEXCIT
      implicit real*8 (a-h,o-z)

c routine to calculate inverse and determinant of matrix a
c assumed to be dimensioned a(nsub,nsub).
c the matrix a is replaced by its inverse.

      parameter (eps=10.d0**(-40))
      dimension a(nsub,nsub)
      dimension ipvt(MELEC),work(MELEC),work2(9),det(2)

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
c      elseif(nsub.eq.3) then
c       work2(1)=a(1,1)
c       work2(2)=a(2,1)
c       work2(3)=a(3,1)
c       work2(4)=a(1,2)
c       work2(5)=a(2,2)
c       work2(6)=a(3,2)
c       work2(7)=a(1,3)
c       work2(8)=a(2,3)
c       work2(9)=a(3,3)
c       determinant= work2(1)*(work2(5)*work2(9)-work2(6)*work2(8))
c    &              -work2(2)*(work2(4)*work2(9)-work2(6)*work2(7))
c    &              +work2(3)*(work2(4)*work2(8)-work2(5)*work2(7))
c       if(dabs(determinant).lt.eps) call fatal_error('MATINV: ndim eq 3')
c       deti=1.d0/determinant
c       a(1,1)= (work2(5)*work2(9)-work2(6)*work2(8))*deti
c       a(2,1)=-(work2(2)*work2(9)-work2(3)*work2(8))*deti
c       a(3,1)= (work2(2)*work2(6)-work2(5)*work2(3))*deti
c       a(1,2)=-(work2(4)*work2(9)-work2(6)*work2(7))*deti
c       a(2,2)= (work2(1)*work2(9)-work2(3)*work2(7))*deti
c       a(3,2)=-(work2(1)*work2(6)-work2(3)*work2(4))*deti
c       a(1,3)= (work2(4)*work2(8)-work2(5)*work2(7))*deti
c       a(2,3)=-(work2(1)*work2(8)-work2(2)*work2(7))*deti
c       a(3,3)= (work2(1)*work2(5)-work2(2)*work2(4))*deti
       else
        call dgetrf(nsub,nsub,a,nsub,ipvt,info)
        if(info.gt.0) then
          write(6,'(''MATINV: u(k,k)=0 with k= '',i5)') info
          call fatal_error('MATINV: info ne 0 in dgetrf')
        endif
  
        det(1) = 1.0d0
        det(2) = 0.0d0
        ten = 10.0d0
        do 50 i = 1, nsub
          if (ipvt(i) .ne. i) det(1) = -det(1)
          det(1) = a(i,i)*det(1)
c        ...exit
          if (det(1) .eq. 0.0d0) go to 60
   10     if (dabs(det(1)) .ge. 1.0d0) go to 20
          det(1) = ten*det(1)
          det(2) = det(2) - 1.0d0
          go to 10
   20     continue
   30     if (dabs(det(1)) .lt. ten) go to 40
            det(1) = det(1)/ten
            det(2) = det(2) + 1.0d0
          go to 30
   40     continue
   50   continue
   60   continue

        determinant = det(1)*10.0**det(2)

        call dgetri(nsub,a,nsub,ipvt,work,MELEC,info)

      endif

      return

      end

