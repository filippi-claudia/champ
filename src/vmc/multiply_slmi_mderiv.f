      subroutine multiply_slmi_mderiv_simple(nel,work_mat,work,slmi,xmat)

      implicit real*8(a-h,o-z)

      dimension slmi(*),work_mat(*),xmat(nel,*),work(*)

      msh=-nel
      do m=1,nel
         msh=msh+nel
         jsh=-nel
         do j=1,nel
            work(j)=0.0d0
            jsh=jsh+nel
            do i=1,nel
               work(j)=work(j)+work_mat(i+jsh)*slmi(i+msh)
            enddo
         enddo
         do n=1,nel
            xmat(m,n)=0.0d0
            jsh=-nel
            do j=1,nel
               jsh=jsh+nel
               xmat(m,n)=xmat(m,n)+work(j)*slmi(n+jsh)
            enddo
         enddo
      enddo
      
      end subroutine
