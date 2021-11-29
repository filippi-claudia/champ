      subroutine multiply_slmi_mderiv_simple(nel,work_mat,work,slmi,xmat)

      use precision_kinds, only: dp
      implicit none

      integer :: i, j, jsh, m, msh
      integer :: n, nel

      real(dp), dimension(*) :: slmi
      real(dp), dimension(*) :: work_mat
      real(dp), dimension(nel, *) :: xmat
      real(dp), dimension(*) :: work


      msh=-nel
      do m=1,nel
        msh=msh+nel
        jsh=-nel
        do j=1,nel
          work(j)=0
          jsh=jsh+nel
          do i=1,nel
            work(j)=work(j)+work_mat(i+jsh)*slmi(i+msh)
          enddo
        enddo
       do n=1,nel
         xmat(m,n)=0
         jsh=-nel
         do j=1,nel
           jsh=jsh+nel
           xmat(m,n)=xmat(m,n)+work(j)*slmi(n+jsh)
         enddo
       enddo
      enddo
          
      return
      end
