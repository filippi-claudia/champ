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
      do 20 m=1,nel
        msh=msh+nel
        jsh=-nel
        do 15 j=1,nel
          work(j)=0
          jsh=jsh+nel
          do 10 i=1,nel
   10       work(j)=work(j)+work_mat(i+jsh)*slmi(i+msh)
   15     continue
       do 20 n=1,nel
         xmat(m,n)=0
         jsh=-nel
         do 20 j=1,nel
           jsh=jsh+nel
           xmat(m,n)=xmat(m,n)+work(j)*slmi(n+jsh)
   20     continue
          
      return
      end
