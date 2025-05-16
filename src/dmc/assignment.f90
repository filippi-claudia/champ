! Linear assignment problem (LAP) solver
! -----------------------------------------------
! The first method (method = 1) is a brute force
! solution that computes all the permutations (n!)
! using Heap's algorithm - time complexity O(n!)
! The second method (method = 2) uses the Munkres
! algorithm (also called the Hungarian algorithm) to
! solve the problem - time complexity O(n^3)
! -----------------------------------------------
! F. Brieuc - March 2017
! Modified by Cyrus Umrigar, Oct 2022
! Modified by Tyler Anderson, Feb 2022
! Modified by Joris van de Nes, Feb 2025

module assignment_mod 
   implicit none

   real*8, dimension(:,:), allocatable :: CC    ! cost matrix (max, min sum)
   integer :: n      ! dimension of CC - assumed to be a (nxn) square matrix

   ! following variables are use only for Munkres algorithm
   integer, dimension(:,:), allocatable :: M    ! mask matrix
   integer, dimension(:), allocatable :: rowCover, colCover !cover row and cols
   integer :: pathRow0 = 0, pathCol0 = 0  ! starting point for path finding part

contains

!-------------------------------------------------------------------------------

    !> Assigns all the electrons to their fragments.                           \
    !> method = 1 uses the O(n!) heap algrithm for the assignments             \
    !> method = 2 uses the O(n^3) Munkres algorithm for the assignment (aka the
    !> Hungarian algoritm).                                                    \
    !> The final asssignment is stored in ifragelec
    subroutine assign_elecs(x,method)
       use control, only: ipr
       use fragments, only: nfrag, ifragcent, ifragelec
       use precision_kinds, only: dp
       use system, only: znuc, iwctype, cent, ncent, nelec
       
       implicit none
       real(dp) :: x(3,*)
       integer :: method
       real*8 :: sumSol, pot ! minimal sum
       integer :: i, j, k, ic, nuclj(nelec), znuc_tot
       real*8 :: sum

       
       znuc_tot=0
       do ic=1,ncent
           znuc_tot=znuc_tot+int(znuc(iwctype(ic)))
       enddo
       if(znuc_tot.NE.nelec) stop ('error : Fragments not implemented for ionic systems.') 
       
       
       ! read/generate CC and initialize
       n=nelec
       if (.NOT.allocated(CC)) allocate(CC(n,n))

       do i=1,n
          j=0
          do ic=1,ncent
             pot=-znuc(iwctype(ic))/norm2(x(:,i)-cent(:,ic))
             do k=1,int(znuc(iwctype(ic)))
               j=j+1
               CC(j,i)=pot
               nuclj(j)=ic
             enddo
          enddo
       enddo

       if(n<=20.AND.(ipr.GE.1)) then
           write(6, '(''LAP Cost Matrix='')')
           call printMatrix(CC,n)
       endif
       
       select case(method)
       case (1)
          ! solve lap by computing all permutations using Heap's algorithm
          call heap(ifragelec)
       case(2)
          ! Munkres' (Hungarian) algorithm
          call munkres(ifragelec)
       case default
          stop ('error : wrong value for parameter method (1 or 2)')
       end select
    
       if (ipr.GE.1) then
          ! print the solution
          sumSol = 0
          do i = 1, n
             sumSol = sumSol + CC(ifragelec(i),i)
          enddo
          write(6, '(a29,f9.3)') 'The minimum possible LAP sum is: ', sumSol
          write(6, '(a38)') 'obtained using the following elements:'
          write(6, '(a)') 'row i, col j, CC(i,j), sum'
          sum = 0
          do i = 1, n
             sum = sum + CC(ifragelec(i),i)
             write(6,'(i4,x,i4,x,f6.3,x,f9.3)') i, ifragelec(i), CC(ifragelec(i),i), sum
          enddo
       endif

       do i=1,nelec
          !After heap/munkres is called, ifragelec(i) is the 'job' index j 
          !assigned to electron i.
           
          !The nucleus corresponding to job j is nuclj(j), so 
          !nuclj(ifragelec(i)) is the nucleus assigned to electron i.
           
          !The fragment corresponding to nucleus ic is ifragnucl(ic), so 
          !ifragnucl(nuclj(ifragelec(i))) gives the fragment assigned to 
          !electron i.
          ifragelec(i)=ifragcent(nuclj(ifragelec(i)))
       enddo
    end subroutine assign_elecs

!-------------------------------------------------------------------------------
    subroutine printMatrix(M,d)
       implicit none
    
       integer, intent(in) :: d
    !  integer, dimension(d,d), intent(in) :: M
       real*8, dimension(d,d), intent(in) :: M
       integer :: i, j
    
       write(6, '(a12)') 'Cost matrix:'
    
       if (d == 1) then
          write(6,'(i4)') M(1,1)
       else
          do i = 1, d
             do j = 1, d - 1
    !           write(6, '(i4)', advance="no") M(j,i)
                write(6, '(f8.3)', advance="no") M(j,i)
             enddo
    !        write(6, '(i4)') M(d,i)
             write(6, '(f8.3)') M(d,i)
          enddo
       endif
    
       write(6,*) ''
    
    end subroutine printMatrix

!-------------------------------------------------------------------------------
    subroutine heap(jSol)
       ! Compute all the possible (n!) permutations for jSol(i)
       ! using Heap's algorithm and find the smallest sum
       ! Heap, B. R., The Computer Journal 6, 3 (1963)
       ! see : https://en.wikipedia.org/wiki/Heap's_algorithm
       implicit none
    
       integer, dimension(n), intent(out) :: jSol ! solution indices
    
       integer, dimension(:), allocatable :: p  ! for permutations
       integer, dimension(:), allocatable :: j  ! choosen indice of row i
    !  integer :: sum, sum_tmp
       real*8 :: sum, sum_tmp
       integer :: i, tmp
    
       allocate(j(n))
       allocate(p(n))
    
       ! *** Compute all the possible (n!) permutations for jj(i) ***
       ! Heap's algorithm
       sum = 0
       do i = 1, n
          jSol(i) = i
          j(i) = i
          sum = sum + CC(i,i)
          p(i) = 1
       enddo
       sum_tmp = sum
    
       i = 1;
       do while (i <= n)
          if (p(i) < i) then
             if (mod(i,2) == 0) then
                ! i odd - swap i <-> p(i)
                sum_tmp = sum_tmp - CC(j(p(i)),p(i)) - CC(j(i),i)
                tmp = j(p(i))
                j(p(i)) = j(i)
                j(i) = tmp
                sum_tmp = sum_tmp + CC(j(p(i)),p(i)) + CC(j(i),i)
             else
                ! i even - swap i <-> 1
                sum_tmp = sum_tmp - CC(j(1),1) - CC(j(i),i)
                tmp = j(1)
                j(1) = j(i)
                j(i) = tmp
                sum_tmp = sum_tmp + CC(j(1),1) + CC(j(i),i)
             endif
             if (sum_tmp < sum) then
                sum = sum_tmp
                jSol(:) = j(:)
             endif
             p(i) = p(i) + 1
             i = 1
          else
             p(i) = 1
             i = i + 1
          endif
       enddo
    
       deallocate(j)
       deallocate(p)
    
    end subroutine heap

!-------------------------------------------------------------------------------
    subroutine munkres(jSol)
       ! Implementation of the Munkres algorithm (also referred to as the Hungarian
       ! algorithm). J. Munkres, Journal of the SIAM 5, 1 (1957)
       ! The following implementation is based on
       ! http://csclab.murraystate.edu/%7Ebob.pilgrim/445/munkres.html
       implicit none
    
       integer, dimension(n), intent(out) :: jSol ! solution indices
    
       integer :: step, i, j, tmp
       logical :: done
    
       done = .false.
       step = 1
       tmp = 0
    
       allocate(M(n,n))      ! mask matrix - contains starred zeros
       allocate(rowCover(n)) ! to keep track of covered rows
       allocate(colCover(n)) ! to keep track of covered columns
    
       do i = 1, n
          M(:,i) = 0
          rowCover(i) = 0
          colCover(i) = 0
       enddo
    
       do while(.not. done)
          select case(step)
          case(1)
             call step1(step)
          case(2)
             call step2(step)
          case(3)
             call step3(step)
          case(4)
             call step4(step)
          case(5)
             call step5(step)
          case(6)
             call step6(step)
          case default ! done
             do i = 1, n
                do j = 1, n
                   if (M(j,i) == 1) jSol(i) = j
                enddo
             enddo
             done = .true.
          end select
       enddo
    
       deallocate(M)
       deallocate(rowCover)
       deallocate(colCover)
    
    end subroutine munkres

!-------------------------------------------------------------------------------
    subroutine step1(step)
       ! row reduction : for each row find the smallest value and substract it from
       ! all elements of that row. Go to step 2.
       implicit none
    
       integer, intent(out) :: step
    
    !  integer :: minVal, i, j
       integer :: i, j
       real*8 :: minVal
    
       do i = 1, n
          minVal = CC(1,i)
          do j = 1, n
             if (CC(j,i) < minVal) minVal = CC(j,i)
          enddo
          CC(:,i) = CC(:,i) - minVal
       enddo
    
       step = 2
    
    end subroutine step1

!-------------------------------------------------------------------------------
    subroutine step2(step)
       ! Search for zeros.
       ! Find a zero (Z) in the matrix. If no zeros has been previously starred in
       ! its row and column then star Z. Go to step 3.
       implicit none
    
       integer, intent(out) :: step
    
       integer :: i, j
    
       do i = 1, n
          do j = 1, n
             if (CC(j,i) == 0 .and. rowCover(i) == 0 .and. colCover(j) == 0) then
                M(j,i) = 1
                rowCover(i) = 1
                colCover(j) = 1
             endif
          enddo
       enddo
       ! uncovers
       do i = 1, n
          rowCover(i) = 0
          colCover(i) = 0
       enddo
    
       step = 3
    
    end subroutine step2

!-------------------------------------------------------------------------------
    subroutine step3(step)
       ! cover each column containing a starred zero. If n column are covered
       ! the starred zero describe an optimal assignment and we are done otherwise
       ! go to step 4.
       implicit none
    
       integer, intent(out) :: step
    
       integer :: colCount, i, j
    
       colCount = 0
       do i = 1, n
          do j = 1, n
             ! if starred and column is uncovered
             if (M(j,i) == 1 .and. colCover(j) == 0) then
                colCover(j) = 1
                colCount = colCount + 1
             endif
          enddo
       enddo
    
       if (colCount == n) then
          step = 0
       else
          step = 4
       endif
    
    end subroutine step3

!-------------------------------------------------------------------------------
    subroutine step4(step)
       ! Find a uncovered zero and prime it. If there is no starred zero in the row
       ! go to step 5. Otherwise, cover the row and uncover the column containing
       ! the starred zero. Continue until no uncovered zeros is left. Go to step 6.
       implicit none
    
       integer, intent(out) :: step
    
       logical :: done, starInRow
       integer :: i, j, row, col
    
       done = .false.
    
       do while (.not. done)
          ! find an uncovered zero
          row = 0; col = 0
          starInRow = .false.
          loop1: do i = 1, n
             loop2: do j = 1, n
                if (CC(j,i) == 0 .and. rowCover(i) == 0 .and. colCover(j) == 0) then
                   row = i
                   col = j
                   exit loop1
                endif
             enddo loop2
          enddo loop1
          if (row == 0) then !no zero uncoverred left
             done = .true.
             step = 6
          else
             M(col,row) = 2 !primed zero
             ! search if there is a starred zero in the same row
             do j = 1, n
                if (M(j,row) == 1) then
                   starInRow = .true.
                   col = j
                endif
             enddo
             if (starInRow) then ! if there is a starred zero in line
                rowCover(row) = 1
                colCover(col) = 0
             else ! if no starred zero in line
                done = .true.
                step = 5
                pathRow0 = row
                pathCol0 = col
             endif
          endif
       enddo
    
    end subroutine step4
!-------------------------------------------------------------------------------
    
    subroutine step5(step)
       ! Augmenting path algorithm: construct a serie of alternating primed and
       ! starred zeros as follows. Let Z0 be the uncoverd primed zero found in
       ! step 4. Let Z1 be the starred zero in the column of Z0 (if any).
       ! Let Z2 be the primed zero in the row of Z1 (there will always be one).
       ! Continue until the series terminates at a primed zero that has no starred
       ! zero in its column. Then unstar each starred zeros of the series, star
       ! each primed zeros of the series, erase all primes and uncover every line
       ! and columns. Return to step 3.
       implicit none
    
       integer, intent(out) :: step
    
       logical :: done
       integer :: i, j
       integer :: row, col
       integer :: pathCount
       integer, dimension(2*n+1,2) :: path
    
       pathCount = 1
    
       path(pathCount,1) = pathRow0
       path(pathCount,2) = pathCol0
    
       done = .false.
    
       do while (.not. done)
          ! search for a starred zero in column
          row = 0
          col = path(pathCount,2)
          do i = 1, n
             if (M(col,i) == 1) row = i
          enddo
          if (row /= 0) then ! update path
             pathCount = pathCount + 1
             path(pathCount,1) = row
             path(pathCount,2) = path(pathCount-1,2)
          else
             done = .true.
          endif
          if (.not. done) then
             ! search for a prime zero in row
             do j = 1, n
                if (M(j,row) == 2) col = j
             enddo
             ! update path
             pathCount = pathCount + 1
             path(pathCount,1) = path(pathCount-1,1)
             path(pathCount,2) = col
          endif
       enddo
    
       ! augment path
       do i = 1, pathCount
          if(M(path(i,2),path(i,1)) == 1) then
             M(path(i,2),path(i,1)) = 0
          else
             M(path(i,2),path(i,1)) = 1
          endif
       enddo
    
       ! clear covers and erase primes
       do i = 1, n
          rowCover(i) = 0
          colCover(i) = 0
          do j = 1, n
             if (M(j,i) == 2) M(j,i) = 0
          enddo
       enddo
    
       step = 3
    
    end subroutine step5

!-------------------------------------------------------------------------------
    subroutine step6(step)
       ! Search for the smallest uncovered value and add it to covered rows
       ! and substract it from uncovered columns. Return to step 4.
       implicit none
    
       integer, intent(out) :: step
    
    !  integer :: i, j, minVal
       integer :: i, j
       real*8 :: minVal
    
       minVal = huge(i)
    
       do i = 1, n
          do j = 1, n
             if (rowCover(i) == 0 .and. colCover(j) == 0 .and. CC(j,i) < minVal) then
                minVal = CC(j,i)
             endif
          enddo
       enddo
       do i = 1, n
          do j = 1, n
             if (rowCover(i) == 1) CC(j,i) = CC(j,i) + minVal
             if (colCover(j) == 0) CC(j,i) = CC(j,i) - minVal
          enddo
       enddo
    
       step = 4
    
    end subroutine step6

end module assignment_mod
