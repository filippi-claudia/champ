
!>
!!
!! This module contains various routines to generate and update the Hessian
!! used in geometry optimizations.
!!
!!
!! @author Jonas Feldt (j.feldt@utwente.nl)
!! @date November 2018
!!
module optgeo_hessian

  implicit none

  contains

  !>
  !! Initializes a diagonal hessian. Does not check if it is square.
  !!
  !! @param h square matrix which will be initialized to a identity matrix.
  !!
  subroutine optgeo_diagonal (h)

    real(kind=8), intent(out) :: h(:,:)

    integer :: i, j

    h = 0d0
    do i = 1, size (h, 1)
      h(i, i) = 1d0
    enddo

  end subroutine optgeo_diagonal



  !>
  !! Initializes a diagonal hessian with the elements scaled according to the
  !! type of the degree of freedom. The scaling factors are:
  !!  Bond     0.5
  !!  Angle    0.2
  !!  Dihedral 0.1
  !! The degrees of freedom are expected to be sorted in the order:
  !!  bonds, angles, dihedrals
  !!
  !! @param h square Hessian matrix which will be initialized to a scaled
  !!          diagonal matrix
  !! @param index_angles index of the first angle
  !! @param index_dihedrals index of the first dihedral
  !!
  subroutine optgeo_diagonal_scaled (h, index_angles, index_dihedrals)

    real(kind=8), intent(out) :: h(:,:)
    integer, intent(in) :: index_angles, index_dihedrals

    integer :: i

    real(kind=8), parameter :: sb = 0.5d0, sa = 0.2d0, sd = 0.1d0

    h = 0d0

    do i = 1, index_angles - 1
      h(i,i) = sb
    enddo

    do i = index_angles, index_dihedrals - 1
      h(i,i) = sa
    enddo

    do i = index_dihedrals, size (h, 1)
      h(i,i) = sd
    enddo

  end subroutine optgeo_diagonal_scaled


  !>
  !! Updates the hessian with the rank-two Broyden-Fletcher-Shanno-Goldfarb
  !! (BFGS) formula. Guarantees a positive definite updated matrix which is
  !! useful for minimizations but cannot be used for saddle-points. The
  !! Hessian, gradient and coordinates have to be in the same coordinate
  !! system.
  !!
  !! @param h Hessian matrix
  !! @param dx change of the coordinates in the last step
  !! @param dg change of the gradient in the last step
  !!
  subroutine optgeo_update_bfgs (h, dx, dg)
    use contrl_file,    only: ounit
    real(kind=8), intent(inout) :: h(:,:)
    real(kind=8), intent(in) :: dx(:), dg(:)

    real(kind=8), allocatable :: u(:,:)

    integer :: num
    real(kind=8) :: alpha

    write(ounit,*) 'Do not use optgeo_update_bfgs. Not tested yet!'
    stop

    ! initializes
    num = size (h, 1)
    allocate (u(num, num), source=0.0_8)
    u = 0d0

    ! second term, eq. 44, u = h.outer(dx, dx).h / inner(dx, h.dx)
    call dger (num, num, alpha, dx, 1, dx, 1, u, num) ! u = outer(dx, dx)
    u = matmul (u, h)                                 ! u = u.h
    u = matmul (h, u)                                 ! u = h.u
    u = u / dot_product (dx, matmul (h, dx))

    ! first term, eq. 44, h = h + outer (dg, dg) / inner (dg, dx)
    alpha = dot_product (dg, dx)
    call dger (num, num, alpha, dg, 1, dg, 1, h, num)

    h = h - u

  end subroutine optgeo_update_bfgs





end module optgeo_hessian








