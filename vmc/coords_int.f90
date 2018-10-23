!>
!!
!! Subroutines related to transforming coordinates and gradients
!! between cartesian and internal coordinate systems.
!!
!! Jonas Feldt (j.feldt@utwente.nl)
!!
module coords_int

  implicit none

  real(kind=8), allocatable :: bmat (:,:)
  real(kind=8), allocatable :: int_gradients (:)
  real(kind=8), allocatable :: cart_gradients (:)
  integer :: num_cart, num_int, num_centers

  

  !private
  !public



  contains


  !>
  !! Initalizes data for the geometry optimization in z-Matrix internal
  !! coordinates based on the provided cartesian coordinates.
  !!
  !! @param mcent number of atoms
  !!
  subroutine init (mcent)

    integer, intent(in) :: mcent

    num_centers = mcent
    num_cart = 3 * num_centers
    if (num_cart.eq.6) then
      num_int = 1
    else
      num_int = num_cart - 6
    endif

    allocate (bmat(num_int, num_cart))
    allocate (int_gradients(num_int))
    allocate (cart_gradients(num_cart))
    
  end subroutine init



  !>
  !! Computes Wilson's B matrix with the dimensions (3*MCENT)x(nint) i.e.
  !! number of cartesian coordinates times number of internal coordinates.
  !!
  !! @param cart_coords current cartesian coordinates
  !! @param connectivities connectivities as specified in a z-Matrix
  !!
  subroutine compute_bmat (cart_coords, connectivities)

    real(kind=8), dimension(:,:), intent(in) :: cart_coords
    integer, dimension(:,:), intent(in) :: connectivities
    
    real(kind=8) :: vec(3)
    integer :: cb, cc, cd
    integer :: iint = 1 ! ID of current internal coordinate
    integer :: i
    integer :: b1, b2 ! indices for cart. coordinates in bmat

    ! computes all bond contributions
    do i = 2, num_centers ! starts at 2 ignoring the irrelevant 1. row in z Matrix

      ! computes normalized bond vector
      cb = connectivities(1, i)
      vec = cart_coords(1:3, i) - cart_coords(1:3, cb)
      vec = vec / norm2 (vec)


      b1 = b (i)
      b2 = b (cb)
      bmat(iint, b1:b1+2) = -vec
      bmat(iint, b2:b2+2) =  vec
      iint = iint + 1

    enddo

  end subroutine compute_bmat

  !>
  !! Computes the starting index of the atom i in the bmatrix.
  !!
  !! @param i ID of atom
  !! 
  !! @return index in the bmatrix
  !!
  integer function b (i)
    integer :: i, ib
    ib = 3 * (i - 1) + 1
    return
  end function



  !>
  !! Transforms the given cartesian gradients with the current
  !! B matrix.
  !!
  !! @param cart_gradient2d gradients of the cartesian coordinates (3,MCENT)
  !!
  subroutine transform_gradients (cart_gradients2d)
    real(kind=8), dimension(:,:), intent(in) :: cart_gradients2d
    real(kind=8), dimension(:,:), allocatable :: bmattinv

    real(kind=8), dimension(:), allocatable :: s
    real(kind=8), dimension(:,:), allocatable :: u
    real(kind=8), dimension(:,:), allocatable :: a
    real(kind=8), dimension(:), allocatable :: work
    integer, dimension(:), allocatable :: iwork
    integer :: i, irank, info
    integer :: lwork = -1, liwork = -1
    real(kind=8) :: rcond = -1d0


    ! trivially reshapes to vector
    cart_gradients = reshape (cart_gradients2d, shape (cart_gradients))

    ! computes the inverse of the transpose of B (= (B^T)^+)

    allocate (s(num_int))
    allocate (u(num_cart, num_cart))
    u=0d0
    do i=1,num_cart
      u(i,i) = 1d0
    enddo
    allocate (work(1))
    allocate (iwork(1))

    ! query for work space size
    call dgelsd (num_int, num_cart, num_int, bmat, num_int, u, num_cart, s, rcond, irank, work, lwork, iwork, info)
    if (info.ne.0) then
      write (*,*) 'transform_gradients: dgelsd() query for workspace failed.'
      stop
    end if  

    ! allocate optimal work and iwork
    lwork = int (work(1))
    deallocate (work)
    allocate (work(lwork))
    liwork = int (iwork(1))
    deallocate (iwork)
    allocate (iwork(liwork))


    ! do something
    call dgelsd (num_int, num_cart, num_int, bmat, num_int, u, num_cart, s, rcond, irank, work, lwork, iwork, info)

    if (info.gt.0) then
      write (*,*) 'transform_gradients: dgelsd() did not converge.'
      stop
    end if  

    print *,"u"
    do i = 1, num_cart
      write (*,'(9f10.5)') u(i, 1:num_int)
    enddo

    allocate (bmattinv(num_int, num_cart))
    bmattinv = transpose(u(1:9, 1:3))

    print *,"bmattinv"
    do i = 1, num_int
      write (*,'(9f10.5)') bmattinv(i, 1:num_cart)
    enddo

    ! transforms the gradients
    int_gradients = matmul (bmattinv, cart_gradients)

  end subroutine transform_gradients


end module coords_int





