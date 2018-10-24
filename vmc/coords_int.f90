!>
!!
!! Subroutines related to transforming coordinates and gradients
!! between cartesian and internal coordinate systems.
!!
!! @author Jonas Feldt (j.feldt@utwente.nl)
!! @date October 2018
!!
module coords_int

  implicit none

  real(kind=8), allocatable :: bmat (:,:)
  real(kind=8), allocatable :: bmatinv (:,:)

  real(kind=8), allocatable :: int_gradients (:)
  real(kind=8), allocatable :: int_step (:)
  real(kind=8), allocatable :: int_coords (:)

  real(kind=8), allocatable :: cart_gradients (:)
  integer :: num_cart, num_int, num_centers

  logical :: initialized

  

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

    if (initialized.eqv..true.) return

    num_centers = mcent
    num_cart = 3 * num_centers
    if (num_cart.eq.6) then
      num_int = 1
    else
      num_int = num_cart - 6
    endif

    allocate (bmat(num_int, num_cart))
    allocate (bmatinv(num_cart, num_int))

    allocate (int_gradients(num_int))
    allocate (int_step(num_int))
    allocate (int_coords(num_int))

    allocate (cart_gradients(num_cart))

    initialized = .true.
    
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
    real(kind=8) :: rcond = -1d0 ! optional conditioner to remove eigenvalues close to zero


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

    ! queries for work space size
    call dgelsd (num_int, num_cart, num_int, bmat, num_int, u, num_cart, s, rcond, irank, work, lwork, iwork, info)
    if (info.ne.0) then
      write (*,*) 'transform_gradients: dgelsd() query for workspace failed.'
      stop
    end if  

    ! allocates optimal work and iwork
    lwork = int (work(1))
    deallocate (work)
    allocate (work(lwork))
    liwork = int (iwork(1))
    deallocate (iwork)
    allocate (iwork(liwork))


    call dgelsd (num_int, num_cart, num_int, bmat, num_int, u, num_cart, s, rcond, irank, work, lwork, iwork, info)

    if (info.gt.0) then
      write (*,*) 'transform_gradients: dgelsd() did not converge.'
      stop
    end if  

    print *,"u"
    do i = 1, num_cart
      write (*,'(9f10.5)') u(i, 1:num_int)
    enddo

    ! extracts pseudo-inverse
    allocate (bmattinv(num_int, num_cart))
    bmatinv = u(1:9, 1:3) ! save for later
    bmattinv = transpose(bmatinv)

    print *,"bmattinv"
    do i = 1, num_int
      write (*,'(9f10.5)') bmattinv(i, 1:num_cart)
    enddo

    ! transforms the gradients
    int_gradients = matmul (bmattinv, cart_gradients)

  end subroutine transform_gradients



  !>
  !! Computes the step in internal coordinates
  !!
  !! For now just a trivial steepest descent step. This can be improved.
  !!
  !! @param alpha scales the gradient and determines the length of the step
  !!
  subroutine compute_step_int (alpha)

    real(kind=8), intent (in) :: alpha

    int_step = -alpha * int_gradients
    
  end subroutine compute_step_int


  !>
  !! Transforms the step from internal to cartesian coordinates
  !! and computes the new geometry.
  !!
  !! @param int_coords2d internal coordinates in z-Matrix representation (3xMCENT)
  !! @param cart_coords2d cartesian coordinates (3xMCENT)
  !! @param connectivities z Matrix connectivity matrix (3xMCENT)
  !!
  !! Slightly more complicated than it has to be for a z Matrix in preparation
  !! for the generalization to other internal coordinates.
  !!
  subroutine do_step (int_coords2d, cart_coords2d, connectivities)

    real(kind=8), dimension(:,:), intent(inout) :: cart_coords2d
    real(kind=8), dimension(:,:), intent(inout) :: int_coords2d
    integer, dimension(:,:), intent(in) :: connectivities

    integer, parameter :: maxit = 25

    real(kind=8), dimension(num_cart) :: cart_coords
    real(kind=8), dimension(num_int) :: int_dnew
    real :: delta
    integer :: iint = 1
    integer :: ic, it

    ! computes new geometry in internal coordinates as reference
    do ic = 2, num_centers ! loop over bonds
      int_coords(iint) = int_coords2d(1, ic) + int_step(iint)
      iint = iint + 1
    enddo

    print *,"Internal reference"
    write (*,'(f10.5)') (int_coords(iint),iint=1,num_int)

    ! computes new geometry in cartesian coordinates
    cart_coords  = reshape (cart_coords2d, shape (cart_coords)) ! 2d->1d
    cart_coords = cart_coords - matmul (bmatinv, int_step) !TODO why do I need here minus instead of plus, check paper again...
    print *
    print *,"old cart"
    write (*,'(9f10.5)') cart_coords2d
    cart_coords2d  = reshape (cart_coords, shape (cart_coords2d)) ! 1d->2d
    print *
    print *,"new cart"
    write (*,'(9f10.5)') cart_coords2d

    ! transforms back to internal coordinates
    call cart2zmat(num_centers, cart_coords2d, connectivities, int_coords2d)
    !int_coords2d (1, 2) = norm2 (cart_coords2d(1:3, 2) - cart_coords2d(1:3, 1))
    !int_coords2d (1, 3) = norm2 (cart_coords2d(1:3, 3) - cart_coords2d(1:3, 2))

    print *
    print *, 'internal new'
    write (*,'(9f10.5)') int_coords2d

    ! computes difference between reference step and actual step
    iint = 1
    do ic = 2, num_centers ! loop over bonds
      int_dnew(iint) = (int_coords(iint) - int_coords2d(1, ic))
      iint = iint + 1
    enddo
    delta = sqrt(sum(int_dnew**2))

    print *
    print *,"int_step"
    write (*,'(3f10.5)') int_step
    print *,"int_dnew"
    write (*,'(3f10.5)') int_dnew
    print *, "Delta", delta

    it = 0
    do while (delta.gt.1d-6.and.it.lt.maxit)
      print *
      print *,"Iteration", it

      cart_coords = cart_coords - matmul (bmatinv, int_dnew) !TODO why do I need here minus instead of plus? check paper again...
      print *
      print *,"old cart"
      write (*,'(9f10.5)') cart_coords2d
      cart_coords2d  = reshape (cart_coords, shape (cart_coords2d)) ! 1d->2d
      print *
      print *,"new cart"
      write (*,'(9f10.5)') cart_coords2d

      ! transforms back to internal coordinates
      call cart2zmat(num_centers, cart_coords2d, connectivities, int_coords2d)
      !int_coords2d (1, 2) = norm2 (cart_coords2d(1:3, 2) - cart_coords2d(1:3, 1))
      !int_coords2d (1, 3) = norm2 (cart_coords2d(1:3, 3) - cart_coords2d(1:3, 2))

      print *
      print *, 'internal new'
      write (*,'(9f10.5)') int_coords2d

      ! computes difference between reference step and actual step
      iint = 1
      do ic = 2, num_centers ! loop over bonds
        int_dnew(iint) = (int_coords(iint) - int_coords2d(1, ic))
        iint = iint + 1
      enddo
      delta = sqrt(sum(int_dnew**2))

      print *
      print *,"int_step"
      write (*,'(3f10.5)') int_step
      print *,"int_dnew"
      write (*,'(3f10.5)') int_dnew
      print *, "Delta", delta

      it = it + 1
    enddo

    if (delta.gt.1d-6) then
      write (*,*) 'do_step: backtransformation did not converge'
      stop
    endif


  end subroutine do_step


end module coords_int

