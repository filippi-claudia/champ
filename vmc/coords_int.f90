!>
!!
!! Subroutines related to transforming coordinates and gradients
!! between cartesian and internal coordinate systems.
!! 
!! Useful reference: J. Chem. Phys, 117, 2002
!!       (I think there is an error in the angle contribution to B in the paper)
!!
!! @author Jonas Feldt (j.feldt@utwente.nl)
!! @date October 2018
!!
module coords_int

  implicit none

  real(kind=8), allocatable :: b (:,:)
  real(kind=8), allocatable :: binv (:,:)
  real(kind=8), allocatable :: p (:,:)
  real(kind=8), allocatable :: h (:,:)

  real(kind=8), allocatable :: int_gradients (:)
  real(kind=8), allocatable :: int_step (:)
  real(kind=8), allocatable :: int_coords (:)

  real(kind=8), allocatable :: cart_gradients (:)
  integer :: num_cart, num_int, num_centers
  integer :: num_bonds, num_angles, num_dihedrals

  !> Uses the projector P=BB^+
  logical, parameter :: project = .false.
  !> Recalculates B in every iteration while transforming the geometry step
  logical, parameter :: recalculate = .false.

  logical :: initialized = .false.
  
  private

  public :: coords_init
  public :: coords_compute_wilson
  public :: coords_transform_gradients
  public :: coords_compute_step
  public :: coords_transform_step



  contains


  !>
  !! Initalizes data for the geometry optimization in z-Matrix internal
  !! coordinates based on the provided cartesian coordinates.
  !!
  !! @param ncent number of atoms
  !!
  subroutine coords_init (ncent)
    use optgeo_hessian

    integer, intent(in) :: ncent

    if (initialized.eqv..true.) return

    num_centers = ncent
    num_cart = 3 * num_centers
    if (num_cart.eq.6) then
      num_int = 1
    else
      num_int = num_cart - 6
    endif

    num_bonds = num_centers - 1
    num_angles = num_centers - 2
    if (num_angles.lt.0) num_angles = 0
    num_dihedrals = num_centers - 3
    if (num_dihedrals.lt.0) num_dihedrals = 0

    allocate (b(num_int, num_cart))
    allocate (binv(num_cart, num_int))
    if (project) allocate (p(num_int, num_int))

    allocate (int_gradients(num_int))
    allocate (int_step(num_int))
    allocate (int_coords(num_int))

    allocate (cart_gradients(num_cart))

    allocate (h(num_int, num_int))
    call optgeo_diagonal_scaled (h, num_centers, 2 * num_centers - 2)

    initialized = .true.
    
  end subroutine coords_init



  !>
  !! Computes Wilson's B matrix with the dimensions (3*MCENT)x(nint) i.e.
  !! number of cartesian coordinates times number of internal coordinates.
  !!
  !! @param cart_coords current cartesian coordinates
  !! @param connectivities connectivities as specified in a z-Matrix
  !!
  subroutine coords_compute_wilson (cart_coords, connectivities)
    use misc_bond_func, only: cross_product

    real(kind=8), dimension(:,:), intent(in) :: cart_coords
    integer, dimension(:,:), intent(in) :: connectivities
    
    real(kind=8) :: uvec(3), vvec(3), wvec(3)
    real(kind=8) :: un, vn, wn, cijk, sijk
    real(kind=8) :: uxw(3), vxw(3)
    real(kind=8) :: cu, cv, susq, svsq
    integer :: i
    integer :: iint           ! ID of current internal coordinate
    integer :: cb, cc, cd     ! indices of bonded neighbours
    integer :: b1, b2, b3, b4 ! indices for cart. coordinates in b

    b = 0d0

    ! computes all bond contributions
    iint = 1
    do i = 2, num_centers ! starts at 2 ignoring the irrelevant 1. row in z Matrix

      ! computes normalized bond vector
      cb = connectivities(1, i)
      uvec = cart_coords(:, i) - cart_coords(:, cb)
      uvec = uvec / norm2 (uvec)

      b1 = 3 * (i - 1) + 1
      b2 = 3 * (cb - 1) + 1
      b(iint, b1:b1+2) =  uvec
      b(iint, b2:b2+2) = -uvec

      iint = iint + 1
    enddo

    ! computes all angle contributions
    do i = 3, num_centers ! first row in z Matrix with angle

      cb = connectivities(1, i)
      cc = connectivities(2, i)
      ! vectors: i->b->c
      uvec = cart_coords(:, cb) - cart_coords(:,i )
      vvec = cart_coords(:, cc) - cart_coords(:,cb)
      un = norm2 (uvec)
      vn = norm2 (vvec)

      ! computes cos and sin of angle
      cijk = dot_product (uvec, vvec) / un / vn
      sijk = sqrt (1d0 - cijk**2)

      b1 = 3 * (i - 1) + 1
      b2 = 3 * (cb - 1) + 1
      b3 = 3 * (cc - 1) + 1

      b(iint, b1:b1+2) = (-cijk * vn * uvec - un * vvec) / (sijk * un**2 * vn)
      b(iint, b2:b2+2) = ( cijk * vn**2 * uvec - un * vn * uvec &
                          - cijk * un**2 * vvec + un * vn * vvec) / (sijk * un**2 * vn**2)
      b(iint, b3:b3+2) = ( cijk * un * vvec + vn * uvec) / (sijk * un * vn**2)

      iint = iint + 1
    enddo


    do i = 4, num_centers ! first row in z Matrix with dihedral
      cb = connectivities(1, i)
      cc = connectivities(2, i)
      cd = connectivities(3, i)
      ! vectors: i->b->c->d
      uvec = cart_coords(:, i ) - cart_coords(:,cb)
      wvec = cart_coords(:, cc) - cart_coords(:,cb)
      vvec = cart_coords(:, cd) - cart_coords(:,cc)
      un = norm2 (uvec)
      vn = norm2 (vvec)
      wn = norm2 (wvec)
      uvec = uvec / un
      vvec = vvec / vn
      wvec = wvec / wn

      ! cross products
      uxw = cross_product (uvec, wvec)
      vxw = cross_product (vvec, wvec)
      cu =  dot_product (uvec, wvec)
      cv = -dot_product (uvec, wvec)
      susq = 1d0 - cu**2
      svsq = 1d0 - cv**2

      b1 = 3 * (i  - 1) + 1
      b2 = 3 * (cb - 1) + 1
      b3 = 3 * (cc - 1) + 1
      b4 = 3 * (cd - 1) + 1

      b(iint, b1:b1+2) =  uxw / (un * susq)
      b(iint, b2:b2+2) = -uxw / (un * susq) + uxw * cu / (susq * wn) + vxw * cv / (svsq * wn)
      b(iint, b3:b3+2) =  vxw / (vn * svsq) - uxw * cu / (susq * wn) - vxw * cv / (svsq * wn)
      b(iint, b4:b4+2) = -vxw / (vn * svsq)

      iint = iint + 1
    enddo

  end subroutine coords_compute_wilson



  !>
  !! Transforms the given cartesian gradients with the current
  !! B matrix.
  !!
  !! @param cart_gradient2d gradients of the cartesian coordinates (3,MCENT)
  !!
  subroutine coords_transform_gradients (cart_gradients2d)

    real(kind=8), dimension(:,:), intent(in) :: cart_gradients2d

    real(kind=8), dimension(:,:), allocatable :: btinv

    ! Everything related to Lapack
    real(kind=8), dimension(:), allocatable :: s
    real(kind=8), dimension(:,:), allocatable :: u
    real(kind=8), dimension(:,:), allocatable :: a
    real(kind=8), dimension(:), allocatable :: work
    integer, dimension(:), allocatable :: iwork
    integer :: i, irank, info
    integer :: lwork, liwork
    real(kind=8), parameter :: rcond = -1d0 ! optional conditioner to remove eigenvalues close to zero

    real(kind=8),  parameter :: PI  = 4 * atan (1d0)

    lwork = -1
    liwork = -1

    ! trivially reshapes to vector
    cart_gradients = reshape (cart_gradients2d, shape (cart_gradients))

    !
    ! computes the inverse of B -> B^+
    !
    allocate (a, source=b)
    allocate (s(num_int))
    allocate (u(num_cart, num_cart))
    u=0d0
    do i=1,num_cart
      u(i,i) = 1d0
    enddo
    allocate (work(1))
    allocate (iwork(1))

    ! queries for work space size
    call dgelsd (num_int, num_cart, num_int, a, num_int, u, num_cart, s, rcond, irank, work, lwork, iwork, info)
    if (info.ne.0) then
      write (6,*) 'transform_gradients: dgelsd() query for workspace failed.'
      stop
    end if  

    ! allocates optimal work and iwork
    lwork = int (work(1))
    deallocate (work)
    allocate (work(lwork))
    liwork = iwork(1)
    deallocate (iwork)
    allocate (iwork(liwork))


    call dgelsd (num_int, num_cart, num_int, a, num_int, u, num_cart, s, rcond, irank, work, lwork, iwork, info)

    if (info.gt.0) then
      write (6,*) 'transform_gradients: dgelsd() did not converge.'
      stop
    end if  

    ! extracts pseudo-inverse
    allocate (btinv(num_int, num_cart))
    binv = u(1:num_cart, 1:num_int) ! save for later
    btinv = transpose(binv)

    ! transforms the gradients
    int_gradients = matmul (btinv, cart_gradients)

    if (project) then ! computes projector P = BB^+
      p = matmul (b, binv)
      int_gradients = matmul (p, int_gradients)
    endif

  end subroutine coords_transform_gradients



  !>
  !! Computes the step in internal coordinates
  !!
  !! For now just a trivial steepest descent step. This can be improved.
  !!
  !! @param alpha scales the gradient and determines the length of the step
  !!
  subroutine coords_compute_step (alpha)

    real(kind=8), intent (in) :: alpha
    integer :: i

    int_step = -alpha * int_gradients

    ! divide by the diagonal elements of the Hessian matrix !TODO for now only works for this diagonal one
    do i = 1, num_int
      int_step(i) = int_step(i) / h(i, i)
    enddo

    ! Hessian diagonal elements according to type x=0.5, 0.2, 0.1
    !do i = 1, num_bonds
    !  int_step(i) = int_step(i) / 0.5d0
    !enddo
    !do i = num_bonds + 1, num_bonds + num_angles
    !  int_step(i) = int_step(i) / 0.2d0
    !enddo
    !do i = num_angles + num_bonds + 1, num_bonds + num_angles + num_dihedrals
    !  int_step(i) = int_step(i) / 0.1d0
    !enddo
    
  end subroutine coords_compute_step


  !>
  !! Transforms the step from internal to cartesian coordinates
  !! and computes the new geometry.
  !!
  !! @param int_coords2d internal coordinates in z-Matrix representation (3xMCENT)
  !! @param cart_coords2d cartesian coordinates (3xMCENT)
  !! @param connectivities z Matrix connectivity matrix (3xMCENT)
  !!
  !!
  subroutine coords_transform_step (int_coords2d, cart_coords2d, connectivities)

    real(kind=8), dimension(:,:), intent(inout) :: cart_coords2d
    real(kind=8), dimension(:,:), intent(inout) :: int_coords2d
    integer, dimension(:,:), intent(in) :: connectivities

    integer, parameter :: maxit = 25

    real(kind=8), dimension(num_cart) :: cart_coords
    real(kind=8), dimension(num_int)  :: int_dnew
    real(kind=8), dimension(num_cart) :: cart_d
    real :: delta
    integer :: iint
    integer :: ic, it
    integer :: i, j

    iint = 1

    ! saves the original internal coordinates as reference (saves q0)
    do ic = 2, num_centers ! loop over bonds
      int_coords(iint) = int_coords2d(1, ic)
      iint = iint + 1
    enddo
    do ic = 3, num_centers ! loop over angles
      int_coords(iint) = int_coords2d(2, ic)
      iint = iint + 1
    enddo
    do ic = 4, num_centers ! loop over dihedrals
      int_coords(iint) = int_coords2d(3, ic)
      iint = iint + 1
    enddo

    ! computes new geometry in cartesian coordinates (x1, eq. 13)
    cart_coords  = reshape (cart_coords2d, shape (cart_coords)) ! 2d->1d
    cart_coords = cart_coords + matmul (binv, int_step)

    ! updates 2D cartesian coordinates
    do i=1,num_centers
      do j = 1,3
        cart_coords2d(j, i) = cart_coords(3*(i-1)+j)
      enddo
    enddo

    ! transforms back to internal coordinates (q1)
    call cart2zmat(num_centers, cart_coords2d, connectivities, int_coords2d)

    call fix_dihedrals (int_coords, int_coords2d)

    ! computes the guess step in internal coordinates (q1-q0)
    iint = 1
    do ic = 2, num_centers ! loop over bonds
      int_dnew(iint) = (int_coords2d(1, ic) - int_coords(iint))
      iint = iint + 1
    enddo
    do ic = 3, num_centers ! loop over angles
      int_dnew(iint) = (int_coords2d(2, ic) - int_coords(iint))
      iint = iint + 1
    enddo
    do ic = 4, num_centers ! loop over dihedrals
      int_dnew(iint) = (int_coords2d(3, ic) - int_coords(iint))
      iint = iint + 1
    enddo

    ! computes difference between the reference and guess step (delta_q1, eq. 14)
    int_dnew = int_step - int_dnew
    cart_d = matmul (binv, int_dnew)
    delta = sqrt(sum(cart_d**2) / num_cart)

    it = 0
    do while (delta.gt.1d-6.and.it.lt.maxit)

      if (recalculate) then
        call coords_compute_wilson (cart_coords2d, connectivities)
      endif

      ! step in cartesian coordinates
      cart_coords = cart_coords + cart_d

      ! updates 2D cartesian coordinates
      do i=1,num_centers
        do j = 1,3
          cart_coords2d(j, i) = cart_coords(3*(i-1)+j)
        enddo
      enddo

      ! transforms back to internal coordinates
      call cart2zmat(num_centers, cart_coords2d, connectivities, int_coords2d)

      call fix_dihedrals (int_coords, int_coords2d)

      ! computes difference between reference step and actual step
      iint = 1
      do ic = 2, num_centers ! loop over bonds
        int_dnew(iint) = (int_coords2d(1, ic) - int_coords(iint))
        iint = iint + 1
      enddo
      do ic = 3, num_centers ! loop over angles
        int_dnew(iint) = (int_coords2d(2, ic) - int_coords(iint))
        iint = iint + 1
      enddo
      do ic = 4, num_centers ! loop over dihedrals
        int_dnew(iint) = (int_coords2d(3, ic) - int_coords(iint))
        iint = iint + 1
      enddo

      ! computes difference between the reference and guess step
      int_dnew = int_step - int_dnew
      cart_d = matmul (binv, int_dnew)
      delta = sqrt(sum(cart_d**2) / num_cart)

      it = it + 1
    enddo

    if (delta.gt.1d-6) then
      write (6,*) 'do_step: backtransformation did not converge'
      stop
    endif

    ! obtain final set of coordinates with dihedrals in the range of (-pi,pi)
    call cart2zmat(num_centers, cart_coords2d, connectivities, int_coords2d)

  end subroutine coords_transform_step


  subroutine fix_dihedrals (int_reference, int_coords2d)

    real(kind=8), intent(in) :: int_reference(:)
    real(kind=8), intent(inout) :: int_coords2d(:,:)
    integer :: iint, ic

    real(kind=8),  parameter :: PI2  = 8 * atan (1d0)

    iint = num_bonds + num_angles + 1
    do ic = 4, num_centers ! loop over dihedrals

      if (int_reference(iint).lt.0d0.and.int_coords2d(3,ic).gt.0d0) then
        int_coords2d(3, ic) = int_coords2d(3, ic) - PI2
        write (6,*) 'FIX Flip of the dihedral detected. - to +.'
      endif
      if (int_reference(iint).gt.0d0.and.int_coords2d(3,ic).lt.0d0) then
        int_coords2d(3, ic) = int_coords2d(3, ic) + PI2
        write (6,*) 'FIX Flip of the dihedral detected. + to -.'
      endif

      iint = iint + 1
    enddo

  end subroutine fix_dihedrals


end module coords_int

