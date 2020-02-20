module array_utils

  use numeric_kinds, only: dp
  use lapack_wrapper, only: lapack_generalized_eigensolver, lapack_matmul, lapack_matrix_vector, &
       lapack_qr, lapack_solver, lapack_sort
  implicit none

  !> \private
  private
  !> \public
  public :: concatenate, diagonal,eye, generate_diagonal_dominant, norm, &
       initialize_subspace, write_matrix, write_vector, check_deallocate_vector, &
       check_deallocate_matrix, modified_gram_schmidt, diag_mat

contains

    pure function eye(m, n, alpha)
    !> Create a matrix with ones in the diagonal and zero everywhere else
    !> \param m: number of rows
    !> \param n: number of colums
    !> \param alpha: optional diagonal value
    !> \return matrix of size n x m
    integer, intent(in) :: n, m
    real(dp), dimension(m, n) :: eye
    real(dp), intent(in), optional :: alpha
    
    !local variable
    integer :: i, j, min_dim
    real(dp) :: x

    ! check optional values
    x = 1.d0
    if (present(alpha)) x = alpha
    
    ! that's sooooo costly
    do i=1, m
       do j=1, n
          if (i /= j) then
             eye(i, j) = 0.d0
          else
             eye(i, i) = x
          end if
       end do
    end do

   ! why not that ?
   !  eye = 0.0_dp
   !  min_dim = min(m,n)
   !  do i=1, min_dim
   !    eye(i,i) = x
   !  end do
    
  end function eye

  pure function diag_mat(vec)
   real(dp), dimension(:), intent(in) :: vec
   real(dp), dimension(size(vec,1),size(vec,1)) :: diag_mat 

   integer :: i
   diag_mat = 0.0_dp
   do i=1,size(vec,1)
      diag_mat(i,i) = vec(i)
   end do

  end function diag_mat

  pure function norm(vector)
    !> compute the norm-2 of a vector
    real(dp), dimension(:), intent(in) :: vector
    real(dp) :: norm

    norm = sqrt(sum(vector ** 2.d0))

  end function norm
  
  subroutine concatenate(arr, brr)

    !> Concatenate two matrices
    !> \param arr: first array
    !> \param brr: second array
    !> \return arr concatenate brr (overwrites arr)

    real(dp), dimension(:, :), intent(inout), allocatable :: arr
    real(dp), dimension(:, :), intent(in) :: brr
    real(dp), dimension(:, :), allocatable :: tmp_array
    integer :: new_dim, dim_cols, dim_rows
    
    ! dimension
    dim_rows = size(arr, 1)
    dim_cols = size(arr, 2)

    ! Number of columns of the new matrix
    new_dim = dim_cols + size(brr, 2)

    ! move to temporal array
    allocate(tmp_array(dim_rows, new_dim))
    tmp_array(:, :dim_cols) = arr
   
    ! Move to new expanded matrix
    deallocate(arr)
    call move_alloc(tmp_array, arr)

    arr(:, dim_cols + 1:) = brr

  end subroutine concatenate
    
  function generate_diagonal_dominant(m, sparsity, diag_val) result(arr)
    !> Generate a diagonal dominant square matrix of dimension m
    !> \param m dimension of the matrix
    !> \param sparsity magnitude order of the off-diagonal values
      
    integer, intent(in) :: m ! size of the square matrix
    real(dp), optional :: diag_val
    integer :: i, j
    real(dp) :: sparsity 
    real(dp), dimension(m, m) :: arr
    call random_number(arr)

    arr = arr * sparsity
    do j=1, m
       do i=1, m
          if (i > j) then
             arr(i, j) = arr(j, i)
          else if(i == j) then
            if (present(diag_val))then
              arr(i,i) = diag_val
            else
             arr(i, i) = i
            end if
          end if
       end do
    end do

  end function generate_diagonal_dominant

  function diagonal(matrix)
    !> return the diagonal of a matrix
    real(dp), dimension(:, :), intent(in) :: matrix
    real(dp), dimension(size(matrix, 1)) :: diagonal

    ! local variables
    integer :: i, j, m

    ! dimension of the matrix
    m = size(matrix, 1)
    
    do i=1,m
       do j=1,m
          if  (i == j) then
             diagonal(i) = matrix(i, j)
          end if
       end do
    end do

  end function diagonal  



  function initialize_subspace(diag, dim_sub, dim_base) result(precond)
    !> \brief generates a diagonal preconditioner for `matrix`.
    !> \return diagonal matrix

    ! input variable
    real(dp), dimension(:), intent(inout) :: diag
    integer, intent(in) :: dim_sub, dim_base

    ! local variables
    real(dp), dimension(dim_base, dim_sub) :: precond
    integer, dimension(size(diag)) :: keys
    integer :: i, k
    
    ! sort diagonal
    keys = lapack_sort('I', diag)
    ! Fill matrix with zeros
    precond = 0.0_dp

    ! Add one depending on the order of the matrix diagonal
    do i=1, dim_sub
       write(6,*) k
       k = search_key(keys, i)
       precond(k, i) = 1.d0
    end do
    
  end function initialize_subspace

  subroutine modified_gram_schmidt(mat, nstart)
   !> \brief use modifed gram-schmidt orthogonalization on mat
   !> \brief nstart is the index of the first vector to orthogonalize
   
   ! input
   real(dp), dimension(:,:), intent(inout) :: mat
   integer, optional, intent(in) :: nstart

   integer :: i
   integer :: nrows, ncols
   integer :: idx_start
   real(dp), dimension(:), allocatable :: tmp_array

   idx_start = 1
   if (present(nstart)) idx_start = nstart

   nrows = size(mat,1)
   ncols = size(mat,2)

   allocate(tmp_array(nrows))

   do i=idx_start, ncols

      tmp_array = 0.0_dp
      tmp_array = lapack_matrix_vector('T', mat(:,:i-1),mat(:,i))
      mat(:,i) = mat(:,i) - lapack_matrix_vector('N',mat(:,:i-1), tmp_array)
      mat(:,i) = mat(:,i) / norm(mat(:,i))

   end do
   
   deallocate(tmp_array)

  end subroutine modified_gram_schmidt

  function search_key(keys, i) result(k)
    !> \brief Search for a given index `i` in a vector `keys`
    !> \param keys Vector of index
    !> \param i Index to search for
    !> \return index of i inside keys

    integer, dimension(:), intent(in) :: keys
    integer, intent(in) :: i
    integer :: j, k
    
    do j=1,size(keys)
       if (keys(j) == i) then
          k = j
          exit
       end if
    end do

  end function search_key

  subroutine write_matrix(path_file, mtx)
    !> Write matrix to path_file
    character(len=*), intent(in) :: path_file
    real(dp), dimension(:, :), intent(in) :: mtx
    integer :: i, j

    open(unit=314, file=path_file, status="REPLACE")
    do i=1,size(mtx, 1)
       do j=1,size(mtx, 2)
          write(314, *) mtx(i, j)
       end do
    end do
    close(314)
    
  end subroutine write_matrix

  subroutine write_vector(path_file, vector)
    !> Write vector to path_file
    character(len=*), intent(in) :: path_file
    real(dp), dimension(:), intent(in) :: vector
    integer :: i

    open(unit=314, file=path_file, status="REPLACE")
    do i=1,size(vector)
       write(314, *) vector(i)
    end do
    close(314)
    
  end subroutine write_vector  

!   subroutine check_deallocate_matrices( mtx_proj, stx_proj, lambda, eigenvectors_sub, &
!                    ritz_vectors, mtxV, stxV)
!     !> deallocate a matrix if allocated
!     real(dp), dimension(:, :), allocatable, intent(inout) ::  mtx_proj
!     real(dp), dimension(:, :), allocatable, intent(inout) ::  stx_proj
!     real(dp), dimension(:, :), allocatable, intent(inout) ::  lambda
!     real(dp), dimension(:, :), allocatable, intent(inout) ::  eigenvectors_sub
!     real(dp), dimension(:, :), allocatable, intent(inout) ::  ritz_vectors
!     real(dp), dimension(:, :), allocatable, intent(inout), optional ::  mtxV
!     real(dp), dimension(:, :), allocatable, intent(inout), optional ::  stxV
!       call check_deallocate_matrix( mtx_proj)
!       call check_deallocate_matrix( stx_proj)
!       call check_deallocate_matrix( lambda)
!       call check_deallocate_matrix( eigenvectors_sub)
!       call check_deallocate_matrix( ritz_vectors)
!       if( present( mtxV))  call check_deallocate_matrix( mtxV)
!       if( present( stxV))  call check_deallocate_matrix( stxV)

!   end subroutine check_deallocate_matrices

  subroutine check_deallocate_matrix(mtx)

    !> deallocate a matrix if allocated
    real(dp), dimension(:, :), allocatable, intent(inout) ::  mtx

    if (allocated(mtx)) then
       deallocate(mtx)
    end if

  end subroutine check_deallocate_matrix



  subroutine check_deallocate_vector(vec)

   !> deallocate a matrix if allocated
   real(dp), dimension(:), allocatable, intent(inout) ::  vec

   if (allocated(vec)) then
      deallocate(vec)
   end if

 end subroutine check_deallocate_vector
  
end module array_utils
