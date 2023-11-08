    subroutine ortho_orbitals(mo_num,Tmat,Umat)

      use precision_kinds, only: dp

        implicit none

        integer, intent(in) :: mo_num
        real(dp), intent(inout) :: Tmat(mo_num,mo_num)
        real(dp), intent(inout) :: Umat(mo_num,mo_num)

        real(dp) :: f
        real(dp), dimension(mo_num,mo_num) :: Tpotmat
        real(dp), dimension(mo_num,mo_num) :: Tpotmat2

        integer :: i, iter
        logical :: converged

        Tpotmat(:,:)=0.d0
        Umat(:,:)   =0.d0
        do i=1,mo_num
          Tpotmat(i,i)=1.d0
          Umat(i,i)   =1.d0
        end do
        iter=0
        converged=.false.
        do while (.not.converged)
          iter=iter+1
          f = 1.d0 / dble(iter)
          Tpotmat2(:,:) = Tpotmat(:,:) * f
          call dgemm('N','N', mo_num,mo_num,mo_num,1.d0,                   &
              Tpotmat2, size(Tpotmat2,1),                                  &
              Tmat, size(Tmat,1), 0.d0,                                    &
              Tpotmat, size(Tpotmat,1))
          Umat(:,:) = Umat(:,:) + Tpotmat(:,:)

          converged = ( sum(abs(Tpotmat(:,:))) < 1.d-6).or.(iter>30)
        end do
        if(.not.converged) write(6,'(''Orthogonalization of orbitals failed'')')

  return
  end
