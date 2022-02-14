subroutine slm(l,rvec,r2,y,dy,ddy,ddy_lap,dlapy,ider)
!> @brief Compute the zero-th, first, second, and third derivatives of a Ylm function
!> @param[in] l The descriptor for angular momentum quantum number
!> @param[in] rvec The vector of radial coordinates
!> @param[in] r2 The square of the radial coordinate
!> @param[in] ider The order of the derivative
!> @param[out] y Value of the Ylm function
!> @param[out] dy Value of the first derivative of the Ylm function
!> @param[out] ddy Value of the second derivative of the Ylm function
!> @param[out] Value of the third derivative of the Ylm function
!> @param[out] ddy_lap Value of the Laplacian of the Ylm function
!> @param[out] dlapy Value of the derivative of the Laplacian of the Ylm function
!> @author Claudia Filippi
!> @author Ravindra Shinde
!> @date 2022-02-11
!> @details
!>   l   |   1  2  3  4    5   6   7   8   9   10   11  12  13  14  15  16  17  18  19  20
!> ------+-------------------------------------------------------------------------------------
!>   y   |   1  x  y  z    xx  xy  xz  yy  yz  zz   xxx xxy xxz xyy xyz xzz yyy yyz yzz zzz
!>

     use precision_kinds, only: dp
     implicit none

     integer :: i, j, l, ider
     real(dp) :: ddy_lap
     real(dp) :: r2, y
     real(dp), dimension(3) :: rvec
     real(dp), dimension(3) :: dy
     real(dp), dimension(3,3) :: ddy
     real(dp), dimension(3) :: dlapy
     real(dp), parameter :: half = 0.5d0

     real(dp), parameter :: cs = 0.28209479177387814d0      ! 1/sqrt(4*pi)
     real(dp), parameter :: cp = 0.48860251190291992d0      ! sqrt(3/(4*pi))
     real(dp), parameter :: cd1 = 0.630783130505040d0       ! sqrt(5/(4*pi))
     real(dp), parameter :: cd2 = 1.092548430592079d0       ! sqrt(15/(4*pi))

     real(dp), parameter :: cf = 0.74635266518023078d0      ! sqrt(7/(4*pi))
     real(dp), parameter :: cf2 = 1.66889529453113635d0     ! sqrt(35/(4*pi))
     real(dp), parameter :: cf3 = 2.89061144264055405d0     ! sqrt(105/(4*pi))

     ! Initialize all the elements of the arrays to zero
     ddy_lap = 0.d0
     dy = 0.d0           ! dy(1:3) = 0.0
     dlapy = 0.d0        ! dlapy(1:3) = 0.0
     ddy = 0.d0          ! ddy(1:3,1:3) = 0.0


     if (ider.ne.3) then

          select case(l)
               case(1)
                    ! S type
                    y = cs
                    return

               case(2)
                    ! Px type
                    y = cp*rvec(1)
                    dy(1) = cp
                    return

               case(3)
                    ! Py type
                    y = cp*rvec(2)
                    dy(2) = cp
                    return

               case(4)
                    ! Pz type
                    y = cp*rvec(3)
                    dy(3) = cp
                    return

               case(5)
                    ! Dxx type
                    y = cd2*rvec(1)*rvec(1)
                    dy(1) = 2.d0*cd2*rvec(1)
                    ddy(1,1) = 2.d0*cd2
                    return

               case(6)
                    ! Dxy type
                    y = cd2*rvec(1)*rvec(2)
                    dy(1) = cd2*rvec(2)
                    dy(2) = cd2*rvec(1)
                    ddy(1,2) = cd2
                    ddy(2,1) = cd2
                    return

               case(7)
                    ! Dxz type
                    y = cd2*rvec(1)*rvec(3)
                    dy(1) = cd2*rvec(3)
                    dy(3) = cd2*rvec(1)
                    ddy(1,3) = cd2
                    ddy(3,1) = cd2
                    return

               case(8)
                    ! Dyy type
                    y = cd2*rvec(2)*rvec(2)
                    dy(2) = 2.d0*cd2*rvec(2)
                    ddy(2,2) = 2.d0*cd2
                    return

               case(9)
                    ! Dyz type
                    y = cd2*rvec(2)*rvec(3)
                    dy(2) = cd2*rvec(3)
                    dy(3) = cd2*rvec(2)
                    ddy(2,3) = cd2
                    ddy(3,2) = cd2
                    return

               case(10)
                    ! Dzz type
                    y = cd2*rvec(3)*rvec(3)
                    dy(3) = 2.d0*cd2*rvec(3)
                    ddy(3,3) = 2.d0*cd2
                    return

               case(11)
                    ! Fxxx type
                    y=cf*rvec(1)*rvec(1)*rvec(1)
                    dy(1)=cf*3.d0*rvec(1)*rvec(1)
                    ddy_lap=cf*6.d0*rvec(1)
                    return

               case(12)
                    ! Fxxy type
                    y=cf2*rvec(1)*rvec(1)*rvec(2)
                    dy(1)=cf2*2.0d0*rvec(1)*rvec(2)
                    dy(2)=cf2*rvec(1)*rvec(1)
                    ddy_lap=cf2*2.0d0*rvec(2)
                    return

               case(13)
                    ! Fxxz type
                    y=cf2*rvec(1)*rvec(1)*rvec(3)
                    dy(1)=cf2*2.0d0*rvec(1)*rvec(3)
                    dy(3)=cf2*rvec(1)*rvec(1)
                    ddy_lap=cf2*2.0d0*rvec(3)
                    return

               case(14)
                    ! Fxyy type
                    y=cf2*rvec(1)*rvec(2)*rvec(2)
                    dy(1)=cf2*rvec(2)*rvec(2)
                    dy(2)=cf2*2.0d0*rvec(1)*rvec(2)
                    ddy_lap=cf2*2.0d0*rvec(2)
                    return

               case(15)
                    ! Fxyz type
                    y=cf3*rvec(1)*rvec(2)*rvec(3)
                    dy(1)=cf3*rvec(2)*rvec(3)
                    dy(2)=cf3*rvec(1)*rvec(3)
                    dy(3)=cf3*rvec(1)*rvec(2)
                    return

               case(16)
                    ! Fxzz type
                    y=cf2*rvec(1)*rvec(3)*rvec(3)
                    dy(1)=cf2*rvec(3)*rvec(3)
                    dy(3)=cf2*2.0d0*rvec(3)*rvec(1)
                    ddy_lap=cf2*2.0d0*rvec(1)
                    return

               case(17)
                    ! Fyyy type
                    y=cf*rvec(2)*rvec(2)*rvec(2)
                    dy(2)=cf*3.0d0*rvec(2)*rvec(2)
                    ddy_lap=cf*6.0d0*rvec(2)
                    return

               case(18)
                    ! Fyyz type
                    y=cf2*rvec(2)*rvec(2)*rvec(3)
                    dy(2)=cf2*2.0d0*rvec(3)*rvec(3)
                    dy(3)=cf2*rvec(2)*rvec(2)
                    ddy_lap=cf2*2.0d0*rvec(3)
                    return

               case(19)
                    ! Fyzz type
                    y=cf2*rvec(2)*rvec(3)*rvec(3)
                    dy(2)=cf2*rvec(3)*rvec(3)
                    dy(3)=cf2*2.0d0*rvec(3)*rvec(2)
                    ddy_lap=cf2*2.0d0*rvec(2)
                    return

               case(20)
                    ! Fzzz type
                    y=cf*rvec(3)*rvec(3)*rvec(3)
                    dy(3)=cf*3.0d0*rvec(3)*rvec(3)
                    ddy_lap=cf*6.0d0*rvec(3)
                    return

               case default
                    stop 'Ylm: AOs of type s, p, d, and f are implemented only'
          end select

     elseif (ider .eq. 3) then
          select case(l)
               case(1)
                    ! S type
                    y = cs
                    return

               case(2)
                    ! Px type
                    y = cp*rvec(1)
                    dy(1) = cp
                    return

               case(3)
                    ! Py type
                    y = cp*rvec(2)
                    dy(2) = cp
                    return

               case(4)
                    ! Pz type
                    y = cp*rvec(3)
                    dy(3) = cp
                    return

               case(5)
                    ! Dxx type
                    y = cd2*rvec(1)*rvec(1)
                    dy(1) = 2.d0*cd2*rvec(1)
                    ddy(1,1) = 2.d0*cd2
                    return

               case(6)
                    ! Dxy type
                    y = cd2*rvec(1)*rvec(2)
                    dy(1) = cd2*rvec(2)
                    dy(2) = cd2*rvec(1)
                    ddy(1,2) = cd2
                    ddy(2,1) = cd2
                    return

               case(7)
                    ! Dxz type
                    y = cd2*rvec(1)*rvec(3)
                    dy(1) = cd2*rvec(3)
                    dy(3) = cd2*rvec(1)
                    ddy(1,3) = cd2
                    ddy(3,1) = cd2
                    return

               case(8)
                    ! Dyy type
                    y = cd2*rvec(2)*rvec(2)
                    dy(2) = 2.d0*cd2*rvec(2)
                    ddy(2,2) = 2.d0*cd2
                    return

               case(9)
                    ! Dyz type
                    y = cd2*rvec(2)*rvec(3)
                    dy(2) = cd2*rvec(3)
                    dy(3) = cd2*rvec(2)
                    ddy(2,3) = cd2
                    ddy(3,2) = cd2
                    return

               case(10)
                    ! Dzz type
                    y = cd2*rvec(3)*rvec(3)
                    dy(3) = 2.d0*cd2*rvec(3)
                    ddy(3,3) = 2.d0*cd2
                    return

               case(11)
                    ! Fxxx type
                    y=cf*rvec(1)*rvec(1)*rvec(1)
                    dy(1)=cf*3.d0*rvec(1)*rvec(1)
                    ddy_lap=cf*6.d0*rvec(1)
                    ddy(1,1)=ddy_lap
                    dlapy(1)=cf*6.0d0
                    return

               case(12)
                    ! Fxxy type
                    y=cf2*rvec(1)*rvec(1)*rvec(2)
                    dy(1)=cf2*2.0d0*rvec(1)*rvec(2)
                    dy(2)=cf2*rvec(1)*rvec(1)
                    ddy_lap=cf2*2.0d0*rvec(2)
                    ddy(1,1)=ddy_lap
                    ddy(1,2)=cf2*2.0d0*rvec(1)
                    ddy(2,1)=ddy(1,2)
                    dlapy(2)=cf2*2.0d0
                    return

               case(13)
                    ! Fxxz type
                    y=cf2*rvec(1)*rvec(1)*rvec(3)
                    dy(1)=cf2*2.0d0*rvec(1)*rvec(3)
                    dy(3)=cf2*rvec(1)*rvec(1)
                    ddy_lap=cf2*2.0d0*rvec(3)
                    ddy(1,1)=ddy_lap
                    ddy(1,3)=cf2*2.0d0*rvec(1)
                    ddy(3,1)=ddy(1,3)
                    dlapy(3)=cf2*2.0d0
                    return

               case(14)
                    ! Fxyy type
                    y=cf2*rvec(1)*rvec(2)*rvec(2)
                    dy(1)=cf2*rvec(2)*rvec(2)
                    dy(2)=cf2*2.0d0*rvec(1)*rvec(2)
                    ddy_lap=cf2*2.0d0*rvec(1)
                    ddy(2,2)=ddy_lap
                    ddy(1,2)=cf2*2.0d0*rvec(2)
                    ddy(2,1)=ddy(1,2)
                    dlapy(1)=cf2*2.0d0
                    return

               case(15)
                    ! Fxyz type
                    y=cf3*rvec(1)*rvec(2)*rvec(3)
                    dy(1)=cf3*rvec(2)*rvec(3)
                    dy(2)=cf3*rvec(1)*rvec(3)
                    dy(3)=cf3*rvec(1)*rvec(2)
                    ddy(1,2)=cf3*rvec(3)
                    ddy(2,1)=ddy(1,2)
                    ddy(1,3)=cf3*rvec(2)
                    ddy(3,1)=ddy(1,3)
                    ddy(2,3)=cf3*rvec(1)
                    ddy(3,2)=ddy(2,3)
                    return

               case(16)
                    ! Fxzz type
                    y=cf2*rvec(1)*rvec(3)*rvec(3)
                    dy(1)=cf2*rvec(3)*rvec(3)
                    dy(3)=cf2*2.0d0*rvec(3)*rvec(1)
                    ddy_lap=cf2*2.0d0*rvec(1)
                    ddy(3,3)=ddy_lap
                    ddy(1,3)=cf2*2.0d0*rvec(3)
                    ddy(3,1)=ddy(1,3)
                    dlapy(1)=cf2*2.0d0
                    return

               case(17)
                    ! Fyyy type
                    y=cf*rvec(2)*rvec(2)*rvec(2)
                    dy(2)=cf*3.0d0*rvec(2)*rvec(2)
                    ddy_lap=cf*6.0d0*rvec(2)
                    ddy(2,2)=ddy_lap
                    dlapy(2)=cf*6.0d0
                    return

               case(18)
                    ! Fyyz type
                    y=cf2*rvec(2)*rvec(2)*rvec(3)
                    dy(2)=cf2*2.0d0*rvec(2)*rvec(3)
                    dy(3)=cf2*rvec(2)*rvec(2)
                    ddy_lap=cf2*2.0d0*rvec(3)
                    ddy(2,2)=ddy_lap
                    ddy(2,3)=cf2*2.0d0*rvec(2)
                    ddy(3,2)=ddy(2,3)
                    dlapy(3)=cf2*2.0d0
                    return

               case(19)
                    ! Fyzz type
                    y=cf2*rvec(2)*rvec(3)*rvec(3)
                    dy(2)=cf2*rvec(3)*rvec(3)
                    dy(3)=cf2*rvec(3)*rvec(2)
                    ddy_lap=cf2*2.0d0*rvec(2)
                    ddy(3,3)=ddy_lap
                    ddy(2,3)=cf2*2.0d0*rvec(3)
                    ddy(3,2)=ddy(2,3)
                    dlapy(2)=cf2*2.0d0
                    return

               case(20)
                    ! Fzzz type
                    y=cf*rvec(3)*rvec(3)*rvec(3)
                    dy(3)=cf*rvec(3)*rvec(3)
                    ddy_lap=cf*6.0d0*rvec(3)
                    ddy(3,3)=ddy_lap
                    dlapy(3)=cf*6.0d0
                    return

               case default
                    stop 'Ylm: AOs of type s, p, d, and f are implemented only'
          end select

     endif

end subroutine
