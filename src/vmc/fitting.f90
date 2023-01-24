    module fitting_methods

        implicit none

        private
        public  :: exp_fit

        contains

        subroutine exp_fit(x, y, n, a, b)
!! This subroutine takes an array of floating point numbers and fits them to
!! an exponential function in the form of y = a*exp(-b*x).
!! @Author: Ravindra Shinde
!! @Date:   Tue May 18 13:54:00 EDT 2022
!> \param[in] x: Array of floating point numbers.
!> \param[in] y: Array of floating point numbers.
!> \param[in] n: Number of elements in x or y.
!> \param[out] a: Fitted parameter a.
!> \param[out] b: Fitted parameter b.

        use error,   only: fatal_error
        implicit none
        integer, parameter :: dp = kind(1.0d0)
        integer, intent(in) :: n
        real(dp), intent(in) :: x(n), y(n)
        real(dp), intent(out) :: a, b
        real(dp) :: x_sum, y_sum, xx_sum, xw_sum, w_sum, denom
        real(dp) :: x_mean, y_mean, w_mean, xw_mean, xx_mean
        real(dp) :: w(n)

        integer :: i

        ! y = a*exp(-b*x)
        ! take logarithms of both sides
        ! ln(y) = w(y) = A + B*x
        ! Then do least square fitting for a and b
        ! The final value of a would be exp(A) and b would be B

        x_sum = 0.0
        y_sum = 0.0
        w_sum = 0.0
        xw_sum = 0.0
        xx_sum = 0.0

        do i = 1, n
            
            w(i) = log(dabs(y(i)))
            x_sum = x_sum + x(i)
            y_sum = y_sum + y(i)
            w_sum = w_sum + w(i)
            xw_sum = xw_sum + x(i)*w(i)
            xx_sum = xx_sum + x(i)*x(i)

            if(i.lt.n.and.y(i+1)*y(i).lt.0.d0) call fatal_error('FIT NUM BASIS: changing sign')


        end do

        x_mean = x_sum/n
        y_mean = y_sum/n
        w_mean = w_sum/n
        xw_mean = xw_sum/n
        xx_mean = xx_sum/n

        denom = xx_sum - n*x_mean*x_mean
        b = (n*x_mean*w_mean - xw_sum )/denom
        a = dexp((w_mean*xx_sum - x_mean * xw_sum)/denom)

        if(y(1).lt.0) a=-a

    end subroutine exp_fit
end module
