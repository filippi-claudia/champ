module constants
    !> Module Constants. All the constants used in the code will be listed here.
    !> @param hb
    !> @param pi
    !> @param twopi
    !> @author Ravindra Shinde
    !> @email r.l.shinde@utwente.nl
    use precision_kinds, only: dp
    implicit none

    real(dp)  :: hb = 0.5d0
    real(dp), parameter :: pi = 4.0d0*datan(1.0d0)
    real(dp), parameter :: twopi = 2.d0*pi

    private
    public    :: hb, pi, twopi
    save
end module constants
