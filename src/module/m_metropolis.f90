module metropolis
    !> Module Metropolis.
    !> @param imetro
    !> @param delta
    !> @param deltai
    !> @param deltar
    !> @param deltat
    !> @param fbias
    !> @author Ravindra Shinde
    use precision_kinds, only: dp
    implicit none

    integer     :: imetro

    real(dp)    :: delta
    real(dp)    :: deltai
    real(dp)    :: deltar
    real(dp)    :: deltat
    real(dp)    :: fbias

    private
    public      :: imetro, delta, deltai, deltar, deltat, fbias
    save
end module
