module const
    !> Arguments: pi, hb, etrial, delta, deltai, fbias, nelec, imetro, ipr
    use precision_kinds, only: dp
#ifdef QMCKL_FOUND
    use qmckl
#endif
    
    implicit none

    real(dp) :: delta
    real(dp) :: deltai
    real(dp) :: etrial
    real(dp) :: fbias
    real(dp) :: hb
    integer  :: imetro
    integer  :: ipr
    integer  :: nelec
    real(dp) :: pi = 4.0d0*datan(1.0d0)

    private
    public   :: pi, hb, etrial, delta, deltai, fbias, nelec, imetro, ipr

#ifdef QMCKL_FOUND
    integer(qmckl_context), public :: qmckl_ctx
    logical, public :: use_qmckl = .True.
#else
    integer, public :: qmckl_ctx = 0
    integer, public :: QMCKL_SUCCESS = 0
    logical, public :: use_qmckl = .False.
#endif
    
    save
end module const

module const2
    !> Arguments: deltar, deltat
    use precision_kinds, only: dp

    implicit none

    real(dp) :: deltar
    real(dp) :: deltat

    private
    public :: deltar, deltat
    save
end module const2

module constant
    !> Arguments: twopi
    use precision_kinds, only: dp

    implicit none

    real(dp) :: twopi

    private
    public :: twopi
    save
end module constant
