!
!    ------------------------------------------
!     QMCkl - Quantum Monte Carlo kernel library
!     ------------------------------------------
!    
!     Documentation : https://trex-coe.github.io/qmckl
!     Issues        : https://github.com/trex-coe/qmckl/issues
!    
!     BSD 3-Clause License
!     
!     Copyright (c) 2020, TREX Center of Excellence
!     All rights reserved.
!     
!     Redistribution and use in source and binary forms, with or without
!     modification, are permitted provided that the following conditions are met:
!     
!     1. Redistributions of source code must retain the above copyright notice, this
!        list of conditions and the following disclaimer.
!     
!     2. Redistributions in binary form must reproduce the above copyright notice,
!        this list of conditions and the following disclaimer in the documentation
!        and/or other materials provided with the distribution.
!     
!     3. Neither the name of the copyright holder nor the names of its
!        contributors may be used to endorse or promote products derived from
!        this software without specific prior written permission.
!     
!     THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
!     AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
!     IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
!     DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
!     FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
!     DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
!     SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
!     CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
!     OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
!     OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!     
!     
!    
!    
!
module qmckl
  use, intrinsic :: iso_c_binding
integer  , parameter :: qmckl_context = c_int64_t
integer*8, parameter :: QMCKL_NULL_CONTEXT = 0
integer  , parameter :: qmckl_exit_code = c_int32_t

integer(qmckl_exit_code), parameter :: QMCKL_SUCCESS                  = 0
integer(qmckl_exit_code), parameter :: QMCKL_INVALID_ARG_1            = 1
integer(qmckl_exit_code), parameter :: QMCKL_INVALID_ARG_2            = 2
integer(qmckl_exit_code), parameter :: QMCKL_INVALID_ARG_3            = 3
integer(qmckl_exit_code), parameter :: QMCKL_INVALID_ARG_4            = 4
integer(qmckl_exit_code), parameter :: QMCKL_INVALID_ARG_5            = 5
integer(qmckl_exit_code), parameter :: QMCKL_INVALID_ARG_6            = 6
integer(qmckl_exit_code), parameter :: QMCKL_INVALID_ARG_7            = 7
integer(qmckl_exit_code), parameter :: QMCKL_INVALID_ARG_8            = 8
integer(qmckl_exit_code), parameter :: QMCKL_INVALID_ARG_9            = 9
integer(qmckl_exit_code), parameter :: QMCKL_INVALID_ARG_10           = 10
integer(qmckl_exit_code), parameter :: QMCKL_INVALID_ARG_11           = 11
integer(qmckl_exit_code), parameter :: QMCKL_INVALID_ARG_12           = 12
integer(qmckl_exit_code), parameter :: QMCKL_INVALID_ARG_13           = 13
integer(qmckl_exit_code), parameter :: QMCKL_INVALID_ARG_14           = 14
integer(qmckl_exit_code), parameter :: QMCKL_INVALID_ARG_15           = 15
integer(qmckl_exit_code), parameter :: QMCKL_INVALID_ARG_16           = 16
integer(qmckl_exit_code), parameter :: QMCKL_INVALID_ARG_17           = 17
integer(qmckl_exit_code), parameter :: QMCKL_INVALID_ARG_18           = 18
integer(qmckl_exit_code), parameter :: QMCKL_INVALID_ARG_19           = 19
integer(qmckl_exit_code), parameter :: QMCKL_INVALID_ARG_20           = 20
integer(qmckl_exit_code), parameter :: QMCKL_FAILURE                  = 101
integer(qmckl_exit_code), parameter :: QMCKL_ERRNO                    = 102
integer(qmckl_exit_code), parameter :: QMCKL_INVALID_CONTEXT          = 103
integer(qmckl_exit_code), parameter :: QMCKL_ALLOCATION_FAILED        = 104
integer(qmckl_exit_code), parameter :: QMCKL_DEALLOCATION_FAILED      = 105
integer(qmckl_exit_code), parameter :: QMCKL_NOT_PROVIDED             = 106
integer(qmckl_exit_code), parameter :: QMCKL_OUT_OF_BOUNDS            = 107
integer(qmckl_exit_code), parameter :: QMCKL_INVALID_EXIT_CODE        = 108
! Fortran interface


interface
  integer(c_int32_t) function qmckl_set_ao_basis_type (context, &
       basis_type) bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none
    integer (c_int64_t) , intent(in)  , value :: context
    character(c_char)   , intent(in)  , value :: basis_type
  end function qmckl_set_ao_basis_type
end interface

interface
  integer(c_int32_t) function qmckl_set_ao_basis_shell_num(context, &
       num) bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none
    integer (c_int64_t) , intent(in)  , value :: context
    integer (c_int64_t) , intent(in)  , value :: num
  end function qmckl_set_ao_basis_shell_num
end interface

interface
  integer(c_int32_t) function qmckl_set_ao_basis_prim_num(context, &
       num) bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none
    integer (c_int64_t) , intent(in)  , value :: context
    integer (c_int64_t) , intent(in)  , value :: num
  end function qmckl_set_ao_basis_prim_num
end interface

interface
  integer(c_int32_t) function qmckl_set_ao_basis_nucleus_index(context, &
       idx, size_max) bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none
    integer (c_int64_t) , intent(in)  , value :: context
    integer (c_int64_t) , intent(in)          :: idx(*)
    integer (c_int64_t) , intent(in)  , value :: size_max
  end function qmckl_set_ao_basis_nucleus_index
end interface

interface
  integer(c_int32_t) function qmckl_set_ao_basis_nucleus_shell_num(context, &
       shell_num, size_max) bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none
    integer (c_int64_t) , intent(in)  , value :: context
    integer (c_int64_t) , intent(in)          :: shell_num(*)
    integer (c_int64_t) , intent(in)  , value :: size_max
  end function qmckl_set_ao_basis_nucleus_shell_num
end interface

interface
  integer(c_int32_t) function qmckl_set_ao_basis_shell_ang_mom(context, &
       shell_ang_mom, size_max) bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none
    integer (c_int64_t) , intent(in)  , value :: context
    integer (c_int32_t) , intent(in)          :: shell_ang_mom(*)
    integer (c_int64_t) , intent(in)  , value :: size_max
  end function qmckl_set_ao_basis_shell_ang_mom
end interface

interface
  integer(c_int32_t) function qmckl_set_ao_basis_shell_prim_num(context, &
       shell_prim_num, size_max) bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none
    integer (c_int64_t) , intent(in)  , value :: context
    integer (c_int64_t) , intent(in)          :: shell_prim_num(*)
    integer (c_int64_t) , intent(in)  , value :: size_max
  end function qmckl_set_ao_basis_shell_prim_num
end interface

interface
  integer(c_int32_t) function qmckl_set_ao_basis_shell_prim_index(context, &
       shell_prim_index, size_max) bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none
    integer (c_int64_t) , intent(in)  , value :: context
    integer (c_int64_t) , intent(in)          :: shell_prim_index(*)
    integer (c_int64_t) , intent(in)  , value :: size_max
  end function qmckl_set_ao_basis_shell_prim_index
end interface

interface
  integer(c_int32_t) function qmckl_set_ao_basis_shell_factor(context, &
       shell_factor, size_max) bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none
    integer (c_int64_t) , intent(in)  , value :: context
    real    (c_double)  , intent(in)          :: shell_factor(*)
    integer (c_int64_t) , intent(in)  , value :: size_max
  end function qmckl_set_ao_basis_shell_factor
end interface

interface
  integer(c_int32_t) function qmckl_set_ao_basis_exponent(context, &
       exponent, size_max) bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none
    integer (c_int64_t) , intent(in)  , value :: context
    real    (c_double)  , intent(in)          :: exponent(*)
    integer (c_int64_t) , intent(in)  , value :: size_max
  end function qmckl_set_ao_basis_exponent
end interface

interface
  integer(c_int32_t) function qmckl_set_ao_basis_coefficient(context, &
       coefficient, size_max)  bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none
    integer (c_int64_t) , intent(in)  , value :: context
    real    (c_double)  , intent(in)          :: coefficient(*)
    integer (c_int64_t) , intent(in)  , value :: size_max
  end function qmckl_set_ao_basis_coefficient
end interface

interface
  integer(c_int32_t) function qmckl_set_ao_basis_prim_factor(context, &
       prim_factor, size_max) bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none
    integer (c_int64_t) , intent(in)  , value :: context
    real    (c_double)  , intent(in)          :: prim_factor(*)
    integer (c_int64_t) , intent(in)  , value :: size_max
  end function qmckl_set_ao_basis_prim_factor
end interface

interface
  integer(c_int32_t) function qmckl_set_ao_basis_ao_num(context, &
       num) bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none
    integer (c_int64_t) , intent(in)  , value :: context
    integer (c_int64_t) , intent(in)  , value :: num
  end function qmckl_set_ao_basis_ao_num
end interface

interface
  integer(c_int32_t) function qmckl_set_ao_basis_cartesian(context, &
       cartesian) bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none
    integer (c_int64_t) , intent(in)  , value :: context
    logical (c_bool)    , intent(in)  , value :: cartesian
  end function qmckl_set_ao_basis_cartesian
end interface

interface
  integer(c_int32_t) function qmckl_set_ao_basis_ao_factor(context, &
       ao_factor, size_max) bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none
    integer (c_int64_t) , intent(in)  , value :: context
    real    (c_double)  , intent(in)          :: ao_factor(*)
    integer (c_int64_t) , intent(in)  , value :: size_max
  end function qmckl_set_ao_basis_ao_factor
end interface

! Fortran interface


interface
  integer(c_int32_t) function qmckl_get_ao_basis_type (context, &
       basis_type) bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none
    integer (c_int64_t) , intent(in)  , value :: context
    character(c_char)   , intent(out)         :: basis_type
  end function qmckl_get_ao_basis_type
end interface

interface
  integer(c_int32_t) function qmckl_get_ao_basis_shell_num(context, &
       num) bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none
    integer (c_int64_t) , intent(in)  , value :: context
    integer (c_int64_t) , intent(out)         :: num
  end function qmckl_get_ao_basis_shell_num
end interface

interface
  integer(c_int32_t) function qmckl_get_ao_basis_prim_num(context, &
       num) bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none
    integer (c_int64_t) , intent(in)  , value :: context
    integer (c_int64_t) , intent(out)         :: num
  end function qmckl_get_ao_basis_prim_num
end interface

interface
  integer(c_int32_t) function qmckl_get_ao_basis_nucleus_shell_num(context, &
       shell_num, size_max) bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none
    integer (c_int64_t) , intent(in)  , value :: context
    integer (c_int64_t) , intent(out)         :: shell_num(*)
    integer (c_int64_t) , intent(in)  , value :: size_max
  end function qmckl_get_ao_basis_nucleus_shell_num
end interface

interface
  integer(c_int32_t) function qmckl_get_ao_basis_nucleus_index(context, &
       idx, size_max) bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none
    integer (c_int64_t) , intent(in)  , value :: context
    integer (c_int64_t) , intent(out)         :: idx(*)
    integer (c_int64_t) , intent(in)  , value :: size_max
  end function qmckl_get_ao_basis_nucleus_index
end interface

interface
  integer(c_int32_t) function qmckl_get_ao_basis_shell_ang_mom(context, &
       shell_ang_mom, size_max) bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none
    integer (c_int64_t) , intent(in)  , value :: context
    integer (c_int32_t) , intent(out)         :: shell_ang_mom(*)
    integer (c_int64_t) , intent(in)  , value :: size_max
  end function qmckl_get_ao_basis_shell_ang_mom
end interface

interface
  integer(c_int32_t) function qmckl_get_ao_basis_shell_prim_num(context, &
       shell_prim_num, size_max) bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none
    integer (c_int64_t) , intent(in)  , value :: context
    integer (c_int64_t) , intent(out)         :: shell_prim_num(*)
    integer (c_int64_t) , intent(in)  , value :: size_max
  end function qmckl_get_ao_basis_shell_prim_num
end interface

interface
  integer(c_int32_t) function qmckl_get_ao_basis_shell_prim_index(context, &
       shell_prim_index, size_max) bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none
    integer (c_int64_t) , intent(in)  , value :: context
    integer (c_int64_t) , intent(out)         :: shell_prim_index(*)
    integer (c_int64_t) , intent(in)  , value :: size_max
  end function qmckl_get_ao_basis_shell_prim_index
end interface

interface
  integer(c_int32_t) function qmckl_get_ao_basis_shell_factor(context, &
       shell_factor, size_max) bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none
    integer (c_int64_t) , intent(in)  , value :: context
    real    (c_double)  , intent(out)         :: shell_factor(*)
    integer (c_int64_t) , intent(in)  , value :: size_max
  end function qmckl_get_ao_basis_shell_factor
end interface

interface
  integer(c_int32_t) function qmckl_get_ao_basis_exponent(context, &
       exponent, size_max) bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none
    integer (c_int64_t) , intent(in)  , value :: context
    real    (c_double)  , intent(out)         :: exponent(*)
    integer (c_int64_t) , intent(in)  , value :: size_max
  end function qmckl_get_ao_basis_exponent
end interface

interface
  integer(c_int32_t) function qmckl_get_ao_basis_coefficient(context, &
       coefficient, size_max) bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none
    integer (c_int64_t) , intent(in)  , value :: context
    real    (c_double)  , intent(out)         :: coefficient(*)
    integer (c_int64_t) , intent(in)  , value :: size_max
  end function qmckl_get_ao_basis_coefficient
end interface

interface
  integer(c_int32_t) function qmckl_get_ao_basis_prim_factor(context, &
       prim_factor, size_max) bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none
    integer (c_int64_t) , intent(in)  , value :: context
    real    (c_double)  , intent(out)         :: prim_factor(*)
    integer (c_int64_t) , intent(in)  , value :: size_max
  end function qmckl_get_ao_basis_prim_factor
end interface

interface
  integer(c_int32_t) function qmckl_get_ao_basis_ao_num(context, &
       num) bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none
    integer (c_int64_t) , intent(in)  , value :: context
    integer (c_int64_t) , intent(out)         :: num
  end function qmckl_get_ao_basis_ao_num
end interface

interface
  integer(c_int32_t) function qmckl_get_ao_basis_cartesian(context, &
       cartesian) bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none
    integer (c_int64_t) , intent(in)  , value :: context
    logical (c_bool)    , intent(out)         :: cartesian
  end function qmckl_get_ao_basis_cartesian
end interface

interface
  integer(c_int32_t) function qmckl_get_ao_basis_ao_factor(context, &
       ao_factor, size_max) bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none
    integer (c_int64_t) , intent(in)  , value :: context
    real    (c_double)  , intent(out)         :: ao_factor(*)
    integer (c_int64_t) , intent(in)  , value :: size_max
  end function qmckl_get_ao_basis_ao_factor
end interface

interface
  integer(c_int32_t) function qmckl_get_ao_basis_primitive_vgl &
       (context, primitive_vgl, size_max) &
       bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none

    integer (c_int64_t) , intent(in)  , value :: context
    double precision,     intent(out)         :: primitive_vgl(*)
    integer (c_int64_t) , intent(in)  , value :: size_max
  end function
end interface

interface
  integer(c_int32_t) function qmckl_get_ao_basis_shell_vgl &
       (context, shell_vgl, size_max) &
       bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none

    integer (c_int64_t) , intent(in)  , value :: context
    double precision,     intent(out)         :: shell_vgl(*)
    integer (c_int64_t) , intent(in)  , value :: size_max
  end function
end interface

interface
   integer(c_int32_t) function qmckl_get_ao_basis_ao_vgl (context, &
        ao_vgl, size_max) bind(C)
     use, intrinsic :: iso_c_binding
     import
     implicit none
     integer (c_int64_t) , intent(in)  , value :: context
     double precision,     intent(out)         :: ao_vgl(*)
     integer (c_int64_t) , intent(in)  , value :: size_max
   end function qmckl_get_ao_basis_ao_vgl
end interface

interface
   integer(c_int32_t) function qmckl_get_ao_basis_ao_vgl_inplace (context, &
        ao_vgl, size_max) bind(C)
     use, intrinsic :: iso_c_binding
     import
     implicit none
     integer (c_int64_t) , intent(in)  , value :: context
     double precision,     intent(out)         :: ao_vgl(*)
     integer (c_int64_t) , intent(in)  , value :: size_max
   end function qmckl_get_ao_basis_ao_vgl_inplace
end interface

interface
   integer(c_int32_t) function qmckl_get_ao_basis_ao_value (context, &
        ao_value, size_max) bind(C)
     use, intrinsic :: iso_c_binding
     import
     implicit none
     integer (c_int64_t) , intent(in)  , value :: context
     double precision,     intent(out)         :: ao_value(*)
     integer (c_int64_t) , intent(in)  , value :: size_max
   end function qmckl_get_ao_basis_ao_value
end interface

interface
   integer(c_int32_t) function qmckl_get_ao_basis_ao_value_inplace (context, &
        ao_value, size_max) bind(C)
     use, intrinsic :: iso_c_binding
     import
     implicit none
     integer (c_int64_t) , intent(in)  , value :: context
     double precision,     intent(out)         :: ao_value(*)
     integer (c_int64_t) , intent(in)  , value :: size_max
   end function qmckl_get_ao_basis_ao_value_inplace
end interface

interface
   integer(c_int32_t) function qmckl_ao_gaussian_vgl(context, &
        X, R, n, A, VGL, ldv) bind(C)
     use, intrinsic :: iso_c_binding
     integer (c_int64_t) , intent(in) , value :: context
     integer (c_int64_t) , intent(in) , value :: ldv
     integer (c_int64_t) , intent(in) , value :: n
     real    (c_double)  , intent(in)         :: X(3), R(3), A(n)
     real    (c_double)  , intent(out)        :: VGL(ldv,5)
   end function qmckl_ao_gaussian_vgl
end interface



! #+CALL: generate_f_interface(table=qmckl_ao_power_args,rettyp=get_value("CRetType"),fname="qmckl_ao_power")

! #+RESULTS:

interface
  integer(c_int32_t) function qmckl_ao_power &
      (context, n, X, LMAX, P, ldp) &
      bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none

    integer (c_int64_t) , intent(in)  , value :: context
    integer (c_int64_t) , intent(in)  , value :: n
    real    (c_double ) , intent(in)          :: X(n)
    integer (c_int32_t) , intent(in)          :: LMAX(n)
    real    (c_double ) , intent(out)         :: P(ldp,n)
    integer (c_int64_t) , intent(in)  , value :: ldp

  end function qmckl_ao_power
end interface



! #+CALL: generate_f_interface(table=qmckl_ao_polynomial_vgl_args,rettyp=get_value("FRetType"),fname="qmckl_ao_polynomial_vgl_doc" )

! #+RESULTS:

interface
  integer(c_int32_t) function qmckl_ao_polynomial_vgl_doc &
      (context, X, R, lmax, n, L, ldl, VGL, ldv) &
      bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none

    integer (c_int64_t) , intent(in)  , value :: context
    real    (c_double ) , intent(in)          :: X(3)
    real    (c_double ) , intent(in)          :: R(3)
    integer (c_int32_t) , intent(in)  , value :: lmax
    integer (c_int64_t) , intent(inout)        :: n
    integer (c_int32_t) , intent(out)         :: L(ldl,n)
    integer (c_int64_t) , intent(in)  , value :: ldl
    real    (c_double ) , intent(out)         :: VGL(ldv,n)
    integer (c_int64_t) , intent(in)  , value :: ldv

  end function qmckl_ao_polynomial_vgl_doc
end interface



! #+CALL: generate_f_interface(table=qmckl_ao_polynomial_vgl_args,rettyp=get_value("FRetType"),fname="qmckl_ao_polynomial_vgl" )

! #+RESULTS:

interface
  integer(c_int32_t) function qmckl_ao_polynomial_vgl &
      (context, X, R, lmax, n, L, ldl, VGL, ldv) &
      bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none

    integer (c_int64_t) , intent(in)  , value :: context
    real    (c_double ) , intent(in)          :: X(3)
    real    (c_double ) , intent(in)          :: R(3)
    integer (c_int32_t) , intent(in)  , value :: lmax
    integer (c_int64_t) , intent(inout)        :: n
    integer (c_int32_t) , intent(out)         :: L(ldl,n)
    integer (c_int64_t) , intent(in)  , value :: ldl
    real    (c_double ) , intent(out)         :: VGL(ldv,n)
    integer (c_int64_t) , intent(in)  , value :: ldv

  end function qmckl_ao_polynomial_vgl
end interface



! #+CALL: generate_f_interface(table=qmckl_ao_polynomial_vgl_args,rettyp=get_value("FRetType"),fname="qmckl_ao_polynomial_transp_vgl_doc")

! #+RESULTS:

interface
  integer(c_int32_t) function qmckl_ao_polynomial_transp_vgl_doc &
      (context, X, R, lmax, n, L, ldl, VGL, ldv) &
      bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none

    integer (c_int64_t) , intent(in)  , value :: context
    real    (c_double ) , intent(in)          :: X(3)
    real    (c_double ) , intent(in)          :: R(3)
    integer (c_int32_t) , intent(in)  , value :: lmax
    integer (c_int64_t) , intent(inout)        :: n
    integer (c_int32_t) , intent(out)         :: L(ldl,n)
    integer (c_int64_t) , intent(in)  , value :: ldl
    real    (c_double ) , intent(out)         :: VGL(ldv,n)
    integer (c_int64_t) , intent(in)  , value :: ldv

  end function qmckl_ao_polynomial_transp_vgl_doc
end interface



! #+CALL: generate_f_interface(table=qmckl_ao_polynomial_vgl_args,rettyp=get_value("FRetType"),fname="qmckl_ao_polynomial_transp_vgl")

! #+RESULTS:

interface
  integer(c_int32_t) function qmckl_ao_polynomial_transp_vgl &
      (context, X, R, lmax, n, L, ldl, VGL, ldv) &
      bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none

    integer (c_int64_t) , intent(in)  , value :: context
    real    (c_double ) , intent(in)          :: X(3)
    real    (c_double ) , intent(in)          :: R(3)
    integer (c_int32_t) , intent(in)  , value :: lmax
    integer (c_int64_t) , intent(inout)        :: n
    integer (c_int32_t) , intent(out)         :: L(ldl,n)
    integer (c_int64_t) , intent(in)  , value :: ldl
    real    (c_double ) , intent(out)         :: VGL(ldv,n)
    integer (c_int64_t) , intent(in)  , value :: ldv

  end function qmckl_ao_polynomial_transp_vgl
end interface



! #+CALL: generate_f_interface(table=qmckl_dgemm_args,rettyp="qmckl_exit_code",fname="qmckl_dgemm")

! #+RESULTS:

interface
  integer(c_int32_t) function qmckl_dgemm &
      (context, TransA, TransB, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc) &
      bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none

    integer (c_int64_t) , intent(in)  , value :: context
    character           , intent(in)  , value :: TransA
    character           , intent(in)  , value :: TransB
    integer (c_int64_t) , intent(in)  , value :: m
    integer (c_int64_t) , intent(in)  , value :: n
    integer (c_int64_t) , intent(in)  , value :: k
    real    (c_double ) , intent(in)  , value :: alpha
    real    (c_double ) , intent(in)          :: A(lda,*)
    integer (c_int64_t) , intent(in)  , value :: lda
    real    (c_double ) , intent(in)          :: B(ldb,*)
    integer (c_int64_t) , intent(in)  , value :: ldb
    real    (c_double ) , intent(in)  , value :: beta
    real    (c_double ) , intent(out)         :: C(ldc,*)
    integer (c_int64_t) , intent(in)  , value :: ldc

  end function qmckl_dgemm
end interface



! #+CALL: generate_f_interface(table=qmckl_adjugate_args,rettyp="qmckl_exit_code",fname="qmckl_adjugate")

! #+RESULTS:

interface
  integer(c_int32_t) function qmckl_adjugate &
      (context, n, A, lda, B, ldb, det_l) &
      bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none

    integer (c_int64_t) , intent(in)  , value :: context
    integer (c_int64_t) , intent(in)  , value :: n
    real    (c_double ) , intent(in)          :: A(lda,*)
    integer (c_int64_t) , intent(in)  , value :: lda
    real    (c_double ) , intent(out)         :: B(ldb,*)
    integer (c_int64_t) , intent(in)  , value :: ldb
    real    (c_double ) , intent(inout)        :: det_l

  end function qmckl_adjugate
end interface
interface
   integer (qmckl_context) function qmckl_context_create() bind(C)
     use, intrinsic :: iso_c_binding
     import
   end function qmckl_context_create
end interface

! interface
!    integer (qmckl_context) function qmckl_context_copy(context) bind(C)
!      use, intrinsic :: iso_c_binding
!      import
!      integer (qmckl_context), intent(in), value :: context
!    end function qmckl_context_copy
! end interface

interface
   integer (qmckl_exit_code) function qmckl_context_destroy(context) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer (qmckl_context), intent(in), value :: context
   end function qmckl_context_destroy
end interface


! #+CALL: generate_f_interface(table=qmckl_distance_sq_args,rettyp=get_value("FRetType"),fname=get_value("Name"))

! #+RESULTS:

interface
  integer(c_int32_t) function qmckl_distance_sq &
      (context, transa, transb, m, n, A, lda, B, ldb, C, ldc) &
      bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none

    integer (c_int64_t) , intent(in)  , value :: context
    character           , intent(in)  , value :: transa
    character           , intent(in)  , value :: transb
    integer (c_int64_t) , intent(in)  , value :: m
    integer (c_int64_t) , intent(in)  , value :: n
    real    (c_double ) , intent(in)          :: A(lda,*)
    integer (c_int64_t) , intent(in)  , value :: lda
    real    (c_double ) , intent(in)          :: B(ldb,*)
    integer (c_int64_t) , intent(in)  , value :: ldb
    real    (c_double ) , intent(out)         :: C(ldc,n)
    integer (c_int64_t) , intent(in)  , value :: ldc

  end function qmckl_distance_sq
end interface



! #+CALL: generate_f_interface(table=qmckl_distance_args,rettyp=get_value("FRetType"),fname=get_value("Name"))

! #+RESULTS:

interface
  integer(c_int32_t) function qmckl_distance &
      (context, transa, transb, m, n, A, lda, B, ldb, C, ldc) &
      bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none

    integer (c_int64_t) , intent(in)  , value :: context
    character           , intent(in)  , value :: transa
    character           , intent(in)  , value :: transb
    integer (c_int64_t) , intent(in)  , value :: m
    integer (c_int64_t) , intent(in)  , value :: n
    real    (c_double ) , intent(in)          :: A(lda,*)
    integer (c_int64_t) , intent(in)  , value :: lda
    real    (c_double ) , intent(in)          :: B(ldb,*)
    integer (c_int64_t) , intent(in)  , value :: ldb
    real    (c_double ) , intent(out)         :: C(ldc,n)
    integer (c_int64_t) , intent(in)  , value :: ldc

  end function qmckl_distance
end interface



! #+CALL: generate_f_interface(table=qmckl_distance_rescaled_args,rettyp=get_value("FRetType"),fname=get_value("Name"))

! #+RESULTS:

interface
  integer(c_int32_t) function qmckl_distance_rescaled &
      (context, transa, transb, m, n, A, lda, B, ldb, C, ldc, rescale_factor_kappa) &
      bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none

    integer (c_int64_t) , intent(in)  , value :: context
    character           , intent(in)  , value :: transa
    character           , intent(in)  , value :: transb
    integer (c_int64_t) , intent(in)  , value :: m
    integer (c_int64_t) , intent(in)  , value :: n
    real    (c_double ) , intent(in)          :: A(lda,*)
    integer (c_int64_t) , intent(in)  , value :: lda
    real    (c_double ) , intent(in)          :: B(ldb,*)
    integer (c_int64_t) , intent(in)  , value :: ldb
    real    (c_double ) , intent(out)         :: C(ldc,n)
    integer (c_int64_t) , intent(in)  , value :: ldc
    real    (c_double ) , intent(in)  , value :: rescale_factor_kappa

  end function qmckl_distance_rescaled
end interface



! #+CALL: generate_f_interface(table=qmckl_distance_rescaled_deriv_e_args,rettyp=get_value("FRetType"),fname=get_value("Name"))

! #+RESULTS:

interface
  integer(c_int32_t) function qmckl_distance_rescaled_deriv_e &
      (context, transa, transb, m, n, A, lda, B, ldb, C, ldc, rescale_factor_kappa) &
      bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none

    integer (c_int64_t) , intent(in)  , value :: context
    character           , intent(in)  , value :: transa
    character           , intent(in)  , value :: transb
    integer (c_int64_t) , intent(in)  , value :: m
    integer (c_int64_t) , intent(in)  , value :: n
    real    (c_double ) , intent(in)          :: A(lda,*)
    integer (c_int64_t) , intent(in)  , value :: lda
    real    (c_double ) , intent(in)          :: B(ldb,*)
    integer (c_int64_t) , intent(in)  , value :: ldb
    real    (c_double ) , intent(out)         :: C(ldc,n,4)
    integer (c_int64_t) , intent(in)  , value :: ldc
    real    (c_double ) , intent(in)  , value :: rescale_factor_kappa

  end function qmckl_distance_rescaled_deriv_e
end interface
interface
  integer(c_int32_t) function qmckl_set_electron_num(context, alpha, beta) bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none

    integer (c_int64_t) , intent(in)  , value :: context
    integer (c_int64_t) , intent(in)  , value :: alpha
    integer (c_int64_t) , intent(in)  , value :: beta
  end function
end interface

interface
  integer(c_int32_t) function qmckl_set_electron_walk_num(context, walk_num) bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none

    integer (c_int64_t) , intent(in)  , value :: context
    integer (c_int64_t) , intent(in)  , value :: walk_num
  end function
end interface

interface
  integer(c_int32_t) function qmckl_set_electron_coord(context, transp, coord, size_max) bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none

    integer (c_int64_t) , intent(in)  , value :: context
    character           , intent(in)  , value :: transp
    double precision    , intent(in)          :: coord(*)
    integer (c_int64_t) , intent(in)  , value :: size_max
  end function
end interface

interface
  integer(c_int32_t) function qmckl_get_electron_ee_distance(context, distance) &
    bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none
    integer (c_int64_t) , intent(in)  , value :: context
    real    (c_double ) , intent(out)         :: distance(*)
  end function
end interface

interface
  integer(c_int32_t) function qmckl_get_electron_en_distance(context, distance) &
    bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none
    integer (c_int64_t) , intent(in)  , value :: context
    real    (c_double ) , intent(out)         :: distance(*)
  end function
end interface
interface
   subroutine qmckl_string_of_error (error, string) bind(C, name='qmckl_string_of_error_f')
     use, intrinsic :: iso_c_binding
     import
     integer (qmckl_exit_code), intent(in), value :: error
     character, intent(out) :: string(128)
   end subroutine qmckl_string_of_error
end interface
! Fortran interfaces


interface
  integer(c_int32_t) function qmckl_get_mo_basis_mo_num (context, &
       mo_num) bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none
    integer (c_int64_t) , intent(in)  , value :: context
    integer (c_int64_t) , intent(out)         :: mo_num
  end function qmckl_get_mo_basis_mo_num
end interface

interface
  integer(c_int32_t) function qmckl_get_mo_basis_coefficient(context, &
       coefficient, size_max) bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none
    integer (c_int64_t) , intent(in)  , value :: context
    double precision, intent(out)             :: coefficient(*)
    integer (c_int64_t) , intent(in), value   :: size_max
  end function qmckl_get_mo_basis_coefficient
end interface

interface
   integer(c_int32_t) function qmckl_get_mo_basis_mo_value (context, &
        mo_value, size_max) bind(C)
     use, intrinsic :: iso_c_binding
     import
     implicit none

     integer (c_int64_t) , intent(in)  , value :: context
     double precision,     intent(out)         :: mo_value(*)
     integer (c_int64_t) , intent(in)  , value :: size_max
   end function qmckl_get_mo_basis_mo_value
end interface

interface
   integer(c_int32_t) function qmckl_get_mo_basis_mo_value_inplace (context, &
        mo_value, size_max) bind(C)
     use, intrinsic :: iso_c_binding
     import
     implicit none
     integer (c_int64_t) , intent(in)  , value :: context
     double precision,     intent(out)         :: mo_value(*)
     integer (c_int64_t) , intent(in)  , value :: size_max
   end function qmckl_get_mo_basis_mo_value_inplace
end interface

interface
   integer(c_int32_t) function qmckl_get_mo_basis_mo_vgl (context, &
        mo_vgl, size_max) bind(C)
     use, intrinsic :: iso_c_binding
     import
     implicit none

     integer (c_int64_t) , intent(in)  , value :: context
     double precision,     intent(out)         :: mo_vgl(*)
     integer (c_int64_t) , intent(in)  , value :: size_max
   end function qmckl_get_mo_basis_mo_vgl
end interface

interface
   integer(c_int32_t) function qmckl_get_mo_basis_mo_vgl_inplace (context, &
        mo_vgl, size_max) bind(C)
     use, intrinsic :: iso_c_binding
     import
     implicit none
     integer (c_int64_t) , intent(in)  , value :: context
     double precision,     intent(out)         :: mo_vgl(*)
     integer (c_int64_t) , intent(in)  , value :: size_max
   end function qmckl_get_mo_basis_mo_vgl_inplace
end interface
interface
  integer(c_int32_t) function qmckl_get_nucleus_num(context, num) &
    bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none
    integer (c_int64_t) , intent(in)  , value :: context
    integer (c_int64_t) , intent(out)         :: num
  end function qmckl_get_nucleus_num
end interface

interface
  integer(c_int32_t) function qmckl_get_nucleus_charge(context, charge, size_max) &
    bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none
    integer (c_int64_t) , intent(in)  , value :: context
    real    (c_double)  , intent(out)         :: charge(*)
    integer (c_int64_t) , intent(in)  , value :: size_max
  end function qmckl_get_nucleus_charge
end interface

interface
  integer(c_int32_t) function qmckl_get_nucleus_rescale_factor(context, kappa) &
    bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none
    integer (c_int64_t) , intent(in)  , value :: context
    real    (c_double)  , intent(out)         :: kappa
  end function qmckl_get_nucleus_rescale_factor
end interface

interface
  integer(c_int32_t) function qmckl_get_nucleus_coord(context, transp, coord, size_max) &
    bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none
    integer (c_int64_t) , intent(in)  , value :: context
    character(c_char)   , intent(in)  , value :: transp
    real    (c_double)  , intent(out)         :: coord(*)
    integer (c_int64_t) , intent(in)  , value :: size_max
  end function qmckl_get_nucleus_coord
end interface

interface
  integer(c_int32_t) function qmckl_set_nucleus_num(context, num) &
    bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none
    integer (c_int64_t) , intent(in)  , value :: context
    integer (c_int64_t) , intent(in)  , value :: num
  end function qmckl_set_nucleus_num
end interface

interface
  integer(c_int32_t) function qmckl_set_nucleus_charge(context, charge, size_max) &
    bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none
    integer (c_int64_t) , intent(in)  , value :: context
    real    (c_double)  , intent(in)          :: charge(*)
    integer (c_int64_t) , intent(in)  , value :: size_max
  end function
end interface

interface
  integer(c_int32_t) function qmckl_set_nucleus_coord(context, transp, coord, size_max) &
    bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none
    integer (c_int64_t) , intent(in)  , value :: context
    character(c_char)   , intent(in)  , value :: transp
    real    (c_double)  , intent(in)          :: coord(*)
    integer (c_int64_t) , intent(in)  , value :: size_max
  end function qmckl_set_nucleus_coord
end interface

interface
  integer(c_int32_t) function qmckl_set_nucleus_rescale_factor(context, kappa) &
    bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none
    integer (c_int64_t) , intent(in)  , value :: context
    real    (c_double)  , intent(in)  , value :: kappa
  end function qmckl_set_nucleus_rescale_factor
end interface

interface
  integer(c_int32_t) function qmckl_get_nucleus_nn_distance(context, distance, size_max) &
    bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none
    integer (c_int64_t) , intent(in)  , value :: context
    real    (c_double ) , intent(out)         :: distance(*)
    integer (c_int64_t) , intent(in)  , value :: size_max
  end function
end interface

interface
  integer(c_int32_t) function qmckl_get_nucleus_nn_distance_rescaled(context, distance_rescaled, size_max) &
    bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none
    integer (c_int64_t) , intent(in)  , value :: context
    real    (c_double ) , intent(out)         :: distance_rescaled(*)
    integer (c_int64_t) , intent(in)  , value :: size_max
  end function
end interface

interface
  integer(c_int32_t) function qmckl_get_nucleus_repulsion(context, energy) &
    bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none
    integer (c_int64_t) , intent(in)  , value :: context
    real    (c_double ) , intent(out)         :: energy
  end function
end interface
integer, parameter :: QMCKL_DEFAULT_PRECISION        = 53
integer, parameter :: QMCKL_DEFAULT_RANGE            = 11

interface
   integer (qmckl_exit_code) function qmckl_set_numprec_precision(context, precision) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer (qmckl_context), intent(in), value :: context
     integer (c_int32_t), intent(in), value :: precision
   end function qmckl_set_numprec_precision
end interface

interface
   integer (qmckl_exit_code) function qmckl_get_numprec_precision(context) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer (qmckl_context), intent(in), value :: context
   end function qmckl_get_numprec_precision
end interface

interface
   integer (qmckl_exit_code) function qmckl_set_numprec_range(context, range) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer (qmckl_context), intent(in), value :: context
     integer (c_int32_t), intent(in), value :: range
   end function qmckl_set_numprec_range
end interface

interface
   integer (qmckl_exit_code) function qmckl_get_numprec_range(context) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer (qmckl_context), intent(in), value :: context
   end function qmckl_get_numprec_range
end interface

interface
   real (c_double) function qmckl_get_numprec_epsilon(context) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer (qmckl_context), intent(in), value :: context
   end function qmckl_get_numprec_epsilon
end interface
interface
  integer(c_int32_t) function qmckl_get_point_num(context, num) bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none

    integer (c_int64_t) , intent(in)  , value :: context
    integer (c_int64_t) , intent(out)         :: num
  end function
end interface

interface
  integer(c_int32_t) function qmckl_get_point(context, transp, coord, size_max) bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none

    integer (c_int64_t) , intent(in)  , value :: context
    character(c_char)   , intent(in)  , value :: transp
    real    (c_double ) , intent(out)         :: coord(*)
    integer (c_int64_t) , intent(in)          :: size_max
  end function
end interface

interface
  integer(c_int32_t) function qmckl_set_point(context, &
       transp, num, coord, size_max) bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none

    integer (c_int64_t) , intent(in)  , value :: context
    character(c_char)   , intent(in)  , value :: transp
    integer (c_int64_t) , intent(in)  , value :: num
    real    (c_double ) , intent(in)          :: coord(*)
    integer (c_int64_t) , intent(in)  , value :: size_max
  end function
end interface
! Fortran interface                                               :noexport:
!    :PROPERTIES:
!    :Name:     qmckl_sherman_morrison
!    :CRetType: qmckl_exit_code
!    :FRetType: qmckl_exit_code
!    :END:

!    #+CALL: generate_f_interface(table=qmckl_sherman_morrison_args,rettyp=get_value("FRetType"),fname=get_value("Name"))

!    #+RESULTS:

interface
  integer(c_int32_t) function qmckl_sherman_morrison &
      (context, LDS, Dim, N_updates, Updates, Updates_index, breakdown, Slater_inv, determinant) &
      bind(C)

    use, intrinsic :: iso_c_binding
    import
    implicit none

    integer (c_int64_t) , intent(in) , value :: context
    integer (c_int64_t) , intent(in) , value :: LDS
    integer (c_int64_t) , intent(in) , value :: Dim
    integer (c_int64_t) , intent(in) , value :: N_updates
    real    (c_double ) , intent(in)         :: Updates(N_updates*Dim)
    integer (c_int64_t) , intent(in)         :: Updates_index(N_updates)
    real    (c_double ) , intent(in) , value :: breakdown
    real    (c_double ) , intent(inout)      :: Slater_inv(LDS*Dim)
    real    (c_double ) , intent(inout)      :: determinant

  end function qmckl_sherman_morrison
end interface

! Fortran interface                                              :noexport:
!    :PROPERTIES:
!    :Name:     qmckl_woodbury_2
!    :CRetType: qmckl_exit_code
!    :FRetType: qmckl_exit_code
!    :END:

!    #+CALL: generate_f_interface(table=qmckl_woodbury_2_args,rettyp=get_value("FRetType"),fname=get_value("Name"))

!    #+RESULTS:

interface
  integer(c_int32_t) function qmckl_woodbury_2 &
      (context, LDS, Dim, Updates, Updates_index, breakdown, Slater_inv, determinant) &
      bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none

    integer (c_int64_t) , intent(in) , value :: context
    integer (c_int64_t) , intent(in) , value :: LDS
    integer (c_int64_t) , intent(in) , value :: Dim
    real    (c_double ) , intent(in)         :: Updates(2*Dim)
    integer (c_int64_t) , intent(in)         :: Updates_index(2)
    real    (c_double ) , intent(in) , value :: breakdown
    real    (c_double ) , intent(inout)      :: Slater_inv(LDS*Dim)
    real    (c_double ) , intent(inout)      :: determinant

  end function qmckl_woodbury_2
end interface

! Fortran interface                                               :noexport:
!    :PROPERTIES:
!    :Name:     qmckl_woodbury_3
!    :CRetType: qmckl_exit_code
!    :FRetType: qmckl_exit_code
!    :END:

!    #+CALL: generate_f_interface(table=qmckl_woodbury_3_args,rettyp=get_value("FRetType"),fname=get_value("Name"))

!    #+RESULTS:

interface
  integer(c_int32_t) function qmckl_woodbury_3 &
      (context, LDS, Dim, Updates, Updates_index, breakdown, Slater_inv, determinant) &
      bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none

    integer (c_int64_t) , intent(in) , value :: context
    integer (c_int64_t) , intent(in) , value :: LDS
    integer (c_int64_t) , intent(in) , value :: Dim
    real    (c_double ) , intent(in)         :: Updates(3*Dim)
    integer (c_int64_t) , intent(in)         :: Updates_index(3)
    real    (c_double ) , intent(in) , value :: breakdown
    real    (c_double ) , intent(inout)      :: Slater_inv(LDS*Dim)
    real    (c_double ) , intent(inout)      :: determinant

  end function qmckl_woodbury_3
end interface

! Fortran interface                                              :noexport:
!    :PROPERTIES:
!    :Name:     qmckl_sherman_morrison_splitting
!    :CRetType: qmckl_exit_code
!    :FRetType: qmckl_exit_code
!    :END:

!    #+CALL: generate_f_interface(table=qmckl_sherman_morrison_splitting_args,rettyp=get_value("FRetType"),fname=get_value("Name"))

!    #+RESULTS:

interface
  integer(c_int32_t) function qmckl_sherman_morrison_splitting &
      (context, LDS, Dim, N_updates, Updates, Updates_index, breakdown, Slater_inv, determinant) &
      bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none

    integer (c_int64_t) , intent(in) , value :: context
    integer (c_int64_t) , intent(in) , value :: LDS
    integer (c_int64_t) , intent(in) , value :: Dim
    integer (c_int64_t) , intent(in) , value :: N_updates
    real    (c_double ) , intent(in)         :: Updates(N_updates*Dim)
    integer (c_int64_t) , intent(in)         :: Updates_index(N_updates)
    real    (c_double ) , intent(in) , value :: breakdown
    real    (c_double ) , intent(inout)      :: Slater_inv(LDS*Dim)
    real    (c_double ) , intent(inout)      :: determinant

  end function qmckl_sherman_morrison_splitting
end interface

! Fortran interface                                               :noexport:
!    :PROPERTIES:
!    :Name:     qmckl_sherman_morrison_smw32s
!    :CRetType: qmckl_exit_code
!    :FRetType: qmckl_exit_code
!    :END:

!    #+CALL: generate_f_interface(table=qmckl_sherman_morrison_smw32s_args,rettyp=get_value("FRetType"),fname=get_value("Name"))

!    #+RESULTS:

interface
  integer(c_int32_t) function qmckl_sherman_morrison_smw32s &
      (context, LDS, Dim, N_updates, Updates, Updates_index, breakdown, Slater_inv, determinant) &
      bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none

    integer (c_int64_t) , intent(in) , value :: context
    integer (c_int64_t) , intent(in) , value :: LDS
    integer (c_int64_t) , intent(in) , value :: Dim
    integer (c_int64_t) , intent(in) , value :: N_updates
    real    (c_double ) , intent(in)         :: Updates(N_updates*Dim)
    integer (c_int64_t) , intent(in)         :: Updates_index(N_updates)
    real    (c_double ) , intent(in) , value :: breakdown
    real    (c_double ) , intent(inout)      :: Slater_inv(LDS*Dim)
    real    (c_double ) , intent(inout)      :: determinant

  end function qmckl_sherman_morrison_smw32s
end interface
interface
  integer(c_int32_t) function qmckl_trexio_read &
      (context, file_name, size_max) &
      bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none

    integer  (c_int64_t) , intent(in)  , value :: context
    integer  (c_int64_t) , intent(in)  , value :: size_max
    character(c_char   ) , intent(in)          :: file_name(size_max)

  end function qmckl_trexio_read
end interface
end module qmckl
