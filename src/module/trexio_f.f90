module trexio

  use, intrinsic :: iso_c_binding
  implicit none

  integer, parameter :: trexio_exit_code  = c_int32_t
  integer, parameter :: trexio_back_end_t = c_int32_t
  integer, parameter :: trexio_t          = c_size_t

  character(kind=c_char), parameter :: TREXIO_DELIM = c_new_line

integer(trexio_exit_code), parameter :: TREXIO_FAILURE                 = -1
integer(trexio_exit_code), parameter :: TREXIO_SUCCESS                 = 0
integer(trexio_exit_code), parameter :: TREXIO_INVALID_ARG_1           = 1
integer(trexio_exit_code), parameter :: TREXIO_INVALID_ARG_2           = 2
integer(trexio_exit_code), parameter :: TREXIO_INVALID_ARG_3           = 3
integer(trexio_exit_code), parameter :: TREXIO_INVALID_ARG_4           = 4
integer(trexio_exit_code), parameter :: TREXIO_INVALID_ARG_5           = 5
integer(trexio_exit_code), parameter :: TREXIO_END                     = 6
integer(trexio_exit_code), parameter :: TREXIO_READONLY                = 7
integer(trexio_exit_code), parameter :: TREXIO_ERRNO                   = 8
integer(trexio_exit_code), parameter :: TREXIO_INVALID_ID              = 9
integer(trexio_exit_code), parameter :: TREXIO_ALLOCATION_FAILED       = 10
integer(trexio_exit_code), parameter :: TREXIO_HAS_NOT                 = 11
integer(trexio_exit_code), parameter :: TREXIO_INVALID_NUM             = 12
integer(trexio_exit_code), parameter :: TREXIO_ATTR_ALREADY_EXISTS     = 13
integer(trexio_exit_code), parameter :: TREXIO_DSET_ALREADY_EXISTS     = 14
integer(trexio_exit_code), parameter :: TREXIO_OPEN_ERROR              = 15
integer(trexio_exit_code), parameter :: TREXIO_LOCK_ERROR              = 16
integer(trexio_exit_code), parameter :: TREXIO_UNLOCK_ERROR            = 17
integer(trexio_exit_code), parameter :: TREXIO_FILE_ERROR              = 18
integer(trexio_exit_code), parameter :: TREXIO_GROUP_READ_ERROR        = 19
integer(trexio_exit_code), parameter :: TREXIO_GROUP_WRITE_ERROR       = 20
integer(trexio_exit_code), parameter :: TREXIO_ELEM_READ_ERROR         = 21
integer(trexio_exit_code), parameter :: TREXIO_ELEM_WRITE_ERROR        = 22
integer(trexio_exit_code), parameter :: TREXIO_UNSAFE_ARRAY_DIM        = 23
integer(trexio_exit_code), parameter :: TREXIO_ATTR_MISSING            = 24
integer(trexio_exit_code), parameter :: TREXIO_DSET_MISSING            = 25
integer(trexio_exit_code), parameter :: TREXIO_BACK_END_MISSING        = 26
integer(trexio_exit_code), parameter :: TREXIO_INVALID_ARG_6           = 27
integer(trexio_exit_code), parameter :: TREXIO_INVALID_ARG_7           = 28
integer(trexio_exit_code), parameter :: TREXIO_INVALID_ARG_8           = 29
integer(trexio_exit_code), parameter :: TREXIO_INVALID_STR_LEN         = 30
integer(trexio_exit_code), parameter :: TREXIO_INT_SIZE_OVERFLOW       = 31
integer(trexio_exit_code), parameter :: TREXIO_SAFE_MODE               = 32
integer(trexio_exit_code), parameter :: TREXIO_INVALID_ELECTRON_NUM    = 33
integer(trexio_exit_code), parameter :: TREXIO_INVALID_DETERMINANT_NUM = 34
integer(trexio_exit_code), parameter :: TREXIO_INVALID_STATE           = 35
integer(trexio_exit_code), parameter :: TREXIO_VERSION_PARSING_ISSUE   = 36

interface
   subroutine trexio_string_of_error (error, string) bind(C, name='trexio_string_of_error_f')
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_exit_code), intent(in), value :: error
     character(kind=c_char), intent(out)          :: string(128)
   end subroutine trexio_string_of_error
end interface

integer(trexio_back_end_t), parameter :: TREXIO_HDF5 = 0
  integer(trexio_back_end_t), parameter :: TREXIO_TEXT = 1
! integer(trexio_back_end_t), parameter :: TREXIO_JSON = 2
  integer(trexio_back_end_t), parameter :: TREXIO_INVALID_BACK_END = 2
  integer(trexio_back_end_t), parameter :: TREXIO_AUTO = TREXIO_INVALID_BACK_END

interface
   logical(c_bool) function trexio_has_back_end (back_end) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_back_end_t), intent(in), value :: back_end
   end function trexio_has_back_end
end interface

interface
   logical(c_bool) function trexio_has_backend (back_end) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_back_end_t), intent(in), value :: back_end
   end function trexio_has_backend
end interface

interface
   integer(trexio_t) function trexio_open_c (filename, mode, back_end, rc_open) bind(C, name="trexio_open")
     use, intrinsic :: iso_c_binding
     import
     character(kind=c_char), dimension(*)          :: filename
     character(kind=c_char), intent(in), value     :: mode
     integer(trexio_back_end_t), intent(in), value :: back_end
     integer(trexio_exit_code), intent(out)        :: rc_open
   end function trexio_open_c
end interface

interface
   integer(trexio_exit_code) function trexio_set_one_based(trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
   end function trexio_set_one_based
end interface

interface
   integer(trexio_exit_code) function trexio_close (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
   end function trexio_close
end interface

interface
   integer(trexio_exit_code) function trexio_flush (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
   end function trexio_flush
end interface

interface
   integer(trexio_exit_code) function trexio_inquire_c (filename) bind(C, name="trexio_inquire")
     use, intrinsic :: iso_c_binding
     import
     character(kind=c_char), dimension(*)       :: filename
   end function trexio_inquire_c
end interface

interface
   integer(trexio_exit_code) function trexio_set_state (trex_file, state) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int32_t), intent(in), value :: state
   end function trexio_set_state
end interface

interface
   integer(trexio_exit_code) function trexio_get_state (trex_file, state) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int32_t), intent(out)       :: state
   end function trexio_get_state
end interface

interface
   integer function trexio_info () bind(C)
     use, intrinsic :: iso_c_binding
   end function trexio_info
end interface

interface
   integer(trexio_exit_code) function trexio_to_bitfield_list_c(list, occupied_num, det_list, N_int) &
   bind(C, name="trexio_to_bitfield_list")
     use, intrinsic :: iso_c_binding
     import
     integer(c_int32_t), intent(in)        :: list(*)
     integer(c_int32_t), intent(in), value :: occupied_num
     integer(c_int64_t), intent(inout)     :: det_list(*)
     integer(c_int32_t), intent(in), value :: N_int
   end function trexio_to_bitfield_list_c
end interface

interface
   integer(trexio_exit_code) function trexio_to_orbital_list_c(N_int, d1, list, occupied_num) &
   bind(C, name="trexio_to_orbital_list")
     use, intrinsic :: iso_c_binding
     import
     integer(c_int32_t), intent(in), value :: N_int
     integer(c_int64_t), intent(in)        :: d1(*)
     integer(c_int32_t), intent(inout)     :: list(*)
     integer(c_int32_t), intent(inout)     :: occupied_num
   end function trexio_to_orbital_list_c
end interface

interface
   integer(trexio_exit_code) function trexio_to_orbital_list_up_dn_c(N_int, d1, list_up, list_dn, occ_num_up, occ_num_dn) &
   bind(C, name="trexio_to_orbital_list_up_dn")
     use, intrinsic :: iso_c_binding
     import
     integer(c_int32_t), intent(in), value :: N_int
     integer(c_int64_t), intent(in)        :: d1(*)
     integer(c_int32_t), intent(inout)     :: list_up(*)
     integer(c_int32_t), intent(inout)     :: list_dn(*)
     integer(c_int32_t), intent(inout)     :: occ_num_up
     integer(c_int32_t), intent(inout)     :: occ_num_dn
   end function trexio_to_orbital_list_up_dn_c
end interface

character(len = 12) :: TREXIO_PACKAGE_VERSION = "2.3.0"
integer :: TREXIO_VERSION_MAJOR = 2
integer :: TREXIO_VERSION_MINOR = 3
integer :: TREXIO_VERSION_PATCH = 0
character(len = 64) :: TREXIO_GIT_HASH = "813577e4a6d0a617ebafadea3306f10a9f1ad237"

interface
   integer(trexio_exit_code) function trexio_delete_metadata (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value  :: trex_file
   end function trexio_delete_metadata
end interface

interface
   integer(trexio_exit_code) function trexio_delete_electron (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value  :: trex_file
   end function trexio_delete_electron
end interface

interface
   integer(trexio_exit_code) function trexio_delete_nucleus (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value  :: trex_file
   end function trexio_delete_nucleus
end interface

interface
   integer(trexio_exit_code) function trexio_delete_ecp (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value  :: trex_file
   end function trexio_delete_ecp
end interface

interface
   integer(trexio_exit_code) function trexio_delete_basis (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value  :: trex_file
   end function trexio_delete_basis
end interface

interface
   integer(trexio_exit_code) function trexio_delete_ao (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value  :: trex_file
   end function trexio_delete_ao
end interface

interface
   integer(trexio_exit_code) function trexio_delete_ao_1e_int (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value  :: trex_file
   end function trexio_delete_ao_1e_int
end interface

interface
   integer(trexio_exit_code) function trexio_delete_ao_2e_int (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value  :: trex_file
   end function trexio_delete_ao_2e_int
end interface

interface
   integer(trexio_exit_code) function trexio_delete_mo (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value  :: trex_file
   end function trexio_delete_mo
end interface

interface
   integer(trexio_exit_code) function trexio_delete_mo_1e_int (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value  :: trex_file
   end function trexio_delete_mo_1e_int
end interface

interface
   integer(trexio_exit_code) function trexio_delete_mo_2e_int (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value  :: trex_file
   end function trexio_delete_mo_2e_int
end interface

interface
   integer(trexio_exit_code) function trexio_delete_determinant (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value  :: trex_file
   end function trexio_delete_determinant
end interface

interface
   integer(trexio_exit_code) function trexio_delete_csf (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value  :: trex_file
   end function trexio_delete_csf
end interface

interface
   integer(trexio_exit_code) function trexio_delete_state (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value  :: trex_file
   end function trexio_delete_state
end interface

interface
   integer(trexio_exit_code) function trexio_delete_rdm (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value  :: trex_file
   end function trexio_delete_rdm
end interface

interface
   integer(trexio_exit_code) function trexio_delete_cell (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value  :: trex_file
   end function trexio_delete_cell
end interface

interface
   integer(trexio_exit_code) function trexio_delete_pbc (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value  :: trex_file
   end function trexio_delete_pbc
end interface

interface
   integer(trexio_exit_code) function trexio_delete_qmc (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value  :: trex_file
   end function trexio_delete_qmc
end interface

interface
   integer(trexio_exit_code) function trexio_has_metadata_code_num (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
   end function trexio_has_metadata_code_num
end interface

interface
   integer(trexio_exit_code) function trexio_has_metadata_author_num (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
   end function trexio_has_metadata_author_num
end interface

interface
   integer(trexio_exit_code) function trexio_has_metadata_unsafe (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
   end function trexio_has_metadata_unsafe
end interface

interface
   integer(trexio_exit_code) function trexio_has_electron_num (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
   end function trexio_has_electron_num
end interface

interface
   integer(trexio_exit_code) function trexio_has_electron_up_num (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
   end function trexio_has_electron_up_num
end interface

interface
   integer(trexio_exit_code) function trexio_has_electron_dn_num (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
   end function trexio_has_electron_dn_num
end interface

interface
   integer(trexio_exit_code) function trexio_has_nucleus_num (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
   end function trexio_has_nucleus_num
end interface

interface
   integer(trexio_exit_code) function trexio_has_nucleus_repulsion (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
   end function trexio_has_nucleus_repulsion
end interface

interface
   integer(trexio_exit_code) function trexio_has_ecp_num (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
   end function trexio_has_ecp_num
end interface

interface
   integer(trexio_exit_code) function trexio_has_basis_prim_num (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
   end function trexio_has_basis_prim_num
end interface

interface
   integer(trexio_exit_code) function trexio_has_basis_shell_num (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
   end function trexio_has_basis_shell_num
end interface

interface
   integer(trexio_exit_code) function trexio_has_ao_cartesian (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
   end function trexio_has_ao_cartesian
end interface

interface
   integer(trexio_exit_code) function trexio_has_ao_num (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
   end function trexio_has_ao_num
end interface

interface
   integer(trexio_exit_code) function trexio_has_mo_num (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
   end function trexio_has_mo_num
end interface

interface
   integer(trexio_exit_code) function trexio_has_determinant_num (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
   end function trexio_has_determinant_num
end interface

interface
   integer(trexio_exit_code) function trexio_has_csf_num (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
   end function trexio_has_csf_num
end interface

interface
   integer(trexio_exit_code) function trexio_has_state_num (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
   end function trexio_has_state_num
end interface

interface
   integer(trexio_exit_code) function trexio_has_state_id (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
   end function trexio_has_state_id
end interface

interface
   integer(trexio_exit_code) function trexio_has_pbc_periodic (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
   end function trexio_has_pbc_periodic
end interface

interface
   integer(trexio_exit_code) function trexio_has_qmc_num (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
   end function trexio_has_qmc_num
end interface

interface
   integer(trexio_exit_code) function trexio_has_metadata_package_version (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value  :: trex_file
   end function trexio_has_metadata_package_version
end interface

interface
   integer(trexio_exit_code) function trexio_has_metadata_description (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value  :: trex_file
   end function trexio_has_metadata_description
end interface

interface
   integer(trexio_exit_code) function trexio_has_nucleus_point_group (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value  :: trex_file
   end function trexio_has_nucleus_point_group
end interface

interface
   integer(trexio_exit_code) function trexio_has_basis_type (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value  :: trex_file
   end function trexio_has_basis_type
end interface

interface
   integer(trexio_exit_code) function trexio_has_mo_type (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value  :: trex_file
   end function trexio_has_mo_type
end interface

interface
   integer(trexio_exit_code) function trexio_has_state_current_label (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value  :: trex_file
   end function trexio_has_state_current_label
end interface

interface
   integer(trexio_exit_code) function trexio_has_determinant_coefficient (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
   end function trexio_has_determinant_coefficient
end interface

interface
   integer(trexio_exit_code) function trexio_has_csf_coefficient (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
   end function trexio_has_csf_coefficient
end interface

interface
   integer(trexio_exit_code) function trexio_has_nucleus_charge (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
   end function trexio_has_nucleus_charge
end interface

interface
   integer(trexio_exit_code) function trexio_has_nucleus_coord (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
   end function trexio_has_nucleus_coord
end interface

interface
   integer(trexio_exit_code) function trexio_has_ecp_max_ang_mom_plus_1 (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
   end function trexio_has_ecp_max_ang_mom_plus_1
end interface

interface
   integer(trexio_exit_code) function trexio_has_ecp_z_core (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
   end function trexio_has_ecp_z_core
end interface

interface
   integer(trexio_exit_code) function trexio_has_ecp_ang_mom (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
   end function trexio_has_ecp_ang_mom
end interface

interface
   integer(trexio_exit_code) function trexio_has_ecp_nucleus_index (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
   end function trexio_has_ecp_nucleus_index
end interface

interface
   integer(trexio_exit_code) function trexio_has_ecp_exponent (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
   end function trexio_has_ecp_exponent
end interface

interface
   integer(trexio_exit_code) function trexio_has_ecp_coefficient (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
   end function trexio_has_ecp_coefficient
end interface

interface
   integer(trexio_exit_code) function trexio_has_ecp_power (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
   end function trexio_has_ecp_power
end interface

interface
   integer(trexio_exit_code) function trexio_has_basis_nucleus_index (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
   end function trexio_has_basis_nucleus_index
end interface

interface
   integer(trexio_exit_code) function trexio_has_basis_shell_ang_mom (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
   end function trexio_has_basis_shell_ang_mom
end interface

interface
   integer(trexio_exit_code) function trexio_has_basis_shell_factor (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
   end function trexio_has_basis_shell_factor
end interface

interface
   integer(trexio_exit_code) function trexio_has_basis_shell_index (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
   end function trexio_has_basis_shell_index
end interface

interface
   integer(trexio_exit_code) function trexio_has_basis_exponent (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
   end function trexio_has_basis_exponent
end interface

interface
   integer(trexio_exit_code) function trexio_has_basis_coefficient (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
   end function trexio_has_basis_coefficient
end interface

interface
   integer(trexio_exit_code) function trexio_has_basis_prim_factor (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
   end function trexio_has_basis_prim_factor
end interface

interface
   integer(trexio_exit_code) function trexio_has_ao_shell (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
   end function trexio_has_ao_shell
end interface

interface
   integer(trexio_exit_code) function trexio_has_ao_normalization (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
   end function trexio_has_ao_normalization
end interface

interface
   integer(trexio_exit_code) function trexio_has_ao_1e_int_overlap (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
   end function trexio_has_ao_1e_int_overlap
end interface

interface
   integer(trexio_exit_code) function trexio_has_ao_1e_int_kinetic (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
   end function trexio_has_ao_1e_int_kinetic
end interface

interface
   integer(trexio_exit_code) function trexio_has_ao_1e_int_potential_n_e (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
   end function trexio_has_ao_1e_int_potential_n_e
end interface

interface
   integer(trexio_exit_code) function trexio_has_ao_1e_int_ecp (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
   end function trexio_has_ao_1e_int_ecp
end interface

interface
   integer(trexio_exit_code) function trexio_has_ao_1e_int_core_hamiltonian (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
   end function trexio_has_ao_1e_int_core_hamiltonian
end interface

interface
   integer(trexio_exit_code) function trexio_has_ao_1e_int_overlap_im (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
   end function trexio_has_ao_1e_int_overlap_im
end interface

interface
   integer(trexio_exit_code) function trexio_has_ao_1e_int_kinetic_im (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
   end function trexio_has_ao_1e_int_kinetic_im
end interface

interface
   integer(trexio_exit_code) function trexio_has_ao_1e_int_potential_n_e_im (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
   end function trexio_has_ao_1e_int_potential_n_e_im
end interface

interface
   integer(trexio_exit_code) function trexio_has_ao_1e_int_ecp_im (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
   end function trexio_has_ao_1e_int_ecp_im
end interface

interface
   integer(trexio_exit_code) function trexio_has_ao_1e_int_core_hamiltonian_im (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
   end function trexio_has_ao_1e_int_core_hamiltonian_im
end interface

interface
   integer(trexio_exit_code) function trexio_has_mo_coefficient (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
   end function trexio_has_mo_coefficient
end interface

interface
   integer(trexio_exit_code) function trexio_has_mo_coefficient_im (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
   end function trexio_has_mo_coefficient_im
end interface

interface
   integer(trexio_exit_code) function trexio_has_mo_occupation (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
   end function trexio_has_mo_occupation
end interface

interface
   integer(trexio_exit_code) function trexio_has_mo_energy (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
   end function trexio_has_mo_energy
end interface

interface
   integer(trexio_exit_code) function trexio_has_mo_spin (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
   end function trexio_has_mo_spin
end interface

interface
   integer(trexio_exit_code) function trexio_has_mo_1e_int_overlap (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
   end function trexio_has_mo_1e_int_overlap
end interface

interface
   integer(trexio_exit_code) function trexio_has_mo_1e_int_kinetic (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
   end function trexio_has_mo_1e_int_kinetic
end interface

interface
   integer(trexio_exit_code) function trexio_has_mo_1e_int_potential_n_e (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
   end function trexio_has_mo_1e_int_potential_n_e
end interface

interface
   integer(trexio_exit_code) function trexio_has_mo_1e_int_ecp (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
   end function trexio_has_mo_1e_int_ecp
end interface

interface
   integer(trexio_exit_code) function trexio_has_mo_1e_int_core_hamiltonian (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
   end function trexio_has_mo_1e_int_core_hamiltonian
end interface

interface
   integer(trexio_exit_code) function trexio_has_mo_1e_int_overlap_im (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
   end function trexio_has_mo_1e_int_overlap_im
end interface

interface
   integer(trexio_exit_code) function trexio_has_mo_1e_int_kinetic_im (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
   end function trexio_has_mo_1e_int_kinetic_im
end interface

interface
   integer(trexio_exit_code) function trexio_has_mo_1e_int_potential_n_e_im (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
   end function trexio_has_mo_1e_int_potential_n_e_im
end interface

interface
   integer(trexio_exit_code) function trexio_has_mo_1e_int_ecp_im (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
   end function trexio_has_mo_1e_int_ecp_im
end interface

interface
   integer(trexio_exit_code) function trexio_has_mo_1e_int_core_hamiltonian_im (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
   end function trexio_has_mo_1e_int_core_hamiltonian_im
end interface

interface
   integer(trexio_exit_code) function trexio_has_rdm_1e (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
   end function trexio_has_rdm_1e
end interface

interface
   integer(trexio_exit_code) function trexio_has_rdm_1e_up (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
   end function trexio_has_rdm_1e_up
end interface

interface
   integer(trexio_exit_code) function trexio_has_rdm_1e_dn (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
   end function trexio_has_rdm_1e_dn
end interface

interface
   integer(trexio_exit_code) function trexio_has_cell_a (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
   end function trexio_has_cell_a
end interface

interface
   integer(trexio_exit_code) function trexio_has_cell_b (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
   end function trexio_has_cell_b
end interface

interface
   integer(trexio_exit_code) function trexio_has_cell_c (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
   end function trexio_has_cell_c
end interface

interface
   integer(trexio_exit_code) function trexio_has_pbc_k_point (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
   end function trexio_has_pbc_k_point
end interface

interface
   integer(trexio_exit_code) function trexio_has_qmc_point (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
   end function trexio_has_qmc_point
end interface

interface
   integer(trexio_exit_code) function trexio_has_qmc_psi (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
   end function trexio_has_qmc_psi
end interface

interface
   integer(trexio_exit_code) function trexio_has_qmc_e_loc (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
   end function trexio_has_qmc_e_loc
end interface

interface
   integer(trexio_exit_code) function trexio_has_ao_2e_int_eri (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
   end function trexio_has_ao_2e_int_eri
end interface

interface
   integer(trexio_exit_code) function trexio_has_ao_2e_int_eri_lr (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
   end function trexio_has_ao_2e_int_eri_lr
end interface

interface
   integer(trexio_exit_code) function trexio_has_mo_2e_int_eri (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
   end function trexio_has_mo_2e_int_eri
end interface

interface
   integer(trexio_exit_code) function trexio_has_mo_2e_int_eri_lr (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
   end function trexio_has_mo_2e_int_eri_lr
end interface

interface
   integer(trexio_exit_code) function trexio_has_csf_det_coefficient (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
   end function trexio_has_csf_det_coefficient
end interface

interface
   integer(trexio_exit_code) function trexio_has_rdm_2e (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
   end function trexio_has_rdm_2e
end interface

interface
   integer(trexio_exit_code) function trexio_has_rdm_2e_upup (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
   end function trexio_has_rdm_2e_upup
end interface

interface
   integer(trexio_exit_code) function trexio_has_rdm_2e_dndn (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
   end function trexio_has_rdm_2e_dndn
end interface

interface
   integer(trexio_exit_code) function trexio_has_rdm_2e_updn (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
   end function trexio_has_rdm_2e_updn
end interface

interface
   integer(trexio_exit_code) function trexio_has_rdm_2e_dnup (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
   end function trexio_has_rdm_2e_dnup
end interface

interface
   integer(trexio_exit_code) function trexio_has_metadata_code (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
   end function trexio_has_metadata_code
end interface

interface
   integer(trexio_exit_code) function trexio_has_metadata_author (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
   end function trexio_has_metadata_author
end interface

interface
   integer(trexio_exit_code) function trexio_has_nucleus_label (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
   end function trexio_has_nucleus_label
end interface

interface
   integer(trexio_exit_code) function trexio_has_mo_class (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
   end function trexio_has_mo_class
end interface

interface
   integer(trexio_exit_code) function trexio_has_mo_symmetry (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
   end function trexio_has_mo_symmetry
end interface

interface
   integer(trexio_exit_code) function trexio_has_state_label (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
   end function trexio_has_state_label
end interface

interface
   integer(trexio_exit_code) function trexio_has_state_file_name (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
   end function trexio_has_state_file_name
end interface

interface
   integer(trexio_exit_code) function trexio_has_metadata (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
   end function trexio_has_metadata
end interface

interface
   integer(trexio_exit_code) function trexio_has_electron (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
   end function trexio_has_electron
end interface

interface
   integer(trexio_exit_code) function trexio_has_nucleus (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
   end function trexio_has_nucleus
end interface

interface
   integer(trexio_exit_code) function trexio_has_ecp (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
   end function trexio_has_ecp
end interface

interface
   integer(trexio_exit_code) function trexio_has_basis (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
   end function trexio_has_basis
end interface

interface
   integer(trexio_exit_code) function trexio_has_ao (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
   end function trexio_has_ao
end interface

interface
   integer(trexio_exit_code) function trexio_has_ao_1e_int (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
   end function trexio_has_ao_1e_int
end interface

interface
   integer(trexio_exit_code) function trexio_has_ao_2e_int (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
   end function trexio_has_ao_2e_int
end interface

interface
   integer(trexio_exit_code) function trexio_has_mo (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
   end function trexio_has_mo
end interface

interface
   integer(trexio_exit_code) function trexio_has_mo_1e_int (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
   end function trexio_has_mo_1e_int
end interface

interface
   integer(trexio_exit_code) function trexio_has_mo_2e_int (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
   end function trexio_has_mo_2e_int
end interface

interface
   integer(trexio_exit_code) function trexio_has_determinant (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
   end function trexio_has_determinant
end interface

interface
   integer(trexio_exit_code) function trexio_has_csf (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
   end function trexio_has_csf
end interface

interface
   integer(trexio_exit_code) function trexio_has_state (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
   end function trexio_has_state
end interface

interface
   integer(trexio_exit_code) function trexio_has_rdm (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
   end function trexio_has_rdm
end interface

interface
   integer(trexio_exit_code) function trexio_has_cell (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
   end function trexio_has_cell
end interface

interface
   integer(trexio_exit_code) function trexio_has_pbc (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
   end function trexio_has_pbc
end interface

interface
   integer(trexio_exit_code) function trexio_has_qmc (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
   end function trexio_has_qmc
end interface

interface
   integer(trexio_exit_code) function trexio_read_metadata_code_num_32 (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int32_t), intent(out) :: num
   end function trexio_read_metadata_code_num_32
end interface

interface
   integer(trexio_exit_code) function trexio_read_metadata_author_num_32 (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int32_t), intent(out) :: num
   end function trexio_read_metadata_author_num_32
end interface

interface
   integer(trexio_exit_code) function trexio_read_metadata_unsafe_32 (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int32_t), intent(out) :: num
   end function trexio_read_metadata_unsafe_32
end interface

interface
   integer(trexio_exit_code) function trexio_read_electron_num_32 (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int32_t), intent(out) :: num
   end function trexio_read_electron_num_32
end interface

interface
   integer(trexio_exit_code) function trexio_read_electron_up_num_32 (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int32_t), intent(out) :: num
   end function trexio_read_electron_up_num_32
end interface

interface
   integer(trexio_exit_code) function trexio_read_electron_dn_num_32 (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int32_t), intent(out) :: num
   end function trexio_read_electron_dn_num_32
end interface

interface
   integer(trexio_exit_code) function trexio_read_nucleus_num_32 (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int32_t), intent(out) :: num
   end function trexio_read_nucleus_num_32
end interface

interface
   integer(trexio_exit_code) function trexio_read_nucleus_repulsion_32 (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_float), intent(out) :: num
   end function trexio_read_nucleus_repulsion_32
end interface

interface
   integer(trexio_exit_code) function trexio_read_ecp_num_32 (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int32_t), intent(out) :: num
   end function trexio_read_ecp_num_32
end interface

interface
   integer(trexio_exit_code) function trexio_read_basis_prim_num_32 (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int32_t), intent(out) :: num
   end function trexio_read_basis_prim_num_32
end interface

interface
   integer(trexio_exit_code) function trexio_read_basis_shell_num_32 (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int32_t), intent(out) :: num
   end function trexio_read_basis_shell_num_32
end interface

interface
   integer(trexio_exit_code) function trexio_read_ao_cartesian_32 (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int32_t), intent(out) :: num
   end function trexio_read_ao_cartesian_32
end interface

interface
   integer(trexio_exit_code) function trexio_read_ao_num_32 (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int32_t), intent(out) :: num
   end function trexio_read_ao_num_32
end interface

interface
   integer(trexio_exit_code) function trexio_read_mo_num_32 (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int32_t), intent(out) :: num
   end function trexio_read_mo_num_32
end interface

interface
   integer(trexio_exit_code) function trexio_read_determinant_num_32 (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int32_t), intent(out) :: num
   end function trexio_read_determinant_num_32
end interface

interface
   integer(trexio_exit_code) function trexio_read_csf_num_32 (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int32_t), intent(out) :: num
   end function trexio_read_csf_num_32
end interface

interface
   integer(trexio_exit_code) function trexio_read_state_num_32 (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int32_t), intent(out) :: num
   end function trexio_read_state_num_32
end interface

interface
   integer(trexio_exit_code) function trexio_read_state_id_32 (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int32_t), intent(out) :: num
   end function trexio_read_state_id_32
end interface

interface
   integer(trexio_exit_code) function trexio_read_pbc_periodic_32 (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int32_t), intent(out) :: num
   end function trexio_read_pbc_periodic_32
end interface

interface
   integer(trexio_exit_code) function trexio_read_qmc_num_32 (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int32_t), intent(out) :: num
   end function trexio_read_qmc_num_32
end interface

interface
   integer(trexio_exit_code) function trexio_read_metadata_code_num_64 (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int64_t), intent(out) :: num
   end function trexio_read_metadata_code_num_64
end interface

interface
   integer(trexio_exit_code) function trexio_read_metadata_author_num_64 (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int64_t), intent(out) :: num
   end function trexio_read_metadata_author_num_64
end interface

interface
   integer(trexio_exit_code) function trexio_read_metadata_unsafe_64 (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int64_t), intent(out) :: num
   end function trexio_read_metadata_unsafe_64
end interface

interface
   integer(trexio_exit_code) function trexio_read_electron_num_64 (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int64_t), intent(out) :: num
   end function trexio_read_electron_num_64
end interface

interface
   integer(trexio_exit_code) function trexio_read_electron_up_num_64 (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int64_t), intent(out) :: num
   end function trexio_read_electron_up_num_64
end interface

interface
   integer(trexio_exit_code) function trexio_read_electron_dn_num_64 (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int64_t), intent(out) :: num
   end function trexio_read_electron_dn_num_64
end interface

interface
   integer(trexio_exit_code) function trexio_read_nucleus_num_64 (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int64_t), intent(out) :: num
   end function trexio_read_nucleus_num_64
end interface

interface
   integer(trexio_exit_code) function trexio_read_nucleus_repulsion_64 (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(out) :: num
   end function trexio_read_nucleus_repulsion_64
end interface

interface
   integer(trexio_exit_code) function trexio_read_ecp_num_64 (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int64_t), intent(out) :: num
   end function trexio_read_ecp_num_64
end interface

interface
   integer(trexio_exit_code) function trexio_read_basis_prim_num_64 (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int64_t), intent(out) :: num
   end function trexio_read_basis_prim_num_64
end interface

interface
   integer(trexio_exit_code) function trexio_read_basis_shell_num_64 (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int64_t), intent(out) :: num
   end function trexio_read_basis_shell_num_64
end interface

interface
   integer(trexio_exit_code) function trexio_read_ao_cartesian_64 (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int64_t), intent(out) :: num
   end function trexio_read_ao_cartesian_64
end interface

interface
   integer(trexio_exit_code) function trexio_read_ao_num_64 (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int64_t), intent(out) :: num
   end function trexio_read_ao_num_64
end interface

interface
   integer(trexio_exit_code) function trexio_read_mo_num_64 (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int64_t), intent(out) :: num
   end function trexio_read_mo_num_64
end interface

interface
   integer(trexio_exit_code) function trexio_read_determinant_num_64 (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int64_t), intent(out) :: num
   end function trexio_read_determinant_num_64
end interface

interface
   integer(trexio_exit_code) function trexio_read_csf_num_64 (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int64_t), intent(out) :: num
   end function trexio_read_csf_num_64
end interface

interface
   integer(trexio_exit_code) function trexio_read_state_num_64 (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int64_t), intent(out) :: num
   end function trexio_read_state_num_64
end interface

interface
   integer(trexio_exit_code) function trexio_read_state_id_64 (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int64_t), intent(out) :: num
   end function trexio_read_state_id_64
end interface

interface
   integer(trexio_exit_code) function trexio_read_pbc_periodic_64 (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int64_t), intent(out) :: num
   end function trexio_read_pbc_periodic_64
end interface

interface
   integer(trexio_exit_code) function trexio_read_qmc_num_64 (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int64_t), intent(out) :: num
   end function trexio_read_qmc_num_64
end interface

interface
   integer(trexio_exit_code) function trexio_read_metadata_code_num (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int32_t), intent(out) :: num
   end function trexio_read_metadata_code_num
end interface

interface
   integer(trexio_exit_code) function trexio_read_metadata_author_num (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int32_t), intent(out) :: num
   end function trexio_read_metadata_author_num
end interface

interface
   integer(trexio_exit_code) function trexio_read_metadata_unsafe (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int32_t), intent(out) :: num
   end function trexio_read_metadata_unsafe
end interface

interface
   integer(trexio_exit_code) function trexio_read_electron_num (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int32_t), intent(out) :: num
   end function trexio_read_electron_num
end interface

interface
   integer(trexio_exit_code) function trexio_read_electron_up_num (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int32_t), intent(out) :: num
   end function trexio_read_electron_up_num
end interface

interface
   integer(trexio_exit_code) function trexio_read_electron_dn_num (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int32_t), intent(out) :: num
   end function trexio_read_electron_dn_num
end interface

interface
   integer(trexio_exit_code) function trexio_read_nucleus_num (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int32_t), intent(out) :: num
   end function trexio_read_nucleus_num
end interface

interface
   integer(trexio_exit_code) function trexio_read_nucleus_repulsion (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(out) :: num
   end function trexio_read_nucleus_repulsion
end interface

interface
   integer(trexio_exit_code) function trexio_read_ecp_num (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int32_t), intent(out) :: num
   end function trexio_read_ecp_num
end interface

interface
   integer(trexio_exit_code) function trexio_read_basis_prim_num (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int32_t), intent(out) :: num
   end function trexio_read_basis_prim_num
end interface

interface
   integer(trexio_exit_code) function trexio_read_basis_shell_num (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int32_t), intent(out) :: num
   end function trexio_read_basis_shell_num
end interface

interface
   integer(trexio_exit_code) function trexio_read_ao_cartesian (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int32_t), intent(out) :: num
   end function trexio_read_ao_cartesian
end interface

interface
   integer(trexio_exit_code) function trexio_read_ao_num (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int32_t), intent(out) :: num
   end function trexio_read_ao_num
end interface

interface
   integer(trexio_exit_code) function trexio_read_mo_num (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int32_t), intent(out) :: num
   end function trexio_read_mo_num
end interface

interface
   integer(trexio_exit_code) function trexio_read_determinant_num (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int32_t), intent(out) :: num
   end function trexio_read_determinant_num
end interface

interface
   integer(trexio_exit_code) function trexio_read_csf_num (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int32_t), intent(out) :: num
   end function trexio_read_csf_num
end interface

interface
   integer(trexio_exit_code) function trexio_read_state_num (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int32_t), intent(out) :: num
   end function trexio_read_state_num
end interface

interface
   integer(trexio_exit_code) function trexio_read_state_id (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int32_t), intent(out) :: num
   end function trexio_read_state_id
end interface

interface
   integer(trexio_exit_code) function trexio_read_pbc_periodic (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int32_t), intent(out) :: num
   end function trexio_read_pbc_periodic
end interface

interface
   integer(trexio_exit_code) function trexio_read_qmc_num (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int32_t), intent(out) :: num
   end function trexio_read_qmc_num
end interface

interface
   integer(trexio_exit_code) function trexio_read_metadata_package_version_c (trex_file, str, max_str_len) &
           bind(C, name="trexio_read_metadata_package_version")
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value  :: trex_file
     character(kind=c_char), intent(out)   :: str(*)
     integer(c_int32_t), intent(in), value :: max_str_len
   end function trexio_read_metadata_package_version_c
end interface

interface
   integer(trexio_exit_code) function trexio_read_metadata_description_c (trex_file, str, max_str_len) &
           bind(C, name="trexio_read_metadata_description")
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value  :: trex_file
     character(kind=c_char), intent(out)   :: str(*)
     integer(c_int32_t), intent(in), value :: max_str_len
   end function trexio_read_metadata_description_c
end interface

interface
   integer(trexio_exit_code) function trexio_read_nucleus_point_group_c (trex_file, str, max_str_len) &
           bind(C, name="trexio_read_nucleus_point_group")
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value  :: trex_file
     character(kind=c_char), intent(out)   :: str(*)
     integer(c_int32_t), intent(in), value :: max_str_len
   end function trexio_read_nucleus_point_group_c
end interface

interface
   integer(trexio_exit_code) function trexio_read_basis_type_c (trex_file, str, max_str_len) &
           bind(C, name="trexio_read_basis_type")
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value  :: trex_file
     character(kind=c_char), intent(out)   :: str(*)
     integer(c_int32_t), intent(in), value :: max_str_len
   end function trexio_read_basis_type_c
end interface

interface
   integer(trexio_exit_code) function trexio_read_mo_type_c (trex_file, str, max_str_len) &
           bind(C, name="trexio_read_mo_type")
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value  :: trex_file
     character(kind=c_char), intent(out)   :: str(*)
     integer(c_int32_t), intent(in), value :: max_str_len
   end function trexio_read_mo_type_c
end interface

interface
   integer(trexio_exit_code) function trexio_read_state_current_label_c (trex_file, str, max_str_len) &
           bind(C, name="trexio_read_state_current_label")
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value  :: trex_file
     character(kind=c_char), intent(out)   :: str(*)
     integer(c_int32_t), intent(in), value :: max_str_len
   end function trexio_read_state_current_label_c
end interface

interface
   integer(trexio_exit_code) function trexio_read_safe_determinant_coefficient (trex_file, &
                                                   offset_file, buffer_size, &
                                                   dset, dset_size) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int64_t), intent(in), value :: offset_file
     integer(c_int64_t), intent(inout)     :: buffer_size
     real(c_double), intent(out)           :: dset(*)
     integer(c_int64_t), intent(in), value :: dset_size
   end function trexio_read_safe_determinant_coefficient
end interface

interface
   integer(trexio_exit_code) function trexio_read_determinant_coefficient(trex_file, &
                                              offset_file, buffer_size, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int64_t), intent(in), value :: offset_file
     integer(c_int64_t), intent(inout)     :: buffer_size
     real(c_double), intent(out)           :: dset(*)
   end function trexio_read_determinant_coefficient
end interface

interface
   integer(trexio_exit_code) function trexio_read_determinant_coefficient_size (trex_file, &
                                                   size_max) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int64_t), intent(out) :: size_max
   end function trexio_read_determinant_coefficient_size
end interface

interface
   integer(trexio_exit_code) function trexio_read_safe_csf_coefficient (trex_file, &
                                                   offset_file, buffer_size, &
                                                   dset, dset_size) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int64_t), intent(in), value :: offset_file
     integer(c_int64_t), intent(inout)     :: buffer_size
     real(c_double), intent(out)           :: dset(*)
     integer(c_int64_t), intent(in), value :: dset_size
   end function trexio_read_safe_csf_coefficient
end interface

interface
   integer(trexio_exit_code) function trexio_read_csf_coefficient(trex_file, &
                                              offset_file, buffer_size, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int64_t), intent(in), value :: offset_file
     integer(c_int64_t), intent(inout)     :: buffer_size
     real(c_double), intent(out)           :: dset(*)
   end function trexio_read_csf_coefficient
end interface

interface
   integer(trexio_exit_code) function trexio_read_csf_coefficient_size (trex_file, &
                                                   size_max) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int64_t), intent(out) :: size_max
   end function trexio_read_csf_coefficient_size
end interface

interface
   integer(trexio_exit_code) function trexio_read_nucleus_charge_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_float), intent(out) :: dset(*)
   end function trexio_read_nucleus_charge_32
end interface

interface
   integer(trexio_exit_code) function trexio_read_nucleus_coord_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_float), intent(out) :: dset(*)
   end function trexio_read_nucleus_coord_32
end interface

interface
   integer(trexio_exit_code) function trexio_read_ecp_max_ang_mom_plus_1_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int32_t), intent(out) :: dset(*)
   end function trexio_read_ecp_max_ang_mom_plus_1_32
end interface

interface
   integer(trexio_exit_code) function trexio_read_ecp_z_core_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int32_t), intent(out) :: dset(*)
   end function trexio_read_ecp_z_core_32
end interface

interface
   integer(trexio_exit_code) function trexio_read_ecp_ang_mom_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int32_t), intent(out) :: dset(*)
   end function trexio_read_ecp_ang_mom_32
end interface

interface
   integer(trexio_exit_code) function trexio_read_ecp_nucleus_index_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int32_t), intent(out) :: dset(*)
   end function trexio_read_ecp_nucleus_index_32
end interface

interface
   integer(trexio_exit_code) function trexio_read_ecp_exponent_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_float), intent(out) :: dset(*)
   end function trexio_read_ecp_exponent_32
end interface

interface
   integer(trexio_exit_code) function trexio_read_ecp_coefficient_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_float), intent(out) :: dset(*)
   end function trexio_read_ecp_coefficient_32
end interface

interface
   integer(trexio_exit_code) function trexio_read_ecp_power_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int32_t), intent(out) :: dset(*)
   end function trexio_read_ecp_power_32
end interface

interface
   integer(trexio_exit_code) function trexio_read_basis_nucleus_index_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int32_t), intent(out) :: dset(*)
   end function trexio_read_basis_nucleus_index_32
end interface

interface
   integer(trexio_exit_code) function trexio_read_basis_shell_ang_mom_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int32_t), intent(out) :: dset(*)
   end function trexio_read_basis_shell_ang_mom_32
end interface

interface
   integer(trexio_exit_code) function trexio_read_basis_shell_factor_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_float), intent(out) :: dset(*)
   end function trexio_read_basis_shell_factor_32
end interface

interface
   integer(trexio_exit_code) function trexio_read_basis_shell_index_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int32_t), intent(out) :: dset(*)
   end function trexio_read_basis_shell_index_32
end interface

interface
   integer(trexio_exit_code) function trexio_read_basis_exponent_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_float), intent(out) :: dset(*)
   end function trexio_read_basis_exponent_32
end interface

interface
   integer(trexio_exit_code) function trexio_read_basis_coefficient_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_float), intent(out) :: dset(*)
   end function trexio_read_basis_coefficient_32
end interface

interface
   integer(trexio_exit_code) function trexio_read_basis_prim_factor_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_float), intent(out) :: dset(*)
   end function trexio_read_basis_prim_factor_32
end interface

interface
   integer(trexio_exit_code) function trexio_read_ao_shell_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int32_t), intent(out) :: dset(*)
   end function trexio_read_ao_shell_32
end interface

interface
   integer(trexio_exit_code) function trexio_read_ao_normalization_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_float), intent(out) :: dset(*)
   end function trexio_read_ao_normalization_32
end interface

interface
   integer(trexio_exit_code) function trexio_read_ao_1e_int_overlap_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_float), intent(out) :: dset(*)
   end function trexio_read_ao_1e_int_overlap_32
end interface

interface
   integer(trexio_exit_code) function trexio_read_ao_1e_int_kinetic_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_float), intent(out) :: dset(*)
   end function trexio_read_ao_1e_int_kinetic_32
end interface

interface
   integer(trexio_exit_code) function trexio_read_ao_1e_int_potential_n_e_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_float), intent(out) :: dset(*)
   end function trexio_read_ao_1e_int_potential_n_e_32
end interface

interface
   integer(trexio_exit_code) function trexio_read_ao_1e_int_ecp_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_float), intent(out) :: dset(*)
   end function trexio_read_ao_1e_int_ecp_32
end interface

interface
   integer(trexio_exit_code) function trexio_read_ao_1e_int_core_hamiltonian_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_float), intent(out) :: dset(*)
   end function trexio_read_ao_1e_int_core_hamiltonian_32
end interface

interface
   integer(trexio_exit_code) function trexio_read_ao_1e_int_overlap_im_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_float), intent(out) :: dset(*)
   end function trexio_read_ao_1e_int_overlap_im_32
end interface

interface
   integer(trexio_exit_code) function trexio_read_ao_1e_int_kinetic_im_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_float), intent(out) :: dset(*)
   end function trexio_read_ao_1e_int_kinetic_im_32
end interface

interface
   integer(trexio_exit_code) function trexio_read_ao_1e_int_potential_n_e_im_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_float), intent(out) :: dset(*)
   end function trexio_read_ao_1e_int_potential_n_e_im_32
end interface

interface
   integer(trexio_exit_code) function trexio_read_ao_1e_int_ecp_im_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_float), intent(out) :: dset(*)
   end function trexio_read_ao_1e_int_ecp_im_32
end interface

interface
   integer(trexio_exit_code) function trexio_read_ao_1e_int_core_hamiltonian_im_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_float), intent(out) :: dset(*)
   end function trexio_read_ao_1e_int_core_hamiltonian_im_32
end interface

interface
   integer(trexio_exit_code) function trexio_read_mo_coefficient_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_float), intent(out) :: dset(*)
   end function trexio_read_mo_coefficient_32
end interface

interface
   integer(trexio_exit_code) function trexio_read_mo_coefficient_im_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_float), intent(out) :: dset(*)
   end function trexio_read_mo_coefficient_im_32
end interface

interface
   integer(trexio_exit_code) function trexio_read_mo_occupation_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_float), intent(out) :: dset(*)
   end function trexio_read_mo_occupation_32
end interface

interface
   integer(trexio_exit_code) function trexio_read_mo_energy_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_float), intent(out) :: dset(*)
   end function trexio_read_mo_energy_32
end interface

interface
   integer(trexio_exit_code) function trexio_read_mo_spin_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int32_t), intent(out) :: dset(*)
   end function trexio_read_mo_spin_32
end interface

interface
   integer(trexio_exit_code) function trexio_read_mo_1e_int_overlap_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_float), intent(out) :: dset(*)
   end function trexio_read_mo_1e_int_overlap_32
end interface

interface
   integer(trexio_exit_code) function trexio_read_mo_1e_int_kinetic_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_float), intent(out) :: dset(*)
   end function trexio_read_mo_1e_int_kinetic_32
end interface

interface
   integer(trexio_exit_code) function trexio_read_mo_1e_int_potential_n_e_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_float), intent(out) :: dset(*)
   end function trexio_read_mo_1e_int_potential_n_e_32
end interface

interface
   integer(trexio_exit_code) function trexio_read_mo_1e_int_ecp_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_float), intent(out) :: dset(*)
   end function trexio_read_mo_1e_int_ecp_32
end interface

interface
   integer(trexio_exit_code) function trexio_read_mo_1e_int_core_hamiltonian_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_float), intent(out) :: dset(*)
   end function trexio_read_mo_1e_int_core_hamiltonian_32
end interface

interface
   integer(trexio_exit_code) function trexio_read_mo_1e_int_overlap_im_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_float), intent(out) :: dset(*)
   end function trexio_read_mo_1e_int_overlap_im_32
end interface

interface
   integer(trexio_exit_code) function trexio_read_mo_1e_int_kinetic_im_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_float), intent(out) :: dset(*)
   end function trexio_read_mo_1e_int_kinetic_im_32
end interface

interface
   integer(trexio_exit_code) function trexio_read_mo_1e_int_potential_n_e_im_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_float), intent(out) :: dset(*)
   end function trexio_read_mo_1e_int_potential_n_e_im_32
end interface

interface
   integer(trexio_exit_code) function trexio_read_mo_1e_int_ecp_im_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_float), intent(out) :: dset(*)
   end function trexio_read_mo_1e_int_ecp_im_32
end interface

interface
   integer(trexio_exit_code) function trexio_read_mo_1e_int_core_hamiltonian_im_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_float), intent(out) :: dset(*)
   end function trexio_read_mo_1e_int_core_hamiltonian_im_32
end interface

interface
   integer(trexio_exit_code) function trexio_read_rdm_1e_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_float), intent(out) :: dset(*)
   end function trexio_read_rdm_1e_32
end interface

interface
   integer(trexio_exit_code) function trexio_read_rdm_1e_up_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_float), intent(out) :: dset(*)
   end function trexio_read_rdm_1e_up_32
end interface

interface
   integer(trexio_exit_code) function trexio_read_rdm_1e_dn_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_float), intent(out) :: dset(*)
   end function trexio_read_rdm_1e_dn_32
end interface

interface
   integer(trexio_exit_code) function trexio_read_cell_a_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_float), intent(out) :: dset(*)
   end function trexio_read_cell_a_32
end interface

interface
   integer(trexio_exit_code) function trexio_read_cell_b_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_float), intent(out) :: dset(*)
   end function trexio_read_cell_b_32
end interface

interface
   integer(trexio_exit_code) function trexio_read_cell_c_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_float), intent(out) :: dset(*)
   end function trexio_read_cell_c_32
end interface

interface
   integer(trexio_exit_code) function trexio_read_pbc_k_point_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_float), intent(out) :: dset(*)
   end function trexio_read_pbc_k_point_32
end interface

interface
   integer(trexio_exit_code) function trexio_read_qmc_point_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_float), intent(out) :: dset(*)
   end function trexio_read_qmc_point_32
end interface

interface
   integer(trexio_exit_code) function trexio_read_qmc_psi_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_float), intent(out) :: dset(*)
   end function trexio_read_qmc_psi_32
end interface

interface
   integer(trexio_exit_code) function trexio_read_qmc_e_loc_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_float), intent(out) :: dset(*)
   end function trexio_read_qmc_e_loc_32
end interface

interface
   integer(trexio_exit_code) function trexio_read_nucleus_charge_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(out) :: dset(*)
   end function trexio_read_nucleus_charge_64
end interface

interface
   integer(trexio_exit_code) function trexio_read_nucleus_coord_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(out) :: dset(*)
   end function trexio_read_nucleus_coord_64
end interface

interface
   integer(trexio_exit_code) function trexio_read_ecp_max_ang_mom_plus_1_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int64_t), intent(out) :: dset(*)
   end function trexio_read_ecp_max_ang_mom_plus_1_64
end interface

interface
   integer(trexio_exit_code) function trexio_read_ecp_z_core_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int64_t), intent(out) :: dset(*)
   end function trexio_read_ecp_z_core_64
end interface

interface
   integer(trexio_exit_code) function trexio_read_ecp_ang_mom_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int64_t), intent(out) :: dset(*)
   end function trexio_read_ecp_ang_mom_64
end interface

interface
   integer(trexio_exit_code) function trexio_read_ecp_nucleus_index_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int64_t), intent(out) :: dset(*)
   end function trexio_read_ecp_nucleus_index_64
end interface

interface
   integer(trexio_exit_code) function trexio_read_ecp_exponent_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(out) :: dset(*)
   end function trexio_read_ecp_exponent_64
end interface

interface
   integer(trexio_exit_code) function trexio_read_ecp_coefficient_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(out) :: dset(*)
   end function trexio_read_ecp_coefficient_64
end interface

interface
   integer(trexio_exit_code) function trexio_read_ecp_power_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int64_t), intent(out) :: dset(*)
   end function trexio_read_ecp_power_64
end interface

interface
   integer(trexio_exit_code) function trexio_read_basis_nucleus_index_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int64_t), intent(out) :: dset(*)
   end function trexio_read_basis_nucleus_index_64
end interface

interface
   integer(trexio_exit_code) function trexio_read_basis_shell_ang_mom_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int64_t), intent(out) :: dset(*)
   end function trexio_read_basis_shell_ang_mom_64
end interface

interface
   integer(trexio_exit_code) function trexio_read_basis_shell_factor_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(out) :: dset(*)
   end function trexio_read_basis_shell_factor_64
end interface

interface
   integer(trexio_exit_code) function trexio_read_basis_shell_index_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int64_t), intent(out) :: dset(*)
   end function trexio_read_basis_shell_index_64
end interface

interface
   integer(trexio_exit_code) function trexio_read_basis_exponent_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(out) :: dset(*)
   end function trexio_read_basis_exponent_64
end interface

interface
   integer(trexio_exit_code) function trexio_read_basis_coefficient_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(out) :: dset(*)
   end function trexio_read_basis_coefficient_64
end interface

interface
   integer(trexio_exit_code) function trexio_read_basis_prim_factor_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(out) :: dset(*)
   end function trexio_read_basis_prim_factor_64
end interface

interface
   integer(trexio_exit_code) function trexio_read_ao_shell_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int64_t), intent(out) :: dset(*)
   end function trexio_read_ao_shell_64
end interface

interface
   integer(trexio_exit_code) function trexio_read_ao_normalization_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(out) :: dset(*)
   end function trexio_read_ao_normalization_64
end interface

interface
   integer(trexio_exit_code) function trexio_read_ao_1e_int_overlap_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(out) :: dset(*)
   end function trexio_read_ao_1e_int_overlap_64
end interface

interface
   integer(trexio_exit_code) function trexio_read_ao_1e_int_kinetic_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(out) :: dset(*)
   end function trexio_read_ao_1e_int_kinetic_64
end interface

interface
   integer(trexio_exit_code) function trexio_read_ao_1e_int_potential_n_e_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(out) :: dset(*)
   end function trexio_read_ao_1e_int_potential_n_e_64
end interface

interface
   integer(trexio_exit_code) function trexio_read_ao_1e_int_ecp_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(out) :: dset(*)
   end function trexio_read_ao_1e_int_ecp_64
end interface

interface
   integer(trexio_exit_code) function trexio_read_ao_1e_int_core_hamiltonian_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(out) :: dset(*)
   end function trexio_read_ao_1e_int_core_hamiltonian_64
end interface

interface
   integer(trexio_exit_code) function trexio_read_ao_1e_int_overlap_im_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(out) :: dset(*)
   end function trexio_read_ao_1e_int_overlap_im_64
end interface

interface
   integer(trexio_exit_code) function trexio_read_ao_1e_int_kinetic_im_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(out) :: dset(*)
   end function trexio_read_ao_1e_int_kinetic_im_64
end interface

interface
   integer(trexio_exit_code) function trexio_read_ao_1e_int_potential_n_e_im_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(out) :: dset(*)
   end function trexio_read_ao_1e_int_potential_n_e_im_64
end interface

interface
   integer(trexio_exit_code) function trexio_read_ao_1e_int_ecp_im_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(out) :: dset(*)
   end function trexio_read_ao_1e_int_ecp_im_64
end interface

interface
   integer(trexio_exit_code) function trexio_read_ao_1e_int_core_hamiltonian_im_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(out) :: dset(*)
   end function trexio_read_ao_1e_int_core_hamiltonian_im_64
end interface

interface
   integer(trexio_exit_code) function trexio_read_mo_coefficient_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(out) :: dset(*)
   end function trexio_read_mo_coefficient_64
end interface

interface
   integer(trexio_exit_code) function trexio_read_mo_coefficient_im_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(out) :: dset(*)
   end function trexio_read_mo_coefficient_im_64
end interface

interface
   integer(trexio_exit_code) function trexio_read_mo_occupation_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(out) :: dset(*)
   end function trexio_read_mo_occupation_64
end interface

interface
   integer(trexio_exit_code) function trexio_read_mo_energy_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(out) :: dset(*)
   end function trexio_read_mo_energy_64
end interface

interface
   integer(trexio_exit_code) function trexio_read_mo_spin_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int64_t), intent(out) :: dset(*)
   end function trexio_read_mo_spin_64
end interface

interface
   integer(trexio_exit_code) function trexio_read_mo_1e_int_overlap_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(out) :: dset(*)
   end function trexio_read_mo_1e_int_overlap_64
end interface

interface
   integer(trexio_exit_code) function trexio_read_mo_1e_int_kinetic_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(out) :: dset(*)
   end function trexio_read_mo_1e_int_kinetic_64
end interface

interface
   integer(trexio_exit_code) function trexio_read_mo_1e_int_potential_n_e_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(out) :: dset(*)
   end function trexio_read_mo_1e_int_potential_n_e_64
end interface

interface
   integer(trexio_exit_code) function trexio_read_mo_1e_int_ecp_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(out) :: dset(*)
   end function trexio_read_mo_1e_int_ecp_64
end interface

interface
   integer(trexio_exit_code) function trexio_read_mo_1e_int_core_hamiltonian_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(out) :: dset(*)
   end function trexio_read_mo_1e_int_core_hamiltonian_64
end interface

interface
   integer(trexio_exit_code) function trexio_read_mo_1e_int_overlap_im_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(out) :: dset(*)
   end function trexio_read_mo_1e_int_overlap_im_64
end interface

interface
   integer(trexio_exit_code) function trexio_read_mo_1e_int_kinetic_im_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(out) :: dset(*)
   end function trexio_read_mo_1e_int_kinetic_im_64
end interface

interface
   integer(trexio_exit_code) function trexio_read_mo_1e_int_potential_n_e_im_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(out) :: dset(*)
   end function trexio_read_mo_1e_int_potential_n_e_im_64
end interface

interface
   integer(trexio_exit_code) function trexio_read_mo_1e_int_ecp_im_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(out) :: dset(*)
   end function trexio_read_mo_1e_int_ecp_im_64
end interface

interface
   integer(trexio_exit_code) function trexio_read_mo_1e_int_core_hamiltonian_im_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(out) :: dset(*)
   end function trexio_read_mo_1e_int_core_hamiltonian_im_64
end interface

interface
   integer(trexio_exit_code) function trexio_read_rdm_1e_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(out) :: dset(*)
   end function trexio_read_rdm_1e_64
end interface

interface
   integer(trexio_exit_code) function trexio_read_rdm_1e_up_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(out) :: dset(*)
   end function trexio_read_rdm_1e_up_64
end interface

interface
   integer(trexio_exit_code) function trexio_read_rdm_1e_dn_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(out) :: dset(*)
   end function trexio_read_rdm_1e_dn_64
end interface

interface
   integer(trexio_exit_code) function trexio_read_cell_a_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(out) :: dset(*)
   end function trexio_read_cell_a_64
end interface

interface
   integer(trexio_exit_code) function trexio_read_cell_b_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(out) :: dset(*)
   end function trexio_read_cell_b_64
end interface

interface
   integer(trexio_exit_code) function trexio_read_cell_c_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(out) :: dset(*)
   end function trexio_read_cell_c_64
end interface

interface
   integer(trexio_exit_code) function trexio_read_pbc_k_point_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(out) :: dset(*)
   end function trexio_read_pbc_k_point_64
end interface

interface
   integer(trexio_exit_code) function trexio_read_qmc_point_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(out) :: dset(*)
   end function trexio_read_qmc_point_64
end interface

interface
   integer(trexio_exit_code) function trexio_read_qmc_psi_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(out) :: dset(*)
   end function trexio_read_qmc_psi_64
end interface

interface
   integer(trexio_exit_code) function trexio_read_qmc_e_loc_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(out) :: dset(*)
   end function trexio_read_qmc_e_loc_64
end interface

interface
   integer(trexio_exit_code) function trexio_read_nucleus_charge (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(out) :: dset(*)
   end function trexio_read_nucleus_charge
end interface

interface
   integer(trexio_exit_code) function trexio_read_nucleus_coord (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(out) :: dset(*)
   end function trexio_read_nucleus_coord
end interface

interface
   integer(trexio_exit_code) function trexio_read_ecp_max_ang_mom_plus_1 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int32_t), intent(out) :: dset(*)
   end function trexio_read_ecp_max_ang_mom_plus_1
end interface

interface
   integer(trexio_exit_code) function trexio_read_ecp_z_core (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int32_t), intent(out) :: dset(*)
   end function trexio_read_ecp_z_core
end interface

interface
   integer(trexio_exit_code) function trexio_read_ecp_ang_mom (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int32_t), intent(out) :: dset(*)
   end function trexio_read_ecp_ang_mom
end interface

interface
   integer(trexio_exit_code) function trexio_read_ecp_nucleus_index (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int32_t), intent(out) :: dset(*)
   end function trexio_read_ecp_nucleus_index
end interface

interface
   integer(trexio_exit_code) function trexio_read_ecp_exponent (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(out) :: dset(*)
   end function trexio_read_ecp_exponent
end interface

interface
   integer(trexio_exit_code) function trexio_read_ecp_coefficient (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(out) :: dset(*)
   end function trexio_read_ecp_coefficient
end interface

interface
   integer(trexio_exit_code) function trexio_read_ecp_power (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int32_t), intent(out) :: dset(*)
   end function trexio_read_ecp_power
end interface

interface
   integer(trexio_exit_code) function trexio_read_basis_nucleus_index (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int32_t), intent(out) :: dset(*)
   end function trexio_read_basis_nucleus_index
end interface

interface
   integer(trexio_exit_code) function trexio_read_basis_shell_ang_mom (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int32_t), intent(out) :: dset(*)
   end function trexio_read_basis_shell_ang_mom
end interface

interface
   integer(trexio_exit_code) function trexio_read_basis_shell_factor (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(out) :: dset(*)
   end function trexio_read_basis_shell_factor
end interface

interface
   integer(trexio_exit_code) function trexio_read_basis_shell_index (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int32_t), intent(out) :: dset(*)
   end function trexio_read_basis_shell_index
end interface

interface
   integer(trexio_exit_code) function trexio_read_basis_exponent (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(out) :: dset(*)
   end function trexio_read_basis_exponent
end interface

interface
   integer(trexio_exit_code) function trexio_read_basis_coefficient (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(out) :: dset(*)
   end function trexio_read_basis_coefficient
end interface

interface
   integer(trexio_exit_code) function trexio_read_basis_prim_factor (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(out) :: dset(*)
   end function trexio_read_basis_prim_factor
end interface

interface
   integer(trexio_exit_code) function trexio_read_ao_shell (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int32_t), intent(out) :: dset(*)
   end function trexio_read_ao_shell
end interface

interface
   integer(trexio_exit_code) function trexio_read_ao_normalization (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(out) :: dset(*)
   end function trexio_read_ao_normalization
end interface

interface
   integer(trexio_exit_code) function trexio_read_ao_1e_int_overlap (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(out) :: dset(*)
   end function trexio_read_ao_1e_int_overlap
end interface

interface
   integer(trexio_exit_code) function trexio_read_ao_1e_int_kinetic (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(out) :: dset(*)
   end function trexio_read_ao_1e_int_kinetic
end interface

interface
   integer(trexio_exit_code) function trexio_read_ao_1e_int_potential_n_e (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(out) :: dset(*)
   end function trexio_read_ao_1e_int_potential_n_e
end interface

interface
   integer(trexio_exit_code) function trexio_read_ao_1e_int_ecp (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(out) :: dset(*)
   end function trexio_read_ao_1e_int_ecp
end interface

interface
   integer(trexio_exit_code) function trexio_read_ao_1e_int_core_hamiltonian (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(out) :: dset(*)
   end function trexio_read_ao_1e_int_core_hamiltonian
end interface

interface
   integer(trexio_exit_code) function trexio_read_ao_1e_int_overlap_im (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(out) :: dset(*)
   end function trexio_read_ao_1e_int_overlap_im
end interface

interface
   integer(trexio_exit_code) function trexio_read_ao_1e_int_kinetic_im (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(out) :: dset(*)
   end function trexio_read_ao_1e_int_kinetic_im
end interface

interface
   integer(trexio_exit_code) function trexio_read_ao_1e_int_potential_n_e_im (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(out) :: dset(*)
   end function trexio_read_ao_1e_int_potential_n_e_im
end interface

interface
   integer(trexio_exit_code) function trexio_read_ao_1e_int_ecp_im (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(out) :: dset(*)
   end function trexio_read_ao_1e_int_ecp_im
end interface

interface
   integer(trexio_exit_code) function trexio_read_ao_1e_int_core_hamiltonian_im (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(out) :: dset(*)
   end function trexio_read_ao_1e_int_core_hamiltonian_im
end interface

interface
   integer(trexio_exit_code) function trexio_read_mo_coefficient (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(out) :: dset(*)
   end function trexio_read_mo_coefficient
end interface

interface
   integer(trexio_exit_code) function trexio_read_mo_coefficient_im (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(out) :: dset(*)
   end function trexio_read_mo_coefficient_im
end interface

interface
   integer(trexio_exit_code) function trexio_read_mo_occupation (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(out) :: dset(*)
   end function trexio_read_mo_occupation
end interface

interface
   integer(trexio_exit_code) function trexio_read_mo_energy (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(out) :: dset(*)
   end function trexio_read_mo_energy
end interface

interface
   integer(trexio_exit_code) function trexio_read_mo_spin (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int32_t), intent(out) :: dset(*)
   end function trexio_read_mo_spin
end interface

interface
   integer(trexio_exit_code) function trexio_read_mo_1e_int_overlap (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(out) :: dset(*)
   end function trexio_read_mo_1e_int_overlap
end interface

interface
   integer(trexio_exit_code) function trexio_read_mo_1e_int_kinetic (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(out) :: dset(*)
   end function trexio_read_mo_1e_int_kinetic
end interface

interface
   integer(trexio_exit_code) function trexio_read_mo_1e_int_potential_n_e (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(out) :: dset(*)
   end function trexio_read_mo_1e_int_potential_n_e
end interface

interface
   integer(trexio_exit_code) function trexio_read_mo_1e_int_ecp (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(out) :: dset(*)
   end function trexio_read_mo_1e_int_ecp
end interface

interface
   integer(trexio_exit_code) function trexio_read_mo_1e_int_core_hamiltonian (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(out) :: dset(*)
   end function trexio_read_mo_1e_int_core_hamiltonian
end interface

interface
   integer(trexio_exit_code) function trexio_read_mo_1e_int_overlap_im (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(out) :: dset(*)
   end function trexio_read_mo_1e_int_overlap_im
end interface

interface
   integer(trexio_exit_code) function trexio_read_mo_1e_int_kinetic_im (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(out) :: dset(*)
   end function trexio_read_mo_1e_int_kinetic_im
end interface

interface
   integer(trexio_exit_code) function trexio_read_mo_1e_int_potential_n_e_im (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(out) :: dset(*)
   end function trexio_read_mo_1e_int_potential_n_e_im
end interface

interface
   integer(trexio_exit_code) function trexio_read_mo_1e_int_ecp_im (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(out) :: dset(*)
   end function trexio_read_mo_1e_int_ecp_im
end interface

interface
   integer(trexio_exit_code) function trexio_read_mo_1e_int_core_hamiltonian_im (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(out) :: dset(*)
   end function trexio_read_mo_1e_int_core_hamiltonian_im
end interface

interface
   integer(trexio_exit_code) function trexio_read_rdm_1e (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(out) :: dset(*)
   end function trexio_read_rdm_1e
end interface

interface
   integer(trexio_exit_code) function trexio_read_rdm_1e_up (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(out) :: dset(*)
   end function trexio_read_rdm_1e_up
end interface

interface
   integer(trexio_exit_code) function trexio_read_rdm_1e_dn (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(out) :: dset(*)
   end function trexio_read_rdm_1e_dn
end interface

interface
   integer(trexio_exit_code) function trexio_read_cell_a (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(out) :: dset(*)
   end function trexio_read_cell_a
end interface

interface
   integer(trexio_exit_code) function trexio_read_cell_b (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(out) :: dset(*)
   end function trexio_read_cell_b
end interface

interface
   integer(trexio_exit_code) function trexio_read_cell_c (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(out) :: dset(*)
   end function trexio_read_cell_c
end interface

interface
   integer(trexio_exit_code) function trexio_read_pbc_k_point (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(out) :: dset(*)
   end function trexio_read_pbc_k_point
end interface

interface
   integer(trexio_exit_code) function trexio_read_qmc_point (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(out) :: dset(*)
   end function trexio_read_qmc_point
end interface

interface
   integer(trexio_exit_code) function trexio_read_qmc_psi (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(out) :: dset(*)
   end function trexio_read_qmc_psi
end interface

interface
   integer(trexio_exit_code) function trexio_read_qmc_e_loc (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(out) :: dset(*)
   end function trexio_read_qmc_e_loc
end interface

interface
   integer(trexio_exit_code) function trexio_read_ao_2e_int_eri (trex_file, &
                                              offset_file, buffer_size, &
                                              index_sparse, value_sparse) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int64_t), intent(in), value :: offset_file
     integer(c_int64_t), intent(inout)     :: buffer_size
     integer(c_int32_t), intent(out)       :: index_sparse(*)
     real(c_double), intent(out)           :: value_sparse(*)
   end function trexio_read_ao_2e_int_eri
end interface

interface
   integer(trexio_exit_code) function trexio_read_safe_ao_2e_int_eri (trex_file, &
                                                   offset_file, buffer_size, &
                                                   index_sparse, index_size, &
                                                   value_sparse, value_size) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int64_t), intent(in), value :: offset_file
     integer(c_int64_t), intent(inout)     :: buffer_size
     integer(c_int32_t), intent(out)       :: index_sparse(*)
     integer(c_int64_t), intent(in), value :: index_size
     real(c_double), intent(out)           :: value_sparse(*)
     integer(c_int64_t), intent(in), value :: value_size
   end function trexio_read_safe_ao_2e_int_eri
end interface

interface
   integer(trexio_exit_code) function trexio_read_ao_2e_int_eri_lr (trex_file, &
                                              offset_file, buffer_size, &
                                              index_sparse, value_sparse) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int64_t), intent(in), value :: offset_file
     integer(c_int64_t), intent(inout)     :: buffer_size
     integer(c_int32_t), intent(out)       :: index_sparse(*)
     real(c_double), intent(out)           :: value_sparse(*)
   end function trexio_read_ao_2e_int_eri_lr
end interface

interface
   integer(trexio_exit_code) function trexio_read_safe_ao_2e_int_eri_lr (trex_file, &
                                                   offset_file, buffer_size, &
                                                   index_sparse, index_size, &
                                                   value_sparse, value_size) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int64_t), intent(in), value :: offset_file
     integer(c_int64_t), intent(inout)     :: buffer_size
     integer(c_int32_t), intent(out)       :: index_sparse(*)
     integer(c_int64_t), intent(in), value :: index_size
     real(c_double), intent(out)           :: value_sparse(*)
     integer(c_int64_t), intent(in), value :: value_size
   end function trexio_read_safe_ao_2e_int_eri_lr
end interface

interface
   integer(trexio_exit_code) function trexio_read_mo_2e_int_eri (trex_file, &
                                              offset_file, buffer_size, &
                                              index_sparse, value_sparse) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int64_t), intent(in), value :: offset_file
     integer(c_int64_t), intent(inout)     :: buffer_size
     integer(c_int32_t), intent(out)       :: index_sparse(*)
     real(c_double), intent(out)           :: value_sparse(*)
   end function trexio_read_mo_2e_int_eri
end interface

interface
   integer(trexio_exit_code) function trexio_read_safe_mo_2e_int_eri (trex_file, &
                                                   offset_file, buffer_size, &
                                                   index_sparse, index_size, &
                                                   value_sparse, value_size) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int64_t), intent(in), value :: offset_file
     integer(c_int64_t), intent(inout)     :: buffer_size
     integer(c_int32_t), intent(out)       :: index_sparse(*)
     integer(c_int64_t), intent(in), value :: index_size
     real(c_double), intent(out)           :: value_sparse(*)
     integer(c_int64_t), intent(in), value :: value_size
   end function trexio_read_safe_mo_2e_int_eri
end interface

interface
   integer(trexio_exit_code) function trexio_read_mo_2e_int_eri_lr (trex_file, &
                                              offset_file, buffer_size, &
                                              index_sparse, value_sparse) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int64_t), intent(in), value :: offset_file
     integer(c_int64_t), intent(inout)     :: buffer_size
     integer(c_int32_t), intent(out)       :: index_sparse(*)
     real(c_double), intent(out)           :: value_sparse(*)
   end function trexio_read_mo_2e_int_eri_lr
end interface

interface
   integer(trexio_exit_code) function trexio_read_safe_mo_2e_int_eri_lr (trex_file, &
                                                   offset_file, buffer_size, &
                                                   index_sparse, index_size, &
                                                   value_sparse, value_size) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int64_t), intent(in), value :: offset_file
     integer(c_int64_t), intent(inout)     :: buffer_size
     integer(c_int32_t), intent(out)       :: index_sparse(*)
     integer(c_int64_t), intent(in), value :: index_size
     real(c_double), intent(out)           :: value_sparse(*)
     integer(c_int64_t), intent(in), value :: value_size
   end function trexio_read_safe_mo_2e_int_eri_lr
end interface

interface
   integer(trexio_exit_code) function trexio_read_csf_det_coefficient (trex_file, &
                                              offset_file, buffer_size, &
                                              index_sparse, value_sparse) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int64_t), intent(in), value :: offset_file
     integer(c_int64_t), intent(inout)     :: buffer_size
     integer(c_int32_t), intent(out)       :: index_sparse(*)
     real(c_double), intent(out)           :: value_sparse(*)
   end function trexio_read_csf_det_coefficient
end interface

interface
   integer(trexio_exit_code) function trexio_read_safe_csf_det_coefficient (trex_file, &
                                                   offset_file, buffer_size, &
                                                   index_sparse, index_size, &
                                                   value_sparse, value_size) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int64_t), intent(in), value :: offset_file
     integer(c_int64_t), intent(inout)     :: buffer_size
     integer(c_int32_t), intent(out)       :: index_sparse(*)
     integer(c_int64_t), intent(in), value :: index_size
     real(c_double), intent(out)           :: value_sparse(*)
     integer(c_int64_t), intent(in), value :: value_size
   end function trexio_read_safe_csf_det_coefficient
end interface

interface
   integer(trexio_exit_code) function trexio_read_rdm_2e (trex_file, &
                                              offset_file, buffer_size, &
                                              index_sparse, value_sparse) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int64_t), intent(in), value :: offset_file
     integer(c_int64_t), intent(inout)     :: buffer_size
     integer(c_int32_t), intent(out)       :: index_sparse(*)
     real(c_double), intent(out)           :: value_sparse(*)
   end function trexio_read_rdm_2e
end interface

interface
   integer(trexio_exit_code) function trexio_read_safe_rdm_2e (trex_file, &
                                                   offset_file, buffer_size, &
                                                   index_sparse, index_size, &
                                                   value_sparse, value_size) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int64_t), intent(in), value :: offset_file
     integer(c_int64_t), intent(inout)     :: buffer_size
     integer(c_int32_t), intent(out)       :: index_sparse(*)
     integer(c_int64_t), intent(in), value :: index_size
     real(c_double), intent(out)           :: value_sparse(*)
     integer(c_int64_t), intent(in), value :: value_size
   end function trexio_read_safe_rdm_2e
end interface

interface
   integer(trexio_exit_code) function trexio_read_rdm_2e_upup (trex_file, &
                                              offset_file, buffer_size, &
                                              index_sparse, value_sparse) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int64_t), intent(in), value :: offset_file
     integer(c_int64_t), intent(inout)     :: buffer_size
     integer(c_int32_t), intent(out)       :: index_sparse(*)
     real(c_double), intent(out)           :: value_sparse(*)
   end function trexio_read_rdm_2e_upup
end interface

interface
   integer(trexio_exit_code) function trexio_read_safe_rdm_2e_upup (trex_file, &
                                                   offset_file, buffer_size, &
                                                   index_sparse, index_size, &
                                                   value_sparse, value_size) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int64_t), intent(in), value :: offset_file
     integer(c_int64_t), intent(inout)     :: buffer_size
     integer(c_int32_t), intent(out)       :: index_sparse(*)
     integer(c_int64_t), intent(in), value :: index_size
     real(c_double), intent(out)           :: value_sparse(*)
     integer(c_int64_t), intent(in), value :: value_size
   end function trexio_read_safe_rdm_2e_upup
end interface

interface
   integer(trexio_exit_code) function trexio_read_rdm_2e_dndn (trex_file, &
                                              offset_file, buffer_size, &
                                              index_sparse, value_sparse) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int64_t), intent(in), value :: offset_file
     integer(c_int64_t), intent(inout)     :: buffer_size
     integer(c_int32_t), intent(out)       :: index_sparse(*)
     real(c_double), intent(out)           :: value_sparse(*)
   end function trexio_read_rdm_2e_dndn
end interface

interface
   integer(trexio_exit_code) function trexio_read_safe_rdm_2e_dndn (trex_file, &
                                                   offset_file, buffer_size, &
                                                   index_sparse, index_size, &
                                                   value_sparse, value_size) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int64_t), intent(in), value :: offset_file
     integer(c_int64_t), intent(inout)     :: buffer_size
     integer(c_int32_t), intent(out)       :: index_sparse(*)
     integer(c_int64_t), intent(in), value :: index_size
     real(c_double), intent(out)           :: value_sparse(*)
     integer(c_int64_t), intent(in), value :: value_size
   end function trexio_read_safe_rdm_2e_dndn
end interface

interface
   integer(trexio_exit_code) function trexio_read_rdm_2e_updn (trex_file, &
                                              offset_file, buffer_size, &
                                              index_sparse, value_sparse) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int64_t), intent(in), value :: offset_file
     integer(c_int64_t), intent(inout)     :: buffer_size
     integer(c_int32_t), intent(out)       :: index_sparse(*)
     real(c_double), intent(out)           :: value_sparse(*)
   end function trexio_read_rdm_2e_updn
end interface

interface
   integer(trexio_exit_code) function trexio_read_safe_rdm_2e_updn (trex_file, &
                                                   offset_file, buffer_size, &
                                                   index_sparse, index_size, &
                                                   value_sparse, value_size) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int64_t), intent(in), value :: offset_file
     integer(c_int64_t), intent(inout)     :: buffer_size
     integer(c_int32_t), intent(out)       :: index_sparse(*)
     integer(c_int64_t), intent(in), value :: index_size
     real(c_double), intent(out)           :: value_sparse(*)
     integer(c_int64_t), intent(in), value :: value_size
   end function trexio_read_safe_rdm_2e_updn
end interface

interface
   integer(trexio_exit_code) function trexio_read_rdm_2e_dnup (trex_file, &
                                              offset_file, buffer_size, &
                                              index_sparse, value_sparse) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int64_t), intent(in), value :: offset_file
     integer(c_int64_t), intent(inout)     :: buffer_size
     integer(c_int32_t), intent(out)       :: index_sparse(*)
     real(c_double), intent(out)           :: value_sparse(*)
   end function trexio_read_rdm_2e_dnup
end interface

interface
   integer(trexio_exit_code) function trexio_read_safe_rdm_2e_dnup (trex_file, &
                                                   offset_file, buffer_size, &
                                                   index_sparse, index_size, &
                                                   value_sparse, value_size) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int64_t), intent(in), value :: offset_file
     integer(c_int64_t), intent(inout)     :: buffer_size
     integer(c_int32_t), intent(out)       :: index_sparse(*)
     integer(c_int64_t), intent(in), value :: index_size
     real(c_double), intent(out)           :: value_sparse(*)
     integer(c_int64_t), intent(in), value :: value_size
   end function trexio_read_safe_rdm_2e_dnup
end interface

interface
   integer(trexio_exit_code) function trexio_read_ao_2e_int_eri_size (trex_file, &
                                                   size_max) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int64_t), intent(out) :: size_max
   end function trexio_read_ao_2e_int_eri_size
end interface

interface
   integer(trexio_exit_code) function trexio_read_ao_2e_int_eri_lr_size (trex_file, &
                                                   size_max) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int64_t), intent(out) :: size_max
   end function trexio_read_ao_2e_int_eri_lr_size
end interface

interface
   integer(trexio_exit_code) function trexio_read_mo_2e_int_eri_size (trex_file, &
                                                   size_max) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int64_t), intent(out) :: size_max
   end function trexio_read_mo_2e_int_eri_size
end interface

interface
   integer(trexio_exit_code) function trexio_read_mo_2e_int_eri_lr_size (trex_file, &
                                                   size_max) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int64_t), intent(out) :: size_max
   end function trexio_read_mo_2e_int_eri_lr_size
end interface

interface
   integer(trexio_exit_code) function trexio_read_csf_det_coefficient_size (trex_file, &
                                                   size_max) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int64_t), intent(out) :: size_max
   end function trexio_read_csf_det_coefficient_size
end interface

interface
   integer(trexio_exit_code) function trexio_read_rdm_2e_size (trex_file, &
                                                   size_max) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int64_t), intent(out) :: size_max
   end function trexio_read_rdm_2e_size
end interface

interface
   integer(trexio_exit_code) function trexio_read_rdm_2e_upup_size (trex_file, &
                                                   size_max) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int64_t), intent(out) :: size_max
   end function trexio_read_rdm_2e_upup_size
end interface

interface
   integer(trexio_exit_code) function trexio_read_rdm_2e_dndn_size (trex_file, &
                                                   size_max) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int64_t), intent(out) :: size_max
   end function trexio_read_rdm_2e_dndn_size
end interface

interface
   integer(trexio_exit_code) function trexio_read_rdm_2e_updn_size (trex_file, &
                                                   size_max) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int64_t), intent(out) :: size_max
   end function trexio_read_rdm_2e_updn_size
end interface

interface
   integer(trexio_exit_code) function trexio_read_rdm_2e_dnup_size (trex_file, &
                                                   size_max) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int64_t), intent(out) :: size_max
   end function trexio_read_rdm_2e_dnup_size
end interface

interface
   integer(trexio_exit_code) function trexio_read_metadata_code_low (trex_file, dset, max_str_len) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(c_int64_t), intent(in), value  :: trex_file
     character(kind=c_char), intent(out)    :: dset(*)
     integer(c_int32_t), intent(in), value  :: max_str_len
   end function trexio_read_metadata_code_low
end interface

interface
   integer(trexio_exit_code) function trexio_read_metadata_author_low (trex_file, dset, max_str_len) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(c_int64_t), intent(in), value  :: trex_file
     character(kind=c_char), intent(out)    :: dset(*)
     integer(c_int32_t), intent(in), value  :: max_str_len
   end function trexio_read_metadata_author_low
end interface

interface
   integer(trexio_exit_code) function trexio_read_nucleus_label_low (trex_file, dset, max_str_len) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(c_int64_t), intent(in), value  :: trex_file
     character(kind=c_char), intent(out)    :: dset(*)
     integer(c_int32_t), intent(in), value  :: max_str_len
   end function trexio_read_nucleus_label_low
end interface

interface
   integer(trexio_exit_code) function trexio_read_mo_class_low (trex_file, dset, max_str_len) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(c_int64_t), intent(in), value  :: trex_file
     character(kind=c_char), intent(out)    :: dset(*)
     integer(c_int32_t), intent(in), value  :: max_str_len
   end function trexio_read_mo_class_low
end interface

interface
   integer(trexio_exit_code) function trexio_read_mo_symmetry_low (trex_file, dset, max_str_len) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(c_int64_t), intent(in), value  :: trex_file
     character(kind=c_char), intent(out)    :: dset(*)
     integer(c_int32_t), intent(in), value  :: max_str_len
   end function trexio_read_mo_symmetry_low
end interface

interface
   integer(trexio_exit_code) function trexio_read_state_label_low (trex_file, dset, max_str_len) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(c_int64_t), intent(in), value  :: trex_file
     character(kind=c_char), intent(out)    :: dset(*)
     integer(c_int32_t), intent(in), value  :: max_str_len
   end function trexio_read_state_label_low
end interface

interface
   integer(trexio_exit_code) function trexio_read_state_file_name_low (trex_file, dset, max_str_len) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(c_int64_t), intent(in), value  :: trex_file
     character(kind=c_char), intent(out)    :: dset(*)
     integer(c_int32_t), intent(in), value  :: max_str_len
   end function trexio_read_state_file_name_low
end interface

interface
   integer(trexio_exit_code) function trexio_write_metadata_code_num_32 (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int32_t), intent(in), value :: num
   end function trexio_write_metadata_code_num_32
end interface

interface
   integer(trexio_exit_code) function trexio_write_metadata_author_num_32 (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int32_t), intent(in), value :: num
   end function trexio_write_metadata_author_num_32
end interface

interface
   integer(trexio_exit_code) function trexio_write_metadata_unsafe_32 (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int32_t), intent(in), value :: num
   end function trexio_write_metadata_unsafe_32
end interface

interface
   integer(trexio_exit_code) function trexio_write_electron_num_32 (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int32_t), intent(in), value :: num
   end function trexio_write_electron_num_32
end interface

interface
   integer(trexio_exit_code) function trexio_write_electron_up_num_32 (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int32_t), intent(in), value :: num
   end function trexio_write_electron_up_num_32
end interface

interface
   integer(trexio_exit_code) function trexio_write_electron_dn_num_32 (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int32_t), intent(in), value :: num
   end function trexio_write_electron_dn_num_32
end interface

interface
   integer(trexio_exit_code) function trexio_write_nucleus_num_32 (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int32_t), intent(in), value :: num
   end function trexio_write_nucleus_num_32
end interface

interface
   integer(trexio_exit_code) function trexio_write_nucleus_repulsion_32 (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_float), intent(in), value :: num
   end function trexio_write_nucleus_repulsion_32
end interface

interface
   integer(trexio_exit_code) function trexio_write_ecp_num_32 (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int32_t), intent(in), value :: num
   end function trexio_write_ecp_num_32
end interface

interface
   integer(trexio_exit_code) function trexio_write_basis_prim_num_32 (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int32_t), intent(in), value :: num
   end function trexio_write_basis_prim_num_32
end interface

interface
   integer(trexio_exit_code) function trexio_write_basis_shell_num_32 (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int32_t), intent(in), value :: num
   end function trexio_write_basis_shell_num_32
end interface

interface
   integer(trexio_exit_code) function trexio_write_ao_cartesian_32 (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int32_t), intent(in), value :: num
   end function trexio_write_ao_cartesian_32
end interface

interface
   integer(trexio_exit_code) function trexio_write_ao_num_32 (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int32_t), intent(in), value :: num
   end function trexio_write_ao_num_32
end interface

interface
   integer(trexio_exit_code) function trexio_write_mo_num_32 (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int32_t), intent(in), value :: num
   end function trexio_write_mo_num_32
end interface

interface
   integer(trexio_exit_code) function trexio_write_state_num_32 (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int32_t), intent(in), value :: num
   end function trexio_write_state_num_32
end interface

interface
   integer(trexio_exit_code) function trexio_write_state_id_32 (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int32_t), intent(in), value :: num
   end function trexio_write_state_id_32
end interface

interface
   integer(trexio_exit_code) function trexio_write_pbc_periodic_32 (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int32_t), intent(in), value :: num
   end function trexio_write_pbc_periodic_32
end interface

interface
   integer(trexio_exit_code) function trexio_write_qmc_num_32 (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int32_t), intent(in), value :: num
   end function trexio_write_qmc_num_32
end interface

interface
   integer(trexio_exit_code) function trexio_write_metadata_code_num_64 (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int64_t), intent(in), value :: num
   end function trexio_write_metadata_code_num_64
end interface

interface
   integer(trexio_exit_code) function trexio_write_metadata_author_num_64 (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int64_t), intent(in), value :: num
   end function trexio_write_metadata_author_num_64
end interface

interface
   integer(trexio_exit_code) function trexio_write_metadata_unsafe_64 (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int64_t), intent(in), value :: num
   end function trexio_write_metadata_unsafe_64
end interface

interface
   integer(trexio_exit_code) function trexio_write_electron_num_64 (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int64_t), intent(in), value :: num
   end function trexio_write_electron_num_64
end interface

interface
   integer(trexio_exit_code) function trexio_write_electron_up_num_64 (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int64_t), intent(in), value :: num
   end function trexio_write_electron_up_num_64
end interface

interface
   integer(trexio_exit_code) function trexio_write_electron_dn_num_64 (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int64_t), intent(in), value :: num
   end function trexio_write_electron_dn_num_64
end interface

interface
   integer(trexio_exit_code) function trexio_write_nucleus_num_64 (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int64_t), intent(in), value :: num
   end function trexio_write_nucleus_num_64
end interface

interface
   integer(trexio_exit_code) function trexio_write_nucleus_repulsion_64 (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(in), value :: num
   end function trexio_write_nucleus_repulsion_64
end interface

interface
   integer(trexio_exit_code) function trexio_write_ecp_num_64 (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int64_t), intent(in), value :: num
   end function trexio_write_ecp_num_64
end interface

interface
   integer(trexio_exit_code) function trexio_write_basis_prim_num_64 (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int64_t), intent(in), value :: num
   end function trexio_write_basis_prim_num_64
end interface

interface
   integer(trexio_exit_code) function trexio_write_basis_shell_num_64 (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int64_t), intent(in), value :: num
   end function trexio_write_basis_shell_num_64
end interface

interface
   integer(trexio_exit_code) function trexio_write_ao_cartesian_64 (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int64_t), intent(in), value :: num
   end function trexio_write_ao_cartesian_64
end interface

interface
   integer(trexio_exit_code) function trexio_write_ao_num_64 (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int64_t), intent(in), value :: num
   end function trexio_write_ao_num_64
end interface

interface
   integer(trexio_exit_code) function trexio_write_mo_num_64 (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int64_t), intent(in), value :: num
   end function trexio_write_mo_num_64
end interface

interface
   integer(trexio_exit_code) function trexio_write_state_num_64 (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int64_t), intent(in), value :: num
   end function trexio_write_state_num_64
end interface

interface
   integer(trexio_exit_code) function trexio_write_state_id_64 (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int64_t), intent(in), value :: num
   end function trexio_write_state_id_64
end interface

interface
   integer(trexio_exit_code) function trexio_write_pbc_periodic_64 (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int64_t), intent(in), value :: num
   end function trexio_write_pbc_periodic_64
end interface

interface
   integer(trexio_exit_code) function trexio_write_qmc_num_64 (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int64_t), intent(in), value :: num
   end function trexio_write_qmc_num_64
end interface

interface
   integer(trexio_exit_code) function trexio_write_metadata_code_num (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int32_t), intent(in), value :: num
   end function trexio_write_metadata_code_num
end interface

interface
   integer(trexio_exit_code) function trexio_write_metadata_author_num (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int32_t), intent(in), value :: num
   end function trexio_write_metadata_author_num
end interface

interface
   integer(trexio_exit_code) function trexio_write_metadata_unsafe (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int32_t), intent(in), value :: num
   end function trexio_write_metadata_unsafe
end interface

interface
   integer(trexio_exit_code) function trexio_write_electron_num (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int32_t), intent(in), value :: num
   end function trexio_write_electron_num
end interface

interface
   integer(trexio_exit_code) function trexio_write_electron_up_num (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int32_t), intent(in), value :: num
   end function trexio_write_electron_up_num
end interface

interface
   integer(trexio_exit_code) function trexio_write_electron_dn_num (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int32_t), intent(in), value :: num
   end function trexio_write_electron_dn_num
end interface

interface
   integer(trexio_exit_code) function trexio_write_nucleus_num (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int32_t), intent(in), value :: num
   end function trexio_write_nucleus_num
end interface

interface
   integer(trexio_exit_code) function trexio_write_nucleus_repulsion (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(in), value :: num
   end function trexio_write_nucleus_repulsion
end interface

interface
   integer(trexio_exit_code) function trexio_write_ecp_num (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int32_t), intent(in), value :: num
   end function trexio_write_ecp_num
end interface

interface
   integer(trexio_exit_code) function trexio_write_basis_prim_num (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int32_t), intent(in), value :: num
   end function trexio_write_basis_prim_num
end interface

interface
   integer(trexio_exit_code) function trexio_write_basis_shell_num (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int32_t), intent(in), value :: num
   end function trexio_write_basis_shell_num
end interface

interface
   integer(trexio_exit_code) function trexio_write_ao_cartesian (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int32_t), intent(in), value :: num
   end function trexio_write_ao_cartesian
end interface

interface
   integer(trexio_exit_code) function trexio_write_ao_num (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int32_t), intent(in), value :: num
   end function trexio_write_ao_num
end interface

interface
   integer(trexio_exit_code) function trexio_write_mo_num (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int32_t), intent(in), value :: num
   end function trexio_write_mo_num
end interface

interface
   integer(trexio_exit_code) function trexio_write_state_num (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int32_t), intent(in), value :: num
   end function trexio_write_state_num
end interface

interface
   integer(trexio_exit_code) function trexio_write_state_id (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int32_t), intent(in), value :: num
   end function trexio_write_state_id
end interface

interface
   integer(trexio_exit_code) function trexio_write_pbc_periodic (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int32_t), intent(in), value :: num
   end function trexio_write_pbc_periodic
end interface

interface
   integer(trexio_exit_code) function trexio_write_qmc_num (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int32_t), intent(in), value :: num
   end function trexio_write_qmc_num
end interface

interface
   integer(trexio_exit_code) function trexio_write_metadata_package_version_c (trex_file, str, max_str_len) &
           bind(C, name="trexio_write_metadata_package_version")
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value  :: trex_file
     character(kind=c_char), intent(in)    :: str(*)
     integer(c_int32_t), intent(in), value :: max_str_len
   end function trexio_write_metadata_package_version_c
end interface

interface
   integer(trexio_exit_code) function trexio_write_metadata_description_c (trex_file, str, max_str_len) &
           bind(C, name="trexio_write_metadata_description")
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value  :: trex_file
     character(kind=c_char), intent(in)    :: str(*)
     integer(c_int32_t), intent(in), value :: max_str_len
   end function trexio_write_metadata_description_c
end interface

interface
   integer(trexio_exit_code) function trexio_write_nucleus_point_group_c (trex_file, str, max_str_len) &
           bind(C, name="trexio_write_nucleus_point_group")
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value  :: trex_file
     character(kind=c_char), intent(in)    :: str(*)
     integer(c_int32_t), intent(in), value :: max_str_len
   end function trexio_write_nucleus_point_group_c
end interface

interface
   integer(trexio_exit_code) function trexio_write_basis_type_c (trex_file, str, max_str_len) &
           bind(C, name="trexio_write_basis_type")
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value  :: trex_file
     character(kind=c_char), intent(in)    :: str(*)
     integer(c_int32_t), intent(in), value :: max_str_len
   end function trexio_write_basis_type_c
end interface

interface
   integer(trexio_exit_code) function trexio_write_mo_type_c (trex_file, str, max_str_len) &
           bind(C, name="trexio_write_mo_type")
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value  :: trex_file
     character(kind=c_char), intent(in)    :: str(*)
     integer(c_int32_t), intent(in), value :: max_str_len
   end function trexio_write_mo_type_c
end interface

interface
   integer(trexio_exit_code) function trexio_write_state_current_label_c (trex_file, str, max_str_len) &
           bind(C, name="trexio_write_state_current_label")
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value  :: trex_file
     character(kind=c_char), intent(in)    :: str(*)
     integer(c_int32_t), intent(in), value :: max_str_len
   end function trexio_write_state_current_label_c
end interface

interface
   integer(trexio_exit_code) function trexio_write_determinant_coefficient(trex_file, &
                                               offset_file, buffer_size, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int64_t), intent(in), value :: offset_file
     integer(c_int64_t), intent(in), value :: buffer_size
     real(c_double), intent(in)            :: dset(*)
   end function trexio_write_determinant_coefficient
end interface

interface
   integer(trexio_exit_code) function trexio_write_safe_determinant_coefficient (trex_file, &
                                                    offset_file, buffer_size, &
                                                    dset, dset_size) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int64_t), intent(in), value :: offset_file
     integer(c_int64_t), intent(in), value :: buffer_size
     real(c_double), intent(in)            :: dset(*)
     integer(c_int64_t), intent(in), value :: dset_size
   end function trexio_write_safe_determinant_coefficient
end interface

interface
   integer(trexio_exit_code) function trexio_write_csf_coefficient(trex_file, &
                                               offset_file, buffer_size, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int64_t), intent(in), value :: offset_file
     integer(c_int64_t), intent(in), value :: buffer_size
     real(c_double), intent(in)            :: dset(*)
   end function trexio_write_csf_coefficient
end interface

interface
   integer(trexio_exit_code) function trexio_write_safe_csf_coefficient (trex_file, &
                                                    offset_file, buffer_size, &
                                                    dset, dset_size) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int64_t), intent(in), value :: offset_file
     integer(c_int64_t), intent(in), value :: buffer_size
     real(c_double), intent(in)            :: dset(*)
     integer(c_int64_t), intent(in), value :: dset_size
   end function trexio_write_safe_csf_coefficient
end interface

interface
   integer(trexio_exit_code) function trexio_write_nucleus_charge_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_float), intent(in) :: dset(*)
   end function trexio_write_nucleus_charge_32
end interface

interface
   integer(trexio_exit_code) function trexio_write_nucleus_coord_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_float), intent(in) :: dset(*)
   end function trexio_write_nucleus_coord_32
end interface

interface
   integer(trexio_exit_code) function trexio_write_ecp_max_ang_mom_plus_1_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int32_t), intent(in) :: dset(*)
   end function trexio_write_ecp_max_ang_mom_plus_1_32
end interface

interface
   integer(trexio_exit_code) function trexio_write_ecp_z_core_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int32_t), intent(in) :: dset(*)
   end function trexio_write_ecp_z_core_32
end interface

interface
   integer(trexio_exit_code) function trexio_write_ecp_ang_mom_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int32_t), intent(in) :: dset(*)
   end function trexio_write_ecp_ang_mom_32
end interface

interface
   integer(trexio_exit_code) function trexio_write_ecp_nucleus_index_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int32_t), intent(in) :: dset(*)
   end function trexio_write_ecp_nucleus_index_32
end interface

interface
   integer(trexio_exit_code) function trexio_write_ecp_exponent_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_float), intent(in) :: dset(*)
   end function trexio_write_ecp_exponent_32
end interface

interface
   integer(trexio_exit_code) function trexio_write_ecp_coefficient_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_float), intent(in) :: dset(*)
   end function trexio_write_ecp_coefficient_32
end interface

interface
   integer(trexio_exit_code) function trexio_write_ecp_power_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int32_t), intent(in) :: dset(*)
   end function trexio_write_ecp_power_32
end interface

interface
   integer(trexio_exit_code) function trexio_write_basis_nucleus_index_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int32_t), intent(in) :: dset(*)
   end function trexio_write_basis_nucleus_index_32
end interface

interface
   integer(trexio_exit_code) function trexio_write_basis_shell_ang_mom_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int32_t), intent(in) :: dset(*)
   end function trexio_write_basis_shell_ang_mom_32
end interface

interface
   integer(trexio_exit_code) function trexio_write_basis_shell_factor_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_float), intent(in) :: dset(*)
   end function trexio_write_basis_shell_factor_32
end interface

interface
   integer(trexio_exit_code) function trexio_write_basis_shell_index_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int32_t), intent(in) :: dset(*)
   end function trexio_write_basis_shell_index_32
end interface

interface
   integer(trexio_exit_code) function trexio_write_basis_exponent_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_float), intent(in) :: dset(*)
   end function trexio_write_basis_exponent_32
end interface

interface
   integer(trexio_exit_code) function trexio_write_basis_coefficient_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_float), intent(in) :: dset(*)
   end function trexio_write_basis_coefficient_32
end interface

interface
   integer(trexio_exit_code) function trexio_write_basis_prim_factor_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_float), intent(in) :: dset(*)
   end function trexio_write_basis_prim_factor_32
end interface

interface
   integer(trexio_exit_code) function trexio_write_ao_shell_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int32_t), intent(in) :: dset(*)
   end function trexio_write_ao_shell_32
end interface

interface
   integer(trexio_exit_code) function trexio_write_ao_normalization_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_float), intent(in) :: dset(*)
   end function trexio_write_ao_normalization_32
end interface

interface
   integer(trexio_exit_code) function trexio_write_ao_1e_int_overlap_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_float), intent(in) :: dset(*)
   end function trexio_write_ao_1e_int_overlap_32
end interface

interface
   integer(trexio_exit_code) function trexio_write_ao_1e_int_kinetic_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_float), intent(in) :: dset(*)
   end function trexio_write_ao_1e_int_kinetic_32
end interface

interface
   integer(trexio_exit_code) function trexio_write_ao_1e_int_potential_n_e_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_float), intent(in) :: dset(*)
   end function trexio_write_ao_1e_int_potential_n_e_32
end interface

interface
   integer(trexio_exit_code) function trexio_write_ao_1e_int_ecp_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_float), intent(in) :: dset(*)
   end function trexio_write_ao_1e_int_ecp_32
end interface

interface
   integer(trexio_exit_code) function trexio_write_ao_1e_int_core_hamiltonian_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_float), intent(in) :: dset(*)
   end function trexio_write_ao_1e_int_core_hamiltonian_32
end interface

interface
   integer(trexio_exit_code) function trexio_write_ao_1e_int_overlap_im_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_float), intent(in) :: dset(*)
   end function trexio_write_ao_1e_int_overlap_im_32
end interface

interface
   integer(trexio_exit_code) function trexio_write_ao_1e_int_kinetic_im_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_float), intent(in) :: dset(*)
   end function trexio_write_ao_1e_int_kinetic_im_32
end interface

interface
   integer(trexio_exit_code) function trexio_write_ao_1e_int_potential_n_e_im_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_float), intent(in) :: dset(*)
   end function trexio_write_ao_1e_int_potential_n_e_im_32
end interface

interface
   integer(trexio_exit_code) function trexio_write_ao_1e_int_ecp_im_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_float), intent(in) :: dset(*)
   end function trexio_write_ao_1e_int_ecp_im_32
end interface

interface
   integer(trexio_exit_code) function trexio_write_ao_1e_int_core_hamiltonian_im_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_float), intent(in) :: dset(*)
   end function trexio_write_ao_1e_int_core_hamiltonian_im_32
end interface

interface
   integer(trexio_exit_code) function trexio_write_mo_coefficient_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_float), intent(in) :: dset(*)
   end function trexio_write_mo_coefficient_32
end interface

interface
   integer(trexio_exit_code) function trexio_write_mo_coefficient_im_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_float), intent(in) :: dset(*)
   end function trexio_write_mo_coefficient_im_32
end interface

interface
   integer(trexio_exit_code) function trexio_write_mo_occupation_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_float), intent(in) :: dset(*)
   end function trexio_write_mo_occupation_32
end interface

interface
   integer(trexio_exit_code) function trexio_write_mo_energy_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_float), intent(in) :: dset(*)
   end function trexio_write_mo_energy_32
end interface

interface
   integer(trexio_exit_code) function trexio_write_mo_spin_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int32_t), intent(in) :: dset(*)
   end function trexio_write_mo_spin_32
end interface

interface
   integer(trexio_exit_code) function trexio_write_mo_1e_int_overlap_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_float), intent(in) :: dset(*)
   end function trexio_write_mo_1e_int_overlap_32
end interface

interface
   integer(trexio_exit_code) function trexio_write_mo_1e_int_kinetic_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_float), intent(in) :: dset(*)
   end function trexio_write_mo_1e_int_kinetic_32
end interface

interface
   integer(trexio_exit_code) function trexio_write_mo_1e_int_potential_n_e_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_float), intent(in) :: dset(*)
   end function trexio_write_mo_1e_int_potential_n_e_32
end interface

interface
   integer(trexio_exit_code) function trexio_write_mo_1e_int_ecp_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_float), intent(in) :: dset(*)
   end function trexio_write_mo_1e_int_ecp_32
end interface

interface
   integer(trexio_exit_code) function trexio_write_mo_1e_int_core_hamiltonian_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_float), intent(in) :: dset(*)
   end function trexio_write_mo_1e_int_core_hamiltonian_32
end interface

interface
   integer(trexio_exit_code) function trexio_write_mo_1e_int_overlap_im_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_float), intent(in) :: dset(*)
   end function trexio_write_mo_1e_int_overlap_im_32
end interface

interface
   integer(trexio_exit_code) function trexio_write_mo_1e_int_kinetic_im_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_float), intent(in) :: dset(*)
   end function trexio_write_mo_1e_int_kinetic_im_32
end interface

interface
   integer(trexio_exit_code) function trexio_write_mo_1e_int_potential_n_e_im_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_float), intent(in) :: dset(*)
   end function trexio_write_mo_1e_int_potential_n_e_im_32
end interface

interface
   integer(trexio_exit_code) function trexio_write_mo_1e_int_ecp_im_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_float), intent(in) :: dset(*)
   end function trexio_write_mo_1e_int_ecp_im_32
end interface

interface
   integer(trexio_exit_code) function trexio_write_mo_1e_int_core_hamiltonian_im_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_float), intent(in) :: dset(*)
   end function trexio_write_mo_1e_int_core_hamiltonian_im_32
end interface

interface
   integer(trexio_exit_code) function trexio_write_rdm_1e_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_float), intent(in) :: dset(*)
   end function trexio_write_rdm_1e_32
end interface

interface
   integer(trexio_exit_code) function trexio_write_rdm_1e_up_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_float), intent(in) :: dset(*)
   end function trexio_write_rdm_1e_up_32
end interface

interface
   integer(trexio_exit_code) function trexio_write_rdm_1e_dn_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_float), intent(in) :: dset(*)
   end function trexio_write_rdm_1e_dn_32
end interface

interface
   integer(trexio_exit_code) function trexio_write_cell_a_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_float), intent(in) :: dset(*)
   end function trexio_write_cell_a_32
end interface

interface
   integer(trexio_exit_code) function trexio_write_cell_b_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_float), intent(in) :: dset(*)
   end function trexio_write_cell_b_32
end interface

interface
   integer(trexio_exit_code) function trexio_write_cell_c_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_float), intent(in) :: dset(*)
   end function trexio_write_cell_c_32
end interface

interface
   integer(trexio_exit_code) function trexio_write_pbc_k_point_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_float), intent(in) :: dset(*)
   end function trexio_write_pbc_k_point_32
end interface

interface
   integer(trexio_exit_code) function trexio_write_qmc_point_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_float), intent(in) :: dset(*)
   end function trexio_write_qmc_point_32
end interface

interface
   integer(trexio_exit_code) function trexio_write_qmc_psi_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_float), intent(in) :: dset(*)
   end function trexio_write_qmc_psi_32
end interface

interface
   integer(trexio_exit_code) function trexio_write_qmc_e_loc_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_float), intent(in) :: dset(*)
   end function trexio_write_qmc_e_loc_32
end interface

interface
   integer(trexio_exit_code) function trexio_write_nucleus_charge_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(in) :: dset(*)
   end function trexio_write_nucleus_charge_64
end interface

interface
   integer(trexio_exit_code) function trexio_write_nucleus_coord_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(in) :: dset(*)
   end function trexio_write_nucleus_coord_64
end interface

interface
   integer(trexio_exit_code) function trexio_write_ecp_max_ang_mom_plus_1_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int64_t), intent(in) :: dset(*)
   end function trexio_write_ecp_max_ang_mom_plus_1_64
end interface

interface
   integer(trexio_exit_code) function trexio_write_ecp_z_core_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int64_t), intent(in) :: dset(*)
   end function trexio_write_ecp_z_core_64
end interface

interface
   integer(trexio_exit_code) function trexio_write_ecp_ang_mom_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int64_t), intent(in) :: dset(*)
   end function trexio_write_ecp_ang_mom_64
end interface

interface
   integer(trexio_exit_code) function trexio_write_ecp_nucleus_index_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int64_t), intent(in) :: dset(*)
   end function trexio_write_ecp_nucleus_index_64
end interface

interface
   integer(trexio_exit_code) function trexio_write_ecp_exponent_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(in) :: dset(*)
   end function trexio_write_ecp_exponent_64
end interface

interface
   integer(trexio_exit_code) function trexio_write_ecp_coefficient_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(in) :: dset(*)
   end function trexio_write_ecp_coefficient_64
end interface

interface
   integer(trexio_exit_code) function trexio_write_ecp_power_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int64_t), intent(in) :: dset(*)
   end function trexio_write_ecp_power_64
end interface

interface
   integer(trexio_exit_code) function trexio_write_basis_nucleus_index_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int64_t), intent(in) :: dset(*)
   end function trexio_write_basis_nucleus_index_64
end interface

interface
   integer(trexio_exit_code) function trexio_write_basis_shell_ang_mom_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int64_t), intent(in) :: dset(*)
   end function trexio_write_basis_shell_ang_mom_64
end interface

interface
   integer(trexio_exit_code) function trexio_write_basis_shell_factor_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(in) :: dset(*)
   end function trexio_write_basis_shell_factor_64
end interface

interface
   integer(trexio_exit_code) function trexio_write_basis_shell_index_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int64_t), intent(in) :: dset(*)
   end function trexio_write_basis_shell_index_64
end interface

interface
   integer(trexio_exit_code) function trexio_write_basis_exponent_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(in) :: dset(*)
   end function trexio_write_basis_exponent_64
end interface

interface
   integer(trexio_exit_code) function trexio_write_basis_coefficient_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(in) :: dset(*)
   end function trexio_write_basis_coefficient_64
end interface

interface
   integer(trexio_exit_code) function trexio_write_basis_prim_factor_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(in) :: dset(*)
   end function trexio_write_basis_prim_factor_64
end interface

interface
   integer(trexio_exit_code) function trexio_write_ao_shell_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int64_t), intent(in) :: dset(*)
   end function trexio_write_ao_shell_64
end interface

interface
   integer(trexio_exit_code) function trexio_write_ao_normalization_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(in) :: dset(*)
   end function trexio_write_ao_normalization_64
end interface

interface
   integer(trexio_exit_code) function trexio_write_ao_1e_int_overlap_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(in) :: dset(*)
   end function trexio_write_ao_1e_int_overlap_64
end interface

interface
   integer(trexio_exit_code) function trexio_write_ao_1e_int_kinetic_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(in) :: dset(*)
   end function trexio_write_ao_1e_int_kinetic_64
end interface

interface
   integer(trexio_exit_code) function trexio_write_ao_1e_int_potential_n_e_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(in) :: dset(*)
   end function trexio_write_ao_1e_int_potential_n_e_64
end interface

interface
   integer(trexio_exit_code) function trexio_write_ao_1e_int_ecp_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(in) :: dset(*)
   end function trexio_write_ao_1e_int_ecp_64
end interface

interface
   integer(trexio_exit_code) function trexio_write_ao_1e_int_core_hamiltonian_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(in) :: dset(*)
   end function trexio_write_ao_1e_int_core_hamiltonian_64
end interface

interface
   integer(trexio_exit_code) function trexio_write_ao_1e_int_overlap_im_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(in) :: dset(*)
   end function trexio_write_ao_1e_int_overlap_im_64
end interface

interface
   integer(trexio_exit_code) function trexio_write_ao_1e_int_kinetic_im_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(in) :: dset(*)
   end function trexio_write_ao_1e_int_kinetic_im_64
end interface

interface
   integer(trexio_exit_code) function trexio_write_ao_1e_int_potential_n_e_im_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(in) :: dset(*)
   end function trexio_write_ao_1e_int_potential_n_e_im_64
end interface

interface
   integer(trexio_exit_code) function trexio_write_ao_1e_int_ecp_im_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(in) :: dset(*)
   end function trexio_write_ao_1e_int_ecp_im_64
end interface

interface
   integer(trexio_exit_code) function trexio_write_ao_1e_int_core_hamiltonian_im_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(in) :: dset(*)
   end function trexio_write_ao_1e_int_core_hamiltonian_im_64
end interface

interface
   integer(trexio_exit_code) function trexio_write_mo_coefficient_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(in) :: dset(*)
   end function trexio_write_mo_coefficient_64
end interface

interface
   integer(trexio_exit_code) function trexio_write_mo_coefficient_im_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(in) :: dset(*)
   end function trexio_write_mo_coefficient_im_64
end interface

interface
   integer(trexio_exit_code) function trexio_write_mo_occupation_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(in) :: dset(*)
   end function trexio_write_mo_occupation_64
end interface

interface
   integer(trexio_exit_code) function trexio_write_mo_energy_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(in) :: dset(*)
   end function trexio_write_mo_energy_64
end interface

interface
   integer(trexio_exit_code) function trexio_write_mo_spin_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int64_t), intent(in) :: dset(*)
   end function trexio_write_mo_spin_64
end interface

interface
   integer(trexio_exit_code) function trexio_write_mo_1e_int_overlap_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(in) :: dset(*)
   end function trexio_write_mo_1e_int_overlap_64
end interface

interface
   integer(trexio_exit_code) function trexio_write_mo_1e_int_kinetic_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(in) :: dset(*)
   end function trexio_write_mo_1e_int_kinetic_64
end interface

interface
   integer(trexio_exit_code) function trexio_write_mo_1e_int_potential_n_e_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(in) :: dset(*)
   end function trexio_write_mo_1e_int_potential_n_e_64
end interface

interface
   integer(trexio_exit_code) function trexio_write_mo_1e_int_ecp_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(in) :: dset(*)
   end function trexio_write_mo_1e_int_ecp_64
end interface

interface
   integer(trexio_exit_code) function trexio_write_mo_1e_int_core_hamiltonian_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(in) :: dset(*)
   end function trexio_write_mo_1e_int_core_hamiltonian_64
end interface

interface
   integer(trexio_exit_code) function trexio_write_mo_1e_int_overlap_im_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(in) :: dset(*)
   end function trexio_write_mo_1e_int_overlap_im_64
end interface

interface
   integer(trexio_exit_code) function trexio_write_mo_1e_int_kinetic_im_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(in) :: dset(*)
   end function trexio_write_mo_1e_int_kinetic_im_64
end interface

interface
   integer(trexio_exit_code) function trexio_write_mo_1e_int_potential_n_e_im_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(in) :: dset(*)
   end function trexio_write_mo_1e_int_potential_n_e_im_64
end interface

interface
   integer(trexio_exit_code) function trexio_write_mo_1e_int_ecp_im_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(in) :: dset(*)
   end function trexio_write_mo_1e_int_ecp_im_64
end interface

interface
   integer(trexio_exit_code) function trexio_write_mo_1e_int_core_hamiltonian_im_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(in) :: dset(*)
   end function trexio_write_mo_1e_int_core_hamiltonian_im_64
end interface

interface
   integer(trexio_exit_code) function trexio_write_rdm_1e_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(in) :: dset(*)
   end function trexio_write_rdm_1e_64
end interface

interface
   integer(trexio_exit_code) function trexio_write_rdm_1e_up_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(in) :: dset(*)
   end function trexio_write_rdm_1e_up_64
end interface

interface
   integer(trexio_exit_code) function trexio_write_rdm_1e_dn_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(in) :: dset(*)
   end function trexio_write_rdm_1e_dn_64
end interface

interface
   integer(trexio_exit_code) function trexio_write_cell_a_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(in) :: dset(*)
   end function trexio_write_cell_a_64
end interface

interface
   integer(trexio_exit_code) function trexio_write_cell_b_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(in) :: dset(*)
   end function trexio_write_cell_b_64
end interface

interface
   integer(trexio_exit_code) function trexio_write_cell_c_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(in) :: dset(*)
   end function trexio_write_cell_c_64
end interface

interface
   integer(trexio_exit_code) function trexio_write_pbc_k_point_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(in) :: dset(*)
   end function trexio_write_pbc_k_point_64
end interface

interface
   integer(trexio_exit_code) function trexio_write_qmc_point_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(in) :: dset(*)
   end function trexio_write_qmc_point_64
end interface

interface
   integer(trexio_exit_code) function trexio_write_qmc_psi_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(in) :: dset(*)
   end function trexio_write_qmc_psi_64
end interface

interface
   integer(trexio_exit_code) function trexio_write_qmc_e_loc_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(in) :: dset(*)
   end function trexio_write_qmc_e_loc_64
end interface

interface
   integer(trexio_exit_code) function trexio_write_nucleus_charge (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(in) :: dset(*)
   end function trexio_write_nucleus_charge
end interface

interface
   integer(trexio_exit_code) function trexio_write_nucleus_coord (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(in) :: dset(*)
   end function trexio_write_nucleus_coord
end interface

interface
   integer(trexio_exit_code) function trexio_write_ecp_max_ang_mom_plus_1 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int32_t), intent(in) :: dset(*)
   end function trexio_write_ecp_max_ang_mom_plus_1
end interface

interface
   integer(trexio_exit_code) function trexio_write_ecp_z_core (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int32_t), intent(in) :: dset(*)
   end function trexio_write_ecp_z_core
end interface

interface
   integer(trexio_exit_code) function trexio_write_ecp_ang_mom (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int32_t), intent(in) :: dset(*)
   end function trexio_write_ecp_ang_mom
end interface

interface
   integer(trexio_exit_code) function trexio_write_ecp_nucleus_index (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int32_t), intent(in) :: dset(*)
   end function trexio_write_ecp_nucleus_index
end interface

interface
   integer(trexio_exit_code) function trexio_write_ecp_exponent (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(in) :: dset(*)
   end function trexio_write_ecp_exponent
end interface

interface
   integer(trexio_exit_code) function trexio_write_ecp_coefficient (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(in) :: dset(*)
   end function trexio_write_ecp_coefficient
end interface

interface
   integer(trexio_exit_code) function trexio_write_ecp_power (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int32_t), intent(in) :: dset(*)
   end function trexio_write_ecp_power
end interface

interface
   integer(trexio_exit_code) function trexio_write_basis_nucleus_index (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int32_t), intent(in) :: dset(*)
   end function trexio_write_basis_nucleus_index
end interface

interface
   integer(trexio_exit_code) function trexio_write_basis_shell_ang_mom (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int32_t), intent(in) :: dset(*)
   end function trexio_write_basis_shell_ang_mom
end interface

interface
   integer(trexio_exit_code) function trexio_write_basis_shell_factor (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(in) :: dset(*)
   end function trexio_write_basis_shell_factor
end interface

interface
   integer(trexio_exit_code) function trexio_write_basis_shell_index (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int32_t), intent(in) :: dset(*)
   end function trexio_write_basis_shell_index
end interface

interface
   integer(trexio_exit_code) function trexio_write_basis_exponent (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(in) :: dset(*)
   end function trexio_write_basis_exponent
end interface

interface
   integer(trexio_exit_code) function trexio_write_basis_coefficient (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(in) :: dset(*)
   end function trexio_write_basis_coefficient
end interface

interface
   integer(trexio_exit_code) function trexio_write_basis_prim_factor (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(in) :: dset(*)
   end function trexio_write_basis_prim_factor
end interface

interface
   integer(trexio_exit_code) function trexio_write_ao_shell (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int32_t), intent(in) :: dset(*)
   end function trexio_write_ao_shell
end interface

interface
   integer(trexio_exit_code) function trexio_write_ao_normalization (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(in) :: dset(*)
   end function trexio_write_ao_normalization
end interface

interface
   integer(trexio_exit_code) function trexio_write_ao_1e_int_overlap (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(in) :: dset(*)
   end function trexio_write_ao_1e_int_overlap
end interface

interface
   integer(trexio_exit_code) function trexio_write_ao_1e_int_kinetic (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(in) :: dset(*)
   end function trexio_write_ao_1e_int_kinetic
end interface

interface
   integer(trexio_exit_code) function trexio_write_ao_1e_int_potential_n_e (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(in) :: dset(*)
   end function trexio_write_ao_1e_int_potential_n_e
end interface

interface
   integer(trexio_exit_code) function trexio_write_ao_1e_int_ecp (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(in) :: dset(*)
   end function trexio_write_ao_1e_int_ecp
end interface

interface
   integer(trexio_exit_code) function trexio_write_ao_1e_int_core_hamiltonian (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(in) :: dset(*)
   end function trexio_write_ao_1e_int_core_hamiltonian
end interface

interface
   integer(trexio_exit_code) function trexio_write_ao_1e_int_overlap_im (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(in) :: dset(*)
   end function trexio_write_ao_1e_int_overlap_im
end interface

interface
   integer(trexio_exit_code) function trexio_write_ao_1e_int_kinetic_im (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(in) :: dset(*)
   end function trexio_write_ao_1e_int_kinetic_im
end interface

interface
   integer(trexio_exit_code) function trexio_write_ao_1e_int_potential_n_e_im (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(in) :: dset(*)
   end function trexio_write_ao_1e_int_potential_n_e_im
end interface

interface
   integer(trexio_exit_code) function trexio_write_ao_1e_int_ecp_im (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(in) :: dset(*)
   end function trexio_write_ao_1e_int_ecp_im
end interface

interface
   integer(trexio_exit_code) function trexio_write_ao_1e_int_core_hamiltonian_im (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(in) :: dset(*)
   end function trexio_write_ao_1e_int_core_hamiltonian_im
end interface

interface
   integer(trexio_exit_code) function trexio_write_mo_coefficient (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(in) :: dset(*)
   end function trexio_write_mo_coefficient
end interface

interface
   integer(trexio_exit_code) function trexio_write_mo_coefficient_im (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(in) :: dset(*)
   end function trexio_write_mo_coefficient_im
end interface

interface
   integer(trexio_exit_code) function trexio_write_mo_occupation (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(in) :: dset(*)
   end function trexio_write_mo_occupation
end interface

interface
   integer(trexio_exit_code) function trexio_write_mo_energy (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(in) :: dset(*)
   end function trexio_write_mo_energy
end interface

interface
   integer(trexio_exit_code) function trexio_write_mo_spin (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int32_t), intent(in) :: dset(*)
   end function trexio_write_mo_spin
end interface

interface
   integer(trexio_exit_code) function trexio_write_mo_1e_int_overlap (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(in) :: dset(*)
   end function trexio_write_mo_1e_int_overlap
end interface

interface
   integer(trexio_exit_code) function trexio_write_mo_1e_int_kinetic (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(in) :: dset(*)
   end function trexio_write_mo_1e_int_kinetic
end interface

interface
   integer(trexio_exit_code) function trexio_write_mo_1e_int_potential_n_e (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(in) :: dset(*)
   end function trexio_write_mo_1e_int_potential_n_e
end interface

interface
   integer(trexio_exit_code) function trexio_write_mo_1e_int_ecp (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(in) :: dset(*)
   end function trexio_write_mo_1e_int_ecp
end interface

interface
   integer(trexio_exit_code) function trexio_write_mo_1e_int_core_hamiltonian (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(in) :: dset(*)
   end function trexio_write_mo_1e_int_core_hamiltonian
end interface

interface
   integer(trexio_exit_code) function trexio_write_mo_1e_int_overlap_im (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(in) :: dset(*)
   end function trexio_write_mo_1e_int_overlap_im
end interface

interface
   integer(trexio_exit_code) function trexio_write_mo_1e_int_kinetic_im (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(in) :: dset(*)
   end function trexio_write_mo_1e_int_kinetic_im
end interface

interface
   integer(trexio_exit_code) function trexio_write_mo_1e_int_potential_n_e_im (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(in) :: dset(*)
   end function trexio_write_mo_1e_int_potential_n_e_im
end interface

interface
   integer(trexio_exit_code) function trexio_write_mo_1e_int_ecp_im (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(in) :: dset(*)
   end function trexio_write_mo_1e_int_ecp_im
end interface

interface
   integer(trexio_exit_code) function trexio_write_mo_1e_int_core_hamiltonian_im (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(in) :: dset(*)
   end function trexio_write_mo_1e_int_core_hamiltonian_im
end interface

interface
   integer(trexio_exit_code) function trexio_write_rdm_1e (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(in) :: dset(*)
   end function trexio_write_rdm_1e
end interface

interface
   integer(trexio_exit_code) function trexio_write_rdm_1e_up (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(in) :: dset(*)
   end function trexio_write_rdm_1e_up
end interface

interface
   integer(trexio_exit_code) function trexio_write_rdm_1e_dn (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(in) :: dset(*)
   end function trexio_write_rdm_1e_dn
end interface

interface
   integer(trexio_exit_code) function trexio_write_cell_a (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(in) :: dset(*)
   end function trexio_write_cell_a
end interface

interface
   integer(trexio_exit_code) function trexio_write_cell_b (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(in) :: dset(*)
   end function trexio_write_cell_b
end interface

interface
   integer(trexio_exit_code) function trexio_write_cell_c (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(in) :: dset(*)
   end function trexio_write_cell_c
end interface

interface
   integer(trexio_exit_code) function trexio_write_pbc_k_point (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(in) :: dset(*)
   end function trexio_write_pbc_k_point
end interface

interface
   integer(trexio_exit_code) function trexio_write_qmc_point (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(in) :: dset(*)
   end function trexio_write_qmc_point
end interface

interface
   integer(trexio_exit_code) function trexio_write_qmc_psi (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(in) :: dset(*)
   end function trexio_write_qmc_psi
end interface

interface
   integer(trexio_exit_code) function trexio_write_qmc_e_loc (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     real(c_double), intent(in) :: dset(*)
   end function trexio_write_qmc_e_loc
end interface

interface
   integer(trexio_exit_code) function trexio_write_ao_2e_int_eri (trex_file, &
                                               offset_file, buffer_size, &
                                               index_sparse, value_sparse) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int64_t), intent(in), value :: offset_file
     integer(c_int64_t), intent(in), value :: buffer_size
     integer(c_int32_t), intent(in)        :: index_sparse(*)
     real(c_double), intent(in)            :: value_sparse(*)
   end function trexio_write_ao_2e_int_eri
end interface

interface
   integer(trexio_exit_code) function trexio_write_safe_ao_2e_int_eri (trex_file, &
                                                    offset_file, buffer_size, &
                                                    index_sparse, index_size, &
                                                    value_sparse, value_size) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int64_t), intent(in), value :: offset_file
     integer(c_int64_t), intent(in), value :: buffer_size
     integer(c_int32_t), intent(in)        :: index_sparse(*)
     integer(c_int64_t), intent(in), value :: index_size
     real(c_double), intent(in)            :: value_sparse(*)
     integer(c_int64_t), intent(in), value :: value_size
   end function trexio_write_safe_ao_2e_int_eri
end interface

interface
   integer(trexio_exit_code) function trexio_write_ao_2e_int_eri_lr (trex_file, &
                                               offset_file, buffer_size, &
                                               index_sparse, value_sparse) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int64_t), intent(in), value :: offset_file
     integer(c_int64_t), intent(in), value :: buffer_size
     integer(c_int32_t), intent(in)        :: index_sparse(*)
     real(c_double), intent(in)            :: value_sparse(*)
   end function trexio_write_ao_2e_int_eri_lr
end interface

interface
   integer(trexio_exit_code) function trexio_write_safe_ao_2e_int_eri_lr (trex_file, &
                                                    offset_file, buffer_size, &
                                                    index_sparse, index_size, &
                                                    value_sparse, value_size) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int64_t), intent(in), value :: offset_file
     integer(c_int64_t), intent(in), value :: buffer_size
     integer(c_int32_t), intent(in)        :: index_sparse(*)
     integer(c_int64_t), intent(in), value :: index_size
     real(c_double), intent(in)            :: value_sparse(*)
     integer(c_int64_t), intent(in), value :: value_size
   end function trexio_write_safe_ao_2e_int_eri_lr
end interface

interface
   integer(trexio_exit_code) function trexio_write_mo_2e_int_eri (trex_file, &
                                               offset_file, buffer_size, &
                                               index_sparse, value_sparse) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int64_t), intent(in), value :: offset_file
     integer(c_int64_t), intent(in), value :: buffer_size
     integer(c_int32_t), intent(in)        :: index_sparse(*)
     real(c_double), intent(in)            :: value_sparse(*)
   end function trexio_write_mo_2e_int_eri
end interface

interface
   integer(trexio_exit_code) function trexio_write_safe_mo_2e_int_eri (trex_file, &
                                                    offset_file, buffer_size, &
                                                    index_sparse, index_size, &
                                                    value_sparse, value_size) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int64_t), intent(in), value :: offset_file
     integer(c_int64_t), intent(in), value :: buffer_size
     integer(c_int32_t), intent(in)        :: index_sparse(*)
     integer(c_int64_t), intent(in), value :: index_size
     real(c_double), intent(in)            :: value_sparse(*)
     integer(c_int64_t), intent(in), value :: value_size
   end function trexio_write_safe_mo_2e_int_eri
end interface

interface
   integer(trexio_exit_code) function trexio_write_mo_2e_int_eri_lr (trex_file, &
                                               offset_file, buffer_size, &
                                               index_sparse, value_sparse) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int64_t), intent(in), value :: offset_file
     integer(c_int64_t), intent(in), value :: buffer_size
     integer(c_int32_t), intent(in)        :: index_sparse(*)
     real(c_double), intent(in)            :: value_sparse(*)
   end function trexio_write_mo_2e_int_eri_lr
end interface

interface
   integer(trexio_exit_code) function trexio_write_safe_mo_2e_int_eri_lr (trex_file, &
                                                    offset_file, buffer_size, &
                                                    index_sparse, index_size, &
                                                    value_sparse, value_size) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int64_t), intent(in), value :: offset_file
     integer(c_int64_t), intent(in), value :: buffer_size
     integer(c_int32_t), intent(in)        :: index_sparse(*)
     integer(c_int64_t), intent(in), value :: index_size
     real(c_double), intent(in)            :: value_sparse(*)
     integer(c_int64_t), intent(in), value :: value_size
   end function trexio_write_safe_mo_2e_int_eri_lr
end interface

interface
   integer(trexio_exit_code) function trexio_write_csf_det_coefficient (trex_file, &
                                               offset_file, buffer_size, &
                                               index_sparse, value_sparse) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int64_t), intent(in), value :: offset_file
     integer(c_int64_t), intent(in), value :: buffer_size
     integer(c_int32_t), intent(in)        :: index_sparse(*)
     real(c_double), intent(in)            :: value_sparse(*)
   end function trexio_write_csf_det_coefficient
end interface

interface
   integer(trexio_exit_code) function trexio_write_safe_csf_det_coefficient (trex_file, &
                                                    offset_file, buffer_size, &
                                                    index_sparse, index_size, &
                                                    value_sparse, value_size) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int64_t), intent(in), value :: offset_file
     integer(c_int64_t), intent(in), value :: buffer_size
     integer(c_int32_t), intent(in)        :: index_sparse(*)
     integer(c_int64_t), intent(in), value :: index_size
     real(c_double), intent(in)            :: value_sparse(*)
     integer(c_int64_t), intent(in), value :: value_size
   end function trexio_write_safe_csf_det_coefficient
end interface

interface
   integer(trexio_exit_code) function trexio_write_rdm_2e (trex_file, &
                                               offset_file, buffer_size, &
                                               index_sparse, value_sparse) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int64_t), intent(in), value :: offset_file
     integer(c_int64_t), intent(in), value :: buffer_size
     integer(c_int32_t), intent(in)        :: index_sparse(*)
     real(c_double), intent(in)            :: value_sparse(*)
   end function trexio_write_rdm_2e
end interface

interface
   integer(trexio_exit_code) function trexio_write_safe_rdm_2e (trex_file, &
                                                    offset_file, buffer_size, &
                                                    index_sparse, index_size, &
                                                    value_sparse, value_size) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int64_t), intent(in), value :: offset_file
     integer(c_int64_t), intent(in), value :: buffer_size
     integer(c_int32_t), intent(in)        :: index_sparse(*)
     integer(c_int64_t), intent(in), value :: index_size
     real(c_double), intent(in)            :: value_sparse(*)
     integer(c_int64_t), intent(in), value :: value_size
   end function trexio_write_safe_rdm_2e
end interface

interface
   integer(trexio_exit_code) function trexio_write_rdm_2e_upup (trex_file, &
                                               offset_file, buffer_size, &
                                               index_sparse, value_sparse) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int64_t), intent(in), value :: offset_file
     integer(c_int64_t), intent(in), value :: buffer_size
     integer(c_int32_t), intent(in)        :: index_sparse(*)
     real(c_double), intent(in)            :: value_sparse(*)
   end function trexio_write_rdm_2e_upup
end interface

interface
   integer(trexio_exit_code) function trexio_write_safe_rdm_2e_upup (trex_file, &
                                                    offset_file, buffer_size, &
                                                    index_sparse, index_size, &
                                                    value_sparse, value_size) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int64_t), intent(in), value :: offset_file
     integer(c_int64_t), intent(in), value :: buffer_size
     integer(c_int32_t), intent(in)        :: index_sparse(*)
     integer(c_int64_t), intent(in), value :: index_size
     real(c_double), intent(in)            :: value_sparse(*)
     integer(c_int64_t), intent(in), value :: value_size
   end function trexio_write_safe_rdm_2e_upup
end interface

interface
   integer(trexio_exit_code) function trexio_write_rdm_2e_dndn (trex_file, &
                                               offset_file, buffer_size, &
                                               index_sparse, value_sparse) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int64_t), intent(in), value :: offset_file
     integer(c_int64_t), intent(in), value :: buffer_size
     integer(c_int32_t), intent(in)        :: index_sparse(*)
     real(c_double), intent(in)            :: value_sparse(*)
   end function trexio_write_rdm_2e_dndn
end interface

interface
   integer(trexio_exit_code) function trexio_write_safe_rdm_2e_dndn (trex_file, &
                                                    offset_file, buffer_size, &
                                                    index_sparse, index_size, &
                                                    value_sparse, value_size) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int64_t), intent(in), value :: offset_file
     integer(c_int64_t), intent(in), value :: buffer_size
     integer(c_int32_t), intent(in)        :: index_sparse(*)
     integer(c_int64_t), intent(in), value :: index_size
     real(c_double), intent(in)            :: value_sparse(*)
     integer(c_int64_t), intent(in), value :: value_size
   end function trexio_write_safe_rdm_2e_dndn
end interface

interface
   integer(trexio_exit_code) function trexio_write_rdm_2e_updn (trex_file, &
                                               offset_file, buffer_size, &
                                               index_sparse, value_sparse) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int64_t), intent(in), value :: offset_file
     integer(c_int64_t), intent(in), value :: buffer_size
     integer(c_int32_t), intent(in)        :: index_sparse(*)
     real(c_double), intent(in)            :: value_sparse(*)
   end function trexio_write_rdm_2e_updn
end interface

interface
   integer(trexio_exit_code) function trexio_write_safe_rdm_2e_updn (trex_file, &
                                                    offset_file, buffer_size, &
                                                    index_sparse, index_size, &
                                                    value_sparse, value_size) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int64_t), intent(in), value :: offset_file
     integer(c_int64_t), intent(in), value :: buffer_size
     integer(c_int32_t), intent(in)        :: index_sparse(*)
     integer(c_int64_t), intent(in), value :: index_size
     real(c_double), intent(in)            :: value_sparse(*)
     integer(c_int64_t), intent(in), value :: value_size
   end function trexio_write_safe_rdm_2e_updn
end interface

interface
   integer(trexio_exit_code) function trexio_write_rdm_2e_dnup (trex_file, &
                                               offset_file, buffer_size, &
                                               index_sparse, value_sparse) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int64_t), intent(in), value :: offset_file
     integer(c_int64_t), intent(in), value :: buffer_size
     integer(c_int32_t), intent(in)        :: index_sparse(*)
     real(c_double), intent(in)            :: value_sparse(*)
   end function trexio_write_rdm_2e_dnup
end interface

interface
   integer(trexio_exit_code) function trexio_write_safe_rdm_2e_dnup (trex_file, &
                                                    offset_file, buffer_size, &
                                                    index_sparse, index_size, &
                                                    value_sparse, value_size) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int64_t), intent(in), value :: offset_file
     integer(c_int64_t), intent(in), value :: buffer_size
     integer(c_int32_t), intent(in)        :: index_sparse(*)
     integer(c_int64_t), intent(in), value :: index_size
     real(c_double), intent(in)            :: value_sparse(*)
     integer(c_int64_t), intent(in), value :: value_size
   end function trexio_write_safe_rdm_2e_dnup
end interface

interface
   integer(trexio_exit_code) function trexio_write_metadata_code_low (trex_file, dset, max_str_len) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     character(kind=c_char), intent(in)    :: dset(*)
     integer(c_int32_t), intent(in), value :: max_str_len
   end function trexio_write_metadata_code_low
end interface

interface
   integer(trexio_exit_code) function trexio_write_metadata_author_low (trex_file, dset, max_str_len) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     character(kind=c_char), intent(in)    :: dset(*)
     integer(c_int32_t), intent(in), value :: max_str_len
   end function trexio_write_metadata_author_low
end interface

interface
   integer(trexio_exit_code) function trexio_write_nucleus_label_low (trex_file, dset, max_str_len) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     character(kind=c_char), intent(in)    :: dset(*)
     integer(c_int32_t), intent(in), value :: max_str_len
   end function trexio_write_nucleus_label_low
end interface

interface
   integer(trexio_exit_code) function trexio_write_mo_class_low (trex_file, dset, max_str_len) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     character(kind=c_char), intent(in)    :: dset(*)
     integer(c_int32_t), intent(in), value :: max_str_len
   end function trexio_write_mo_class_low
end interface

interface
   integer(trexio_exit_code) function trexio_write_mo_symmetry_low (trex_file, dset, max_str_len) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     character(kind=c_char), intent(in)    :: dset(*)
     integer(c_int32_t), intent(in), value :: max_str_len
   end function trexio_write_mo_symmetry_low
end interface

interface
   integer(trexio_exit_code) function trexio_write_state_label_low (trex_file, dset, max_str_len) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     character(kind=c_char), intent(in)    :: dset(*)
     integer(c_int32_t), intent(in), value :: max_str_len
   end function trexio_write_state_label_low
end interface

interface
   integer(trexio_exit_code) function trexio_write_state_file_name_low (trex_file, dset, max_str_len) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     character(kind=c_char), intent(in)    :: dset(*)
     integer(c_int32_t), intent(in), value :: max_str_len
   end function trexio_write_state_file_name_low
end interface

interface
   integer(trexio_exit_code) function trexio_has_determinant_list (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
   end function trexio_has_determinant_list
end interface

interface
   integer(trexio_exit_code) function trexio_get_int64_num (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int32_t), intent(out) :: num
   end function trexio_get_int64_num
end interface
interface
   integer(trexio_exit_code) function trexio_read_determinant_list(trex_file, &
                                              offset_file, buffer_size, list) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int64_t), intent(in), value :: offset_file
     integer(c_int64_t), intent(inout)     :: buffer_size
     integer(c_int64_t), intent(out)       :: list(*)
   end function trexio_read_determinant_list
end interface

interface
   integer(trexio_exit_code) function trexio_read_safe_determinant_list (trex_file, &
                                                   offset_file, buffer_size, &
                                                   list, list_size) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int64_t), intent(in), value :: offset_file
     integer(c_int64_t), intent(inout)     :: buffer_size
     integer(c_int64_t), intent(out)       :: list(*)
     integer(c_int64_t), intent(in), value :: list_size
   end function trexio_read_safe_determinant_list
end interface
interface
   integer(trexio_exit_code) function trexio_write_determinant_list (trex_file, &
                                               offset_file, buffer_size, list) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int64_t), intent(in), value :: offset_file
     integer(c_int64_t), intent(in), value :: buffer_size
     integer(c_int64_t), intent(in)        :: list(*)
   end function trexio_write_determinant_list
end interface

interface
   integer(trexio_exit_code) function trexio_write_safe_determinant_list (trex_file, &
                                                    offset_file, buffer_size, &
                                                    list, list_size) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer(trexio_t), intent(in), value :: trex_file
     integer(c_int64_t), intent(in), value :: offset_file
     integer(c_int64_t), intent(in), value :: buffer_size
     integer(c_int64_t), intent(in)        :: list(*)
     integer(c_int64_t), intent(in), value :: list_size
   end function trexio_write_safe_determinant_list
end interface
contains
   integer(trexio_t) function trexio_open (filename, mode, back_end, rc_open)
     use, intrinsic :: iso_c_binding, only : c_null_char
     implicit none
     character(len=*), intent(in)                    :: filename
     character, intent(in), value                    :: mode
     integer(trexio_back_end_t), intent(in), value   :: back_end
     integer(trexio_exit_code), intent(out)          :: rc_open
     character(len=len_trim(filename)+1)             :: filename_c
     integer(trexio_exit_code) :: rc

     filename_c = trim(filename) // c_null_char
     trexio_open = trexio_open_c(filename_c, mode, back_end, rc_open)
     if (trexio_open == 0_8 .or. rc_open /= TREXIO_SUCCESS) then
       return
     endif
     rc = trexio_set_one_based(trexio_open)
     if (rc /= TREXIO_SUCCESS) then
        rc = trexio_close(trexio_open)
        trexio_open = 0_8
     endif
   end function trexio_open

integer(trexio_exit_code) function trexio_inquire (filename)
  use, intrinsic :: iso_c_binding
  implicit none
  character(len=*), intent(in)        :: filename
  character(len=len_trim(filename)+1) :: filename_c

  filename_c = trim(filename) // c_null_char
  trexio_inquire = trexio_inquire_c(filename_c)
end function trexio_inquire

integer(trexio_exit_code) function trexio_to_bitfield_list(list, occupied_num, det_list, N_int)
  use, intrinsic :: iso_c_binding
  implicit none

  integer(c_int32_t), intent(in)        :: list(*)
  integer(c_int32_t), intent(in), value :: occupied_num
  integer(c_int64_t), intent(out)       :: det_list(*)
  integer(c_int32_t), intent(in), value :: N_int
  integer(c_int32_t)                    :: list_0based(occupied_num)

  integer :: i
  do i = 1,occupied_num
    list_0based(i) = list(i) - 1
  enddo

  trexio_to_bitfield_list = trexio_to_bitfield_list_c(list_0based, occupied_num, det_list, N_int)
  if (trexio_to_bitfield_list /= TREXIO_SUCCESS) then
    return
  endif

end function trexio_to_bitfield_list


integer(trexio_exit_code) function trexio_to_orbital_list(N_int, d1, list, occupied_num)
  use, intrinsic :: iso_c_binding
  implicit none

  integer(c_int32_t), intent(in), value :: N_int
  integer(c_int64_t), intent(in)        :: d1(*)
  integer(c_int32_t), intent(out)       :: list(*)
  integer(c_int32_t), intent(out)       :: occupied_num

  integer :: i

  trexio_to_orbital_list = trexio_to_orbital_list_c(N_int, d1, list, occupied_num)
  if (trexio_to_orbital_list /= TREXIO_SUCCESS) then
    return
  endif

  do i = 1,occupied_num
    list(i) = list(i) + 1
  enddo
end function trexio_to_orbital_list


integer(trexio_exit_code) function trexio_to_orbital_list_up_dn(N_int, d1, list_up, list_dn, occ_num_up, occ_num_dn)
  use, intrinsic :: iso_c_binding
  implicit none
  integer(c_int32_t), intent(in), value :: N_int
  integer(c_int64_t), intent(in)        :: d1(*)
  integer(c_int32_t), intent(out)       :: list_up(*)
  integer(c_int32_t), intent(out)       :: list_dn(*)
  integer(c_int32_t), intent(out)       :: occ_num_up
  integer(c_int32_t), intent(out)       :: occ_num_dn

  integer :: i

  trexio_to_orbital_list_up_dn = trexio_to_orbital_list_up_dn_c(N_int, d1, list_up, list_dn, occ_num_up, occ_num_dn)
  if (trexio_to_orbital_list_up_dn /= TREXIO_SUCCESS) then
    return
  endif

  do i = 1,occ_num_up
    list_up(i) = list_up(i) + 1
  enddo
  do i = 1,occ_num_dn
    list_dn(i) = list_dn(i) + 1
  enddo
end function trexio_to_orbital_list_up_dn

subroutine trexio_strarray2str(str_array, max_num_str, max_len_str, str_res)
  use, intrinsic :: iso_c_binding, only : c_null_char
  implicit none

  integer(c_int64_t), intent(in), value   :: max_num_str  ! number of elements in strign array
  integer, intent(in), value   :: max_len_str  ! maximum length of a string in an array
  character(len=*), intent(in)  :: str_array(*)
  character(len=:), allocatable, intent(out) :: str_res
  integer(c_int64_t) :: i

  str_res = ''
  do i = 1, max_num_str
    str_res = str_res // trim(str_array(i)) // TREXIO_DELIM
  enddo
  str_res = str_res // c_null_char

end subroutine trexio_strarray2str

subroutine trexio_str2strarray(str_flat, max_num_str, max_len_str, str_array)
  implicit none

  integer(c_int64_t), intent(in), value   :: max_num_str  ! number of elements in strign array
  integer, intent(in), value              :: max_len_str  ! maximum length of a string in an array
  character(kind=c_char), intent(in)      :: str_flat(*)
  character(len=*), intent(inout)         :: str_array(*)

  character(len=max_len_str)  :: tmp_str
  integer(c_int64_t) :: i, j, k, ind, len_flat

  len_flat = (max_len_str+1)*max_num_str + 1

  ind=1
  do i=1,max_num_str
    k = 1
    tmp_str=''
    do j=ind,len_flat
      if (str_flat(j) == TREXIO_DELIM) then
        ind=j+1
        exit
      endif
      tmp_str(k:k) = str_flat(j)
      k = k + 1
    enddo
    str_array(i)=tmp_str
  enddo

end subroutine trexio_str2strarray

subroutine trexio_assert(trexio_rc, check_rc, success_message)
  implicit none

  integer, intent(in), value :: trexio_rc
  integer, intent(in), value :: check_rc
  character(len=*), intent(in), optional  :: success_message

  character*(128) :: str

  if (trexio_rc == check_rc) then
    if (present(success_message)) write(*,*) success_message
  else
    call trexio_string_of_error(trexio_rc, str)
    print *, trim(str)
    stop 1
  endif

end subroutine trexio_assert
integer(trexio_exit_code) function trexio_read_metadata_package_version (trex_file, str, max_str_len)
  implicit none
  integer(trexio_t), intent(in), value  :: trex_file
  integer(c_int32_t), intent(in), value :: max_str_len
  character, intent(out) :: str(*)

  trexio_read_metadata_package_version = trexio_read_metadata_package_version_c(trex_file, str, max_str_len)

end function trexio_read_metadata_package_version

integer(trexio_exit_code) function trexio_read_metadata_description (trex_file, str, max_str_len)
  implicit none
  integer(trexio_t), intent(in), value  :: trex_file
  integer(c_int32_t), intent(in), value :: max_str_len
  character, intent(out) :: str(*)

  trexio_read_metadata_description = trexio_read_metadata_description_c(trex_file, str, max_str_len)

end function trexio_read_metadata_description

integer(trexio_exit_code) function trexio_read_nucleus_point_group (trex_file, str, max_str_len)
  implicit none
  integer(trexio_t), intent(in), value  :: trex_file
  integer(c_int32_t), intent(in), value :: max_str_len
  character, intent(out) :: str(*)

  trexio_read_nucleus_point_group = trexio_read_nucleus_point_group_c(trex_file, str, max_str_len)

end function trexio_read_nucleus_point_group

integer(trexio_exit_code) function trexio_read_basis_type (trex_file, str, max_str_len)
  implicit none
  integer(trexio_t), intent(in), value  :: trex_file
  integer(c_int32_t), intent(in), value :: max_str_len
  character, intent(out) :: str(*)

  trexio_read_basis_type = trexio_read_basis_type_c(trex_file, str, max_str_len)

end function trexio_read_basis_type

integer(trexio_exit_code) function trexio_read_mo_type (trex_file, str, max_str_len)
  implicit none
  integer(trexio_t), intent(in), value  :: trex_file
  integer(c_int32_t), intent(in), value :: max_str_len
  character, intent(out) :: str(*)

  trexio_read_mo_type = trexio_read_mo_type_c(trex_file, str, max_str_len)

end function trexio_read_mo_type

integer(trexio_exit_code) function trexio_read_state_current_label (trex_file, str, max_str_len)
  implicit none
  integer(trexio_t), intent(in), value  :: trex_file
  integer(c_int32_t), intent(in), value :: max_str_len
  character, intent(out) :: str(*)

  trexio_read_state_current_label = trexio_read_state_current_label_c(trex_file, str, max_str_len)

end function trexio_read_state_current_label

integer(trexio_exit_code) function trexio_read_metadata_code (trex_file, dset, max_str_len)
  implicit none
  integer(trexio_t), intent(in), value :: trex_file
  integer(c_int32_t), intent(in), value :: max_str_len
  character(len=*), intent(inout) :: dset(*)

  character, allocatable :: str_compiled(:)
  integer(c_int64_t) :: metadata_code_num
  integer(trexio_exit_code) :: rc

  rc = trexio_read_metadata_code_num_64(trex_file, metadata_code_num)
  if (rc /= TREXIO_SUCCESS) trexio_read_metadata_code = rc

  allocate(str_compiled(metadata_code_num*(max_str_len+1)+1))

  rc = trexio_read_metadata_code_low(trex_file, str_compiled, max_str_len)
  if (rc /= TREXIO_SUCCESS) then
    deallocate(str_compiled)
    trexio_read_metadata_code = rc
  else
    call trexio_str2strarray(str_compiled, metadata_code_num, max_str_len, dset)
    deallocate(str_compiled)
    trexio_read_metadata_code = TREXIO_SUCCESS
  endif

end function trexio_read_metadata_code

integer(trexio_exit_code) function trexio_read_metadata_author (trex_file, dset, max_str_len)
  implicit none
  integer(trexio_t), intent(in), value :: trex_file
  integer(c_int32_t), intent(in), value :: max_str_len
  character(len=*), intent(inout) :: dset(*)

  character, allocatable :: str_compiled(:)
  integer(c_int64_t) :: metadata_author_num
  integer(trexio_exit_code) :: rc

  rc = trexio_read_metadata_author_num_64(trex_file, metadata_author_num)
  if (rc /= TREXIO_SUCCESS) trexio_read_metadata_author = rc

  allocate(str_compiled(metadata_author_num*(max_str_len+1)+1))

  rc = trexio_read_metadata_author_low(trex_file, str_compiled, max_str_len)
  if (rc /= TREXIO_SUCCESS) then
    deallocate(str_compiled)
    trexio_read_metadata_author = rc
  else
    call trexio_str2strarray(str_compiled, metadata_author_num, max_str_len, dset)
    deallocate(str_compiled)
    trexio_read_metadata_author = TREXIO_SUCCESS
  endif

end function trexio_read_metadata_author

integer(trexio_exit_code) function trexio_read_nucleus_label (trex_file, dset, max_str_len)
  implicit none
  integer(trexio_t), intent(in), value :: trex_file
  integer(c_int32_t), intent(in), value :: max_str_len
  character(len=*), intent(inout) :: dset(*)

  character, allocatable :: str_compiled(:)
  integer(c_int64_t) :: nucleus_num
  integer(trexio_exit_code) :: rc

  rc = trexio_read_nucleus_num_64(trex_file, nucleus_num)
  if (rc /= TREXIO_SUCCESS) trexio_read_nucleus_label = rc

  allocate(str_compiled(nucleus_num*(max_str_len+1)+1))

  rc = trexio_read_nucleus_label_low(trex_file, str_compiled, max_str_len)
  if (rc /= TREXIO_SUCCESS) then
    deallocate(str_compiled)
    trexio_read_nucleus_label = rc
  else
    call trexio_str2strarray(str_compiled, nucleus_num, max_str_len, dset)
    deallocate(str_compiled)
    trexio_read_nucleus_label = TREXIO_SUCCESS
  endif

end function trexio_read_nucleus_label

integer(trexio_exit_code) function trexio_read_mo_class (trex_file, dset, max_str_len)
  implicit none
  integer(trexio_t), intent(in), value :: trex_file
  integer(c_int32_t), intent(in), value :: max_str_len
  character(len=*), intent(inout) :: dset(*)

  character, allocatable :: str_compiled(:)
  integer(c_int64_t) :: mo_num
  integer(trexio_exit_code) :: rc

  rc = trexio_read_mo_num_64(trex_file, mo_num)
  if (rc /= TREXIO_SUCCESS) trexio_read_mo_class = rc

  allocate(str_compiled(mo_num*(max_str_len+1)+1))

  rc = trexio_read_mo_class_low(trex_file, str_compiled, max_str_len)
  if (rc /= TREXIO_SUCCESS) then
    deallocate(str_compiled)
    trexio_read_mo_class = rc
  else
    call trexio_str2strarray(str_compiled, mo_num, max_str_len, dset)
    deallocate(str_compiled)
    trexio_read_mo_class = TREXIO_SUCCESS
  endif

end function trexio_read_mo_class

integer(trexio_exit_code) function trexio_read_mo_symmetry (trex_file, dset, max_str_len)
  implicit none
  integer(trexio_t), intent(in), value :: trex_file
  integer(c_int32_t), intent(in), value :: max_str_len
  character(len=*), intent(inout) :: dset(*)

  character, allocatable :: str_compiled(:)
  integer(c_int64_t) :: mo_num
  integer(trexio_exit_code) :: rc

  rc = trexio_read_mo_num_64(trex_file, mo_num)
  if (rc /= TREXIO_SUCCESS) trexio_read_mo_symmetry = rc

  allocate(str_compiled(mo_num*(max_str_len+1)+1))

  rc = trexio_read_mo_symmetry_low(trex_file, str_compiled, max_str_len)
  if (rc /= TREXIO_SUCCESS) then
    deallocate(str_compiled)
    trexio_read_mo_symmetry = rc
  else
    call trexio_str2strarray(str_compiled, mo_num, max_str_len, dset)
    deallocate(str_compiled)
    trexio_read_mo_symmetry = TREXIO_SUCCESS
  endif

end function trexio_read_mo_symmetry

integer(trexio_exit_code) function trexio_read_state_label (trex_file, dset, max_str_len)
  implicit none
  integer(trexio_t), intent(in), value :: trex_file
  integer(c_int32_t), intent(in), value :: max_str_len
  character(len=*), intent(inout) :: dset(*)

  character, allocatable :: str_compiled(:)
  integer(c_int64_t) :: state_num
  integer(trexio_exit_code) :: rc

  rc = trexio_read_state_num_64(trex_file, state_num)
  if (rc /= TREXIO_SUCCESS) trexio_read_state_label = rc

  allocate(str_compiled(state_num*(max_str_len+1)+1))

  rc = trexio_read_state_label_low(trex_file, str_compiled, max_str_len)
  if (rc /= TREXIO_SUCCESS) then
    deallocate(str_compiled)
    trexio_read_state_label = rc
  else
    call trexio_str2strarray(str_compiled, state_num, max_str_len, dset)
    deallocate(str_compiled)
    trexio_read_state_label = TREXIO_SUCCESS
  endif

end function trexio_read_state_label

integer(trexio_exit_code) function trexio_read_state_file_name (trex_file, dset, max_str_len)
  implicit none
  integer(trexio_t), intent(in), value :: trex_file
  integer(c_int32_t), intent(in), value :: max_str_len
  character(len=*), intent(inout) :: dset(*)

  character, allocatable :: str_compiled(:)
  integer(c_int64_t) :: state_num
  integer(trexio_exit_code) :: rc

  rc = trexio_read_state_num_64(trex_file, state_num)
  if (rc /= TREXIO_SUCCESS) trexio_read_state_file_name = rc

  allocate(str_compiled(state_num*(max_str_len+1)+1))

  rc = trexio_read_state_file_name_low(trex_file, str_compiled, max_str_len)
  if (rc /= TREXIO_SUCCESS) then
    deallocate(str_compiled)
    trexio_read_state_file_name = rc
  else
    call trexio_str2strarray(str_compiled, state_num, max_str_len, dset)
    deallocate(str_compiled)
    trexio_read_state_file_name = TREXIO_SUCCESS
  endif

end function trexio_read_state_file_name

integer(trexio_exit_code) function trexio_write_metadata_package_version (trex_file, str, max_str_len)
  use, intrinsic :: iso_c_binding, only : c_null_char
  implicit none
  integer(trexio_t), intent(in), value  :: trex_file
  integer(c_int32_t), intent(in), value :: max_str_len
  character(len=*), intent(in) :: str

  character(len=len_trim(str)+1) :: str_c

  str_c = trim(str) // c_null_char

  trexio_write_metadata_package_version = trexio_write_metadata_package_version_c(trex_file, str_c, max_str_len)

end function trexio_write_metadata_package_version

integer(trexio_exit_code) function trexio_write_metadata_description (trex_file, str, max_str_len)
  use, intrinsic :: iso_c_binding, only : c_null_char
  implicit none
  integer(trexio_t), intent(in), value  :: trex_file
  integer(c_int32_t), intent(in), value :: max_str_len
  character(len=*), intent(in) :: str

  character(len=len_trim(str)+1) :: str_c

  str_c = trim(str) // c_null_char

  trexio_write_metadata_description = trexio_write_metadata_description_c(trex_file, str_c, max_str_len)

end function trexio_write_metadata_description

integer(trexio_exit_code) function trexio_write_nucleus_point_group (trex_file, str, max_str_len)
  use, intrinsic :: iso_c_binding, only : c_null_char
  implicit none
  integer(trexio_t), intent(in), value  :: trex_file
  integer(c_int32_t), intent(in), value :: max_str_len
  character(len=*), intent(in) :: str

  character(len=len_trim(str)+1) :: str_c

  str_c = trim(str) // c_null_char

  trexio_write_nucleus_point_group = trexio_write_nucleus_point_group_c(trex_file, str_c, max_str_len)

end function trexio_write_nucleus_point_group

integer(trexio_exit_code) function trexio_write_basis_type (trex_file, str, max_str_len)
  use, intrinsic :: iso_c_binding, only : c_null_char
  implicit none
  integer(trexio_t), intent(in), value  :: trex_file
  integer(c_int32_t), intent(in), value :: max_str_len
  character(len=*), intent(in) :: str

  character(len=len_trim(str)+1) :: str_c

  str_c = trim(str) // c_null_char

  trexio_write_basis_type = trexio_write_basis_type_c(trex_file, str_c, max_str_len)

end function trexio_write_basis_type

integer(trexio_exit_code) function trexio_write_mo_type (trex_file, str, max_str_len)
  use, intrinsic :: iso_c_binding, only : c_null_char
  implicit none
  integer(trexio_t), intent(in), value  :: trex_file
  integer(c_int32_t), intent(in), value :: max_str_len
  character(len=*), intent(in) :: str

  character(len=len_trim(str)+1) :: str_c

  str_c = trim(str) // c_null_char

  trexio_write_mo_type = trexio_write_mo_type_c(trex_file, str_c, max_str_len)

end function trexio_write_mo_type

integer(trexio_exit_code) function trexio_write_state_current_label (trex_file, str, max_str_len)
  use, intrinsic :: iso_c_binding, only : c_null_char
  implicit none
  integer(trexio_t), intent(in), value  :: trex_file
  integer(c_int32_t), intent(in), value :: max_str_len
  character(len=*), intent(in) :: str

  character(len=len_trim(str)+1) :: str_c

  str_c = trim(str) // c_null_char

  trexio_write_state_current_label = trexio_write_state_current_label_c(trex_file, str_c, max_str_len)

end function trexio_write_state_current_label

integer(trexio_exit_code) function trexio_write_metadata_code (trex_file, dset, max_str_len)
  implicit none
  integer(trexio_t), intent(in), value :: trex_file
  integer(c_int32_t), intent(in), value :: max_str_len
  character(len=*), intent(in) :: dset(*)

  character(len=:), allocatable :: str_compiled
  integer(c_int64_t) :: metadata_code_num
  integer(trexio_exit_code) :: rc

  rc = trexio_read_metadata_code_num_64(trex_file, metadata_code_num)
  if (rc /= TREXIO_SUCCESS) then
    trexio_write_metadata_code = rc
  else
    call trexio_strarray2str(dset, metadata_code_num, max_str_len, str_compiled)
    trexio_write_metadata_code = trexio_write_metadata_code_low(trex_file, str_compiled, max_str_len)
  endif

end function trexio_write_metadata_code

integer(trexio_exit_code) function trexio_write_metadata_author (trex_file, dset, max_str_len)
  implicit none
  integer(trexio_t), intent(in), value :: trex_file
  integer(c_int32_t), intent(in), value :: max_str_len
  character(len=*), intent(in) :: dset(*)

  character(len=:), allocatable :: str_compiled
  integer(c_int64_t) :: metadata_author_num
  integer(trexio_exit_code) :: rc

  rc = trexio_read_metadata_author_num_64(trex_file, metadata_author_num)
  if (rc /= TREXIO_SUCCESS) then
    trexio_write_metadata_author = rc
  else
    call trexio_strarray2str(dset, metadata_author_num, max_str_len, str_compiled)
    trexio_write_metadata_author = trexio_write_metadata_author_low(trex_file, str_compiled, max_str_len)
  endif

end function trexio_write_metadata_author

integer(trexio_exit_code) function trexio_write_nucleus_label (trex_file, dset, max_str_len)
  implicit none
  integer(trexio_t), intent(in), value :: trex_file
  integer(c_int32_t), intent(in), value :: max_str_len
  character(len=*), intent(in) :: dset(*)

  character(len=:), allocatable :: str_compiled
  integer(c_int64_t) :: nucleus_num
  integer(trexio_exit_code) :: rc

  rc = trexio_read_nucleus_num_64(trex_file, nucleus_num)
  if (rc /= TREXIO_SUCCESS) then
    trexio_write_nucleus_label = rc
  else
    call trexio_strarray2str(dset, nucleus_num, max_str_len, str_compiled)
    trexio_write_nucleus_label = trexio_write_nucleus_label_low(trex_file, str_compiled, max_str_len)
  endif

end function trexio_write_nucleus_label

integer(trexio_exit_code) function trexio_write_mo_class (trex_file, dset, max_str_len)
  implicit none
  integer(trexio_t), intent(in), value :: trex_file
  integer(c_int32_t), intent(in), value :: max_str_len
  character(len=*), intent(in) :: dset(*)

  character(len=:), allocatable :: str_compiled
  integer(c_int64_t) :: mo_num
  integer(trexio_exit_code) :: rc

  rc = trexio_read_mo_num_64(trex_file, mo_num)
  if (rc /= TREXIO_SUCCESS) then
    trexio_write_mo_class = rc
  else
    call trexio_strarray2str(dset, mo_num, max_str_len, str_compiled)
    trexio_write_mo_class = trexio_write_mo_class_low(trex_file, str_compiled, max_str_len)
  endif

end function trexio_write_mo_class

integer(trexio_exit_code) function trexio_write_mo_symmetry (trex_file, dset, max_str_len)
  implicit none
  integer(trexio_t), intent(in), value :: trex_file
  integer(c_int32_t), intent(in), value :: max_str_len
  character(len=*), intent(in) :: dset(*)

  character(len=:), allocatable :: str_compiled
  integer(c_int64_t) :: mo_num
  integer(trexio_exit_code) :: rc

  rc = trexio_read_mo_num_64(trex_file, mo_num)
  if (rc /= TREXIO_SUCCESS) then
    trexio_write_mo_symmetry = rc
  else
    call trexio_strarray2str(dset, mo_num, max_str_len, str_compiled)
    trexio_write_mo_symmetry = trexio_write_mo_symmetry_low(trex_file, str_compiled, max_str_len)
  endif

end function trexio_write_mo_symmetry

integer(trexio_exit_code) function trexio_write_state_label (trex_file, dset, max_str_len)
  implicit none
  integer(trexio_t), intent(in), value :: trex_file
  integer(c_int32_t), intent(in), value :: max_str_len
  character(len=*), intent(in) :: dset(*)

  character(len=:), allocatable :: str_compiled
  integer(c_int64_t) :: state_num
  integer(trexio_exit_code) :: rc

  rc = trexio_read_state_num_64(trex_file, state_num)
  if (rc /= TREXIO_SUCCESS) then
    trexio_write_state_label = rc
  else
    call trexio_strarray2str(dset, state_num, max_str_len, str_compiled)
    trexio_write_state_label = trexio_write_state_label_low(trex_file, str_compiled, max_str_len)
  endif

end function trexio_write_state_label

integer(trexio_exit_code) function trexio_write_state_file_name (trex_file, dset, max_str_len)
  implicit none
  integer(trexio_t), intent(in), value :: trex_file
  integer(c_int32_t), intent(in), value :: max_str_len
  character(len=*), intent(in) :: dset(*)

  character(len=:), allocatable :: str_compiled
  integer(c_int64_t) :: state_num
  integer(trexio_exit_code) :: rc

  rc = trexio_read_state_num_64(trex_file, state_num)
  if (rc /= TREXIO_SUCCESS) then
    trexio_write_state_file_name = rc
  else
    call trexio_strarray2str(dset, state_num, max_str_len, str_compiled)
    trexio_write_state_file_name = trexio_write_state_file_name_low(trex_file, str_compiled, max_str_len)
  endif

end function trexio_write_state_file_name

end module trexio
