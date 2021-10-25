module trexio

  use, intrinsic :: iso_c_binding
  implicit none

  integer, parameter :: trexio_exit_code = 4
  integer, parameter :: trexio_backend = 4

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
integer(trexio_exit_code), parameter :: TREXIO_INVALID_STR_LEN         = 30

interface
   subroutine trexio_string_of_error (error, string) bind(C, name='trexio_string_of_error_f')
     use, intrinsic :: iso_c_binding
     import
     integer (trexio_exit_code), intent(in), value :: error
     character, intent(out) :: string(128)
   end subroutine trexio_string_of_error
end interface

integer(trexio_backend), parameter :: TREXIO_HDF5 = 0
  integer(trexio_backend), parameter :: TREXIO_TEXT = 1
! integer(trexio_backend), parameter :: TREXIO_JSON = 2
  integer(trexio_backend), parameter :: TREXIO_INVALID_BACK_END = 2

interface
   integer(8) function trexio_open_c (filename, mode, backend, rc_open) bind(C, name="trexio_open")
     use, intrinsic :: iso_c_binding
     import
     character(kind=c_char), dimension(*)       :: filename
     character, intent(in), value               :: mode
     integer(trexio_backend), intent(in), value :: backend
     integer(trexio_exit_code), intent(out)     :: rc_open
   end function trexio_open_c
end interface

interface
   integer function trexio_set_one_based(trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
   end function trexio_set_one_based
end interface

interface
   integer function trexio_close (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
   end function trexio_close
end interface

character(len = 12) :: TREXIO_PACKAGE_VERSION = "2.0.0"
integer(4) :: TREXIO_VERSION_MAJOR = 2
integer(4) :: TREXIO_VERSION_MINOR = 0
integer(4) :: TREXIO_VERSION_PATCH = 0

interface
   integer function trexio_has_metadata_code_num (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
   end function trexio_has_metadata_code_num
end interface

interface
   integer function trexio_has_metadata_author_num (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
   end function trexio_has_metadata_author_num
end interface

interface
   integer function trexio_has_electron_up_num (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
   end function trexio_has_electron_up_num
end interface

interface
   integer function trexio_has_electron_dn_num (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
   end function trexio_has_electron_dn_num
end interface

interface
   integer function trexio_has_nucleus_num (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
   end function trexio_has_nucleus_num
end interface

interface
   integer function trexio_has_ecp_num (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
   end function trexio_has_ecp_num
end interface

interface
   integer function trexio_has_basis_prim_num (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
   end function trexio_has_basis_prim_num
end interface

interface
   integer function trexio_has_basis_shell_num (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
   end function trexio_has_basis_shell_num
end interface

interface
   integer function trexio_has_ao_cartesian (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
   end function trexio_has_ao_cartesian
end interface

interface
   integer function trexio_has_ao_num (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
   end function trexio_has_ao_num
end interface

interface
   integer function trexio_has_mo_num (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
   end function trexio_has_mo_num
end interface

interface
   integer function trexio_has_metadata_package_version (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
   end function trexio_has_metadata_package_version
end interface

interface
   integer function trexio_has_metadata_description (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
   end function trexio_has_metadata_description
end interface

interface
   integer function trexio_has_nucleus_point_group (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
   end function trexio_has_nucleus_point_group
end interface

interface
   integer function trexio_has_basis_type (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
   end function trexio_has_basis_type
end interface

interface
   integer function trexio_has_mo_type (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
   end function trexio_has_mo_type
end interface

interface
   integer function trexio_has_nucleus_charge (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
   end function trexio_has_nucleus_charge
end interface

interface
   integer function trexio_has_nucleus_coord (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
   end function trexio_has_nucleus_coord
end interface

interface
   integer function trexio_has_ecp_max_ang_mom_plus_1 (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
   end function trexio_has_ecp_max_ang_mom_plus_1
end interface

interface
   integer function trexio_has_ecp_z_core (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
   end function trexio_has_ecp_z_core
end interface

interface
   integer function trexio_has_ecp_ang_mom (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
   end function trexio_has_ecp_ang_mom
end interface

interface
   integer function trexio_has_ecp_nucleus_index (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
   end function trexio_has_ecp_nucleus_index
end interface

interface
   integer function trexio_has_ecp_exponent (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
   end function trexio_has_ecp_exponent
end interface

interface
   integer function trexio_has_ecp_coefficient (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
   end function trexio_has_ecp_coefficient
end interface

interface
   integer function trexio_has_ecp_power (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
   end function trexio_has_ecp_power
end interface

interface
   integer function trexio_has_basis_nucleus_index (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
   end function trexio_has_basis_nucleus_index
end interface

interface
   integer function trexio_has_basis_ang_mom (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
   end function trexio_has_basis_ang_mom
end interface

interface
   integer function trexio_has_basis_shell_factor (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
   end function trexio_has_basis_shell_factor
end interface

interface
   integer function trexio_has_basis_shell_index (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
   end function trexio_has_basis_shell_index
end interface

interface
   integer function trexio_has_basis_exponent (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
   end function trexio_has_basis_exponent
end interface

interface
   integer function trexio_has_basis_coefficient (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
   end function trexio_has_basis_coefficient
end interface

interface
   integer function trexio_has_basis_prim_factor (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
   end function trexio_has_basis_prim_factor
end interface

interface
   integer function trexio_has_ao_shell (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
   end function trexio_has_ao_shell
end interface

interface
   integer function trexio_has_ao_normalization (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
   end function trexio_has_ao_normalization
end interface

interface
   integer function trexio_has_ao_1e_int_overlap (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
   end function trexio_has_ao_1e_int_overlap
end interface

interface
   integer function trexio_has_ao_1e_int_kinetic (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
   end function trexio_has_ao_1e_int_kinetic
end interface

interface
   integer function trexio_has_ao_1e_int_potential_n_e (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
   end function trexio_has_ao_1e_int_potential_n_e
end interface

interface
   integer function trexio_has_ao_1e_int_ecp_local (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
   end function trexio_has_ao_1e_int_ecp_local
end interface

interface
   integer function trexio_has_ao_1e_int_ecp_non_local (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
   end function trexio_has_ao_1e_int_ecp_non_local
end interface

interface
   integer function trexio_has_ao_1e_int_core_hamiltonian (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
   end function trexio_has_ao_1e_int_core_hamiltonian
end interface

interface
   integer function trexio_has_ao_2e_int_eri (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
   end function trexio_has_ao_2e_int_eri
end interface

interface
   integer function trexio_has_ao_2e_int_eri_lr (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
   end function trexio_has_ao_2e_int_eri_lr
end interface

interface
   integer function trexio_has_mo_coefficient (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
   end function trexio_has_mo_coefficient
end interface

interface
   integer function trexio_has_mo_occupation (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
   end function trexio_has_mo_occupation
end interface

interface
   integer function trexio_has_mo_1e_int_overlap (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
   end function trexio_has_mo_1e_int_overlap
end interface

interface
   integer function trexio_has_mo_1e_int_kinetic (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
   end function trexio_has_mo_1e_int_kinetic
end interface

interface
   integer function trexio_has_mo_1e_int_potential_n_e (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
   end function trexio_has_mo_1e_int_potential_n_e
end interface

interface
   integer function trexio_has_mo_1e_int_ecp_local (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
   end function trexio_has_mo_1e_int_ecp_local
end interface

interface
   integer function trexio_has_mo_1e_int_ecp_non_local (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
   end function trexio_has_mo_1e_int_ecp_non_local
end interface

interface
   integer function trexio_has_mo_1e_int_core_hamiltonian (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
   end function trexio_has_mo_1e_int_core_hamiltonian
end interface

interface
   integer function trexio_has_mo_2e_int_eri (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
   end function trexio_has_mo_2e_int_eri
end interface

interface
   integer function trexio_has_mo_2e_int_eri_lr (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
   end function trexio_has_mo_2e_int_eri_lr
end interface

interface
   integer function trexio_has_metadata_code (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
   end function trexio_has_metadata_code
end interface

interface
   integer function trexio_has_metadata_author (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
   end function trexio_has_metadata_author
end interface

interface
   integer function trexio_has_nucleus_label (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
   end function trexio_has_nucleus_label
end interface

interface
   integer function trexio_has_mo_class (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
   end function trexio_has_mo_class
end interface

interface
   integer function trexio_has_mo_symmetry (trex_file) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
   end function trexio_has_mo_symmetry
end interface

interface
   integer function trexio_read_metadata_code_num_32 (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     integer(4), intent(out) :: num
   end function trexio_read_metadata_code_num_32
end interface

interface
   integer function trexio_read_metadata_author_num_32 (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     integer(4), intent(out) :: num
   end function trexio_read_metadata_author_num_32
end interface

interface
   integer function trexio_read_electron_up_num_32 (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     integer(4), intent(out) :: num
   end function trexio_read_electron_up_num_32
end interface

interface
   integer function trexio_read_electron_dn_num_32 (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     integer(4), intent(out) :: num
   end function trexio_read_electron_dn_num_32
end interface

interface
   integer function trexio_read_nucleus_num_32 (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     integer(4), intent(out) :: num
   end function trexio_read_nucleus_num_32
end interface

interface
   integer function trexio_read_ecp_num_32 (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     integer(4), intent(out) :: num
   end function trexio_read_ecp_num_32
end interface

interface
   integer function trexio_read_basis_prim_num_32 (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     integer(4), intent(out) :: num
   end function trexio_read_basis_prim_num_32
end interface

interface
   integer function trexio_read_basis_shell_num_32 (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     integer(4), intent(out) :: num
   end function trexio_read_basis_shell_num_32
end interface

interface
   integer function trexio_read_ao_cartesian_32 (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     integer(4), intent(out) :: num
   end function trexio_read_ao_cartesian_32
end interface

interface
   integer function trexio_read_ao_num_32 (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     integer(4), intent(out) :: num
   end function trexio_read_ao_num_32
end interface

interface
   integer function trexio_read_mo_num_32 (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     integer(4), intent(out) :: num
   end function trexio_read_mo_num_32
end interface

interface
   integer function trexio_read_metadata_code_num_64 (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     integer(8), intent(out) :: num
   end function trexio_read_metadata_code_num_64
end interface

interface
   integer function trexio_read_metadata_author_num_64 (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     integer(8), intent(out) :: num
   end function trexio_read_metadata_author_num_64
end interface

interface
   integer function trexio_read_electron_up_num_64 (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     integer(8), intent(out) :: num
   end function trexio_read_electron_up_num_64
end interface

interface
   integer function trexio_read_electron_dn_num_64 (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     integer(8), intent(out) :: num
   end function trexio_read_electron_dn_num_64
end interface

interface
   integer function trexio_read_nucleus_num_64 (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     integer(8), intent(out) :: num
   end function trexio_read_nucleus_num_64
end interface

interface
   integer function trexio_read_ecp_num_64 (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     integer(8), intent(out) :: num
   end function trexio_read_ecp_num_64
end interface

interface
   integer function trexio_read_basis_prim_num_64 (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     integer(8), intent(out) :: num
   end function trexio_read_basis_prim_num_64
end interface

interface
   integer function trexio_read_basis_shell_num_64 (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     integer(8), intent(out) :: num
   end function trexio_read_basis_shell_num_64
end interface

interface
   integer function trexio_read_ao_cartesian_64 (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     integer(8), intent(out) :: num
   end function trexio_read_ao_cartesian_64
end interface

interface
   integer function trexio_read_ao_num_64 (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     integer(8), intent(out) :: num
   end function trexio_read_ao_num_64
end interface

interface
   integer function trexio_read_mo_num_64 (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     integer(8), intent(out) :: num
   end function trexio_read_mo_num_64
end interface

interface
   integer function trexio_read_metadata_code_num (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     integer(4), intent(out) :: num
   end function trexio_read_metadata_code_num
end interface

interface
   integer function trexio_read_metadata_author_num (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     integer(4), intent(out) :: num
   end function trexio_read_metadata_author_num
end interface

interface
   integer function trexio_read_electron_up_num (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     integer(4), intent(out) :: num
   end function trexio_read_electron_up_num
end interface

interface
   integer function trexio_read_electron_dn_num (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     integer(4), intent(out) :: num
   end function trexio_read_electron_dn_num
end interface

interface
   integer function trexio_read_nucleus_num (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     integer(4), intent(out) :: num
   end function trexio_read_nucleus_num
end interface

interface
   integer function trexio_read_ecp_num (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     integer(4), intent(out) :: num
   end function trexio_read_ecp_num
end interface

interface
   integer function trexio_read_basis_prim_num (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     integer(4), intent(out) :: num
   end function trexio_read_basis_prim_num
end interface

interface
   integer function trexio_read_basis_shell_num (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     integer(4), intent(out) :: num
   end function trexio_read_basis_shell_num
end interface

interface
   integer function trexio_read_ao_cartesian (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     integer(4), intent(out) :: num
   end function trexio_read_ao_cartesian
end interface

interface
   integer function trexio_read_ao_num (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     integer(4), intent(out) :: num
   end function trexio_read_ao_num
end interface

interface
   integer function trexio_read_mo_num (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     integer(4), intent(out) :: num
   end function trexio_read_mo_num
end interface

interface
   integer function trexio_read_metadata_package_version_c (trex_file, str, max_str_len) &
           bind(C, name="trexio_read_metadata_package_version")
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     character, intent(out) :: str(*)
     integer(4), intent(in), value :: max_str_len
   end function trexio_read_metadata_package_version_c
end interface

interface
   integer function trexio_read_metadata_description_c (trex_file, str, max_str_len) &
           bind(C, name="trexio_read_metadata_description")
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     character, intent(out) :: str(*)
     integer(4), intent(in), value :: max_str_len
   end function trexio_read_metadata_description_c
end interface

interface
   integer function trexio_read_nucleus_point_group_c (trex_file, str, max_str_len) &
           bind(C, name="trexio_read_nucleus_point_group")
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     character, intent(out) :: str(*)
     integer(4), intent(in), value :: max_str_len
   end function trexio_read_nucleus_point_group_c
end interface

interface
   integer function trexio_read_basis_type_c (trex_file, str, max_str_len) &
           bind(C, name="trexio_read_basis_type")
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     character, intent(out) :: str(*)
     integer(4), intent(in), value :: max_str_len
   end function trexio_read_basis_type_c
end interface

interface
   integer function trexio_read_mo_type_c (trex_file, str, max_str_len) &
           bind(C, name="trexio_read_mo_type")
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     character, intent(out) :: str(*)
     integer(4), intent(in), value :: max_str_len
   end function trexio_read_mo_type_c
end interface

interface
   integer function trexio_read_nucleus_charge_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(4), intent(out) :: dset(*)
   end function trexio_read_nucleus_charge_32
end interface

interface
   integer function trexio_read_nucleus_coord_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(4), intent(out) :: dset(*)
   end function trexio_read_nucleus_coord_32
end interface

interface
   integer function trexio_read_ecp_max_ang_mom_plus_1_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     integer(4), intent(out) :: dset(*)
   end function trexio_read_ecp_max_ang_mom_plus_1_32
end interface

interface
   integer function trexio_read_ecp_z_core_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     integer(4), intent(out) :: dset(*)
   end function trexio_read_ecp_z_core_32
end interface

interface
   integer function trexio_read_ecp_ang_mom_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     integer(4), intent(out) :: dset(*)
   end function trexio_read_ecp_ang_mom_32
end interface

interface
   integer function trexio_read_ecp_nucleus_index_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     integer(4), intent(out) :: dset(*)
   end function trexio_read_ecp_nucleus_index_32
end interface

interface
   integer function trexio_read_ecp_exponent_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(4), intent(out) :: dset(*)
   end function trexio_read_ecp_exponent_32
end interface

interface
   integer function trexio_read_ecp_coefficient_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(4), intent(out) :: dset(*)
   end function trexio_read_ecp_coefficient_32
end interface

interface
   integer function trexio_read_ecp_power_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     integer(4), intent(out) :: dset(*)
   end function trexio_read_ecp_power_32
end interface

interface
   integer function trexio_read_basis_nucleus_index_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     integer(4), intent(out) :: dset(*)
   end function trexio_read_basis_nucleus_index_32
end interface

interface
   integer function trexio_read_basis_ang_mom_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     integer(4), intent(out) :: dset(*)
   end function trexio_read_basis_ang_mom_32
end interface

interface
   integer function trexio_read_basis_shell_factor_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(4), intent(out) :: dset(*)
   end function trexio_read_basis_shell_factor_32
end interface

interface
   integer function trexio_read_basis_shell_index_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     integer(4), intent(out) :: dset(*)
   end function trexio_read_basis_shell_index_32
end interface

interface
   integer function trexio_read_basis_exponent_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(4), intent(out) :: dset(*)
   end function trexio_read_basis_exponent_32
end interface

interface
   integer function trexio_read_basis_coefficient_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(4), intent(out) :: dset(*)
   end function trexio_read_basis_coefficient_32
end interface

interface
   integer function trexio_read_basis_prim_factor_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(4), intent(out) :: dset(*)
   end function trexio_read_basis_prim_factor_32
end interface

interface
   integer function trexio_read_ao_shell_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     integer(4), intent(out) :: dset(*)
   end function trexio_read_ao_shell_32
end interface

interface
   integer function trexio_read_ao_normalization_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(4), intent(out) :: dset(*)
   end function trexio_read_ao_normalization_32
end interface

interface
   integer function trexio_read_ao_1e_int_overlap_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(4), intent(out) :: dset(*)
   end function trexio_read_ao_1e_int_overlap_32
end interface

interface
   integer function trexio_read_ao_1e_int_kinetic_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(4), intent(out) :: dset(*)
   end function trexio_read_ao_1e_int_kinetic_32
end interface

interface
   integer function trexio_read_ao_1e_int_potential_n_e_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(4), intent(out) :: dset(*)
   end function trexio_read_ao_1e_int_potential_n_e_32
end interface

interface
   integer function trexio_read_ao_1e_int_ecp_local_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(4), intent(out) :: dset(*)
   end function trexio_read_ao_1e_int_ecp_local_32
end interface

interface
   integer function trexio_read_ao_1e_int_ecp_non_local_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(4), intent(out) :: dset(*)
   end function trexio_read_ao_1e_int_ecp_non_local_32
end interface

interface
   integer function trexio_read_ao_1e_int_core_hamiltonian_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(4), intent(out) :: dset(*)
   end function trexio_read_ao_1e_int_core_hamiltonian_32
end interface

interface
   integer function trexio_read_ao_2e_int_eri_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(4), intent(out) :: dset(*)
   end function trexio_read_ao_2e_int_eri_32
end interface

interface
   integer function trexio_read_ao_2e_int_eri_lr_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(4), intent(out) :: dset(*)
   end function trexio_read_ao_2e_int_eri_lr_32
end interface

interface
   integer function trexio_read_mo_coefficient_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(4), intent(out) :: dset(*)
   end function trexio_read_mo_coefficient_32
end interface

interface
   integer function trexio_read_mo_occupation_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(4), intent(out) :: dset(*)
   end function trexio_read_mo_occupation_32
end interface

interface
   integer function trexio_read_mo_1e_int_overlap_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(4), intent(out) :: dset(*)
   end function trexio_read_mo_1e_int_overlap_32
end interface

interface
   integer function trexio_read_mo_1e_int_kinetic_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(4), intent(out) :: dset(*)
   end function trexio_read_mo_1e_int_kinetic_32
end interface

interface
   integer function trexio_read_mo_1e_int_potential_n_e_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(4), intent(out) :: dset(*)
   end function trexio_read_mo_1e_int_potential_n_e_32
end interface

interface
   integer function trexio_read_mo_1e_int_ecp_local_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(4), intent(out) :: dset(*)
   end function trexio_read_mo_1e_int_ecp_local_32
end interface

interface
   integer function trexio_read_mo_1e_int_ecp_non_local_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(4), intent(out) :: dset(*)
   end function trexio_read_mo_1e_int_ecp_non_local_32
end interface

interface
   integer function trexio_read_mo_1e_int_core_hamiltonian_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(4), intent(out) :: dset(*)
   end function trexio_read_mo_1e_int_core_hamiltonian_32
end interface

interface
   integer function trexio_read_mo_2e_int_eri_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(4), intent(out) :: dset(*)
   end function trexio_read_mo_2e_int_eri_32
end interface

interface
   integer function trexio_read_mo_2e_int_eri_lr_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(4), intent(out) :: dset(*)
   end function trexio_read_mo_2e_int_eri_lr_32
end interface

interface
   integer function trexio_read_nucleus_charge_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(8), intent(out) :: dset(*)
   end function trexio_read_nucleus_charge_64
end interface

interface
   integer function trexio_read_nucleus_coord_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(8), intent(out) :: dset(*)
   end function trexio_read_nucleus_coord_64
end interface

interface
   integer function trexio_read_ecp_max_ang_mom_plus_1_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     integer(8), intent(out) :: dset(*)
   end function trexio_read_ecp_max_ang_mom_plus_1_64
end interface

interface
   integer function trexio_read_ecp_z_core_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     integer(8), intent(out) :: dset(*)
   end function trexio_read_ecp_z_core_64
end interface

interface
   integer function trexio_read_ecp_ang_mom_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     integer(8), intent(out) :: dset(*)
   end function trexio_read_ecp_ang_mom_64
end interface

interface
   integer function trexio_read_ecp_nucleus_index_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     integer(8), intent(out) :: dset(*)
   end function trexio_read_ecp_nucleus_index_64
end interface

interface
   integer function trexio_read_ecp_exponent_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(8), intent(out) :: dset(*)
   end function trexio_read_ecp_exponent_64
end interface

interface
   integer function trexio_read_ecp_coefficient_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(8), intent(out) :: dset(*)
   end function trexio_read_ecp_coefficient_64
end interface

interface
   integer function trexio_read_ecp_power_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     integer(8), intent(out) :: dset(*)
   end function trexio_read_ecp_power_64
end interface

interface
   integer function trexio_read_basis_nucleus_index_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     integer(8), intent(out) :: dset(*)
   end function trexio_read_basis_nucleus_index_64
end interface

interface
   integer function trexio_read_basis_ang_mom_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     integer(8), intent(out) :: dset(*)
   end function trexio_read_basis_ang_mom_64
end interface

interface
   integer function trexio_read_basis_shell_factor_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(8), intent(out) :: dset(*)
   end function trexio_read_basis_shell_factor_64
end interface

interface
   integer function trexio_read_basis_shell_index_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     integer(8), intent(out) :: dset(*)
   end function trexio_read_basis_shell_index_64
end interface

interface
   integer function trexio_read_basis_exponent_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(8), intent(out) :: dset(*)
   end function trexio_read_basis_exponent_64
end interface

interface
   integer function trexio_read_basis_coefficient_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(8), intent(out) :: dset(*)
   end function trexio_read_basis_coefficient_64
end interface

interface
   integer function trexio_read_basis_prim_factor_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(8), intent(out) :: dset(*)
   end function trexio_read_basis_prim_factor_64
end interface

interface
   integer function trexio_read_ao_shell_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     integer(8), intent(out) :: dset(*)
   end function trexio_read_ao_shell_64
end interface

interface
   integer function trexio_read_ao_normalization_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(8), intent(out) :: dset(*)
   end function trexio_read_ao_normalization_64
end interface

interface
   integer function trexio_read_ao_1e_int_overlap_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(8), intent(out) :: dset(*)
   end function trexio_read_ao_1e_int_overlap_64
end interface

interface
   integer function trexio_read_ao_1e_int_kinetic_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(8), intent(out) :: dset(*)
   end function trexio_read_ao_1e_int_kinetic_64
end interface

interface
   integer function trexio_read_ao_1e_int_potential_n_e_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(8), intent(out) :: dset(*)
   end function trexio_read_ao_1e_int_potential_n_e_64
end interface

interface
   integer function trexio_read_ao_1e_int_ecp_local_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(8), intent(out) :: dset(*)
   end function trexio_read_ao_1e_int_ecp_local_64
end interface

interface
   integer function trexio_read_ao_1e_int_ecp_non_local_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(8), intent(out) :: dset(*)
   end function trexio_read_ao_1e_int_ecp_non_local_64
end interface

interface
   integer function trexio_read_ao_1e_int_core_hamiltonian_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(8), intent(out) :: dset(*)
   end function trexio_read_ao_1e_int_core_hamiltonian_64
end interface

interface
   integer function trexio_read_ao_2e_int_eri_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(8), intent(out) :: dset(*)
   end function trexio_read_ao_2e_int_eri_64
end interface

interface
   integer function trexio_read_ao_2e_int_eri_lr_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(8), intent(out) :: dset(*)
   end function trexio_read_ao_2e_int_eri_lr_64
end interface

interface
   integer function trexio_read_mo_coefficient_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(8), intent(out) :: dset(*)
   end function trexio_read_mo_coefficient_64
end interface

interface
   integer function trexio_read_mo_occupation_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(8), intent(out) :: dset(*)
   end function trexio_read_mo_occupation_64
end interface

interface
   integer function trexio_read_mo_1e_int_overlap_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(8), intent(out) :: dset(*)
   end function trexio_read_mo_1e_int_overlap_64
end interface

interface
   integer function trexio_read_mo_1e_int_kinetic_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(8), intent(out) :: dset(*)
   end function trexio_read_mo_1e_int_kinetic_64
end interface

interface
   integer function trexio_read_mo_1e_int_potential_n_e_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(8), intent(out) :: dset(*)
   end function trexio_read_mo_1e_int_potential_n_e_64
end interface

interface
   integer function trexio_read_mo_1e_int_ecp_local_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(8), intent(out) :: dset(*)
   end function trexio_read_mo_1e_int_ecp_local_64
end interface

interface
   integer function trexio_read_mo_1e_int_ecp_non_local_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(8), intent(out) :: dset(*)
   end function trexio_read_mo_1e_int_ecp_non_local_64
end interface

interface
   integer function trexio_read_mo_1e_int_core_hamiltonian_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(8), intent(out) :: dset(*)
   end function trexio_read_mo_1e_int_core_hamiltonian_64
end interface

interface
   integer function trexio_read_mo_2e_int_eri_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(8), intent(out) :: dset(*)
   end function trexio_read_mo_2e_int_eri_64
end interface

interface
   integer function trexio_read_mo_2e_int_eri_lr_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(8), intent(out) :: dset(*)
   end function trexio_read_mo_2e_int_eri_lr_64
end interface

interface
   integer function trexio_read_nucleus_charge (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(8), intent(out) :: dset(*)
   end function trexio_read_nucleus_charge
end interface

interface
   integer function trexio_read_nucleus_coord (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(8), intent(out) :: dset(*)
   end function trexio_read_nucleus_coord
end interface

interface
   integer function trexio_read_ecp_max_ang_mom_plus_1 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     integer(4), intent(out) :: dset(*)
   end function trexio_read_ecp_max_ang_mom_plus_1
end interface

interface
   integer function trexio_read_ecp_z_core (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     integer(4), intent(out) :: dset(*)
   end function trexio_read_ecp_z_core
end interface

interface
   integer function trexio_read_ecp_ang_mom (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     integer(4), intent(out) :: dset(*)
   end function trexio_read_ecp_ang_mom
end interface

interface
   integer function trexio_read_ecp_nucleus_index (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     integer(4), intent(out) :: dset(*)
   end function trexio_read_ecp_nucleus_index
end interface

interface
   integer function trexio_read_ecp_exponent (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(8), intent(out) :: dset(*)
   end function trexio_read_ecp_exponent
end interface

interface
   integer function trexio_read_ecp_coefficient (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(8), intent(out) :: dset(*)
   end function trexio_read_ecp_coefficient
end interface

interface
   integer function trexio_read_ecp_power (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     integer(4), intent(out) :: dset(*)
   end function trexio_read_ecp_power
end interface

interface
   integer function trexio_read_basis_nucleus_index (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     integer(4), intent(out) :: dset(*)
   end function trexio_read_basis_nucleus_index
end interface

interface
   integer function trexio_read_basis_ang_mom (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     integer(4), intent(out) :: dset(*)
   end function trexio_read_basis_ang_mom
end interface

interface
   integer function trexio_read_basis_shell_factor (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(8), intent(out) :: dset(*)
   end function trexio_read_basis_shell_factor
end interface

interface
   integer function trexio_read_basis_shell_index (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     integer(4), intent(out) :: dset(*)
   end function trexio_read_basis_shell_index
end interface

interface
   integer function trexio_read_basis_exponent (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(8), intent(out) :: dset(*)
   end function trexio_read_basis_exponent
end interface

interface
   integer function trexio_read_basis_coefficient (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(8), intent(out) :: dset(*)
   end function trexio_read_basis_coefficient
end interface

interface
   integer function trexio_read_basis_prim_factor (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(8), intent(out) :: dset(*)
   end function trexio_read_basis_prim_factor
end interface

interface
   integer function trexio_read_ao_shell (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     integer(4), intent(out) :: dset(*)
   end function trexio_read_ao_shell
end interface

interface
   integer function trexio_read_ao_normalization (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(8), intent(out) :: dset(*)
   end function trexio_read_ao_normalization
end interface

interface
   integer function trexio_read_ao_1e_int_overlap (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(8), intent(out) :: dset(*)
   end function trexio_read_ao_1e_int_overlap
end interface

interface
   integer function trexio_read_ao_1e_int_kinetic (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(8), intent(out) :: dset(*)
   end function trexio_read_ao_1e_int_kinetic
end interface

interface
   integer function trexio_read_ao_1e_int_potential_n_e (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(8), intent(out) :: dset(*)
   end function trexio_read_ao_1e_int_potential_n_e
end interface

interface
   integer function trexio_read_ao_1e_int_ecp_local (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(8), intent(out) :: dset(*)
   end function trexio_read_ao_1e_int_ecp_local
end interface

interface
   integer function trexio_read_ao_1e_int_ecp_non_local (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(8), intent(out) :: dset(*)
   end function trexio_read_ao_1e_int_ecp_non_local
end interface

interface
   integer function trexio_read_ao_1e_int_core_hamiltonian (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(8), intent(out) :: dset(*)
   end function trexio_read_ao_1e_int_core_hamiltonian
end interface

interface
   integer function trexio_read_ao_2e_int_eri (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(8), intent(out) :: dset(*)
   end function trexio_read_ao_2e_int_eri
end interface

interface
   integer function trexio_read_ao_2e_int_eri_lr (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(8), intent(out) :: dset(*)
   end function trexio_read_ao_2e_int_eri_lr
end interface

interface
   integer function trexio_read_mo_coefficient (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(8), intent(out) :: dset(*)
   end function trexio_read_mo_coefficient
end interface

interface
   integer function trexio_read_mo_occupation (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(8), intent(out) :: dset(*)
   end function trexio_read_mo_occupation
end interface

interface
   integer function trexio_read_mo_1e_int_overlap (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(8), intent(out) :: dset(*)
   end function trexio_read_mo_1e_int_overlap
end interface

interface
   integer function trexio_read_mo_1e_int_kinetic (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(8), intent(out) :: dset(*)
   end function trexio_read_mo_1e_int_kinetic
end interface

interface
   integer function trexio_read_mo_1e_int_potential_n_e (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(8), intent(out) :: dset(*)
   end function trexio_read_mo_1e_int_potential_n_e
end interface

interface
   integer function trexio_read_mo_1e_int_ecp_local (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(8), intent(out) :: dset(*)
   end function trexio_read_mo_1e_int_ecp_local
end interface

interface
   integer function trexio_read_mo_1e_int_ecp_non_local (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(8), intent(out) :: dset(*)
   end function trexio_read_mo_1e_int_ecp_non_local
end interface

interface
   integer function trexio_read_mo_1e_int_core_hamiltonian (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(8), intent(out) :: dset(*)
   end function trexio_read_mo_1e_int_core_hamiltonian
end interface

interface
   integer function trexio_read_mo_2e_int_eri (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(8), intent(out) :: dset(*)
   end function trexio_read_mo_2e_int_eri
end interface

interface
   integer function trexio_read_mo_2e_int_eri_lr (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(8), intent(out) :: dset(*)
   end function trexio_read_mo_2e_int_eri_lr
end interface

interface
   integer function trexio_read_metadata_code_low (trex_file, dset, max_str_len) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     character, intent(out) :: dset(*)
     integer(4), intent(in), value :: max_str_len
   end function trexio_read_metadata_code_low
end interface

interface
   integer function trexio_read_metadata_author_low (trex_file, dset, max_str_len) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     character, intent(out) :: dset(*)
     integer(4), intent(in), value :: max_str_len
   end function trexio_read_metadata_author_low
end interface

interface
   integer function trexio_read_nucleus_label_low (trex_file, dset, max_str_len) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     character, intent(out) :: dset(*)
     integer(4), intent(in), value :: max_str_len
   end function trexio_read_nucleus_label_low
end interface

interface
   integer function trexio_read_mo_class_low (trex_file, dset, max_str_len) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     character, intent(out) :: dset(*)
     integer(4), intent(in), value :: max_str_len
   end function trexio_read_mo_class_low
end interface

interface
   integer function trexio_read_mo_symmetry_low (trex_file, dset, max_str_len) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     character, intent(out) :: dset(*)
     integer(4), intent(in), value :: max_str_len
   end function trexio_read_mo_symmetry_low
end interface

interface
   integer function trexio_write_metadata_code_num_32 (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     integer(4), intent(in), value :: num
   end function trexio_write_metadata_code_num_32
end interface

interface
   integer function trexio_write_metadata_author_num_32 (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     integer(4), intent(in), value :: num
   end function trexio_write_metadata_author_num_32
end interface

interface
   integer function trexio_write_electron_up_num_32 (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     integer(4), intent(in), value :: num
   end function trexio_write_electron_up_num_32
end interface

interface
   integer function trexio_write_electron_dn_num_32 (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     integer(4), intent(in), value :: num
   end function trexio_write_electron_dn_num_32
end interface

interface
   integer function trexio_write_nucleus_num_32 (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     integer(4), intent(in), value :: num
   end function trexio_write_nucleus_num_32
end interface

interface
   integer function trexio_write_ecp_num_32 (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     integer(4), intent(in), value :: num
   end function trexio_write_ecp_num_32
end interface

interface
   integer function trexio_write_basis_prim_num_32 (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     integer(4), intent(in), value :: num
   end function trexio_write_basis_prim_num_32
end interface

interface
   integer function trexio_write_basis_shell_num_32 (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     integer(4), intent(in), value :: num
   end function trexio_write_basis_shell_num_32
end interface

interface
   integer function trexio_write_ao_cartesian_32 (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     integer(4), intent(in), value :: num
   end function trexio_write_ao_cartesian_32
end interface

interface
   integer function trexio_write_ao_num_32 (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     integer(4), intent(in), value :: num
   end function trexio_write_ao_num_32
end interface

interface
   integer function trexio_write_mo_num_32 (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     integer(4), intent(in), value :: num
   end function trexio_write_mo_num_32
end interface

interface
   integer function trexio_write_metadata_code_num_64 (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     integer(8), intent(in), value :: num
   end function trexio_write_metadata_code_num_64
end interface

interface
   integer function trexio_write_metadata_author_num_64 (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     integer(8), intent(in), value :: num
   end function trexio_write_metadata_author_num_64
end interface

interface
   integer function trexio_write_electron_up_num_64 (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     integer(8), intent(in), value :: num
   end function trexio_write_electron_up_num_64
end interface

interface
   integer function trexio_write_electron_dn_num_64 (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     integer(8), intent(in), value :: num
   end function trexio_write_electron_dn_num_64
end interface

interface
   integer function trexio_write_nucleus_num_64 (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     integer(8), intent(in), value :: num
   end function trexio_write_nucleus_num_64
end interface

interface
   integer function trexio_write_ecp_num_64 (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     integer(8), intent(in), value :: num
   end function trexio_write_ecp_num_64
end interface

interface
   integer function trexio_write_basis_prim_num_64 (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     integer(8), intent(in), value :: num
   end function trexio_write_basis_prim_num_64
end interface

interface
   integer function trexio_write_basis_shell_num_64 (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     integer(8), intent(in), value :: num
   end function trexio_write_basis_shell_num_64
end interface

interface
   integer function trexio_write_ao_cartesian_64 (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     integer(8), intent(in), value :: num
   end function trexio_write_ao_cartesian_64
end interface

interface
   integer function trexio_write_ao_num_64 (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     integer(8), intent(in), value :: num
   end function trexio_write_ao_num_64
end interface

interface
   integer function trexio_write_mo_num_64 (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     integer(8), intent(in), value :: num
   end function trexio_write_mo_num_64
end interface

interface
   integer function trexio_write_metadata_code_num (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     integer(4), intent(in), value :: num
   end function trexio_write_metadata_code_num
end interface

interface
   integer function trexio_write_metadata_author_num (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     integer(4), intent(in), value :: num
   end function trexio_write_metadata_author_num
end interface

interface
   integer function trexio_write_electron_up_num (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     integer(4), intent(in), value :: num
   end function trexio_write_electron_up_num
end interface

interface
   integer function trexio_write_electron_dn_num (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     integer(4), intent(in), value :: num
   end function trexio_write_electron_dn_num
end interface

interface
   integer function trexio_write_nucleus_num (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     integer(4), intent(in), value :: num
   end function trexio_write_nucleus_num
end interface

interface
   integer function trexio_write_ecp_num (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     integer(4), intent(in), value :: num
   end function trexio_write_ecp_num
end interface

interface
   integer function trexio_write_basis_prim_num (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     integer(4), intent(in), value :: num
   end function trexio_write_basis_prim_num
end interface

interface
   integer function trexio_write_basis_shell_num (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     integer(4), intent(in), value :: num
   end function trexio_write_basis_shell_num
end interface

interface
   integer function trexio_write_ao_cartesian (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     integer(4), intent(in), value :: num
   end function trexio_write_ao_cartesian
end interface

interface
   integer function trexio_write_ao_num (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     integer(4), intent(in), value :: num
   end function trexio_write_ao_num
end interface

interface
   integer function trexio_write_mo_num (trex_file, num) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     integer(4), intent(in), value :: num
   end function trexio_write_mo_num
end interface

interface
   integer function trexio_write_metadata_package_version_c (trex_file, str, max_str_len) &
           bind(C, name="trexio_write_metadata_package_version")
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     character, intent(in) :: str(*)
     integer(4), intent(in), value :: max_str_len
   end function trexio_write_metadata_package_version_c
end interface

interface
   integer function trexio_write_metadata_description_c (trex_file, str, max_str_len) &
           bind(C, name="trexio_write_metadata_description")
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     character, intent(in) :: str(*)
     integer(4), intent(in), value :: max_str_len
   end function trexio_write_metadata_description_c
end interface

interface
   integer function trexio_write_nucleus_point_group_c (trex_file, str, max_str_len) &
           bind(C, name="trexio_write_nucleus_point_group")
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     character, intent(in) :: str(*)
     integer(4), intent(in), value :: max_str_len
   end function trexio_write_nucleus_point_group_c
end interface

interface
   integer function trexio_write_basis_type_c (trex_file, str, max_str_len) &
           bind(C, name="trexio_write_basis_type")
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     character, intent(in) :: str(*)
     integer(4), intent(in), value :: max_str_len
   end function trexio_write_basis_type_c
end interface

interface
   integer function trexio_write_mo_type_c (trex_file, str, max_str_len) &
           bind(C, name="trexio_write_mo_type")
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     character, intent(in) :: str(*)
     integer(4), intent(in), value :: max_str_len
   end function trexio_write_mo_type_c
end interface

interface
   integer function trexio_write_nucleus_charge_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(4), intent(in) :: dset(*)
   end function trexio_write_nucleus_charge_32
end interface

interface
   integer function trexio_write_nucleus_coord_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(4), intent(in) :: dset(*)
   end function trexio_write_nucleus_coord_32
end interface

interface
   integer function trexio_write_ecp_max_ang_mom_plus_1_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     integer(4), intent(in) :: dset(*)
   end function trexio_write_ecp_max_ang_mom_plus_1_32
end interface

interface
   integer function trexio_write_ecp_z_core_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     integer(4), intent(in) :: dset(*)
   end function trexio_write_ecp_z_core_32
end interface

interface
   integer function trexio_write_ecp_ang_mom_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     integer(4), intent(in) :: dset(*)
   end function trexio_write_ecp_ang_mom_32
end interface

interface
   integer function trexio_write_ecp_nucleus_index_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     integer(4), intent(in) :: dset(*)
   end function trexio_write_ecp_nucleus_index_32
end interface

interface
   integer function trexio_write_ecp_exponent_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(4), intent(in) :: dset(*)
   end function trexio_write_ecp_exponent_32
end interface

interface
   integer function trexio_write_ecp_coefficient_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(4), intent(in) :: dset(*)
   end function trexio_write_ecp_coefficient_32
end interface

interface
   integer function trexio_write_ecp_power_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     integer(4), intent(in) :: dset(*)
   end function trexio_write_ecp_power_32
end interface

interface
   integer function trexio_write_basis_nucleus_index_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     integer(4), intent(in) :: dset(*)
   end function trexio_write_basis_nucleus_index_32
end interface

interface
   integer function trexio_write_basis_ang_mom_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     integer(4), intent(in) :: dset(*)
   end function trexio_write_basis_ang_mom_32
end interface

interface
   integer function trexio_write_basis_shell_factor_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(4), intent(in) :: dset(*)
   end function trexio_write_basis_shell_factor_32
end interface

interface
   integer function trexio_write_basis_shell_index_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     integer(4), intent(in) :: dset(*)
   end function trexio_write_basis_shell_index_32
end interface

interface
   integer function trexio_write_basis_exponent_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(4), intent(in) :: dset(*)
   end function trexio_write_basis_exponent_32
end interface

interface
   integer function trexio_write_basis_coefficient_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(4), intent(in) :: dset(*)
   end function trexio_write_basis_coefficient_32
end interface

interface
   integer function trexio_write_basis_prim_factor_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(4), intent(in) :: dset(*)
   end function trexio_write_basis_prim_factor_32
end interface

interface
   integer function trexio_write_ao_shell_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     integer(4), intent(in) :: dset(*)
   end function trexio_write_ao_shell_32
end interface

interface
   integer function trexio_write_ao_normalization_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(4), intent(in) :: dset(*)
   end function trexio_write_ao_normalization_32
end interface

interface
   integer function trexio_write_ao_1e_int_overlap_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(4), intent(in) :: dset(*)
   end function trexio_write_ao_1e_int_overlap_32
end interface

interface
   integer function trexio_write_ao_1e_int_kinetic_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(4), intent(in) :: dset(*)
   end function trexio_write_ao_1e_int_kinetic_32
end interface

interface
   integer function trexio_write_ao_1e_int_potential_n_e_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(4), intent(in) :: dset(*)
   end function trexio_write_ao_1e_int_potential_n_e_32
end interface

interface
   integer function trexio_write_ao_1e_int_ecp_local_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(4), intent(in) :: dset(*)
   end function trexio_write_ao_1e_int_ecp_local_32
end interface

interface
   integer function trexio_write_ao_1e_int_ecp_non_local_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(4), intent(in) :: dset(*)
   end function trexio_write_ao_1e_int_ecp_non_local_32
end interface

interface
   integer function trexio_write_ao_1e_int_core_hamiltonian_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(4), intent(in) :: dset(*)
   end function trexio_write_ao_1e_int_core_hamiltonian_32
end interface

interface
   integer function trexio_write_ao_2e_int_eri_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(4), intent(in) :: dset(*)
   end function trexio_write_ao_2e_int_eri_32
end interface

interface
   integer function trexio_write_ao_2e_int_eri_lr_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(4), intent(in) :: dset(*)
   end function trexio_write_ao_2e_int_eri_lr_32
end interface

interface
   integer function trexio_write_mo_coefficient_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(4), intent(in) :: dset(*)
   end function trexio_write_mo_coefficient_32
end interface

interface
   integer function trexio_write_mo_occupation_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(4), intent(in) :: dset(*)
   end function trexio_write_mo_occupation_32
end interface

interface
   integer function trexio_write_mo_1e_int_overlap_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(4), intent(in) :: dset(*)
   end function trexio_write_mo_1e_int_overlap_32
end interface

interface
   integer function trexio_write_mo_1e_int_kinetic_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(4), intent(in) :: dset(*)
   end function trexio_write_mo_1e_int_kinetic_32
end interface

interface
   integer function trexio_write_mo_1e_int_potential_n_e_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(4), intent(in) :: dset(*)
   end function trexio_write_mo_1e_int_potential_n_e_32
end interface

interface
   integer function trexio_write_mo_1e_int_ecp_local_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(4), intent(in) :: dset(*)
   end function trexio_write_mo_1e_int_ecp_local_32
end interface

interface
   integer function trexio_write_mo_1e_int_ecp_non_local_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(4), intent(in) :: dset(*)
   end function trexio_write_mo_1e_int_ecp_non_local_32
end interface

interface
   integer function trexio_write_mo_1e_int_core_hamiltonian_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(4), intent(in) :: dset(*)
   end function trexio_write_mo_1e_int_core_hamiltonian_32
end interface

interface
   integer function trexio_write_mo_2e_int_eri_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(4), intent(in) :: dset(*)
   end function trexio_write_mo_2e_int_eri_32
end interface

interface
   integer function trexio_write_mo_2e_int_eri_lr_32 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(4), intent(in) :: dset(*)
   end function trexio_write_mo_2e_int_eri_lr_32
end interface

interface
   integer function trexio_write_nucleus_charge_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(8), intent(in) :: dset(*)
   end function trexio_write_nucleus_charge_64
end interface

interface
   integer function trexio_write_nucleus_coord_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(8), intent(in) :: dset(*)
   end function trexio_write_nucleus_coord_64
end interface

interface
   integer function trexio_write_ecp_max_ang_mom_plus_1_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     integer(8), intent(in) :: dset(*)
   end function trexio_write_ecp_max_ang_mom_plus_1_64
end interface

interface
   integer function trexio_write_ecp_z_core_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     integer(8), intent(in) :: dset(*)
   end function trexio_write_ecp_z_core_64
end interface

interface
   integer function trexio_write_ecp_ang_mom_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     integer(8), intent(in) :: dset(*)
   end function trexio_write_ecp_ang_mom_64
end interface

interface
   integer function trexio_write_ecp_nucleus_index_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     integer(8), intent(in) :: dset(*)
   end function trexio_write_ecp_nucleus_index_64
end interface

interface
   integer function trexio_write_ecp_exponent_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(8), intent(in) :: dset(*)
   end function trexio_write_ecp_exponent_64
end interface

interface
   integer function trexio_write_ecp_coefficient_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(8), intent(in) :: dset(*)
   end function trexio_write_ecp_coefficient_64
end interface

interface
   integer function trexio_write_ecp_power_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     integer(8), intent(in) :: dset(*)
   end function trexio_write_ecp_power_64
end interface

interface
   integer function trexio_write_basis_nucleus_index_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     integer(8), intent(in) :: dset(*)
   end function trexio_write_basis_nucleus_index_64
end interface

interface
   integer function trexio_write_basis_ang_mom_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     integer(8), intent(in) :: dset(*)
   end function trexio_write_basis_ang_mom_64
end interface

interface
   integer function trexio_write_basis_shell_factor_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(8), intent(in) :: dset(*)
   end function trexio_write_basis_shell_factor_64
end interface

interface
   integer function trexio_write_basis_shell_index_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     integer(8), intent(in) :: dset(*)
   end function trexio_write_basis_shell_index_64
end interface

interface
   integer function trexio_write_basis_exponent_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(8), intent(in) :: dset(*)
   end function trexio_write_basis_exponent_64
end interface

interface
   integer function trexio_write_basis_coefficient_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(8), intent(in) :: dset(*)
   end function trexio_write_basis_coefficient_64
end interface

interface
   integer function trexio_write_basis_prim_factor_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(8), intent(in) :: dset(*)
   end function trexio_write_basis_prim_factor_64
end interface

interface
   integer function trexio_write_ao_shell_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     integer(8), intent(in) :: dset(*)
   end function trexio_write_ao_shell_64
end interface

interface
   integer function trexio_write_ao_normalization_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(8), intent(in) :: dset(*)
   end function trexio_write_ao_normalization_64
end interface

interface
   integer function trexio_write_ao_1e_int_overlap_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(8), intent(in) :: dset(*)
   end function trexio_write_ao_1e_int_overlap_64
end interface

interface
   integer function trexio_write_ao_1e_int_kinetic_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(8), intent(in) :: dset(*)
   end function trexio_write_ao_1e_int_kinetic_64
end interface

interface
   integer function trexio_write_ao_1e_int_potential_n_e_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(8), intent(in) :: dset(*)
   end function trexio_write_ao_1e_int_potential_n_e_64
end interface

interface
   integer function trexio_write_ao_1e_int_ecp_local_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(8), intent(in) :: dset(*)
   end function trexio_write_ao_1e_int_ecp_local_64
end interface

interface
   integer function trexio_write_ao_1e_int_ecp_non_local_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(8), intent(in) :: dset(*)
   end function trexio_write_ao_1e_int_ecp_non_local_64
end interface

interface
   integer function trexio_write_ao_1e_int_core_hamiltonian_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(8), intent(in) :: dset(*)
   end function trexio_write_ao_1e_int_core_hamiltonian_64
end interface

interface
   integer function trexio_write_ao_2e_int_eri_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(8), intent(in) :: dset(*)
   end function trexio_write_ao_2e_int_eri_64
end interface

interface
   integer function trexio_write_ao_2e_int_eri_lr_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(8), intent(in) :: dset(*)
   end function trexio_write_ao_2e_int_eri_lr_64
end interface

interface
   integer function trexio_write_mo_coefficient_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(8), intent(in) :: dset(*)
   end function trexio_write_mo_coefficient_64
end interface

interface
   integer function trexio_write_mo_occupation_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(8), intent(in) :: dset(*)
   end function trexio_write_mo_occupation_64
end interface

interface
   integer function trexio_write_mo_1e_int_overlap_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(8), intent(in) :: dset(*)
   end function trexio_write_mo_1e_int_overlap_64
end interface

interface
   integer function trexio_write_mo_1e_int_kinetic_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(8), intent(in) :: dset(*)
   end function trexio_write_mo_1e_int_kinetic_64
end interface

interface
   integer function trexio_write_mo_1e_int_potential_n_e_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(8), intent(in) :: dset(*)
   end function trexio_write_mo_1e_int_potential_n_e_64
end interface

interface
   integer function trexio_write_mo_1e_int_ecp_local_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(8), intent(in) :: dset(*)
   end function trexio_write_mo_1e_int_ecp_local_64
end interface

interface
   integer function trexio_write_mo_1e_int_ecp_non_local_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(8), intent(in) :: dset(*)
   end function trexio_write_mo_1e_int_ecp_non_local_64
end interface

interface
   integer function trexio_write_mo_1e_int_core_hamiltonian_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(8), intent(in) :: dset(*)
   end function trexio_write_mo_1e_int_core_hamiltonian_64
end interface

interface
   integer function trexio_write_mo_2e_int_eri_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(8), intent(in) :: dset(*)
   end function trexio_write_mo_2e_int_eri_64
end interface

interface
   integer function trexio_write_mo_2e_int_eri_lr_64 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(8), intent(in) :: dset(*)
   end function trexio_write_mo_2e_int_eri_lr_64
end interface

interface
   integer function trexio_write_nucleus_charge (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(8), intent(in) :: dset(*)
   end function trexio_write_nucleus_charge
end interface

interface
   integer function trexio_write_nucleus_coord (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(8), intent(in) :: dset(*)
   end function trexio_write_nucleus_coord
end interface

interface
   integer function trexio_write_ecp_max_ang_mom_plus_1 (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     integer(4), intent(in) :: dset(*)
   end function trexio_write_ecp_max_ang_mom_plus_1
end interface

interface
   integer function trexio_write_ecp_z_core (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     integer(4), intent(in) :: dset(*)
   end function trexio_write_ecp_z_core
end interface

interface
   integer function trexio_write_ecp_ang_mom (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     integer(4), intent(in) :: dset(*)
   end function trexio_write_ecp_ang_mom
end interface

interface
   integer function trexio_write_ecp_nucleus_index (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     integer(4), intent(in) :: dset(*)
   end function trexio_write_ecp_nucleus_index
end interface

interface
   integer function trexio_write_ecp_exponent (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(8), intent(in) :: dset(*)
   end function trexio_write_ecp_exponent
end interface

interface
   integer function trexio_write_ecp_coefficient (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(8), intent(in) :: dset(*)
   end function trexio_write_ecp_coefficient
end interface

interface
   integer function trexio_write_ecp_power (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     integer(4), intent(in) :: dset(*)
   end function trexio_write_ecp_power
end interface

interface
   integer function trexio_write_basis_nucleus_index (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     integer(4), intent(in) :: dset(*)
   end function trexio_write_basis_nucleus_index
end interface

interface
   integer function trexio_write_basis_ang_mom (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     integer(4), intent(in) :: dset(*)
   end function trexio_write_basis_ang_mom
end interface

interface
   integer function trexio_write_basis_shell_factor (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(8), intent(in) :: dset(*)
   end function trexio_write_basis_shell_factor
end interface

interface
   integer function trexio_write_basis_shell_index (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     integer(4), intent(in) :: dset(*)
   end function trexio_write_basis_shell_index
end interface

interface
   integer function trexio_write_basis_exponent (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(8), intent(in) :: dset(*)
   end function trexio_write_basis_exponent
end interface

interface
   integer function trexio_write_basis_coefficient (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(8), intent(in) :: dset(*)
   end function trexio_write_basis_coefficient
end interface

interface
   integer function trexio_write_basis_prim_factor (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(8), intent(in) :: dset(*)
   end function trexio_write_basis_prim_factor
end interface

interface
   integer function trexio_write_ao_shell (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     integer(4), intent(in) :: dset(*)
   end function trexio_write_ao_shell
end interface

interface
   integer function trexio_write_ao_normalization (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(8), intent(in) :: dset(*)
   end function trexio_write_ao_normalization
end interface

interface
   integer function trexio_write_ao_1e_int_overlap (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(8), intent(in) :: dset(*)
   end function trexio_write_ao_1e_int_overlap
end interface

interface
   integer function trexio_write_ao_1e_int_kinetic (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(8), intent(in) :: dset(*)
   end function trexio_write_ao_1e_int_kinetic
end interface

interface
   integer function trexio_write_ao_1e_int_potential_n_e (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(8), intent(in) :: dset(*)
   end function trexio_write_ao_1e_int_potential_n_e
end interface

interface
   integer function trexio_write_ao_1e_int_ecp_local (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(8), intent(in) :: dset(*)
   end function trexio_write_ao_1e_int_ecp_local
end interface

interface
   integer function trexio_write_ao_1e_int_ecp_non_local (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(8), intent(in) :: dset(*)
   end function trexio_write_ao_1e_int_ecp_non_local
end interface

interface
   integer function trexio_write_ao_1e_int_core_hamiltonian (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(8), intent(in) :: dset(*)
   end function trexio_write_ao_1e_int_core_hamiltonian
end interface

interface
   integer function trexio_write_ao_2e_int_eri (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(8), intent(in) :: dset(*)
   end function trexio_write_ao_2e_int_eri
end interface

interface
   integer function trexio_write_ao_2e_int_eri_lr (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(8), intent(in) :: dset(*)
   end function trexio_write_ao_2e_int_eri_lr
end interface

interface
   integer function trexio_write_mo_coefficient (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(8), intent(in) :: dset(*)
   end function trexio_write_mo_coefficient
end interface

interface
   integer function trexio_write_mo_occupation (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(8), intent(in) :: dset(*)
   end function trexio_write_mo_occupation
end interface

interface
   integer function trexio_write_mo_1e_int_overlap (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(8), intent(in) :: dset(*)
   end function trexio_write_mo_1e_int_overlap
end interface

interface
   integer function trexio_write_mo_1e_int_kinetic (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(8), intent(in) :: dset(*)
   end function trexio_write_mo_1e_int_kinetic
end interface

interface
   integer function trexio_write_mo_1e_int_potential_n_e (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(8), intent(in) :: dset(*)
   end function trexio_write_mo_1e_int_potential_n_e
end interface

interface
   integer function trexio_write_mo_1e_int_ecp_local (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(8), intent(in) :: dset(*)
   end function trexio_write_mo_1e_int_ecp_local
end interface

interface
   integer function trexio_write_mo_1e_int_ecp_non_local (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(8), intent(in) :: dset(*)
   end function trexio_write_mo_1e_int_ecp_non_local
end interface

interface
   integer function trexio_write_mo_1e_int_core_hamiltonian (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(8), intent(in) :: dset(*)
   end function trexio_write_mo_1e_int_core_hamiltonian
end interface

interface
   integer function trexio_write_mo_2e_int_eri (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(8), intent(in) :: dset(*)
   end function trexio_write_mo_2e_int_eri
end interface

interface
   integer function trexio_write_mo_2e_int_eri_lr (trex_file, dset) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     real(8), intent(in) :: dset(*)
   end function trexio_write_mo_2e_int_eri_lr
end interface

interface
   integer function trexio_write_metadata_code_low (trex_file, dset, max_str_len) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     character, intent(in) :: dset(*)
     integer(4), intent(in), value :: max_str_len
   end function trexio_write_metadata_code_low
end interface

interface
   integer function trexio_write_metadata_author_low (trex_file, dset, max_str_len) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     character, intent(in) :: dset(*)
     integer(4), intent(in), value :: max_str_len
   end function trexio_write_metadata_author_low
end interface

interface
   integer function trexio_write_nucleus_label_low (trex_file, dset, max_str_len) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     character, intent(in) :: dset(*)
     integer(4), intent(in), value :: max_str_len
   end function trexio_write_nucleus_label_low
end interface

interface
   integer function trexio_write_mo_class_low (trex_file, dset, max_str_len) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     character, intent(in) :: dset(*)
     integer(4), intent(in), value :: max_str_len
   end function trexio_write_mo_class_low
end interface

interface
   integer function trexio_write_mo_symmetry_low (trex_file, dset, max_str_len) bind(C)
     use, intrinsic :: iso_c_binding
     integer(8), intent(in), value :: trex_file
     character, intent(in) :: dset(*)
     integer(4), intent(in), value :: max_str_len
   end function trexio_write_mo_symmetry_low
end interface

contains
   integer(8) function trexio_open (filename, mode, backend, rc_open)
     use, intrinsic :: iso_c_binding, only : c_null_char
     implicit none
     character(len=*), intent(in)                 :: filename
     character, intent(in), value                 :: mode
     integer(trexio_backend), intent(in), value   :: backend
     integer(trexio_exit_code), intent(out)       :: rc_open
     character(len=len_trim(filename)+1) :: filename_c
     integer(trexio_exit_code) :: rc

     filename_c = trim(filename) // c_null_char
     trexio_open = trexio_open_c(filename_c, mode, backend, rc_open)
     if (trexio_open == 0_8 .or. rc_open /= TREXIO_SUCCESS) then
       return
     endif
     rc = trexio_set_one_based(trexio_open)
     if (rc /= TREXIO_SUCCESS) then
        rc = trexio_close(trexio_open)
        trexio_open = 0_8
     endif
   end function trexio_open

subroutine trexio_strarray2str(str_array, max_num_str, max_len_str, str_res)
  use, intrinsic :: iso_c_binding, only : c_null_char
  implicit none

  integer(8), intent(in), value   :: max_num_str  ! number of elements in strign array
  integer, intent(in), value   :: max_len_str  ! maximum length of a string in an array
  character(len=*), intent(in)  :: str_array(*)
  character(len=:), allocatable, intent(out) :: str_res
  integer :: i

  str_res = ''
  do i = 1, max_num_str
    str_res = str_res // trim(str_array(i)) // TREXIO_DELIM
  enddo
  str_res = str_res // c_null_char

end subroutine trexio_strarray2str

subroutine trexio_str2strarray(str_flat, max_num_str, max_len_str, str_array)
  implicit none

  integer(8), intent(in), value   :: max_num_str  ! number of elements in strign array
  integer, intent(in), value   :: max_len_str  ! maximum length of a string in an array
  character, intent(in) :: str_flat(*)
  character(len=*), intent(inout)  :: str_array(*)

  character(len=max_len_str)  :: tmp_str
  integer :: i, j, k, ind, offset
  integer(8) :: len_flat

  len_flat = (max_len_str+1)*max_num_str + 1

  ind=1
  offset=1
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
    offset=ind
  enddo

end subroutine trexio_str2strarray

subroutine trexio_assert(trexio_rc, check_rc, success_message)
  implicit none

  integer, intent(in), value   :: trexio_rc
  integer, intent(in), value   :: check_rc
  character(len=*), intent(in), optional  :: success_message

  character*(128) :: str

  if (trexio_rc == check_rc) then
    if (present(success_message)) write(*,*) success_message
  else
    call trexio_string_of_error(trexio_rc, str)
    print *, trim(str)
    call exit(1)
  endif

end subroutine trexio_assert
integer function trexio_read_metadata_package_version (trex_file, str, max_str_len)
  implicit none
  integer(8), intent(in), value :: trex_file
  integer(4), intent(in), value :: max_str_len
  character, intent(out) :: str(*)

  trexio_read_metadata_package_version = trexio_read_metadata_package_version_c(trex_file, str, max_str_len)

end function trexio_read_metadata_package_version

integer function trexio_read_metadata_description (trex_file, str, max_str_len)
  implicit none
  integer(8), intent(in), value :: trex_file
  integer(4), intent(in), value :: max_str_len
  character, intent(out) :: str(*)

  trexio_read_metadata_description = trexio_read_metadata_description_c(trex_file, str, max_str_len)

end function trexio_read_metadata_description

integer function trexio_read_nucleus_point_group (trex_file, str, max_str_len)
  implicit none
  integer(8), intent(in), value :: trex_file
  integer(4), intent(in), value :: max_str_len
  character, intent(out) :: str(*)

  trexio_read_nucleus_point_group = trexio_read_nucleus_point_group_c(trex_file, str, max_str_len)

end function trexio_read_nucleus_point_group

integer function trexio_read_basis_type (trex_file, str, max_str_len)
  implicit none
  integer(8), intent(in), value :: trex_file
  integer(4), intent(in), value :: max_str_len
  character, intent(out) :: str(*)

  trexio_read_basis_type = trexio_read_basis_type_c(trex_file, str, max_str_len)

end function trexio_read_basis_type

integer function trexio_read_mo_type (trex_file, str, max_str_len)
  implicit none
  integer(8), intent(in), value :: trex_file
  integer(4), intent(in), value :: max_str_len
  character, intent(out) :: str(*)

  trexio_read_mo_type = trexio_read_mo_type_c(trex_file, str, max_str_len)

end function trexio_read_mo_type

integer function trexio_read_metadata_code (trex_file, dset, max_str_len)
  implicit none
  integer(8), intent(in), value :: trex_file
  integer(4), intent(in), value :: max_str_len
  character(len=*), intent(inout) :: dset(*)

  character, allocatable :: str_compiled(:)
  integer(8) :: metadata_code_num
  integer :: rc

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

integer function trexio_read_metadata_author (trex_file, dset, max_str_len)
  implicit none
  integer(8), intent(in), value :: trex_file
  integer(4), intent(in), value :: max_str_len
  character(len=*), intent(inout) :: dset(*)

  character, allocatable :: str_compiled(:)
  integer(8) :: metadata_author_num
  integer :: rc

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

integer function trexio_read_nucleus_label (trex_file, dset, max_str_len)
  implicit none
  integer(8), intent(in), value :: trex_file
  integer(4), intent(in), value :: max_str_len
  character(len=*), intent(inout) :: dset(*)

  character, allocatable :: str_compiled(:)
  integer(8) :: nucleus_num
  integer :: rc

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

integer function trexio_read_mo_class (trex_file, dset, max_str_len)
  implicit none
  integer(8), intent(in), value :: trex_file
  integer(4), intent(in), value :: max_str_len
  character(len=*), intent(inout) :: dset(*)

  character, allocatable :: str_compiled(:)
  integer(8) :: mo_num
  integer :: rc

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

integer function trexio_read_mo_symmetry (trex_file, dset, max_str_len)
  implicit none
  integer(8), intent(in), value :: trex_file
  integer(4), intent(in), value :: max_str_len
  character(len=*), intent(inout) :: dset(*)

  character, allocatable :: str_compiled(:)
  integer(8) :: mo_num
  integer :: rc

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

integer function trexio_write_metadata_package_version (trex_file, str, max_str_len)
  use, intrinsic :: iso_c_binding, only : c_null_char
  implicit none
  integer(8), intent(in), value :: trex_file
  integer(4), intent(in), value :: max_str_len
  character(len=*), intent(in) :: str

  character(len=len_trim(str)+1) :: str_c

  str_c = trim(str) // c_null_char

  trexio_write_metadata_package_version = trexio_write_metadata_package_version_c(trex_file, str_c, max_str_len)

end function trexio_write_metadata_package_version

integer function trexio_write_metadata_description (trex_file, str, max_str_len)
  use, intrinsic :: iso_c_binding, only : c_null_char
  implicit none
  integer(8), intent(in), value :: trex_file
  integer(4), intent(in), value :: max_str_len
  character(len=*), intent(in) :: str

  character(len=len_trim(str)+1) :: str_c

  str_c = trim(str) // c_null_char

  trexio_write_metadata_description = trexio_write_metadata_description_c(trex_file, str_c, max_str_len)

end function trexio_write_metadata_description

integer function trexio_write_nucleus_point_group (trex_file, str, max_str_len)
  use, intrinsic :: iso_c_binding, only : c_null_char
  implicit none
  integer(8), intent(in), value :: trex_file
  integer(4), intent(in), value :: max_str_len
  character(len=*), intent(in) :: str

  character(len=len_trim(str)+1) :: str_c

  str_c = trim(str) // c_null_char

  trexio_write_nucleus_point_group = trexio_write_nucleus_point_group_c(trex_file, str_c, max_str_len)

end function trexio_write_nucleus_point_group

integer function trexio_write_basis_type (trex_file, str, max_str_len)
  use, intrinsic :: iso_c_binding, only : c_null_char
  implicit none
  integer(8), intent(in), value :: trex_file
  integer(4), intent(in), value :: max_str_len
  character(len=*), intent(in) :: str

  character(len=len_trim(str)+1) :: str_c

  str_c = trim(str) // c_null_char

  trexio_write_basis_type = trexio_write_basis_type_c(trex_file, str_c, max_str_len)

end function trexio_write_basis_type

integer function trexio_write_mo_type (trex_file, str, max_str_len)
  use, intrinsic :: iso_c_binding, only : c_null_char
  implicit none
  integer(8), intent(in), value :: trex_file
  integer(4), intent(in), value :: max_str_len
  character(len=*), intent(in) :: str

  character(len=len_trim(str)+1) :: str_c

  str_c = trim(str) // c_null_char

  trexio_write_mo_type = trexio_write_mo_type_c(trex_file, str_c, max_str_len)

end function trexio_write_mo_type

integer function trexio_write_metadata_code (trex_file, dset, max_str_len)
  implicit none
  integer(8), intent(in), value :: trex_file
  integer(4), intent(in), value :: max_str_len
  character(len=*), intent(in) :: dset(*)

  character(len=:), allocatable :: str_compiled
  integer(8) :: metadata_code_num
  integer :: rc

  rc = trexio_read_metadata_code_num_64(trex_file, metadata_code_num)
  if (rc /= TREXIO_SUCCESS) then
    trexio_write_metadata_code = rc
  else
    call trexio_strarray2str(dset, metadata_code_num, max_str_len, str_compiled)
    trexio_write_metadata_code = trexio_write_metadata_code_low(trex_file, str_compiled, max_str_len)
  endif

end function trexio_write_metadata_code

integer function trexio_write_metadata_author (trex_file, dset, max_str_len)
  implicit none
  integer(8), intent(in), value :: trex_file
  integer(4), intent(in), value :: max_str_len
  character(len=*), intent(in) :: dset(*)

  character(len=:), allocatable :: str_compiled
  integer(8) :: metadata_author_num
  integer :: rc

  rc = trexio_read_metadata_author_num_64(trex_file, metadata_author_num)
  if (rc /= TREXIO_SUCCESS) then
    trexio_write_metadata_author = rc
  else
    call trexio_strarray2str(dset, metadata_author_num, max_str_len, str_compiled)
    trexio_write_metadata_author = trexio_write_metadata_author_low(trex_file, str_compiled, max_str_len)
  endif

end function trexio_write_metadata_author

integer function trexio_write_nucleus_label (trex_file, dset, max_str_len)
  implicit none
  integer(8), intent(in), value :: trex_file
  integer(4), intent(in), value :: max_str_len
  character(len=*), intent(in) :: dset(*)

  character(len=:), allocatable :: str_compiled
  integer(8) :: nucleus_num
  integer :: rc

  rc = trexio_read_nucleus_num_64(trex_file, nucleus_num)
  if (rc /= TREXIO_SUCCESS) then
    trexio_write_nucleus_label = rc
  else
    call trexio_strarray2str(dset, nucleus_num, max_str_len, str_compiled)
    trexio_write_nucleus_label = trexio_write_nucleus_label_low(trex_file, str_compiled, max_str_len)
  endif

end function trexio_write_nucleus_label

integer function trexio_write_mo_class (trex_file, dset, max_str_len)
  implicit none
  integer(8), intent(in), value :: trex_file
  integer(4), intent(in), value :: max_str_len
  character(len=*), intent(in) :: dset(*)

  character(len=:), allocatable :: str_compiled
  integer(8) :: mo_num
  integer :: rc

  rc = trexio_read_mo_num_64(trex_file, mo_num)
  if (rc /= TREXIO_SUCCESS) then
    trexio_write_mo_class = rc
  else
    call trexio_strarray2str(dset, mo_num, max_str_len, str_compiled)
    trexio_write_mo_class = trexio_write_mo_class_low(trex_file, str_compiled, max_str_len)
  endif

end function trexio_write_mo_class

integer function trexio_write_mo_symmetry (trex_file, dset, max_str_len)
  implicit none
  integer(8), intent(in), value :: trex_file
  integer(4), intent(in), value :: max_str_len
  character(len=*), intent(in) :: dset(*)

  character(len=:), allocatable :: str_compiled
  integer(8) :: mo_num
  integer :: rc

  rc = trexio_read_mo_num_64(trex_file, mo_num)
  if (rc /= TREXIO_SUCCESS) then
    trexio_write_mo_symmetry = rc
  else
    call trexio_strarray2str(dset, mo_num, max_str_len, str_compiled)
    trexio_write_mo_symmetry = trexio_write_mo_symmetry_low(trex_file, str_compiled, max_str_len)
  endif

end function trexio_write_mo_symmetry

end module trexio
