!> @brief TREXIO library interface wrapper.
!> @author CHAMP developers
!> @date 2025
!>
!> @details This file conditionally includes the TREXIO Fortran interface
!> when the TREXIO library is found during compilation. TREXIO (TREX I/O)
!> is a file format and library for storing wave function data from quantum
!> chemistry calculations.
!>
!> The TREXIO library provides:
!> - Standardized format for quantum chemistry data exchange
!> - Support for basis sets, molecular orbitals, and electron integrals
!> - Efficient I/O operations for large quantum chemistry datasets
!>
!> @note This file requires the TREXIO_FOUND preprocessor flag to be defined
!> during compilation to include the TREXIO Fortran interface.
!>
!> @see https://github.com/TREX-CoE/trexio for more information about TREXIO.

#ifdef TREXIO_FOUND
#include "trexio_f.f90"
#endif
