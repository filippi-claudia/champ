---
tags:
    - dependencies
---

# Dependencies

## Required Libraries / Tools

CHAMP requires the following libraries to compile and run:

-  [x] **CMake** >= 3.17
   - Required for building CHAMP

-  [x] **Fortran/C Compiler**
   - gfortran/gcc >= 9.3.0 or Intel Fortran 2020 onwards (ifort/ifx)

-  [x] **MPI** - Message Passing Interface library
   - OpenMPI >= 3.0 or Intel MPI
   - Required for parallel execution

-  [x] **BLAS/LAPACK** - Linear algebra libraries
   - BLAS/LAPACK or Intel MKL
   - Required for matrix operations

## Included Libraries

CHAMP ships with the following libraries:

-  [x] **libFDF** 
   - Input file parser library based on [mpi-libfdf-parser](https://github.com/neelravi/mpi-libfdf-parser)
   - Handles keyword-value pair parsing in input files
   - No separate installation required

## Optional Libraries / Tools

The following libraries are optional and enable additional features:

-  [ ] **TREXIO** >= 2.0.0
   - For reading wavefunction data from TREXIO HDF5 files
   - See [TREXIO installation guide](trexio.md)

-  [ ] **QMCkl** >= 0.2.1
   - High-performance library for quantum Monte Carlo calculation kernels
   - See [QMCkl installation guide](qmckl.md)

-  [ ] **Doxygen**
   - For generating code developer's documentation

## Starting Wavefunction Requirements

CHAMP requires a starting wavefunction and system information to begin QMC calculations. These can be obtained from various quantum chemistry programs:

- **TREXIO files** (recommended) - from Quantum Package, PySCF, or other TREXIO-ready codes
- **GAMESS** - via conversion tools in CHAMP's tools directory
- **Quantum Package** - via `qp run export_trexio`

Python tools are provided in CHAMP's `tools/` directory to convert and extract wavefunction data from these sources into CHAMP-compatible formats.

