---
tags:
    - dependencies
---

# CHAMP Structure and Dependencies

CHAMP relies on various other program packages:

1. [Parser](https://github.com/neelravi/mpi-libfdf-parser):
   An easy-to-use and easy-to-extend keyword-value pair based input file parser written in Fortran 2008.  This parser uses a heavily modified libFDF library and is written by [Ravindra Shinde](https://github.com/neelravi). It can parse keyword-value pairs, blocks of data, and general variables with different physical units in an order-independent manner. Our implementation can handle multiple data types and file formats. The parser is kept as a library in the code, however, it can be easily adapted by any other Fortran-based code.

2. GAMESS:
   For finite systems the starting wavefunction is obtained from the
   quantum chemistry program GAMESS, written by Mike Schmidt and
   collaborators at Iowa State University.

3. GAMESS_Interface:
   The wavefunction produced by GAMESS has to be cast in a form
   suitable for input to CHAMP.  This is a lot more work than first meets
   the eye. We provide a python package inside the CHAMP's tool directory to extract all the necessary information needed from a GAMESS calculation. The tool can also extract information from a TREXIO file in the hdf5 file format. This utility is written by [Ravindra Shinde](https://github.com/neelravi).

4. MOLCAS Interface:
   A python package qc2champ can be used to convert a    MOLCAS or an openMOLCAS calculation into the input files needed by CHAMP. This package is written by [Ravindra Shinde](https://github.com/neelravi). An independent version of the convertor script was added thanks to Csaba Daday and Monika Dash.

5. TREXIO (optional):
    The TREXIO library is a C library for reading and writing the
    TREXIO file format. The TREXIO file format is a HDF5 file format
    for storing the electronic wavefunctions.

6. QMCkl (optional):
    QMCkl is a high-performance library for executing common quantum Monte Carlo calculations kernels.

# Requirements
1. cmake >= 3.17
2. gfortran/gcc >= 9.3.0 or Intel Fortran 2020 onwards
3. BLAS/LAPACK or Intel MKL
4. openMPI >= 3.0 or Intel MPI
5. [Optional] TREXIO library >= 2.0.0
6. [Optional] QMCkl library >= 0.2.1
7. [Optional] doxygen (for documentation)

