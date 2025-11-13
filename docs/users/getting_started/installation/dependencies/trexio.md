---
title: TREXIO
tags:
    - TREXIO
---


!!! warning
    The version of the the library must be greater than 2.0.0.


## Minimal requirements for the installation of trexio (for users):

- Autotools             (autoconf >= 2.69, automake >= 1.11, libtool >= 2.2) or CMake (>= 3.16)
- C compiler            (gcc/icc/clang)
- Fortran compiler      (gfortran/ifort)
- HDF5 library          (>= 1.8) [optional, recommended for high performance]


## Installation procedure from the tarball (for users):

1. Download the `trexio-<version>.tar.gz` file
2. `gzip -cd trexio-<version>.tar.gz | tar xvf -`
3. `cd trexio-<version>`
4. `./configure`
5. `make`
6. `make check`
7. `sudo make install`



## Installation procedure (for CMAKE users):

[CMake](https://cmake.org) users can achieve the same with the following steps (an example of out-of-source build):

1. `cmake -S. -Bbuild`
2. `cd build`
3. `make`
4. `ctest` (or `make test`)
5. `sudo make install`
