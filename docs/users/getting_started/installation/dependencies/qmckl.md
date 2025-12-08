---
title: "QMCkl"
tags:
    - QMCKL
    - optional
---

# QMCkl Library

QMCkl is a high-performance library that provides optimized implementations of common Quantum Monte Carlo calculation kernels. When enabled, CHAMP can take benefit of QMCkl's hardware-optimized routines for improved performance.

!!! warning
    **Minimum version required:** 1.0.0 or higher

## What is QMCkl?

QMCkl exposes the main QMC algorithms through a standard API with portable implementations and hardware-specific optimizations. It provides:

- High-performance kernel implementations
- Hardware acceleration (CPU vectorization, GPU support)
- Standardized API for QMC operations
- Extensive test suite for validation

## Requirements

Before installing QMCkl, ensure you have:

- **C compiler** (gcc, icc, or clang)
- **Fortran compiler** (gfortran, ifort, or compatible)
- **Autotools** (autoconf, automake, libtool)
- **TREXIO library** (QMCkl requires TREXIO to be installed first)

**Note:** TREXIO is a prerequisite for QMCkl. Install TREXIO before proceeding (see [TREXIO installation guide](trexio.md)).

## Installation Methods

### Method 1: From Release Tarball

**Download the latest release from** [QMCkl GitHub Releases](https://github.com/TREX-CoE/qmckl/releases):

```bash
wget https://github.com/TREX-CoE/qmckl/releases/download/vX.Y.Z/qmckl-X.Y.Z.tar.gz
tar -xzvf qmckl-X.Y.Z.tar.gz
cd qmckl-X.Y.Z
```

**Basic installation:**

```bash
./configure --prefix=$HOME/.local
make -j$(nproc)
make check
make install
export QMCKL_DIR=$HOME/.local
```

### Method 2: From Source Repository

Clone the latest development version:

```bash
git clone https://github.com/TREX-CoE/qmckl.git
cd qmckl
```

**Configure and build:**

For documentation version or human readable version:
```bash
./configure --prefix=$HOME/.local
make -j$(nproc)
make check
make install
```

### Optimized Builds

For production use, build with compiler optimizations:

#### Using GNU Compilers

```bash
./configure \
  --prefix=$HOME/.local \
  --enable-hpc \
  CC=gcc \
  CFLAGS="-O3 -march=native -flto -fno-trapping-math -fno-math-errno -ftree-vectorize" \
  FC=gfortran \
  FCFLAGS="-O3 -march=native -flto -ftree-vectorize"
make -j$(nproc)
make check
make install
```

#### Using Intel Compilers

```bash
./configure \
  --prefix=$HOME/.local \
  --enable-hpc \
  --with-icx \
  --with-ifx \
  --with-openmp \
  --disable-python \
  --disable-doc
make -j$(nproc)
make check
make install
```

#### Using LLVM/Clang

```bash
./configure \
  --prefix=$HOME/.local \
  --enable-hpc \
  CC=clang \
  CFLAGS="-O3 -march=native -flto" \
  FC=gfortran \
  FCFLAGS="-O3 -march=native -flto"
make -j$(nproc)
make check
make install
```

### Set Environment Variable

Add QMCkl to your environment:

```
export QMCKL_DIR=$HOME/.local
```

Make it persistent by adding to `~/.bashrc`:
```
echo "export QMCKL_DIR=$HOME/.local" >> ~/.bashrc
```

## Building CHAMP with QMCkl Support

Once both TREXIO and QMCkl are installed, build CHAMP with QMCkl enabled:

```bash
export TREXIO_DIR=/path/to/trexio
export QMCKL_DIR=/path/to/qmckl
cmake -S. -Bbuild \
    -DCMAKE_Fortran_COMPILER=mpiifx \
    -DENABLE_TREXIO=ON \
    -DTREXIO_INCLUDE_DIR=${TREXIO_DIR}/include \
    -DTREXIO_LIBRARY=${TREXIO_DIR}/lib/libtrexio.so \
    -DENABLE_QMCKL=ON \
    -DQMCKL_INCLUDE_DIR=${QMCKL_DIR}/include \
    -DQMCKL_LIBRARY=${QMCKL_DIR}/lib/libqmckl.so
cmake --build build -j
```

CMake will automatically detect QMCkl and report:

```
Looking for QMCkl library:
-- Is QMCkl library found           :: TRUE
-- QMCKL Library include dirs       :: /path/to/include
-- QMCKL Library lib dirs           :: /path/to/lib/libqmckl.so
```

## Verifying Installation

Check QMCkl installation:

```
ls $QMCKL_DIR/lib/libqmckl.*
ls $QMCKL_DIR/include/qmckl.h
```

Run QMCkl tests:

```
cd qmckl-build-directory
make check
```

Test with CHAMP:
```
cmake -S. -Bbuild -DENABLE_TREXIO=ON -DENABLE_QMCKL=ON
```

Look for confirmation messages in CMake output.

## HPC Systems

On HPC systems, QMCkl might be available as a module:

```
module avail qmckl
module load qmckl
```

If not available, install to your home directory following the instructions above.

## Additional Resources


[:fontawesome-solid-book: QMCkl Documentation   :fontawesome-solid-arrow-up-right-from-square:](https://trex-coe.github.io/qmckl/)

[:fontawesome-solid-book: QMCkl API Reference :fontawesome-solid-arrow-up-right-from-square:](https://trex-coe.github.io/qmckl/html/)

[:fontawesome-brands-github: QMCkl GitHub Repository :fontawesome-solid-arrow-up-right-from-square:](https://github.com/TREX-CoE/qmckl)

[:fontawesome-solid-globe: TREX-CoE Project :fontawesome-solid-arrow-up-right-from-square:](https://trex-coe.eu/)
