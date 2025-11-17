---
title: Building from Source
tags:
    - installation
    - source
    - cmake
---

# Building CHAMP from Source

This guide provides detailed instructions for compiling CHAMP from source code using different compilers and configurations. CHAMP uses CMake as its build system, which automatically detects your compiler and libraries.

## Getting the Source Code

### Option 1: Clone from GitHub (Recommended)

```bash
git clone https://github.com/filippi-claudia/champ.git
cd champ
```

### Option 2: Download a Release

```bash
wget https://github.com/filippi-claudia/champ/archive/refs/tags/vX.Y.Z.tar.gz
tar -xzvf vX.Y.Z.tar.gz
cd champ-X.Y.Z
```

## Basic Build Process

The build process consists of two steps:

1. **Configure** - CMake detects your system and generates build files (run once)
2. **Build** - Compile the source code

```bash
# Configure
cmake -S. -Bbuild -DCMAKE_Fortran_COMPILER=<compiler>

# Build
cmake --build build -j$(nproc)
```

The compiled executables will be in `bin/`:

- `bin/vmc.mov1` - Variational Monte Carlo
- `bin/dmc.mov1` - Diffusion Monte Carlo

## Compiler-Specific Build Recipes

### GNU Compilers (gfortran/gcc)

**Basic build with OpenBLAS:**

```bash
cmake -S. -Bbuild \
  -DCMAKE_Fortran_COMPILER=mpif90 \
  -DCMAKE_C_COMPILER=mpicc

cmake --build build -j$(nproc)
```

**With explicit BLAS/LAPACK paths:**

```bash
cmake -S. -Bbuild \
  -DCMAKE_Fortran_COMPILER=mpif90 \
  -DCMAKE_C_COMPILER=mpicc \
  -DBLAS_LIBRARIES="/usr/lib/x86_64-linux-gnu/libopenblas.so" \
  -DLAPACK_LIBRARIES="/usr/lib/x86_64-linux-gnu/liblapack.so"

cmake --build build -j$(nproc)
```

### Intel Compilers (Classic ifort)

**With Intel MPI and MKL (automatic detection):**

```bash
source /opt/intel/oneapi/setvars.sh

cmake -S. -Bbuild \
  -DCMAKE_Fortran_COMPILER=mpiifort \
  -DCMAKE_C_COMPILER=mpiicc \
  -DBLA_VENDOR=Intel10_64lp

cmake --build build -j$(nproc)
```

**With explicit MKL linking:**

```bash
cmake -S. -Bbuild \
  -DCMAKE_Fortran_COMPILER=mpiifort \
  -DCMAKE_C_COMPILER=mpiicc \
  -DBLAS_LIBRARIES="-qmkl=parallel"

cmake --build build -j$(nproc)
```

### Intel Compilers (LLVM-based ifx)

**Using ifx with MKL:**

```bash
source /opt/intel/oneapi/setvars.sh

cmake -S. -Bbuild \
  -DCMAKE_Fortran_COMPILER=mpiifx \
  -DCMAKE_C_COMPILER=mpiicx \
  -DBLAS_LIBRARIES="-qmkl=parallel

cmake --build build -j$(nproc)
```

### LLVM Flang

**Using Flang with OpenBLAS:**

```bash
cmake -S. -Bbuild \
  -DCMAKE_Fortran_COMPILER=mpiflang \
  -DCMAKE_C_COMPILER=mpicc \
  -DBLA_VENDOR=OpenBLAS

cmake --build build -j$(nproc)
```

### NVIDIA HPC SDK (nvfortran)

**For GPU-accelerated builds:**

```bash
module load nvhpc
module load openmpi

cmake -S. -Bbuild \
  -DCMAKE_Fortran_COMPILER=mpif90 \
  -DCMAKE_C_COMPILER=mpicc \
  -DENABLE_GPU=ON \
  -DNVFORTRAN_PATH=/path/to/nvhpc/compilers/bin

cmake --build build -j$(nproc)
```

### Fujitsu Compiler (frt)

**For ARM-based Fujitsu systems:**

```bash
module load fujitu-mpi

cmake -S. -Bbuild \
  -DCMAKE_Fortran_COMPILER=mpifrt \
  -DCMAKE_C_COMPILER=mpifcc

cmake --build build -j$(nproc)
```

## Building with Optional Libraries

### With TREXIO Support

TREXIO enables reading wavefunction data from TREXIO HDF5 files.

**If TREXIO is in standard location:**

```bash
export TREXIO_DIR=/path/to/trexio/installation

cmake -S. -Bbuild \
  -DCMAKE_Fortran_COMPILER=mpif90 \
  -DENABLE_TREXIO=ON

cmake --build build -j$(nproc)
```

**With explicit TREXIO paths:**

```bash
cmake -S. -Bbuild \
  -DCMAKE_Fortran_COMPILER=mpif90 \
  -DENABLE_TREXIO=ON \
  -DTREXIO_LIBRARY=/path/to/lib/libtrexio.so \
  -DTREXIO_INCLUDE_DIR=/path/to/include

cmake --build build -j$(nproc)
```

**With HDF5 in non-standard location:**

```bash
cmake -S. -Bbuild \
  -DCMAKE_Fortran_COMPILER=mpif90 \
  -DENABLE_TREXIO=ON \
  -DHDF5_LIBRARIES=/path/to/lib/libhdf5.so \
  -DHDF5_INCLUDE_DIRS=/path/to/include

cmake --build build -j$(nproc)
```

### With QMCkl Support

QMCkl provides high-performance implementations of QMC kernels.

!!! note 
    TREXIO must be installed and enabled to use QMCkl.

```bash
export TREXIO_DIR=/path/to/trexio
export QMCKL_DIR=/path/to/qmckl

cmake -S. -Bbuild \
  -DCMAKE_Fortran_COMPILER=mpif90 \
  -DENABLE_TREXIO=ON \
  -DENABLE_QMCKL=ON

cmake --build build -j$(nproc)
```

**With explicit paths:**

```bash
cmake -S. -Bbuild \
  -DCMAKE_Fortran_COMPILER=mpif90 \
  -DENABLE_TREXIO=ON \
  -DENABLE_QMCKL=ON \
  -DTREXIO_LIBRARY=/path/to/lib/libtrexio.so \
  -DTREXIO_INCLUDE_DIR=/path/to/trexio/include \
  -DQMCKL_LIBRARY=/path/to/lib/libqmckl.so \
  -DQMCKL_INCLUDE_DIR=/path/to/qmckl/include

cmake --build build -j$(nproc)
```

### Complete Build with All Features

```
export TREXIO_DIR=$HOME/.local
export QMCKL_DIR=$HOME/.local

cmake -S. -Bbuild \
  -DCMAKE_Fortran_COMPILER=mpiifort \
  -DCMAKE_C_COMPILER=mpiicc \
  -DBLA_VENDOR=Intel10_64lp \
  -DENABLE_TREXIO=ON \
  -DENABLE_QMCKL=ON \
  -DCMAKE_BUILD_TYPE=Release

cmake --build build -j$(nproc)
```

## CMake Configuration Options

### Common Options

| Option | Values | Description |
|--------|--------|-------------|
| `CMAKE_Fortran_COMPILER` | compiler path | MPI-wrapped Fortran compiler |
| `CMAKE_C_COMPILER` | compiler path | MPI-wrapped C compiler |
| `CMAKE_BUILD_TYPE` | `Release`, `Debug` | Build type (default: Release) |
| `ENABLE_MPI` | `ON`, `OFF` | Enable MPI support (default: ON) |

### BLAS/LAPACK Options

| Option | Values | Description |
|--------|--------|-------------|
| `BLA_VENDOR` | `Intel10_64lp`, `OpenBLAS`, etc. | BLAS vendor to use |
| `BLA_STATIC` | `ON`, `OFF` | Use static BLAS libraries |
| `BLAS_LIBRARIES` | library paths | Explicit BLAS libraries |
| `LAPACK_LIBRARIES` | library paths | Explicit LAPACK libraries |

### Optional Features

| Option | Values | Description |
|--------|--------|-------------|
| `ENABLE_TREXIO` | `ON`, `OFF` | Enable TREXIO support (default: OFF) |
| `ENABLE_QMCKL` | `ON`, `OFF` | Enable QMCkl support (default: OFF) |
| `ENABLE_GPU` | `ON`, `OFF` | Enable GPU acceleration (default: OFF) |
| `ENABLE_QMMM` | `ON`, `OFF` | Enable QM/MM support (default: OFF) |

### Library Paths

| Option | Type | Description |
|--------|------|-------------|
| `TREXIO_DIR` | directory | TREXIO installation prefix |
| `QMCKL_DIR` | directory | QMCkl installation prefix |
| `TREXIO_LIBRARY` | file path | Path to libtrexio.so |
| `TREXIO_INCLUDE_DIR` | directory | TREXIO header files |
| `QMCKL_LIBRARY` | file path | Path to libqmckl.so |
| `QMCKL_INCLUDE_DIR` | directory | QMCkl header files |
| `HDF5_LIBRARIES` | file path | Path to libhdf5.so |
| `HDF5_INCLUDE_DIRS` | directory | HDF5 header files |

## Build Commands

### Standard Build

```
cmake --build build -j$(nproc)
```

### Clean and Rebuild

```
cmake --build build --clean-first -j$(nproc)
```

### Parallel Build (specify thread count)

```
cmake --build build -j8
```

### Verbose Build (see all commands)

```
cmake --build build -j$(nproc) --verbose
```

## Verifying the Build

### Run Tests

```
cd build
ctest
```

Or with verbose output:
```
ctest --verbose
```

### Check Compiled Executables

```
ls -lh bin/
```

You should see:
- `vmc.mov1`
- `dmc.mov1`
- Other utility programs

### Verify Library Linking

```
ldd bin/vmc.mov1
```

Check that BLAS/LAPACK, MPI, and optional libraries (TREXIO, QMCkl) are properly linked.

## Build Types

### Release Build (Default)

Optimized for performance:
```
cmake -S. -Bbuild -DCMAKE_BUILD_TYPE=Release
```

### Debug Build

With debugging symbols and runtime checks:
```
cmake -S. -Bbuild -DCMAKE_BUILD_TYPE=Debug
```

Useful for development and troubleshooting.

## Reconfiguring

If you need to change configuration options:

**Option 1: Remove build directory and reconfigure**
```bash
rm -rf build bin
cmake -S. -Bbuild [new options]
```

**Option 2: Use ccmake (interactive)**
```
ccmake build
```

## Installation

To install CHAMP to a specific location:

```bash
cmake -S. -Bbuild -DCMAKE_INSTALL_PREFIX=$HOME/.local
cmake --build build -j$(nproc)
cmake --install build
```

Executables will be in `$HOME/.local/bin/`.


## Next Steps

After successful compilation:

1. [Run tests](testsuite.md) to verify installation
2. Learn about the [Command-Line Interface](../cli/index.md)
3. Follow [Tutorials](../../tutorials/index.md) for example calculations
4. Prepare [Input Files](../../preparation/index.md) for your system

## Platform-Specific Guides

For optimized builds on specific supercomputers, see:

- [LUMI](lumi.md)
- [Fugaku](fugaku.md)
- [Snellius](snellius.md)
- [Ubuntu Desktop](desktop.md)

