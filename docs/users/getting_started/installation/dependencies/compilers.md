---
title: Compilers
tags:
    - compilers
    - GNU
    - Intel
    - LLVM
    - NVHPC
    - Fujitsu
---

# Compiler Support

CHAMP supports multiple Fortran and C compilers. The build system automatically detects your compiler and applies appropriate optimization flags.

## Supported Compilers

### GNU Compiler Collection (GCC)

**Minimum version:** GCC 9.3.0 or higher

**Installation:**

Ubuntu/Debian:
```bash
sudo apt-get install -y gfortran gcc
```

Fedora/RHEL:
```bash
sudo dnf install gcc-gfortran gcc
```

macOS (Homebrew):
```bash
brew install gcc
```

**Compiler flags auto-applied:**

- `-O2` - Optimization level 2
- `-cpp` - Enable preprocessor
- `-mcmodel=large` - Large memory model support
- `-ffree-line-length-none` - No limit on free-form line length
- `-D_MPI_` and `-DCLUSTER` - MPI and cluster support

**Debug mode flags:**

- `-fcheck=all` - Runtime checks
- `-fbacktrace` - Backtrace on errors
- `-g` - Debug symbols
- `-Wall -Wextra` - Enable warnings

### Intel Fortran Compiler (Classic)

**Minimum version:** Intel Fortran 2020 or higher

**Installation:**

The Intel Fortran compiler is part of Intel oneAPI HPC Toolkit.

Without sudo access - download installers from [Intel oneAPI Downloads](https://www.intel.com/content/www/us/en/developer/tools/oneapi/toolkits.html):

```bash
chmod +x ./l_BaseKit_*.sh
./l_BaseKit_*.sh

chmod +x ./l_HPCKit_*.sh
./l_HPCKit_*.sh
```

With sudo access:

```
wget https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB
sudo apt-key add GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB
sudo add-apt-repository "deb https://apt.repos.intel.com/oneapi all main"
sudo apt-get update
sudo apt-get install intel-oneapi-compiler-fortran intel-oneapi-mkl intel-oneapi-mpi
```

**Loading the compiler:**
```
source /opt/intel/oneapi/setvars.sh
```

**Compiler flags auto-applied:**
- `-O2` - Optimization level 2
- `-ipo` - Interprocedural optimization
- `-fpp` - Fortran preprocessor
- `-mcmodel=small` - Memory model
- `-shared-intel` - Dynamic linking for Intel libraries
- Includes MKL optimizations for Intel10_64lp or Intel10_64_dyn

### Intel Fortran Compiler (LLVM-based ifx)

**Available in:** Intel oneAPI 2021.1 and later

The newer LLVM-based Intel Fortran compiler (`ifx`) is also supported with similar flags to the classic compiler but without `-ipo` flag.

**Installation:** Same as Intel Classic (part of oneAPI HPC Toolkit)

Without sudo access:

```bash
chmod +x ./l_BaseKit_*.sh
./l_BaseKit_*.sh

chmod +x ./l_HPCKit_*.sh
./l_HPCKit_*.sh
```

With sudo access:

```
wget https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB
sudo apt-key add GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB
sudo add-apt-repository "deb https://apt.repos.intel.com/oneapi all main"
sudo apt-get update
sudo apt-get install intel-oneapi-compiler-fortran intel-oneapi-mkl intel-oneapi-mpi
```

**Loading the compiler:**
```
source /opt/intel/oneapi/setvars.sh
```

### LLVM Flang

**Compilers supported:**

- Classic Flang (`flang`)
- AMD Optimized Flang (AOCC `aocc-flang`)

**Compiler flags applied:**

- `-O2` - Optimization level 2
- `-cpp` - Preprocessor support
- `-mcmodel=large` - Large memory model

**Note:** Do not confuse with `flang-new` (LLVMFlang in CMake), which has different compatibility.

### NVIDIA HPC SDK (NVHPC)

**Use for:** GPU-accelerated calculations with OpenACC

**Installation:** Download from [NVIDIA HPC SDK website](https://developer.nvidia.com/hpc-sdk)

**Compiler flags applied:**

- `-O2` - Optimization
- `-mp=multicore,gpu` - OpenMP for multicore and GPU
- `-acc` - OpenACC support
- `-gpu=cc80` - Target GPU compute capability
- `-Mcudalib=cublas` - Use cuBLAS library
- CUDA libraries: cuFFT, cuBLAS, cuBLASLt, cuDART

**GPU support requirements:**

- CUDA toolkit
- Set `ENABLE_GPU=ON` in CMake
- Optionally set `NVFORTRAN_PATH` to CUDA executable path

### Fujitsu Compiler

**Use for:** Fugaku supercomputer and other Fujitsu systems

**Compiler flags auto-applied:**

- `-O3` - High optimization
- `-Cpp` - Preprocessor
- `-Kopenmp` - OpenMP support
- `-Kparallel -Kfast` - Auto-parallelization and fast math
- `-mcmodel=large` - Large memory model

**Note:** Available on Fugaku supercomputer and ARM-based Fujitsu systems

## Build Type Configuration

CHAMP supports different build types:

**Release (default):**

- Optimized for performance
- Minimal runtime checks

**Debug:**

- Enabled with `-DCMAKE_BUILD_TYPE=DEBUG`
- Additional runtime checks and error tracing
- Debug symbols for debugging with gdb/lldb

## Verifying Your Compiler

Check your compiler version:

```
gfortran --version    # GNU
ifort --version       # Intel Classic
ifx --version         # Intel LLVM
nvfortran --version   # NVIDIA
frt --version         # Fujitsu
flang --version       # LLVM Flang
```

## HPC Systems

On HPC systems, compilers are typically available through module systems:

```
module avail gcc
module avail intel-compilers
module load gcc/11.2.0
```

Refer to your system's documentation for specific module names. Specific recipes
for a select supercomputers are provided here in the documentation.

