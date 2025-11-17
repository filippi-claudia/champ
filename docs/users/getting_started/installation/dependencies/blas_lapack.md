---
title: Linear Algebra Libraries (BLAS/LAPACK)
tags:
    - BLAS
    - LAPACK
    - MKL
    - OpenBLAS
    - BLIS
---

# BLAS and LAPACK Libraries

CHAMP requires BLAS (Basic Linear Algebra Subprograms) and LAPACK (Linear Algebra PACKage) libraries for all linear algebra operations. These are **mandatory dependencies** for building and running CHAMP.

## Overview

CHAMP's build system automatically detects available BLAS/LAPACK implementations. Several vendor-optimized and open-source options are available, offering different performance characteristics depending on your hardware.

## Installation Options

### 1. Intel Math Kernel Library (MKL)

**Recommended for:** Intel and AMD CPUs (best overall performance)

Intel MKL provides highly optimized BLAS/LAPACK implementations with automatic CPU dispatch and threading support.

#### Installation via Intel oneAPI

**Without sudo access:**

Download Intel oneAPI Base Toolkit from the [Intel oneAPI Downloads](https://www.intel.com/content/www/us/en/developer/tools/oneapi/base-toolkit-download.html) page:

```
chmod +x ./l_BaseKit_*.sh
./l_BaseKit_*.sh
```

Follow the installer prompts. After installation, add to your shell configuration:

```
echo "source /path/to/intel/oneapi/setvars.sh" >> ~/.bashrc
source ~/.bashrc
```

**With sudo access (Ubuntu/Debian):**

```
wget https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB
sudo apt-key add GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB
sudo add-apt-repository "deb https://apt.repos.intel.com/oneapi all main"
sudo apt-get update
sudo apt-get install intel-oneapi-mkl intel-oneapi-mkl-devel
```

**Loading MKL:**

```
source /opt/intel/oneapi/setvars.sh
```

**CMake configuration:**

!!! info
    CHAMP automatically links to the Intel MKL libraries if Intel oneAPI compilers are used. See the provided `compile-intel.sh` or `compile-intel.qmckl` scripts.


When using Intel compilers, CMake automatically selects MKL variants:
- Static linking: `BLA_VENDOR=Intel10_64lp`
- Dynamic linking: `BLA_VENDOR=Intel10_64_dyn`

You can explicitly specify:

```
cmake -S. -Bbuild -DBLA_VENDOR=Intel10_64lp
```

### 2. OpenBLAS

**Recommended for:** General purpose, open-source option with good performance

OpenBLAS is a high-performance, open-source BLAS/LAPACK implementation with threading support.

**Ubuntu/Debian:**

```
sudo apt-get install libopenblas-dev liblapack-dev
```

**Fedora/RHEL:**

```
sudo dnf install openblas-devel lapack-devel
```

**macOS (Homebrew):**

```
brew install openblas lapack
```

**Building from source:**

```
git clone https://github.com/OpenMathLib/OpenBLAS.git
cd OpenBLAS
make -j$(nproc)
make PREFIX=/path/to/install install
```

Set environment variables:

```
export OpenBLAS_HOME=/path/to/install
export LD_LIBRARY_PATH=$OpenBLAS_HOME/lib:$LD_LIBRARY_PATH
```

### 3. BLIS + libFLAME

**Recommended for:** AMD CPUs (optimized for AMD architectures)

BLIS provides optimized BLAS operations, while libFLAME provides LAPACK.

**Ubuntu/Debian:**

```
sudo apt-get install libblis-dev libflame-dev
```

**Building from source:**

BLIS:
```
git clone https://github.com/amd/blis.git
cd blis
./configure --prefix=/path/to/install auto
make -j$(nproc)
make install
```

libFLAME:
```
git clone https://github.com/amd/libflame.git
cd libflame
./configure --prefix=/path/to/install --enable-dynamic-build
make -j$(nproc)
make install
```

### 4. Reference BLAS/LAPACK (Netlib)

**Use for:** Testing or when vendor libraries are unavailable (not recommended for production)

**Ubuntu/Debian:**

```
sudo apt-get install libblas-dev liblapack-dev
```

**Fedora/RHEL:**

```
sudo dnf install blas-devel lapack-devel
```

**Note:** Reference implementations are not optimized and will result in poor performance.

### 5. NVIDIA cuBLAS (GPU) 

!!! danger "Caution"
    This part of the code is under development

**For GPU-accelerated calculations with NVHPC compiler:**

When building with NVIDIA HPC SDK and GPU support enabled:

```bash
cmake -S. -Bbuild -DENABLE_GPU=ON -DCMAKE_Fortran_COMPILER=nvfortran
```

CHAMP will automatically link against cuBLAS libraries provided by CUDA toolkit.

## Configuring CHAMP Build

### Automatic Detection

By default, CMake automatically finds available BLAS/LAPACK libraries:

```bash
cmake -S. -Bbuild
```

### Specifying a Vendor

To explicitly select a BLAS/LAPACK vendor:

```bash
cmake -S. -Bbuild -DBLA_VENDOR=OpenBLAS
```

Supported vendor names:

- `Intel10_64lp` (MKL static)
- `Intel10_64_dyn` (MKL dynamic)
- `OpenBLAS`
- `FLAME`
- `ACML`
- `AOCL` (CMake > 3.27)
- `Generic` (Netlib reference)

### Static vs Dynamic Linking

For static linking (recommended for HPC):

```bash
cmake -S. -Bbuild -DBLA_STATIC=ON
```

## Verifying Installation

After configuration, CMake will report the detected libraries:

```bash
-- Using BLAS and LAPACK for the linear algebra calculations!
-- BLAS and LAPACK LIBRARIES:
                                    :: /path/to/libopenblas.so
                                    :: /path/to/liblapack.so
```

## Performance Considerations

!!! tip
    The `QMCKL` library offers granular control of threadings while running CHAMP. 
    More about that in the following sections.


**Threading:** Most modern BLAS libraries support multithreading. Control thread count via:

```bash
export OMP_NUM_THREADS=8          # OpenBLAS, MKL
export OPENBLAS_NUM_THREADS=8     # OpenBLAS specific
export MKL_NUM_THREADS=8          # MKL specific
```

**MPI + Threading:** When running with MPI, typically use 1 thread per MPI rank or hybrid decomposition:

```bash
# Pure MPI
mpirun -np 128 champ/bin/vmc.mov1 -i input.inp -o output.out

# Hybrid MPI + OpenMP (only via MKL)
export MKL_NUM_THREADS=8          # MKL specific
mpirun -np 32 champ/bin/vmc.mov1 -i input.inp -o output.out
```

## HPC Systems

On HPC systems, optimized BLAS/LAPACK are usually available via modules:

```bash
module avail mkl
module avail openblas
module load mkl
```

Consult your system documentation for specific module names and recommended configurations.