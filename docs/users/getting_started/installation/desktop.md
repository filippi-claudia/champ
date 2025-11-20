---
title: Ubuntu Desktop Installation
tags:
    - installation
    - desktop
    - ubuntu
---

# Installing CHAMP on Ubuntu Desktop

This guide provides step-by-step instructions for installing CHAMP on Ubuntu desktop systems. Both GNU and Intel compiler options are covered.

## Prerequisites

- **Ubuntu** 20.04 LTS or newer
- **Internet connection** for downloading packages
- **Sudo privileges** for installing system packages
- **At least 2 GB free disk space**

## Quick Start

For users who want to get started immediately:

```bash
# Install dependencies
sudo apt-get update
sudo apt-get install -y gfortran gcc \
    cmake git openmpi-bin libopenmpi-dev \
    libopenblas-dev liblapack-dev

# Clone CHAMP
git clone https://github.com/filippi-claudia/champ.git
cd champ

# Build
cmake -S. -Bbuild -DCMAKE_Fortran_COMPILER=mpif90
cmake --build build -j$(nproc)

# Test
./bin/vmc.mov1 --version
```

## Method 1: Using GNU Compilers (Recommended for Desktop)

### Step 1: Install Build Dependencies

Update package lists and install required packages:

```bash
sudo apt-get update
sudo apt-get install -y \
    build-essential \
    gfortran \
    gcc \
    g++ \
    cmake \
    git \
    gawk
```

### Step 2: Install MPI

Install OpenMPI for parallel execution:

```bash
sudo apt-get install -y \
    openmpi-bin \
    libopenmpi-dev
```

Verify MPI installation:

```bash
mpif90 --version
mpirun --version
which mpif90
```

### Step 3: Install BLAS/LAPACK

Install linear algebra libraries:

```bash
sudo apt-get install -y \
    libopenblas-dev \
    liblapack-dev \
    libscalapack-openmpi-dev
```

### Step 4: Install HDF5 (Optional, for TREXIO)

If you plan to use TREXIO files:

```bash
sudo apt-get install -y \
    libhdf5-openmpi-dev \
    libhdf5-dev
```

### Step 5: Install TREXIO (Optional)

**Option A: From Ubuntu packages (if available):**

```bash
sudo apt-get install -y libtrexio-dev
```

**Option B: Build from source:**

```bash
# Download TREXIO
wget https://github.com/TREX-CoE/trexio/releases/download/v2.4.2/trexio-2.4.2.tar.gz
tar -xzvf trexio-2.4.2.tar.gz
cd trexio-2.4.2

# Build and install
cmake -S. -Bbuild -DCMAKE_INSTALL_PREFIX=$HOME/.local
cd build
make -j$(nproc)
make install

# Add to environment
export TREXIO_DIR=$HOME/.local
echo "export TREXIO_DIR=$HOME/.local" >> ~/.bashrc

cd ../..
```

### Step 6: Get CHAMP Source Code

Clone the repository:

```bash
git clone https://github.com/filippi-claudia/champ.git
cd champ
```

Or download a specific release:

```bash
wget https://github.com/filippi-claudia/champ/archive/refs/tags/v2.3.0.tar.gz
tar -xzvf v2.3.0.tar.gz
cd champ-2.3.0
```

### Step 7: Configure CHAMP

**Basic build (without TREXIO):**

```bash
cmake -S. -Bbuild \
    -DCMAKE_Fortran_COMPILER=mpif90 \
    -DCMAKE_C_COMPILER=mpicc
```

**With TREXIO support:**

```bash
export TREXIO_DIR=$HOME/.local

cmake -S. -Bbuild \
    -DCMAKE_Fortran_COMPILER=mpif90 \
    -DCMAKE_C_COMPILER=mpicc \
    -DENABLE_TREXIO=ON
```

### Step 8: Build CHAMP

```bash
cmake --build build -j$(nproc)
```

This compiles CHAMP using all available CPU cores. The build takes 5-15 minutes.

### Step 9: Verify Installation

```bash
# Check executables
ls -lh bin/

# Test version
./bin/vmc.mov1 --version

# Check MPI
mpirun -np 2 ./bin/vmc.mov1 --version
```

### Step 10: Run Tests (Optional)

```bash
cd build
ctest
cd ..
```

### Running Calculations

**Serial execution:**

```bash
./bin/vmc.mov1 -i input.inp -o output.out -e error
```

**Parallel execution:**

```bash
mpirun -np 4 ./bin/vmc.mov1 -i input.inp -o output.out -e error
```

## Method 2: Using Intel oneAPI Compilers

Intel compilers can provide better performance on Intel CPUs.

### Step 1: Install Intel oneAPI

Add Intel repository:

```bash
wget https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB
sudo apt-key add GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB
sudo add-apt-repository "deb https://apt.repos.intel.com/oneapi all main"
sudo apt-get update
```

Install Intel oneAPI Base Toolkit and HPC Toolkit:

```bash
sudo apt-get install -y \
    intel-oneapi-compiler-fortran \
    intel-oneapi-compiler-dpcpp-cpp \
    intel-oneapi-mkl \
    intel-oneapi-mkl-devel \
    intel-oneapi-mpi \
    intel-oneapi-mpi-devel
```

### Step 2: Install Additional Tools

```bash
sudo apt-get install -y cmake git gawk
```

### Step 3: Load Intel Environment

```bash
source /opt/intel/oneapi/setvars.sh
```

Add to `~/.bashrc` for automatic loading:

```bash
echo "source /opt/intel/oneapi/setvars.sh > /dev/null 2>&1" >> ~/.bashrc
```

### Step 4: Verify Intel Compilers

```bash
ifort --version
mpiifort --version
```

### Step 5: Get CHAMP Source Code

```bash
git clone https://github.com/filippi-claudia/champ.git
cd champ
```

### Step 6: Configure with Intel Compilers

**Basic build with Intel MKL:**

```bash
cmake -S. -Bbuild \
    -DCMAKE_Fortran_COMPILER=mpiifort \
    -DCMAKE_C_COMPILER=mpiicc \
    -DBLA_VENDOR=Intel10_64lp
```

**With TREXIO support:**

```bash
cmake -S. -Bbuild \
    -DCMAKE_Fortran_COMPILER=mpiifort \
    -DCMAKE_C_COMPILER=mpiicc \
    -DBLA_VENDOR=Intel10_64lp \
    -DENABLE_TREXIO=ON
```

### Step 7: Build CHAMP

```bash
cmake --build build -j$(nproc)
```

### Step 8: Run Calculations

```bash
# Load Intel environment (if not in .bashrc)
source /opt/intel/oneapi/setvars.sh

# Run with Intel MPI
mpirun -np 8 ./bin/vmc.mov1 -i input.inp -o output.out -e error
```

## Adding CHAMP to Your PATH

To use CHAMP from anywhere:

```bash
# Add to PATH
export PATH=/path/to/champ/bin:$PATH

# Make permanent
echo "export PATH=/path/to/champ/bin:$PATH" >> ~/.bashrc
source ~/.bashrc
```

Then you can run:

```bash
vmc.mov1 -i input.inp -o output.out -e error
```

## Optimizing for Your System

### Detect CPU Architecture

```bash
# Check CPU info
lscpu | grep "Model name"
gcc -march=native -Q --help=target | grep march
```

### Optimize Build

Add architecture-specific flags during configuration:

```bash
cmake -S. -Bbuild \
    -DCMAKE_Fortran_COMPILER=mpif90 \
    -DCMAKE_Fortran_FLAGS="-O3 -march=native" \
    -DCMAKE_C_FLAGS="-O3 -march=native"
```

## Common Desktop Use Cases

### Setup for Learning/Development

Minimal setup for testing and learning:

```bash
sudo apt-get install gfortran openmpi-bin libopenmpi-dev libopenblas-dev cmake git
git clone https://github.com/filippi-claudia/champ.git
cd champ
cmake -S. -Bbuild -DCMAKE_Fortran_COMPILER=mpif90
cmake --build build -j4
```

### Setup for Production Calculations

Full setup with all features:

```bash
# Install all dependencies including TREXIO
sudo apt-get install gfortran openmpi-bin libopenmpi-dev libopenblas-dev \
    libhdf5-openmpi-dev cmake git

# Install TREXIO
wget https://github.com/TREX-CoE/trexio/releases/download/v2.4.2/trexio-2.4.2.tar.gz
tar -xzvf trexio-2.4.2.tar.gz
cd trexio-2.4.2
cmake -S. -Bbuild -DCMAKE_INSTALL_PREFIX=$HOME/.local
cd build && make -j$(nproc) && make install
cd ../..

# Build CHAMP with all features
export TREXIO_DIR=$HOME/.local
git clone https://github.com/filippi-claudia/champ.git
cd champ
cmake -S. -Bbuild -DCMAKE_Fortran_COMPILER=mpif90 -DENABLE_TREXIO=ON
cmake --build build -j$(nproc)
```

## Performance Tips for Desktop

### CPU Affinity

For better performance, pin MPI ranks to cores:

```
mpirun -np 4 --bind-to core ./bin/vmc.mov1 -i input.inp -o output.out
```

### Thread Count

Control OpenMP threads (if using threaded BLAS):

```
export OMP_NUM_THREADS=1
mpirun -np $(nproc) ./bin/vmc.mov1 -i input.inp -o output.out
```

### Memory Usage

Monitor memory usage:

```
/usr/bin/time -v mpirun -np 4 ./bin/vmc.mov1 -i input.inp -o output.out
```

## Uninstalling

To remove CHAMP:

```
cd champ
rm -rf build bin
```

To remove dependencies:

```
# GNU compiler stack
sudo apt-get remove gfortran openmpi-bin libopenmpi-dev libopenblas-dev

# Intel oneAPI
sudo apt-get remove intel-oneapi-*
```

## Additional Resources

- [Building from Source (detailed)](from_source.md)
- [Dependencies Guide](dependencies/index.md)
- [Docker Installation](using_docker.md) (alternative method)
