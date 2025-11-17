---
title: TREXIO
tags:
    - TREXIO
    - optional
---

# TREXIO Library

TREXIO is an optional library that enables CHAMP to read wavefunction data from HDF5 files in the TREXIO format. This is the **recommended method** for importing wavefunctions from quantum chemistry codes like Quantum Package, PySCF, and others.

!!! note 
    **Minimum version required:** 2.0.0 or higher

## What is TREXIO?

TREXIO is a standardized file format and library for storing electronic structure data (basis sets, molecular orbitals, integrals, etc.) in an efficient, portable manner. It supports both text-based and HDF5-based storage.

## Requirements

Before installing TREXIO, ensure you have:

- **C compiler** (gcc, icc, or clang)
- **Fortran compiler** (gfortran, ifort, or compatible)
- **CMake** >= 3.16 or **Autotools** (autoconf >= 2.69, automake >= 1.11, libtool >= 2.2)
- **HDF5 library** >= 1.8 (optional but strongly recommended for performance)

### Installing HDF5

**Ubuntu/Debian:**
```bash
sudo apt-get install libhdf5-dev
```

**Fedora/RHEL:**
```bash
sudo dnf install hdf5-devel
```

**macOS (Homebrew):**
```bash
brew install hdf5
```

For HPC systems, HDF5 is typically available via modules:
```bash
module load hdf5
```

## Installation Methods

### Method 1: Using CMake (Recommended)

**From release tarball:**

```bash
wget https://github.com/TREX-CoE/trexio/releases/download/vX.Y.Z/trexio-X.Y.Z.tar.gz
tar -xzvf trexio-X.Y.Z.tar.gz
cd trexio-X.Y.Z
```

**Configure and build:**

```bash
cmake -S. -Bbuild -DCMAKE_INSTALL_PREFIX=$HOME/.local
cd build
make -j$(nproc)
ctest
make install
```

**Set environment variable:**

```bash
export TREXIO_DIR=$HOME/.local
```

Add this to your `~/.bashrc` for persistence:
```bash
echo "export TREXIO_DIR=$HOME/.local" >> ~/.bashrc
```

### Method 2: From Release Tarball (Autotools)

**Download the latest release from** [TREXIO GitHub Releases](https://github.com/TREX-CoE/trexio/releases):

```bash
wget https://github.com/TREX-CoE/trexio/releases/download/vX.Y.Z/trexio-X.Y.Z.tar.gz
tar -xzvf trexio-X.Y.Z.tar.gz
cd trexio-X.Y.Z
```

**Configure and build:**

With default installation to `/usr/local`:
```bash
./configure --with-hdf5=/path/to/hdf5/ CC=icx FC=ifx CFLAGS=-O2 FCFLAGS=-O2
make -j$(nproc)
make check
sudo make install
```

For custom installation directory (no sudo required):
```bash
./configure --prefix=$HOME/.local
make -j$(nproc)
make check
make install
export TREXIO_DIR=$HOME/.local
```


### Method 3: From Source Repo (devel version)

For the latest development version:

```bash
git clone https://github.com/TREX-CoE/trexio.git
cd trexio
cmake -S. -Bbuild -DCMAKE_INSTALL_PREFIX=$HOME/.local
cd build
make -j$(nproc)
ctest
make install
export TREXIO_DIR=$HOME/.local
```

## Building CHAMP with TREXIO Support

Once TREXIO is installed, build CHAMP with TREXIO enabled:

```bash
export TREXIO_DIR=/path/to/trexio

cmake -S. -Bbuild -DCMAKE_Fortran_COMPILER=mpiifx \
    -DENABLE_TREXIO=ON \
    -DTREXIO_INCLUDE_DIR=$TREXIO_DIR/include \
    -DTREXIO_LIBRARY=$TREXIO_DIR/lib/libtrexio.so 
cmake --build build -j
```

CMake will automatically detect TREXIO and report:

```
Looking for TREXIO library:
-- Is TREXIO library found          :: TRUE
-- TREXIO Library include dirs      :: /path/to/include
-- TREXIO Library lib dirs          :: /path/to/lib/libtrexio.so
```

## Using TREXIO Files with CHAMP

With TREXIO support enabled, you can directly use TREXIO files as input:

```
champ/bin/vmc.mov1 -i vmc_trexio.inp -o vmc.out -e error
```

In your input file, specify the TREXIO file:
```perl
load trexio 'benzene.hdf5'
```

## Verifying Installation

Check TREXIO installation:

```
ls $TREXIO_DIR/lib/libtrexio.*
ls $TREXIO_DIR/include/trexio.h
```

Test with CHAMP:
```
cmake -S. -Bbuild -DENABLE_TREXIO=ON
```

Look for the confirmation message in CMake output.

## HPC Systems

On HPC systems, TREXIO might be available as a module:

```
module avail trexio
module load trexio
```

If not available, install to your home directory following the CMake method above.

## Additional Resources


[TREXIO Repository](https://github.com/TREX-CoE/trexio){ .md-button .md-button--primary}

[TREXIO User Guide](https://trex-coe.github.io/trexio/trex.html){ .md-button .md-button--primary}

