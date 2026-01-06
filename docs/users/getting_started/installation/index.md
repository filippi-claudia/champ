---
title: Installation
tags:
    - Installation
---

# Installing CHAMP

CHAMP is a high-performance QMC code that can be installed on various platforms, from desktop workstations to supercomputers. This section provides comprehensive installation guides for different environments and use cases.

## Quick Start

For most users, the basic installation process involves:

- [x] **Installing dependencies** - BLAS/LAPACK, MPI, and TREXIO/QMCkl (Done in previous section of [Dependencies](dependencies/index.md))
- [ ] **Obtaining the source code** - Clone from GitHub or download a release
- [ ] **Configuring with CMake** - Set compiler and library options
- [ ] **Building** - Compile CHAMP
- [ ] **Testing** - Verify the installation

## Installation Methods

### From Source

The recommended method for most users. Provides full control over compilation options and optimizations.

- [**Build from Source**](from_source.md) - General instructions for compiling CHAMP with CMake
- [**Dependencies**](dependencies/index.md) - Required and optional libraries

### Platform-Specific Guides

Installation instructions tailored for specific systems:

#### Desktop/Workstation
- [**Ubuntu Desktop**](desktop.md) - Installation on Ubuntu Linux systems with Intel or GNU compilers

#### Supercomputers
- [**LUMI**](supercomputers/lumi.md) - Finland's LUMI supercomputer (lumi.csc.fi)
- [**Fugaku**](supercomputers/fugaku.md) - Japan's Fugaku supercomputer at RIKEN
- [**Snellius**](supercomputers/snellius.md) - Netherlands' Snellius supercomputer
- [**CCPHead**](supercomputers/ccphead.md) - University of Twente's CCPHead computing cluster

Each guide includes:

- System-specific module loading
- Compiler configuration
- Sample job submission scripts
- Performance optimization tips

## Prerequisites (see [Dependencies](dependencies/index.md))

Before installing CHAMP, ensure you have:

### Required
- **CMake** >= 3.17
- **Fortran/C Compiler** (GCC >= 9.3, Intel Fortran >= 2020, or compatible)
- **MPI Library** (OpenMPI >= 3.0 or Intel MPI)
- **BLAS/LAPACK** (or Intel MKL, OpenBLAS, etc.)

### Optional
- **TREXIO** >= 2.0.0 - For reading wavefunction data from TREXIO files (recommended)
- **QMCkl** >= 1.0.0 - For high-performance QMC kernels
- **HDF5** >= 1.8 - Required by TREXIO

See the [Dependencies Guide](dependencies/index.md) for detailed installation instructions.

## Choosing Your Installation Method

**Use the source installation if:**

- [x] You need custom compiler optimizations
- [x] You're installing on an HPC system
- [x] You want the latest development features
- [x] You need specific library versions or configurations

**Use platform-specific guides if:**

- [x] You're working on a known supercomputer
- [x] You want tested module combinations
- [x] You need job submission examples

## Getting the Source Code

### Latest Release (Stable)

Download the latest stable release from [GitHub Releases](https://github.com/filippi-claudia/champ/releases):

```
wget https://github.com/filippi-claudia/champ/archive/refs/tags/vX.Y.Z.tar.gz
tar -xzvf vX.Y.Z.tar.gz
cd champ-X.Y.Z
```

### Development Version

Clone the repository for the latest features:

```
git clone https://github.com/filippi-claudia/champ.git
cd champ
```

## Need Help?

- Check the [Troubleshooting Guide](../../troubleshooting/index.md)
- Review platform-specific guides for your system
- Consult dependency installation guides for library issues
- Visit the [GitHub repository](https://github.com/filippi-claudia/champ) for issues and discussions

## Contributing

If you've successfully installed CHAMP on a new platform or have improvements to these guides, please consider [contributing](../../contributing/tools.md) to the documentation.