---
title: Getting Started
tags:
    - getting started
    - installation
    - CLI
---

# Getting Started with CHAMP

This section guides you through the essential steps to get CHAMP up and running on your system. Whether you're setting up CHAMP on a desktop workstation or a high-performance computing cluster, these guides will help you install, configure, and begin using CHAMP for quantum Monte Carlo calculations.

## What You'll Find Here

### [Installation](installation/index.md)

Complete installation guides for CHAMP covering multiple platforms and configurations:

- **[Dependencies](installation/dependencies/index.md)** - Installation guides for all required and optional libraries:
    - **[CMake](installation/dependencies/cmake.md)** - Build system setup
    - **[Compilers](installation/dependencies/compilers.md)** - Fortran and C compiler installation for all supported platforms
    - **[BLAS/LAPACK](installation/dependencies/blas_lapack.md)** - Linear algebra libraries from various vendors (Intel MKL, OpenBLAS, BLIS, etc.)
    - **[TREXIO](installation/dependencies/trexio.md)** - Optional library for reading wavefunction data (recommended)
    - **[QMCkl](installation/dependencies/qmckl.md)** - Optional high-performance QMC kernel library
- **[Building from Source](installation/from_source.md)** - Detailed instructions for compiling CHAMP with various compilers (GNU, Intel, LLVM, NVIDIA HPC, Fujitsu) and optional libraries (TREXIO, QMCkl)


### [Command-Line Interface](installation/supercomputers/cli.md)

Learn how to use CHAMP from the command line:

- Command-line options and flags
- Input and output file specifications
- Execution modes (VMC and DMC)
- MPI execution examples

- Common usage patterns

### [Recipes for Supercomputers](installation/supercomputers/index.md)

Specific installation guides and job scripts for major supercomputing centers:

- **[LUMI](installation/supercomputers/lumi.md)** - Recipes for the pre-exascale EuroHPC LUMI supercomputer
- **[Fugaku](installation/supercomputers/fugaku.md)** - Recipes for the Fugaku supercomputer
- **[Snellius](installation/supercomputers/snellius.md)** - Recipes for the Dutch national supercomputer
- **[CCPHead](installation/supercomputers/ccphead.md)** - Recipes for the CCPHead cluster

## Quick Navigation

**New users should follow this path:**

1. Start with **[Installation → Dependencies](installation/dependencies/index.md)** to install required libraries
2. Follow **[Installation → Building from Source](installation/from_source.md)** to compile CHAMP
3. Learn the **[Command-Line Interface](installation/supercomputers/cli.md)** to run calculations
4. Proceed to **[Tutorials](../tutorials/index.md)** for hands-on examples

**If you're on a supercomputer:**

1. Check if your system has a dedicated guide in the Installation section
2. Follow the platform-specific instructions for module loading and compilation
3. Use the provided job script templates for your first runs

## Prerequisites

Before you begin, you should have:

- Basic familiarity with Linux/Unix command line
- Access to a system with a Fortran compiler and MPI
- Understanding of your system's module environment (for HPC systems)

## Support

If you encounter issues during installation or setup:

- Consult the **[Troubleshooting Guide](../troubleshooting/index.md)**
- Review the platform-specific guide for your system
- Check the dependency installation guides for library-specific issues
- Visit the [GitHub repository](https://github.com/filippi-claudia/champ) for community support
