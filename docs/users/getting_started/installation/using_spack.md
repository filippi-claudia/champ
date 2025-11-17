---
title: Installing with Spack
tags:
    - installation
    - spack
    - package manager
---

# Installing CHAMP with Spack

[Spack](https://spack.io/) is a flexible package manager designed for scientific software on HPC systems. It can automatically handle dependencies, compiler selection, and optimizations for your specific hardware. This guide covers two scenarios: using Spack on supercomputers where it's already available, and installing Spack from scratch.

## What is Spack?

Spack simplifies the installation of complex scientific software by:

- Automatically resolving and installing dependencies
- Supporting multiple versions and configurations side-by-side
- Providing reproducible builds
- Optimizing for different compilers and architectures
- Managing module files for easy environment loading

## Scenario 1: Using Spack on Supercomputers

Many HPC systems provide Spack pre-installed or as a loadable module.

### Check if Spack is Available

```
# Check for spack command
which spack

# Or check for spack module
module avail spack
```

### Load Spack (if available as module)

```
module load spack
```

Or source the setup script if provided by your system:

```
source /path/to/spack/share/spack/setup-env.sh
```

### Verify Spack Installation

```
spack --version
spack find
```

### Install CHAMP with Spack

**Basic installation:**

```
spack install champ
```

**With specific compiler:**

```
# List available compilers
spack compilers

# Install with specific compiler (e.g., gcc@11.2.0)
spack install champ %gcc@11.2.0
```

**With TREXIO support:**

```
spack install champ +trexio
```

**With QMCkl support:**

```
spack install champ +trexio +qmckl
```

**With specific MPI implementation:**

```
spack install champ ^openmpi@4.1.4
# or
spack install champ ^intel-mpi
```

**Complete optimized build:**

```
spack install champ %gcc@11.2.0 +trexio +qmckl ^openmpi@4.1.4 ^intel-mkl
```

### Load CHAMP

After installation, load CHAMP into your environment:

```
spack load champ
```

Or create a module file:

```
spack module tcl refresh champ
module load champ
```

### Verify Installation

```
which vmc.mov1
vmc.mov1 --version
```

### Example: Fugaku Supercomputer

On Fugaku, Spack is available and provides optimized packages for the ARM-based Fujitsu A64FX processor:

```
# Load Spack environment
. /vol0004/apps/oss/spack/share/spack/setup-env.sh

# Load required dependencies via Spack
spack load cmake@3.24.3%fj@4.8.1/p5qsrqc
spack load fujitsu-mpi@head%fj@4.8.1
spack load hdf5@1.12.2%fj@4.8.1/tpglq6h
spack load fujitsu-ssl2@head%fj@4.8.1/nndozbk

# Install CHAMP with Fujitsu compiler and dependencies
spack install champ %fj@4.8.1 +trexio ^fujitsu-mpi ^fujitsu-ssl2 ^hdf5@1.12.2

# Load CHAMP
spack load champ

# Verify installation
which vmc.mov1
vmc.mov1 --version
```

**Note:** The Fujitsu compiler (`fj`) and Fujitsu-SSL2 library provide highly optimized BLAS/LAPACK routines specifically tuned for the A64FX architecture, delivering excellent performance on Fugaku.

## Scenario 2: Installing Spack from Scratch

If Spack is not available on your system, you can install it in your home directory.

### Prerequisites

Ensure you have:

- **Git** - To clone Spack repository
- **Python** >= 3.6
- **C/C++ compiler** - For building packages
- **make**, **patch**, **tar**, **gzip**, **bzip2**, **unzip**

### Install Spack

**Clone the Spack repository:**

```
cd ~
git clone -c feature.manyFiles=true https://github.com/spack/spack.git
cd spack
```

**Choose a version (optional but recommended):**

```
# Use latest release
git checkout releases/latest

# Or specific version
git checkout v0.21.0
```

**Set up Spack environment:**

```
# For bash/zsh
. ~/spack/share/spack/setup-env.sh

# Add to ~/.bashrc for automatic loading
echo ". ~/spack/share/spack/setup-env.sh" >> ~/.bashrc
```

**For tcsh/csh:**

```
source ~/spack/share/spack/setup-env.csh
echo "source ~/spack/share/spack/setup-env.csh" >> ~/.tcshrc
```

### Configure Spack

**Detect available compilers:**

```
spack compiler find
spack compilers
```

**Add compilers manually (if needed):**

```
spack compiler add /path/to/compiler/bin
```

**Configure external packages (recommended):**

Create or edit `~/.spack/packages.yaml` to use system libraries:

```
packages:
  mpi:
    buildable: false
    externals:
    - spec: openmpi@4.1.4
      prefix: /usr/local/openmpi
  
  cmake:
    buildable: false
    externals:
    - spec: cmake@3.24.0
      prefix: /usr
  
  blas:
    buildable: false
    externals:
    - spec: openblas@0.3.20
      prefix: /usr
  
  lapack:
    buildable: false
    externals:
    - spec: openblas@0.3.20
      prefix: /usr

  all:
    compiler: [gcc@11.2.0]
    providers:
      mpi: [openmpi]
      blas: [openblas]
      lapack: [openblas]
```

This prevents Spack from building system packages from source.

### Install CHAMP Dependencies

**Install required dependencies:**

```
# Install MPI (if not using system MPI)
spack install openmpi

# Install BLAS/LAPACK (if not using system libraries)
spack install openblas
```

**Install optional dependencies:**

```
# Install HDF5
spack install hdf5 +mpi

# Install TREXIO
spack install trexio

# Install QMCkl
spack install qmckl
```

### Install CHAMP

**Basic installation:**

```
spack install champ
```

**With all features:**

```
spack install champ %gcc@11.2.0 +trexio +qmckl ^openmpi ^openblas
```

**Monitor installation progress:**

Spack will show detailed output during compilation. This may take 30-60 minutes depending on your system and whether dependencies need to be built.

### Load and Use CHAMP

```
# Load CHAMP
spack load champ

# Verify
which vmc.mov1
vmc.mov1 --version

# Check dependencies
spack find -dl champ
```

## Spack Variants and Options

CHAMP Spack package supports several variants (options):

| Variant | Default | Description |
|---------|---------|-------------|
| `+trexio` | OFF | Enable TREXIO support |
| `+qmckl` | OFF | Enable QMCkl support |
| `+mpi` | ON | Enable MPI support |

**Example configurations:**

```
# Minimal build (no optional features)
spack install champ~trexio~qmckl

# TREXIO only
spack install champ+trexio~qmckl

# All features
spack install champ+trexio+qmckl
```

## Advanced Spack Usage

### Specifying Dependencies

**Use specific versions:**

```
spack install champ ^openmpi@4.1.4 ^hdf5@1.12.2
```

**Use specific BLAS implementation:**

```
# With Intel MKL
spack install champ ^intel-mkl

# With OpenBLAS
spack install champ ^openblas

# With BLIS
spack install champ ^blis
```

### Creating Environments

Spack environments allow managing multiple configurations:

```
# Create a new environment
spack env create champ-env

# Activate environment
spack env activate champ-env

# Add packages to environment
spack add champ+trexio+qmckl

# Install all packages in environment
spack install

# Deactivate
spack env deactivate
```

### Viewing Installation Details

```
# Show installed packages
spack find

# Show dependencies
spack find -dl champ

# Show installation path
spack location -i champ

# Show build log
spack cd -b champ
cat spack-build-out.txt
```

### Uninstalling

```
# Uninstall CHAMP
spack uninstall champ

# Uninstall with all dependencies (careful!)
spack uninstall --all champ
```

## Using CHAMP Installed with Spack

### In Interactive Sessions

```
# Load CHAMP
spack load champ

# Run calculation
mpirun -np 8 vmc.mov1 -i input.inp -o output.out -e error
```

### In Job Scripts

Add Spack loading to your job script:

```bash
#!/bin/bash
#SBATCH --job-name=champ-job
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=64
#SBATCH --time=01:00:00

# Load Spack
source /path/to/spack/share/spack/setup-env.sh

# Load CHAMP
spack load champ

# Run calculation
srun vmc.mov1 -i input.inp -o output.out -e error
```

### Using Module Files

Generate module files for easier loading:

```
# Generate modules
spack module tcl refresh

# Add to module path
module use $(spack location -i champ)/share/spack/modules/tcl

# Load as module
module load champ
```

## Troubleshooting

### Spack command not found

Ensure you've sourced the setup script:
```
. ~/spack/share/spack/setup-env.sh
```

### Build fails with compiler errors

Check compiler compatibility:
```
spack compilers
spack install champ %gcc@11.2.0
```

### Dependency conflicts

Try using a fresh environment:
```
spack env create clean-env
spack env activate clean-env
spack install champ
```

### Long build times

Use system packages when possible by configuring `packages.yaml` to avoid rebuilding MPI, BLAS, etc.

### Installation fails

View detailed error log:
```
spack cd -b champ
less spack-build-out.txt
```

## Advantages of Using Spack

- **Automated dependency management** - No manual library installation
- **Multiple versions** - Install different CHAMP versions side-by-side
- **Optimized builds** - Compiler and architecture-specific optimizations
- **Reproducibility** - Spack specs ensure consistent builds
- **Easy updates** - `spack install champ@new-version`

## Best Practices

1. **Use system packages when available** - Configure `packages.yaml` to use system MPI, BLAS, etc.
2. **Create environments** - Use Spack environments for different projects
3. **Document your spec** - Save the full Spack spec for reproducibility
4. **Keep Spack updated** - Regularly update Spack for bug fixes and new packages
5. **Test installation** - Run CHAMP test suite after installation

## Additional Resources

- [Spack Documentation](https://spack.readthedocs.io/)
- [Spack GitHub Repository](https://github.com/spack/spack)
- [Spack Tutorial](https://spack-tutorial.readthedocs.io/)
- [CHAMP Spack Package](https://github.com/spack/spack/tree/develop/var/spack/repos/builtin/packages/champ)

## Next Steps

After installing CHAMP with Spack:

1. Verify the installation with `vmc.mov1 --version`
2. Run the test suite if available
3. Learn about the [Command-Line Interface](../cli/index.md)
4. Follow the [Tutorials](../../tutorials/index.md) for example calculations

