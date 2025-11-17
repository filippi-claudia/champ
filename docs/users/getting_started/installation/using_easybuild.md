---
title: Installing with EasyBuild
tags:
    - installation
    - easybuild
    - package manager
---

# Installing CHAMP with EasyBuild

[EasyBuild](https://easybuild.io/) is a software build and installation framework designed for HPC systems. It automates the process of building and installing scientific software with proper dependency resolution and toolchain management. This guide covers using EasyBuild on supercomputers where it's available, and installing EasyBuild from scratch.

## What is EasyBuild?

EasyBuild simplifies the installation of complex scientific software by:

- Using **toolchains** - standardized combinations of compilers and libraries
- Automatically resolving and building dependencies
- Providing **easyconfig** files - recipes for building software
- Supporting multiple versions and configurations side-by-side
- Generating module files automatically for easy environment loading
- Ensuring reproducible builds across HPC systems

## Scenario 1: Using EasyBuild on Supercomputers

Many HPC systems provide EasyBuild pre-installed or as a loadable module.

### Check if EasyBuild is Available

```
# Check for eb command
which eb

# Or check for EasyBuild module
module avail EasyBuild
```

### Load EasyBuild (if available as module)

```
module load EasyBuild
```

### Verify EasyBuild Installation

```
eb --version
eb --show-system-info
```

### Understanding Toolchains

EasyBuild uses **toolchains** - standardized combinations of compilers and libraries. Common toolchains include:

| Toolchain | Compilers | MPI | Math Libraries |
|-----------|-----------|-----|----------------|
| `foss` | GCC | OpenMPI | OpenBLAS, FFTW, ScaLAPACK |
| `intel` | Intel | Intel MPI | Intel MKL |
| `gompi` | GCC | OpenMPI | - |
| `iimpi` | Intel | Intel MPI | - |

### Search for CHAMP Easyconfig

Check if CHAMP easyconfig exists:

```
# Search for CHAMP easyconfig files
eb -S CHAMP

# Show details of a specific easyconfig
eb CHAMP-2.4.0-foss-2024a.eb --dry-run
```

### Install CHAMP with EasyBuild

**Basic installation with foss toolchain:**

```
eb CHAMP-2.4.0-foss-2024a.eb --robot
```

The `--robot` flag automatically resolves and installs dependencies.

**With Intel toolchain:**

```
eb CHAMP-2.4.0-intel-2024a.eb --robot
```

**Install with TREXIO support:**

If an easyconfig with TREXIO exists:

```
eb CHAMP-2.4.0-foss-2024a-TREXIO.eb --robot
```

**Dry run (see what will be installed):**

```
eb CHAMP-2.4.0-foss-2024a.eb --dry-run
```

### Load CHAMP

After installation, load CHAMP via the generated module:

```
module load CHAMP/2.4.0-foss-2024a
```

### Verify Installation

```
which vmc.mov1
vmc.mov1 --version
module list
```

### Example: Using EasyBuild on HPC Systems

Most HPC centers provide EasyBuild with site-specific configurations:

```
# Load EasyBuild
module load EasyBuild

# List available easyconfigs
eb -S CHAMP

# Install CHAMP
eb CHAMP-2.4.0-foss-2024a.eb --robot

# Load CHAMP
module load CHAMP/2.4.0-foss-2024a

# Verify
which vmc.mov1
vmc.mov1 --version
```

## Scenario 2: Installing EasyBuild from Scratch

If EasyBuild is not available on your system, you can install it in your home directory.

### Prerequisites

Ensure you have:

- **Python** >= 3.6
- **Environment modules** or **Lmod** (for module loading)
- **Internet connection** (for downloading sources)
- **Basic build tools** (gcc, make, etc.)

### Install EasyBuild

**Method 1: Bootstrap script (Recommended)**

```
# Download bootstrap script
wget https://raw.githubusercontent.com/easybuilders/easybuild-framework/develop/easybuild/scripts/bootstrap_eb.py

# Run bootstrap (installs to $HOME/.local/easybuild)
python3 bootstrap_eb.py $HOME/.local/easybuild

# Set up environment
export MODULEPATH=$HOME/.local/easybuild/modules/all:$MODULEPATH
module load EasyBuild
```

**Method 2: Using pip**

```
# Install EasyBuild
pip install --user easybuild

# Verify
eb --version
```

### Configure EasyBuild

**Set up build and install paths:**

Create `~/.easybuild/config.cfg`:

```
[config]
buildpath=/tmp/$USER/easybuild/build
sourcepath=$HOME/.local/easybuild/sources
installpath=$HOME/.local/easybuild
repositorypath=$HOME/.local/easybuild/ebfiles_repo
```

**Or use environment variables:**

```
export EASYBUILD_PREFIX=$HOME/.local/easybuild
export EASYBUILD_BUILDPATH=/tmp/$USER/easybuild/build
export EASYBUILD_INSTALLPATH=$HOME/.local/easybuild
```

Add to `~/.bashrc` for persistence:

```
# EasyBuild configuration
export EASYBUILD_PREFIX=$HOME/.local/easybuild
export MODULEPATH=$HOME/.local/easybuild/modules/all:$MODULEPATH
```

### Install Toolchain

Before installing CHAMP, install a toolchain:

**Install foss toolchain (GCC + OpenMPI + OpenBLAS):**

```
eb foss-2024a.eb --robot
```

This may take 1-2 hours as it builds all toolchain components.

**Or install a minimal toolchain:**

```
# Just GCC and OpenMPI
eb gompi-2022a.eb --robot
```

### Create CHAMP Easyconfig

If CHAMP easyconfig doesn't exist, create one. Save as `CHAMP-2.4.0-foss-2024a.eb`:

```python
easyblock = 'CMakeMake'

name = 'CHAMP'
version = '2.4.0'

homepage = 'https://github.com/filippi-claudia/champ'
description = """CHAMP is a quantum Monte Carlo suite of programs for
electronic structure calculations."""

toolchain = {'name': 'foss', 'version': '2022a'}
toolchainopts = {'usempi': True, 'pic': True}

source_urls = ['https://github.com/filippi-claudia/champ/archive/']
sources = ['v%(version)s.tar.gz']

builddependencies = [
    ('CMake', '3.24.3'),
]

dependencies = [
    ('HDF5', '1.12.2'),
]

configopts = '-DCMAKE_Fortran_COMPILER=mpif90 -DCMAKE_C_COMPILER=mpicc '
configopts += '-DENABLE_MPI=ON '

sanity_check_paths = {
    'files': ['bin/vmc.mov1', 'bin/dmc.mov1'],
    'dirs': ['bin'],
}

moduleclass = 'chem'
```

**For CHAMP with TREXIO support:**

Add to the easyconfig:

```python
dependencies = [
    ('HDF5', '1.12.2'),
    ('TREXIO', '2.4.0'),
]

configopts += '-DENABLE_TREXIO=ON '
```

### Install CHAMP

```
eb CHAMP-2.4.0-foss-2024a.eb --robot
```

This will:
1. Download CHAMP source
2. Build with the foss toolchain
3. Install to `$HOME/.local/easybuild`
4. Generate a module file

### Load and Use CHAMP

```
# Load CHAMP module
module load CHAMP/2.4.0-foss-2024a

# Verify
which vmc.mov1
vmc.mov1 --version

# See loaded dependencies
module list
```

## EasyBuild Options and Flags

Useful EasyBuild command-line options:

| Option | Description |
|--------|-------------|
| `--robot` | Automatically resolve and install dependencies |
| `--dry-run` | Show what would be installed without actually building |
| `--force` | Force rebuild even if already installed |
| `--rebuild` | Rebuild only specified software (not dependencies) |
| `--parallel=N` | Use N parallel processes for building |
| `--trace` | Show detailed trace of build process |

**Example usage:**

```
# Dry run to see what will be built
eb CHAMP-2.4.0-foss-2024a.eb --dry-run

# Parallel build with 8 cores
eb CHAMP-2.4.0-foss-2024a.eb --robot --parallel=8

# Force rebuild
eb CHAMP-2.4.0-foss-2024a.eb --force
```

## Advanced EasyBuild Usage

### Using Different Toolchains

**List available toolchains:**

```
eb --list-toolchains
```

**Install with different toolchain versions:**

```
# With newer foss toolchain
eb CHAMP-2.4.0-foss-2024a.eb --robot

# With Intel toolchain
eb CHAMP-2.4.0-intel-2024a.eb --robot
```

### Customizing Easyconfigs

Create a custom easyconfig by modifying an existing one:

```
# Copy existing easyconfig
cp CHAMP-2.4.0-foss-2024a.eb CHAMP-2.4.0-foss-2024a-custom.eb

# Edit to add options (e.g., enable TREXIO, QMCkl)
# Modify configopts, dependencies, etc.

# Build custom version
eb CHAMP-2.4.0-foss-2024a-custom.eb --robot
```

### Viewing Installation Details

```
# List all installed modules
module avail

# Show module contents
module show CHAMP/2.4.0-foss-2024a

# Find installation directory
eb --search CHAMP
```

### Managing Installations

```
# List installed software
eb --list-installed-software

# Remove a module
rm -rf $EASYBUILD_INSTALLPATH/software/CHAMP/2.4.0-foss-2024a
rm -rf $EASYBUILD_INSTALLPATH/modules/all/CHAMP/2.4.0-foss-2024a.lua
```

## Using CHAMP Installed with EasyBuild

### In Interactive Sessions

```
# Load CHAMP module
module load CHAMP/2.4.0-foss-2024a

# Run calculation
mpirun -np 8 vmc.mov1 -i input.inp -o output.out -e error
```

### In Job Scripts

Add module loading to your job script:

```bash
#!/bin/bash
#SBATCH --job-name=champ-job
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=64
#SBATCH --time=01:00:00

# Load CHAMP module (automatically loads dependencies)
module load CHAMP/2.4.0-foss-2024a

# Run calculation
srun vmc.mov1 -i input.inp -o output.out -e error
```

### Checking Dependencies

```
# See what modules are loaded
module list

# Show module dependencies
module show CHAMP/2.4.0-foss-2024a
```

## Troubleshooting

### EasyBuild command not found

Ensure you've loaded the EasyBuild module or set up your environment:
```
module load EasyBuild
# or
export PATH=$HOME/.local/easybuild/bin:$PATH
```

### Build fails with missing dependencies

Use `--robot` flag to automatically resolve dependencies:
```
eb CHAMP-2.4.0-foss-2024a.eb --robot
```

### Toolchain not found

Install the toolchain first:
```
eb foss-2024a.eb --robot
```

### Module command not found

Install environment modules or Lmod:

**Ubuntu/Debian:**
```
sudo apt-get install environment-modules
```

**Fedora/RHEL:**
```
sudo dnf install environment-modules
```

Then source the initialization script:
```
source /usr/share/modules/init/bash
```

### Build directory full

Clean old build directories:
```
rm -rf /tmp/$USER/easybuild/build/*
```

Or configure a different build path:
```
export EASYBUILD_BUILDPATH=/scratch/$USER/easybuild/build
```

### Installation fails

View detailed log:
```
eb CHAMP-2.4.0-foss-2024a.eb --robot --trace
```

Check log files in:
```
$EASYBUILD_INSTALLPATH/software/CHAMP/2.4.0-foss-2024a/easybuild/*.log
```

## Advantages of Using EasyBuild

- **Standardized toolchains** - Tested compiler and library combinations
- **Automated module generation** - Easy environment management
- **Reproducible builds** - Easyconfig files ensure consistency
- **Wide software support** - Large repository of easyconfigs
- **HPC-focused** - Designed specifically for supercomputing environments
- **Multiple versions** - Install and maintain different versions side-by-side

## Best Practices

1. **Use standard toolchains** - Prefer `foss` or `intel` for compatibility
2. **Keep easyconfigs** - Save custom easyconfigs for reproducibility
3. **Document toolchain versions** - Record which toolchain was used
4. **Test installations** - Run test suite after building
5. **Use robot flag** - Let EasyBuild handle dependencies automatically
6. **Clean build directories** - Regularly clean `/tmp` build files
7. **Share easyconfigs** - Contribute working easyconfigs back to the community

## Creating Easyconfig for Development Version

To build the latest CHAMP from Git:

```python
easyblock = 'CMakeMake'

name = 'CHAMP'
version = 'develop'

homepage = 'https://github.com/filippi-claudia/champ'
description = """CHAMP quantum Monte Carlo suite"""

toolchain = {'name': 'foss', 'version': '2024a'}
toolchainopts = {'usempi': True, 'pic': True}

source_urls = ['https://github.com/filippi-claudia/champ/archive/']
sources = [{
    'download_filename': 'develop.tar.gz',
    'filename': 'CHAMP-develop.tar.gz',
}]

builddependencies = [
    ('CMake', '3.24.3'),
]

dependencies = [
    ('HDF5', '1.12.2'),
    ('TREXIO', '2.4.0'),
]

configopts = '-DCMAKE_Fortran_COMPILER=mpif90 -DCMAKE_C_COMPILER=mpicc '
configopts += '-DENABLE_MPI=ON -DENABLE_TREXIO=ON '

sanity_check_paths = {
    'files': ['bin/vmc.mov1', 'bin/dmc.mov1'],
    'dirs': ['bin'],
}

moduleclass = 'chem'
```

## Additional Resources

- [EasyBuild Documentation](https://easybuild.io/en/latest/)
- [EasyBuild Tutorial](https://tutorial.easybuild.io/)
- [EasyBuild GitHub Repository](https://github.com/easybuilders/easybuild)
- [Easyconfig Repository](https://github.com/easybuilders/easybuild-easyconfigs)

## Next Steps

After installing CHAMP with EasyBuild:

1. Load the module: `module load CHAMP/2.4.0-foss-2024a`
2. Verify with `vmc.mov1 --version`
3. Run the test suite if available
4. Learn about the [Command-Line Interface](../cli/index.md)
5. Follow the [Tutorials](../../tutorials/index.md) for example calculations

