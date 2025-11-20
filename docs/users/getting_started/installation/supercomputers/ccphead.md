---
title: Running on CCPHead
tags:
    - CCPHead
    - University of Twente
    - cluster
---

# Running CHAMP on CCPHead Cluster

[CCPHead](https://www.utwente.nl/en/tnw/ccp/infrastructure/) is the computing cluster of the Computational Chemical Physics (CCP) group at the University of Twente, Netherlands. This cluster provides computational resources for quantum chemistry and physics calculations, including quantum Monte Carlo simulations with CHAMP.

## System Overview

**Location**: Faculty of Science and Technology, University of Twente, Enschede, Netherlands  
**Website**: [www.utwente.nl/en/tnw/ccp/infrastructure](https://www.utwente.nl/en/tnw/ccp/infrastructure/)  
**Hostname**: `ccphead.tnw.utwente.nl`  
**Architecture**: x86_64  
**Available Compilers**:

- Intel Compiler Suite (ifort, icc) with Intel MKL
- GNU Compiler Collection (gfortran, gcc)

**Available MPI**: Intel MPI, OpenMPI  
**Scheduler**: SLURM

## Module Environment Setup

CCPHead uses a module system to manage software environments. Two main compilation options are available: Intel compilers with Intel MKL, or GNU compilers with system BLAS/LAPACK.

### Option 1: Intel Compilers (Recommended)

The Intel compiler toolchain provides optimized performance with Intel MKL for BLAS/LAPACK operations.

#### Load Intel Modules

```bash
# Load Intel compiler
module load compiler

# Load compiler runtime libraries
module load compiler-rt

# Load Intel Math Kernel Library (MKL)
module load mkl

# Load Intel MPI
module load mpi
```

#### Optional Modules

```bash
# Load TREXIO (if available)
module load trexio/2.3.0-intel

# Load Python (for analysis scripts)
module load python3
```

#### Verify Loaded Modules

```bash
module list
```

You should see:

- `compiler`
- `compiler-rt`
- `mkl`
- `mpi`
- Optional: `trexio`, `python3`

### Option 2: GNU Compilers

The GNU compiler toolchain uses system BLAS/LAPACK libraries from the Ubuntu repository.

#### Using System GNU Compilers

No modules are required for basic GNU compilation, as system compilers are available by default:

```bash
# Verify compiler availability
which gfortran  # Should show /usr/bin/gfortran
which mpif90    # Should show /usr/bin/mpif90
```

**Note**: Combining gfortran with Intel MKL requires special handling with the `-mcmodel=large` compiler flag and may cause compatibility issues. It's recommended to use either full Intel toolchain or full GNU toolchain.

## Building CHAMP

### Clone CHAMP Repository

```bash
cd $HOME
git clone https://github.com/filippi-claudia/champ.git
cd champ
```

### Option 1: Build with Intel Compilers

**Without TREXIO:**

```bash
cmake -S. -Bbuild \
  -DCMAKE_Fortran_COMPILER=mpiifort \
  -DCMAKE_C_COMPILER=mpiicc \
  -DCMAKE_BUILD_TYPE=Release
```

**With TREXIO:**

```bash
cmake -S. -Bbuild \
  -DCMAKE_Fortran_COMPILER=mpiifort \
  -DCMAKE_C_COMPILER=mpiicc \
  -DCMAKE_BUILD_TYPE=Release \
  -DENABLE_TREXIO=ON
```


!!! note
    If TREXIO module is loaded, CMake should automatically detect it. Otherwise, specify:
    ```bash
    -DTREXIO_DIR=/path/to/trexio
    ```

### Option 2: Build with GNU Compilers

```bash
cmake -S. -Bbuild \
  -DCMAKE_Fortran_COMPILER=/usr/bin/mpif90 \
  -DCMAKE_C_COMPILER=/usr/bin/mpicc \
  -DCMAKE_BUILD_TYPE=Release
```

CMake will automatically find LAPACK and BLAS from the Ubuntu repository if Intel MKL environment variables are not set.

### Compile CHAMP

```bash
cmake --build build -j8 --clean-first
```

Executables will be in `build/bin/`:

- `vmc.mov1` - VMC executable
- `dmc.mov1` - DMC executable

### Verify Installation

```bash
ls -lh build/bin/
```

## Running CHAMP on CCPHead

CCPHead supports both interactive execution and job submission via SLURM scheduler (if configured).

### Interactive Execution

For small test calculations or debugging, you can run CHAMP interactively.

#### Using Intel MPI

```bash
# Load required modules
module load compiler compiler-rt mkl mpi

# Run on single node (adjust -np for number of cores)
mpirun -np 32 $HOME/champ/bin/vmc.mov1 -i vmc.inp -o vmc.out -e error
```

#### Using OpenMPI with Machine File

For multi-node execution with OpenMPI:

```bash
# Create machine file listing available nodes/cores
# Example machinefile content:
# node1 slots=24
# node2 slots=24
# node3 slots=12

mpirun -np 60 -machinefile machinefile \
  $HOME/champ/bin/vmc.mov1 -i vmc.inp -o vmc.out -e error
```

### SLURM Job Submission

If CCPHead uses SLURM, submit jobs as batch scripts.

#### Sample VMC Job Script

Create `vmc_job.sh`:

```bash
#!/bin/bash
#SBATCH --job-name=champ_vmc      # Job name
#SBATCH --output=vmc.%j.out       # Standard output file
#SBATCH --error=vmc.%j.err        # Standard error file
#SBATCH --nodes=1                 # Number of nodes
#SBATCH --ntasks=32               # Total number of MPI tasks
#SBATCH --ntasks-per-node=32      # MPI tasks per node
#SBATCH --time=12:00:00           # Wall time (HH:MM:SS)

# Load Intel modules
module load compiler
module load compiler-rt
module load mkl
module load mpi

# Optional: Load TREXIO
# module load trexio/2.3.0-intel

# Set paths
CHAMP_BIN=$HOME/champ/build/bin
INPUT_FILE=vmc.inp
OUTPUT_FILE=vmc.out

# Launch VMC calculation
mpirun -np $SLURM_NTASKS $CHAMP_BIN/vmc.mov1 -i $INPUT_FILE -o $OUTPUT_FILE -e error
```

#### Sample DMC Job Script

Create `dmc_job.sh`:

```bash
#!/bin/bash
#SBATCH --job-name=champ_dmc      # Job name
#SBATCH --output=dmc.%j.out       # Standard output file
#SBATCH --error=dmc.%j.err        # Standard error file
#SBATCH --nodes=2                 # Number of nodes
#SBATCH --ntasks=64               # Total number of MPI tasks
#SBATCH --ntasks-per-node=32      # MPI tasks per node
#SBATCH --time=48:00:00           # Wall time (48 hours)

# Load Intel modules
module load compiler
module load compiler-rt
module load mkl
module load mpi

# Optional: Load TREXIO
# module load trexio/2.3.0-intel

# Set paths
CHAMP_BIN=$HOME/champ/bin
VMC_INPUT=vmc.inp
VMC_OUTPUT=vmc.out
DMC_INPUT=dmc.inp
DMC_OUTPUT=dmc.out

# Step 1: Run VMC to generate configurations
echo "Starting VMC calculation..."
mpirun -np $SLURM_NTASKS $CHAMP_BIN/vmc.mov1 -i $VMC_INPUT -o $VMC_OUTPUT -e error

# Step 2: Concatenate configuration files
echo "Merging configuration files..."
cat mc_configs_new* >> mc_configs
rm mc_configs_new*

# Step 3: Run DMC
echo "Starting DMC calculation..."
mpirun -np $SLURM_NTASKS $CHAMP_BIN/dmc.mov1 -i $DMC_INPUT -o $DMC_OUTPUT -e error

echo "DMC calculation completed."
```

#### Submit Jobs

```bash
# Submit VMC job
sbatch vmc_job.sh

# Submit DMC job
sbatch dmc_job.sh

# Check job status
squeue -u $USER

# Cancel a job
scancel <job_id>
```

### Direct Execution (Non-SLURM)

If SLURM is not configured, execute directly with `mpirun`:

```bash
# Single node (24 cores)
mpirun -np 24 ./bin/vmc.mov1 -i vmc.inp -o vmc.out -e error

# Multi-node with machine file
mpirun -np 60 -machinefile machinefile \
  ./bin/dmc.mov1 -i dmc.inp -o dmc.out -e error
```

## Performance Optimization

### Intel MKL Threading

Control MKL thread count for hybrid MPI+OpenMP (if applicable):

```bash
export MKL_NUM_THREADS=1  # Disable MKL threading for pure MPI
export OMP_NUM_THREADS=1  # Disable OpenMP threading
```

### Process Binding

For better performance, bind MPI processes to cores:

```bash
# Intel MPI
export I_MPI_PIN=1
export I_MPI_PIN_DOMAIN=auto

# OpenMPI
mpirun --bind-to core -np 24 ./bin/vmc.mov1 -i vmc.inp -o vmc.out -e error
```

## Additional Resources

- [CCP Infrastructure Page](https://www.utwente.nl/en/tnw/ccp/infrastructure/)
- [Intel MPI Documentation](https://www.intel.com/content/www/us/en/docs/mpi-library/developer-guide-linux/current/overview.html)
- [OpenMPI Documentation](https://www.open-mpi.org/doc/)

## Getting Help

- **CHAMP usage**: Consult the [CHAMP documentation](../../../index.md)
- **Bug reports**: Open an issue on [GitHub](https://github.com/filippi-claudia/champ)

## Next Steps

- Review [Command-Line Interface](cli.md) for execution options
- Explore [Input Preparation](../../../preparation/index.md) for setting up calculations
- Try [Tutorials](../../../tutorials/index.md) for practical examples
- Learn about [Analysis Tools](../../../analysis/index.md) for processing results

