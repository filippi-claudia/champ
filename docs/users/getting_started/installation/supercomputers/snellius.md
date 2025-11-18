---
title: Running on Snellius
icon: material/server
tags:
    - Snellius
    - supercomputer
    - Intel
---

# Running CHAMP on Snellius Supercomputer

[Snellius](https://www.surf.nl/en/services/snellius-the-national-supercomputer) is the Dutch national supercomputer operated by SURF, providing high-performance computing resources for academic research in the Netherlands. The system features a hybrid architecture with both Intel and AMD processors, making it suitable for a wide range of computational chemistry and physics applications.

## System Overview

**Location**: SURF, Amsterdam, Netherlands  
**Website**: [snellius.surf.nl](https://servicedesk.surf.nl/wiki/display/WIKI/Snellius)  
**Architecture**: Hybrid Intel/AMD x86_64  
**Compute Nodes**:

| Node Type | Processor | Cores/Node | Memory |
|-----------|-----------|------------|--------|
| Thin Rome | 2× AMD Rome 7H12 | 128 | 256 GB |
| Thin Genoa| 2× AMD Genoa 9654 | 192 | 384 GB |
| Fat Rome  | 2× AMD Rome 7H12 | 128 | 1 TB |
| Fat Genoa | 2× AMD Genoa 9654 | 192 | 1.5 TB |

**Interconnect**: Mellanox HDR InfiniBand (200 Gbit/s)  
**Scheduler**: SLURM  
**Software Stack**: EasyBuild-based module system

## Module Environment Setup

Snellius uses an EasyBuild-based module system organized by year and toolchain. For CHAMP, we'll use the Intel toolchain 2024a.

### Load Required Modules for Compilation

```bash
# Clear any existing modules
module purge

# Load 2024a software stack
module load 2024a

# Load Intel toolchain (includes Intel compilers and Intel MPI)
module load intel/2024a

# Load CMake
module load CMake/3.29.3-GCCcore-13.3.0

# Load HDF5 with Intel MPI (required for TREXIO support)
module load HDF5/1.14.3-iimpi-2024a
```

### Available Compilers

The `intel/2024a` toolchain includes:

- **Intel Fortran**: `ifort` (Intel Fortran Compiler Classic)
- **Intel MPI wrapper**: `mpiifort` (MPI Fortran wrapper)
- **Intel C/C++**: `icc`/`icpc` (Intel C/C++ Compilers)
- **Intel MKL**: Math Kernel Library (BLAS/LAPACK)

## Installing TREXIO and QMCkl (Optional)

If you need TREXIO and QMCkl support, install them in your user space.

### Install TREXIO

```bash
cd $HOME
wget https://github.com/TREX-CoE/trexio/releases/download/v2.6.0/trexio-2.6.0.tar.gz
tar -xzvf trexio-2.6.0.tar.gz
cd trexio-2.6.0

# Configure with Intel compilers
./configure --prefix=$HOME/trexio FC=mpiifort CC=mpiicc
make -j16
make install

export TREXIO_DIR=$HOME/trexio
```

### Install QMCkl

```bash
cd $HOME
wget https://github.com/TREX-CoE/qmckl/releases/download/v1.0.0/qmckl-1.0.0.tar.gz
tar -xzvf qmckl-1.0.0.tar.gz
cd qmckl-1.0.0

# Configure with Intel compilers
./configure --prefix=$HOME/qmckl --enable-hdf5 FC=mpiifort CC=mpiicc
make -j16
make install

export QMCKL_DIR=$HOME/qmckl
```

Add to your `~/.bashrc`:

```bash
export TREXIO_DIR=$HOME/trexio
export QMCKL_DIR=$HOME/qmckl
export PATH=$TREXIO_DIR/bin:$QMCKL_DIR/bin:$PATH
export LD_LIBRARY_PATH=$TREXIO_DIR/lib:$QMCKL_DIR/lib:$LD_LIBRARY_PATH
```

## Building CHAMP

### Clone CHAMP Repository

```bash
cd $HOME
git clone https://github.com/filippi-claudia/champ.git
cd champ
```

### Configure the Build

**Without TREXIO/QMCkl:**

```bash
cmake -S. -Bbuild \
  -DCMAKE_Fortran_COMPILER=mpiifort \
  -DCMAKE_C_COMPILER=mpiicc \
  -DCMAKE_BUILD_TYPE=Release
```

**With TREXIO/QMCkl:**

```bash
cmake -S. -Bbuild \
  -DCMAKE_Fortran_COMPILER=mpiifort \
  -DCMAKE_C_COMPILER=mpiicc \
  -DCMAKE_BUILD_TYPE=Release \
  -DENABLE_TREXIO=ON \
  -DENABLE_QMCKL=ON \
  -DTREXIO_DIR=$TREXIO_DIR \
  -DQMCKL_DIR=$QMCKL_DIR
```

### Compile CHAMP

```bash
cmake --build build -j16 --clean-first
```

Executables will be in `bin/`:

- `vmc.mov1` - VMC executable
- `dmc.mov1` - DMC executable

### Verify Installation

```bash
ls -lh bin/
```

## Running CHAMP on Snellius

Snellius uses SLURM for job scheduling. All production calculations must be submitted as batch jobs.

### Important Environment Variables

Set this in job scripts for proper Intel MPI operation:

```bash
export I_MPI_PMI_LIBRARY=/usr/lib64/libpmi2.so
```

This ensures Intel MPI uses the correct PMI library for SLURM integration.

### Sample VMC Job Script

Create `vmc_job.sh`:

```bash
#!/bin/bash
#SBATCH --job-name=champ_vmc      # Job name
#SBATCH --output=vmc.%j.out       # Standard output file (%j = job ID)
#SBATCH --error=vmc.%j.err        # Standard error file
#SBATCH --partition=rome          # Partition (queue)
#SBATCH --nodes=1                 # Number of nodes
#SBATCH --ntasks=128              # Total number of MPI tasks
#SBATCH --ntasks-per-node=128     # MPI tasks per node
#SBATCH --cpus-per-task=1         # CPUs per task
#SBATCH --time=0-12:00:00         # Wall time (days-hours:min:sec)
#SBATCH --exclusive               # Exclusive node access
#SBATCH --account=PROJECTID       # Replace with your project ID

# Load modules (same as compilation)
module purge
module load 2024a
module load intel/2024a
module load HDF5/1.14.3-iimpi-2024a

# Set Intel MPI environment
export I_MPI_PMI_LIBRARY=/usr/lib64/libpmi2.so

# Optional: Load TREXIO/QMCkl if compiled with support
# export TREXIO_DIR=$HOME/trexio
# export QMCKL_DIR=$HOME/qmckl
# export LD_LIBRARY_PATH=$TREXIO_DIR/lib:$QMCKL_DIR/lib:$LD_LIBRARY_PATH

# Set paths
CHAMP_BIN=$HOME/champ/bin
INPUT_FILE=vmc.inp
OUTPUT_FILE=vmc.out

# Change to working directory
cd $SLURM_SUBMIT_DIR

# Launch VMC calculation
srun $CHAMP_BIN/vmc.mov1 -i $INPUT_FILE -o $OUTPUT_FILE -e error
```

### Sample DMC Job Script

Create `dmc_job.sh`:

```bash
#!/bin/bash
#SBATCH --job-name=champ_dmc      # Job name
#SBATCH --output=dmc.%j.out       # Standard output file (%j = job ID)
#SBATCH --error=dmc.%j.err        # Standard error file
#SBATCH --partition=rome          # Partition (queue)
#SBATCH --nodes=2                 # Number of nodes
#SBATCH --ntasks=256              # Total number of MPI tasks (2×128)
#SBATCH --ntasks-per-node=128     # MPI tasks per node
#SBATCH --cpus-per-task=1         # CPUs per task
#SBATCH --time=5-00:00:00         # Wall time (5 days)
#SBATCH --exclusive               # Exclusive node access
#SBATCH --account=PROJECTID       # Replace with your project ID

# Load modules
module purge
module load 2024a
module load intel/2024a
module load HDF5/1.14.3-iimpi-2024a

# Set Intel MPI environment
export I_MPI_PMI_LIBRARY=/usr/lib64/libpmi2.so

# Optional: Load TREXIO/QMCkl
# export TREXIO_DIR=$HOME/trexio
# export QMCKL_DIR=$HOME/qmckl
# export LD_LIBRARY_PATH=$TREXIO_DIR/lib:$QMCKL_DIR/lib:$LD_LIBRARY_PATH

# Set paths
CHAMP_BIN=$HOME/champ/bin
VMC_INPUT=vmc.inp
VMC_OUTPUT=vmc.out
DMC_INPUT=dmc.inp
DMC_OUTPUT=dmc.out

# Change to working directory
cd $SLURM_SUBMIT_DIR

# Step 1: Run VMC to generate configurations
echo "Starting VMC calculation..."
srun $CHAMP_BIN/vmc.mov1 -i $VMC_INPUT -o $VMC_OUTPUT -e error

# Step 2: Concatenate configuration files
echo "Merging configuration files..."
cat mc_configs_new* >> mc_configs
rm mc_configs_new*

# Step 3: Run DMC
echo "Starting DMC calculation..."
srun $CHAMP_BIN/dmc.mov1 -i $DMC_INPUT -o $DMC_OUTPUT -e error

echo "DMC calculation completed."
```

### Submit Jobs

```bash
# Submit VMC job
sbatch vmc_job.sh

# Submit DMC job
sbatch dmc_job.sh

# Check job status
squeue -u $USER

# Check detailed job information
scontrol show job <job_id>

# Cancel a job
scancel <job_id>
```

### Monitor Job Progress

```bash
# View output in real-time
tail -f vmc.<jobid>.out

# Check job efficiency after completion
seff <job_id>

# View detailed accounting information
sacct -j <job_id> --format=JobID,JobName,Partition,State,Elapsed,MaxRSS,MaxVMSize
```

## Performance Optimization

### Intel MPI Tuning

For optimal performance with Intel MPI on Snellius:

```bash
# Add to job script
export I_MPI_PIN=1                     # Enable process pinning
export I_MPI_PIN_DOMAIN=auto           # Automatic domain detection
export I_MPI_FABRICS=shm:ofi           # Use shared memory + OFI
```

### CPU Binding

Use SLURM's CPU binding for better performance:

```bash
srun --cpu-bind=cores $CHAMP_BIN/vmc.mov1 -i vmc.inp -o vmc.out -e error
```

### File System Considerations

Snellius uses different file systems:

- `$HOME` - Home directory (limited space, backed up)
- `$TMPDIR` - Local scratch per job (fast, deleted after job)
- `/scratch-shared` - Shared scratch (fast, not backed up, periodic cleanup)

Use scratch for temporary files:

```bash
# In job script
export SCRATCH=/scratch-shared/$USER
mkdir -p $SCRATCH
cd $SCRATCH
# Copy input files, run calculations, copy results back
```

## Troubleshooting

### Intel MPI Initialization Errors

**Problem**: MPI fails to initialize

**Solution**: Set PMI library:

```bash
export I_MPI_PMI_LIBRARY=/usr/lib64/libpmi2.so
```

## Additional Resources

- [Snellius Documentation](https://servicedesk.surf.nl/wiki/display/WIKI/Snellius)
- [SURF User Support](https://servicedesk.surf.nl/wiki/display/WIKI/SURF+Research+IT+Documentation)
- [SLURM Documentation](https://slurm.schedmd.com/)
- [Intel MPI Documentation](https://www.intel.com/content/www/us/en/docs/mpi-library/developer-guide-linux/2021-8/overview.html)
- [EasyBuild Documentation](https://docs.easybuild.io/)

## Getting Help

- **Snellius-specific issues**: Contact SURF support via [servicedesk.surf.nl](https://servicedesk.surf.nl/)
- **CHAMP usage**: Consult the [CHAMP documentation](../../../index.md)
- **Bug reports**: Open an issue on [GitHub](https://github.com/filippi-claudia/champ)

## Next Steps

- Review [Command-Line Interface](cli.md) for execution options
- Explore [Input Preparation](../../../preparation/index.md) for setting up calculations
- Try [Tutorials](../../../tutorials/index.md) for practical examples
- Learn about [Analysis Tools](../../../analysis/index.md) for processing results