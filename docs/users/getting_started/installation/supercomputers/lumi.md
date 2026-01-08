---
title: Running on LUMI
tags:
    - LUMI
    - supercomputer
    - HPE Cray EX
---

# Running CHAMP on LUMI Supercomputer

[LUMI](https://www.lumi-supercomputer.eu/) (Large Unified Modern Infrastructure) is one of Europe's most powerful pre-exascale supercomputers, hosted by CSC in Finland. It features HPE Cray EX architecture with AMD EPYC CPUs and AMD MI250X GPUs, making it ideal for large-scale quantum Monte Carlo calculations.

## System Overview

**Access**: [lumi.csc.fi](https://www.lumi-supercomputer.eu/)  
**Architecture**: HPE Cray EX  
**Compute Nodes**:

- **LUMI-C** (CPU partition): AMD EPYC 7763 (Milan), 128 cores/node
- **LUMI-G** (GPU partition): AMD EPYC 7A53 + 4× AMD MI250X GPUs

**Scheduler**: SLURM  
**Programming Environment**: Cray Programming Environment (CPE)  
**File Systems**: Lustre parallel file system

## Module Environment Setup

LUMI uses the Cray Programming Environment with multiple compiler options. For CHAMP, we recommend using the GNU programming environment.

### Load Required Modules

```bash
# Switch from default Cray to GNU programming environment
module swap PrgEnv-cray PrgEnv-gnu

# Load LUMI software stack
module load LUMI

# Load parallel HDF5 (required for TREXIO)
module load cray-hdf5-parallel

# Load build tools (includes CMake)
module load buildtools

# Target LUMI-C partition
module load craype-x86-milan
```

!!! note
    Module versions may change. Check available versions with:

### Verify Loaded Modules

```bash
module list
```

You should see:

- `PrgEnv-gnu`
- `LUMI`
- `cray-hdf5-parallel`
- `buildtools`
- `craype-x86-milan`

## Installing TREXIO and QMCkl

If you need TREXIO and QMCkl support, install them in your project scratch space:

### Install TREXIO

```bash
cd /scratch/projectX
wget https://github.com/TREX-CoE/trexio/releases/download/v2.6.0/trexio-2.6.0.tar.gz
tar -xzvf trexio-2.6.0.tar.gz
cd trexio-2.6.0

# Optimization flags
OPT_FLAGS="-O3 -march=native -flto -fno-trapping-math -fno-math-errno -funroll-loops"

./configure --prefix=/scratch/projectX/trexio \
    FC=ftn \
    CC=cc \
    CFLAGS="${OPT_FLAGS}" \
    FCFLAGS="${OPT_FLAGS}"

make -j 64
make check
make install

export TREXIO_DIR=/scratch/projectX/trexio
```

### Install QMCkl

```bash
cd /scratch/projectX
wget https://github.com/TREX-CoE/qmckl/releases/download/v1.0.0/qmckl-1.0.0.tar.gz
tar -xzvf qmckl-1.0.0.tar.gz
cd qmckl-1.0.0

# Let QMCkl find trexio library
export TREXIO_DIR=/scratch/projectX/trexio
export CPPFLAGS="-I${TREXIO_DIR}/include"
export LDFLAGS="-L${TREXIO_DIR}/lib"

# Note: HDF5 and MPI libraries are handled implicitly by 'cc'/'ftn' wrappers
export TREXIO_LIBS="-L${TREXIO_DIR}/lib -ltrexio"
export TREXIO_CFLAGS="-I${TREXIO_DIR}/include"

./configure --prefix=/scratch/projectX/qmckl \
    --enable-hpc \
    FC=ftn \
    CC=cc \
    CFLAGS="${OPT_FLAGS}" \
    FCFLAGS="${OPT_FLAGS}" \
    LDFLAGS="${LDFLAGS} -O3 -flto"

make -j 64
make install

export QMCKL_DIR=/scratch/projectX/qmckl
```

Add to your `~/.bashrc`:

```bash
export TREXIO_DIR=/scratch/projectX/trexio
export QMCKL_DIR=/scratch/projectX/qmckl
export PATH=$TREXIO_DIR/bin:$QMCKL_DIR/bin:$PATH
export LD_LIBRARY_PATH=$TREXIO_DIR/lib:$QMCKL_DIR/lib:$LD_LIBRARY_PATH
```

## Building CHAMP

### Clone CHAMP Repository

```bash
cd /scratch/projectX
git clone https://github.com/filippi-claudia/champ.git
cd champ
```

### Configure the Build

Use Cray compiler wrappers `ftn` (Fortran) and `cc` (C), which automatically handle MPI and library linking:

**Without TREXIO/QMCkl:**

```bash

# Optimization Flags (Global)
OPT_FLAGS="-O3 -march=native -flto -fno-trapping-math -fno-math-errno -funroll-loops"

FORTRAN_FLAGS="${OPT_FLAGS} -fallow-argument-mismatch -Wno-argument-mismatch -w"

cmake -S. -Bbuild \
  -DCMAKE_Fortran_COMPILER=ftn \
  -DCMAKE_C_COMPILER=cc \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_Fortran_FLAGS="${FORTRAN_FLAGS}" \
  -DCMAKE_C_FLAGS="${OPT_FLAGS}" \
  -DCMAKE_EXE_LINKER_FLAGS="-flto"
```

**With TREXIO/QMCkl:**

```bash
# Optimization Flags (Global)
OPT_FLAGS="-O3 -march=native -flto -fno-trapping-math -fno-math-errno -funroll-loops"

FORTRAN_FLAGS="${OPT_FLAGS} -fallow-argument-mismatch -Wno-argument-mismatch -w"

cmake -S. -Bbuild \
    -DCMAKE_Fortran_COMPILER=ftn \
    -DCMAKE_C_COMPILER=cc \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_Fortran_FLAGS="${FORTRAN_FLAGS}" \
    -DCMAKE_C_FLAGS="${OPT_FLAGS}" \
    -DCMAKE_EXE_LINKER_FLAGS="-flto" \
    -DENABLE_TREXIO=ON \
    -DTREXIO_INCLUDE_DIR="${TREXIO_DIR}/include" \
    -DTREXIO_LIBRARY="${TREXIO_DIR}/lib/libtrexio.so" \
    -DENABLE_QMCKL=ON \
    -DQMCKL_INCLUDE_DIR="${QMCKL_DIR}/include" \
    -DQMCKL_LIBRARY="${QMCKL_DIR}/lib/libqmckl.so"
```

### Compile CHAMP

```bash
cmake --build build -j 64 --clean-first
```

Executables will be in `bin/`:

- `vmc.mov1` - VMC executable
- `dmc.mov1` - DMC executable

## Running CHAMP on LUMI

LUMI uses SLURM for job scheduling. All production calculations must be submitted as batch jobs.

### SLURM Partitions

Common partitions on LUMI-C:

- `small` - Up to 128 cores (1 node)
- `standard` - 129 to 512 nodes
- `large` - 513+ nodes
- `dev` - Development/testing (limited time)

Check partition details:

```bash
sinfo -s
```

### Environment Variables

Set this before running CHAMP to avoid MPI initialization issues:

```bash
export PMI_NO_PREINITIALIZE=y
```

### Sample VMC Job Script

Create `vmc_job.sh`:

```bash
#!/bin/bash -l
#SBATCH --job-name=champ_vmc      # Job name
#SBATCH --output=vmc.o%j          # Standard output file
#SBATCH --error=vmc.e%j           # Standard error file
#SBATCH --partition=standard      # Use small/standard partition
#SBATCH --nodes=10                # Number of nodes
#SBATCH --ntasks-per-node=128     # MPI tasks per node
#SBATCH --cpus-per-task=1         # CPUs per task
#SBATCH --time=00:30:00           # Wall time (HH:MM:SS)
#SBATCH --account=project_XXXXXX  # Replace with your project ID
#SBATCH --mail-type=END,FAIL      # Email notifications

# Load modules (same as compilation)
module purge
module load LUMI/24.03
module load buildtools/24.03
module load PrgEnv-gnu/8.5.0
module load cray-hdf5-parallel/1.12.2.11
module load craype-x86-milan # Ensure CPU architecture is set to Milan

# Optional: Load TREXIO/QMCkl if compiled with support
# export TREXIO_DIR=/scratch/projectX/trexio
# export QMCKL_DIR=/scratch/projectX/qmckl
# export LD_LIBRARY_PATH=$TREXIO_DIR/lib:$QMCKL_DIR/lib:$LD_LIBRARY_PATH

# Set paths
CHAMP_BIN=/scratch/projectX/champ/bin
INPUT_FILE=vmc.inp
OUTPUT_FILE=vmc.out

# Launch VMC calculation
srun $CHAMP_BIN/vmc.mov1 -i $INPUT_FILE -o $OUTPUT_FILE -e error
```

### Sample DMC Job Script

Create `dmc_job.sh`:

```bash
#!/bin/bash -l
#SBATCH --job-name=champ_dmc      # Job name
#SBATCH --output=dmc.o%j          # Standard output file
#SBATCH --error=dmc.e%j           # Standard error file
#SBATCH --partition=standard      # Use standard partition
#SBATCH --nodes=4                 # Number of nodes
#SBATCH --ntasks-per-node=128     # MPI tasks per node (512 total)
#SBATCH --cpus-per-task=1         # CPUs per task
#SBATCH --time=02:00:00           # Wall time (HH:MM:SS)
#SBATCH --account=project_XXXXXX  # Replace with your project ID
#SBATCH --mail-type=END,FAIL      # Email notifications

# Load modules
module swap PrgEnv-cray PrgEnv-gnu
module load LUMI
module load cray-hdf5-parallel
module load buildtools

# Optional: Load TREXIO/QMCkl
# export TREXIO_DIR=/scratch/projectX/trexio
# export QMCKL_DIR=/scratch/projectX/qmckl
# export LD_LIBRARY_PATH=$TREXIO_DIR/lib:$QMCKL_DIR/lib:$LD_LIBRARY_PATH

# Set paths
CHAMP_BIN=/scratch/projectX/champ/bin
VMC_INPUT=vmc.inp
VMC_OUTPUT=vmc.out
DMC_INPUT=dmc.inp
DMC_OUTPUT=dmc.out

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

# Cancel a job
scancel <job_id>
```

### Monitor Job Progress

```bash
# View output in real-time
tail -f vmc.o<jobid>

# Check job efficiency after completion
seff <jobid>

# View detailed job information
sacct -j <jobid> --format=JobID,JobName,Partition,State,Elapsed,MaxRSS
```

## Performance Optimization

### CPU Binding

For optimal performance, bind MPI processes to CPU cores:

```bash
srun --cpu-bind=cores $CHAMP_BIN/vmc.mov1 -i vmc.inp -o vmc.out -e error
```

### File System Considerations

- Use `/scratch` for temporary files and job output (fast, not backed up)
- Use `/project` for long-term storage (backed up, slower)
- Avoid excessive small file I/O operations

```bash
# Set scratch directory in job script
export SCRATCH=/scratch/project_XXXXXX/$USER
cd $SCRATCH
```

## Additional Resources

[:fontawesome-solid-server: SLURM Documentation   :fontawesome-solid-arrow-up-right-from-square:](https://slurm.schedmd.com/)

[:fontawesome-solid-book: LUMI Documentation   :fontawesome-solid-arrow-up-right-from-square:](https://docs.lumi-supercomputer.eu/firststeps/)


## Getting Help

- **LUMI-specific issues**: Contact LUMI support at [lumi-helpdesk@lumi-supercomputer.eu](mailto:lumi-helpdesk@lumi-supercomputer.eu)
- **Bug reports**: Open an issue on [GitHub](https://github.com/filippi-claudia/champ)

## Next Steps

- Review [Command-Line Interface](cli.md) for execution options
- Explore [Input Preparation](../../../preparation/index.md) for setting up calculations
- Try [Tutorials](../../../tutorials/index.md) for practical examples
- Learn about [Analysis Tools](../../../analysis/index.md) for processing results