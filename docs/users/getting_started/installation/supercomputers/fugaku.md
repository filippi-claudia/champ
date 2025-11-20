---
title: Running on Fugaku
tags:
    - Fugaku
    - supercomputer
    - ARM A64FX
---

# Running CHAMP on Fugaku Supercomputer

[Fugaku](https://www.r-ccs.riken.jp/en/fugaku/) is Japan's flagship supercomputer located at RIKEN Center for Computational Science. Built on ARM architecture with Fujitsu A64FX processors, Fugaku held the #1 position on the TOP500 list and is optimized for high-performance scientific computing.

## System Overview

**Location**: RIKEN Center for Computational Science (R-CCS), Kobe, Japan  
**Website**: [fugaku.r-ccs.riken.jp](https://www.r-ccs.riken.jp/en/fugaku/)  
**Architecture**: ARM-based (Fujitsu A64FX)  
**Processor**: Fujitsu A64FX (48 cores + 4 assistant cores per node)  
**Interconnect**: Tofu Interconnect D (6D mesh/torus)  
**Memory**: 32 GB HBM2 per node  
**Scheduler**: Fujitsu Technical Computing Suite (PJM - Parallel Job Manager)  
**Software Management**: Spack package manager

## Module Environment Setup

Fugaku uses Spack for software package management. The required modules are available through the centralized Spack installation.

### Initialize Spack Environment

```bash
# Source Spack setup (required in every session)
. /vol0004/apps/oss/spack/share/spack/setup-env.sh
```

Add this to your `~/.bashrc` for automatic loading:

```bash
echo ". /vol0004/apps/oss/spack/share/spack/setup-env.sh" >> ~/.bashrc
```

### Load Required Modules for Compilation

Load the specific versions of CMake, Fujitsu MPI, HDF5, and linear algebra libraries:

```bash
# Load CMake
spack load cmake@3.24.3%fj@4.8.1/p5qsrqc

# Load Fujitsu MPI
spack load fujitsu-mpi@head%fj@4.8.1

# Load HDF5 (required for TREXIO support)
spack load hdf5@1.12.2%fj@4.8.1/tpglq6h

# Load Fujitsu Scientific Subroutine Library (SSL2) for BLAS/LAPACK
spack load fujitsu-ssl2@head%fj@4.8.1/nndozbk
```

## Installing TREXIO and QMCkl (Optional)

If you need TREXIO and QMCkl support, install them in your user space.

### Install TREXIO

```bash
cd $HOME
wget https://github.com/TREX-CoE/trexio/releases/download/v2.6.0/trexio-2.6.0.tar.gz
tar -xzvf trexio-2.6.0.tar.gz
cd trexio-2.6.0

# Configure with Fujitsu compilers
./configure --prefix=$HOME/trexio FC=mpifrt CC=mpifcc
make -j8
make install

export TREXIO_DIR=$HOME/trexio
```

### Install QMCkl

```bash
cd $HOME
wget https://github.com/TREX-CoE/qmckl/releases/download/v1.0.0/qmckl-1.0.0.tar.gz
tar -xzvf qmckl-1.0.0.tar.gz
cd qmckl-1.0.0

# Configure with Fujitsu compilers
./configure --prefix=$HOME/qmckl --enable-hdf5 FC=mpifrt CC=mpifcc
make -j8
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

Set up environment variables for Fujitsu MPI compiler wrappers:

```bash
export MPIFC='mpifrt'    # Fujitsu Fortran MPI wrapper
export MPICC='mpifcc'    # Fujitsu C MPI wrapper
```

**Without TREXIO/QMCkl:**

```bash
cmake -S. -Bbuild \
  -DCMAKE_Fortran_COMPILER=${MPIFC} \
  -DCMAKE_C_COMPILER=${MPICC} \
  -DCMAKE_BUILD_TYPE=Release
```

**With TREXIO/QMCkl:**

```bash
cmake -S. -Bbuild \
  -DCMAKE_Fortran_COMPILER=${MPIFC} \
  -DCMAKE_C_COMPILER=${MPICC} \
  -DCMAKE_BUILD_TYPE=Release \
  -DENABLE_TREXIO=ON \
  -DENABLE_QMCKL=ON \
  -DTREXIO_DIR=$TREXIO_DIR \
  -DQMCKL_DIR=$QMCKL_DIR
```

### Compile CHAMP

```bash
cmake --build build -j --clean-first
```

**Note**: On Fugaku, parallel builds (`-j`) automatically use all available cores on the login node.

Executables will be in `bin/`:

- `vmc.mov1` - VMC executable
- `dmc.mov1` - DMC executable

## Running CHAMP on Fugaku

Fugaku uses the Fujitsu PJM (Parallel Job Manager) for job scheduling. All calculations must be submitted as batch jobs using `pjsub`.

### PJM Resource Groups

Common resource groups on Fugaku:

- `small` - Small-scale jobs (up to 384 nodes)
- `large` - Large-scale jobs (385+ nodes)
- `huge` - Very large-scale jobs (55,297+ nodes)
- `small-s5` - Small-scale short jobs (priority queue)

Check available resource groups:

```bash
pjstat --rsc
```

### Basic PJM Directives

Key PJM directives for job scripts:

| Directive | Description | Example |
|-----------|-------------|---------|
| `#PJM -L node=N` | Number of nodes | `#PJM -L node=4` |
| `#PJM -L rscgrp=GROUP` | Resource group | `#PJM -L rscgrp=small` |
| `#PJM -L elapse=HH:MM:SS` | Wall time limit | `#PJM -L elapse=01:00:00` |
| `#PJM -g GROUP` | Project group ID | `#PJM -g hp230349` |
| `#PJM -x PJM_LLIO_GFSCACHE=PATH` | I/O optimization | `#PJM -x PJM_LLIO_GFSCACHE=/vol0004` |

### Sample VMC Job Script

Create `vmc_job.sh`:

```bash
#!/bin/sh
#PJM -L node=1
#PJM -L rscgrp=small-s5
#PJM -L elapse=00:30:00
#PJM -g hp230349
#PJM -x PJM_LLIO_GFSCACHE=/vol0004
#PJM -j
#PJM --mpi proc=48

# Load Spack environment
. /vol0004/apps/oss/spack/share/spack/setup-env.sh

# Load required modules (use module load instead of spack load in job scripts)
module load fujitsu-mpi/head-fj-4.8.1-gncxc6a
module load fujitsu-ssl2/head-fj-4.8.1-r3hdjbl
module load hdf5/1.12.2-fj-4.8.1-tpglq6h

# Optional: Load TREXIO/QMCkl if compiled with support
# export TREXIO_DIR=$HOME/trexio
# export QMCKL_DIR=$HOME/qmckl
# export LD_LIBRARY_PATH=$TREXIO_DIR/lib:$QMCKL_DIR/lib:$LD_LIBRARY_PATH

# Set paths
CHAMP_BIN=$HOME/champ/bin
INPUT_FILE=vmc.inp
OUTPUT_FILE=vmc.out

# Launch VMC calculation
mpiexec $CHAMP_BIN/vmc.mov1 -i $INPUT_FILE -o $OUTPUT_FILE -e error
```

**Key points**:

- `#PJM --mpi proc=48` - Run 48 MPI processes (full node utilization on A64FX)
- `#PJM -j` - Merge stdout and stderr
- `#PJM -x PJM_LLIO_GFSCACHE=/vol0004` - Enable layered I/O for performance

### Sample DMC Job Script

Create `dmc_job.sh`:

```bash
#!/bin/sh
#PJM -L node=4
#PJM -L rscgrp=small
#PJM -L elapse=02:00:00
#PJM -g hp230349
#PJM -x PJM_LLIO_GFSCACHE=/vol0004
#PJM -j
#PJM --mpi proc=192

# Load Spack environment
. /vol0004/apps/oss/spack/share/spack/setup-env.sh

# Load required modules
module load fujitsu-mpi/head-fj-4.8.1-gncxc6a
module load fujitsu-ssl2/head-fj-4.8.1-r3hdjbl
module load hdf5/1.12.2-fj-4.8.1-tpglq6h

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

# Step 1: Run VMC to generate configurations
echo "Starting VMC calculation..."
mpiexec $CHAMP_BIN/vmc.mov1 -i $VMC_INPUT -o $VMC_OUTPUT -e error

# Step 2: Concatenate configuration files
echo "Merging configuration files..."
cat mc_configs_new* >> mc_configs
rm mc_configs_new*

# Step 3: Run DMC
echo "Starting DMC calculation..."
mpiexec $CHAMP_BIN/dmc.mov1 -i $DMC_INPUT -o $DMC_OUTPUT -e error

echo "DMC calculation completed."
```

**Multi-node configuration**:

- 4 nodes × 48 cores/node = 192 MPI processes
- `#PJM --mpi proc=192` for full utilization

### Submit Jobs

```bash
# Submit VMC job
pjsub vmc_job.sh

# Submit DMC job
pjsub dmc_job.sh

# Check job status
pjstat

# Check specific job
pjstat -v <job_id>

# Cancel a job
pjdel <job_id>
```

### Monitor Job Progress

```bash
# View job statistics
pjstat -v <job_id>

# Check job history
pjstat --history -u $USER

# View output file during execution
tail -f vmc_job.sh.o<job_id>
```

## Performance Optimization

### ARM SVE Optimization

The Fujitsu A64FX processor supports ARM Scalable Vector Extension (SVE). Ensure the compiler uses appropriate optimization flags (typically handled automatically by the Fujitsu compiler).

### Process Binding

For optimal performance, MPI processes should be bound to cores:

```bash
# PJM automatically handles process binding on Fugaku
# Verify with:
export PLE_MPI_STD_EMPTYFILE=off
```

### File System Optimization

Use the layered I/O cache for improved performance:

```bash
#PJM -x PJM_LLIO_GFSCACHE=/vol0004
```

This enables caching for file system operations, reducing I/O overhead.

## Troubleshooting Memory Issues

**Problem**: Out of memory errors

**Solution**: Fugaku nodes have 32 GB HBM2. Reduce walker population or use more nodes:

```bash
#PJM -L node=8  # Increase number of nodes
```

## Additional Resources

- [Fugaku User Guide](https://www.fugaku.r-ccs.riken.jp/doc_root/en/user_guides/)
- [PJM User Guide](https://www.fugaku.r-ccs.riken.jp/doc_root/en/manuals/)
- [Fujitsu Compiler Manual](https://www.fugaku.r-ccs.riken.jp/doc_root/en/programming_guides/)
- [Spack Documentation](https://spack.readthedocs.io/)

## Getting Help

- **Fugaku-specific issues**: Contact R-CCS support through the user portal
- **CHAMP usage**: Consult the [CHAMP documentation](../../../index.md)
- **Bug reports**: Open an issue on [GitHub](https://github.com/filippi-claudia/champ)

## Next Steps

- Review [Command-Line Interface](cli.md) for execution options
- Explore [Input Preparation](../../../preparation/index.md) for setting up calculations
- Try [Tutorials](../../../tutorials/index.md) for practical examples
- Learn about [Analysis Tools](../../../analysis/index.md) for processing results
