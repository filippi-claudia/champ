---
title: Running CHAMP
tags:
    - usage
    - execution
---

# Running CHAMP

Once CHAMP is successfully installed, you're ready to perform quantum Monte Carlo calculations. This section provides comprehensive guidance on running CHAMP across different computing environments, from interactive sessions on desktop workstations to large-scale parallel jobs on supercomputers.

## Overview

Running CHAMP involves:

1. **Understanding the command-line interface** - CHAMP executables and their flags
2. **Preparing input files** - Configuration files for VMC/DMC calculations
3. **Executing calculations** - Running CHAMP interactively or via job schedulers
4. **Platform-specific considerations** - Module loading, compiler environments, and job submission

## Running on Different Platforms

### Supercomputers & HPC Clusters

Large-scale production calculations typically run on high-performance computing systems using job schedulers (SLURM, PBS, PJM). Each system has specific requirements for:

- Module environment setup
- Compiler toolchains
- MPI configurations
- Job submission syntax
- Resource allocation

#### Supported Systems

We provide detailed guides for the following supercomputing facilities:

| System | Location | Scheduler | Guide |
|--------|----------|-----------|-------|
| **LUMI** | CSC, Finland | SLURM | [LUMI Guide](lumi.md) |
| **Fugaku** | RIKEN, Japan | PJM | [Fugaku Guide](fugaku.md) |
| **Snellius** | SURF, Netherlands | SLURM | [Snellius Guide](snellius.md) |
| **CCPHead** | University of Twente | SLURM | [CCPHead Guide](ccphead.md) |

Each guide includes:

- ✓ Required module loading commands
- ✓ Environment configuration
- ✓ Compilation instructions
- ✓ Sample job submission scripts
- ✓ Platform-specific optimizations

## Next Steps

After setting up CHAMP execution:

1. **Learn the CLI** - Review [command-line options](cli.md) and execution modes
2. **Choose your platform** - Follow the appropriate system-specific guide

## Getting Help

- **Platform-specific questions**: Contact your HPC support team
- **CHAMP usage**: Consult the [documentation](../../index.md)
- **Bug reports**: Open an issue on [GitHub](https://github.com/filippi-claudia/champ)
- **General questions**: Check the [FAQ](../../troubleshooting/faq.md)
