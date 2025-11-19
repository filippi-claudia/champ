---
title: Using TREXIO Files
icon: material/file-code
tags:
    - TREXIO
    - HDF5
    - input format
---

# Using TREXIO Files as Input

[TREXIO](https://trex-coe.github.io/trexio/) (TREX Input/Output) is a modern, standardized file format for storing quantum chemistry system and wavefunction data. CHAMP supports reading wavefunction information directly from TREXIO files in either HDF5 or text backend format, simplifying input preparation, reducing errors, and massive improvements in I/O overheads.

## Overview

TREXIO files can contain all wavefunction data needed for CHAMP calculations except Jastrow factors:

**Included in TREXIO:**

- [x] Molecular geometry (atomic coordinates and charges)
- [x] Basis set definitions (Gaussian primitives and contractions)
- [x] Basis set pointers (atom-to-basis mapping)
- [x] Molecular orbital coefficients (LCAO)
- [x] Determinants or CSFs (wavefunction expansion)
- [x] Effective core potentials (ECPs/pseudopotentials) or All electron calculations
- [x] Orbital symmetries and eigenvalues

**Not included in TREXIO:**

- [ ] Jastrow parameters (must be provided separately)
- [ ] Jastrow derivatives (optional, for optimization)

## Advantages of TREXIO Format

- **Standardized**: Compatible across multiple QMC codes (CHAMP, QMCPack, TurboRVB)
- **Portable**: HDF5 backend is binary, cross-platform, and self-describing
- **Efficient**: Fast I/O with HDF5 library
- **Validated**: Built-in data integrity checks
- **Comprehensive**: Single file replaces many individual text files
- **Supported**: Export available from GAMESS, Quantum Package, PySCF, and other codes

## Basic CHAMP Input with TREXIO

A minimal CHAMP input file using TREXIO looks like:

```perl
%module general
    title           'VMC calculation with TREXIO'
    mode            'vmc'
%endmodule

load trexio          molecule.hdf5

load jastrow         jastrow.jas
load jastrow_der     jastrow.der

%module blocking_vmc
    vmc_nstep     20
    vmc_nblk      100
    vmc_nblkeq    1
    vmc_nconf_new 0
%endmodule
```

**Key points**:

- `load trexio molecule.hdf5` - Loads geometry, basis, orbitals, determinants from single file
- `load jastrow jastrow.jas` - Jastrow factors still need separate file
- `load jastrow_der jastrow.der` - Optional derivatives for optimization

## Creating TREXIO Files

### Method 1: Using trexio_tools

The `trexio_tools` package provides command-line utilities for converting quantum chemistry outputs to TREXIO format.

#### Installation

```bash
pip install trexio_tools
```

This installs the `trexio` command-line tool.

#### Verify Installation

```bash
trexio --version
```

#### Converting from GAMESS Output

```bash
trexio convert-from \
  --type gamess \
  --input gamess.out \
  --motype "RHF" \
  --back_end=HDF5 \
  molecule.hdf5
```

**Parameters**:

- `--type` - Input file format (gamess, pyscf, qp)
- `--input` - Path to quantum chemistry output file
- `--motype` - Molecular orbital type
- `--back_end` - Storage format (HDF5 or TEXT)
- Last argument - Output TREXIO filename

**Supported MO types**:

- `RHF` - Restricted Hartree-Fock
- `ROHF` - Restricted open-shell Hartree-Fock
- `MCSCF` - Multi-configurational SCF
- `NATURAL` - Natural orbitals
- `GUGA` - GUGA orbitals
- Check `trexio convert-from --help` for complete list


#### Converting from PySCF

```python
from pyscf import gto, scf
from trexio_tools import pyscf_to_trexio

# Perform calculation
mol = gto.M(atom='H 0 0 0; H 0 0 0.74', basis='cc-pvdz')
mf = scf.RHF(mol).run()

# Export to TREXIO
pyscf_to_trexio(mf, 'molecule.hdf5', backend='hdf5')
```

#### Help and Options

```bash
# General help
trexio --help

# Conversion help
trexio convert-from --help

# Inspection help
trexio inspect --help
```

### Method 2: Using Quantum Chemistry Code Plugins

Some quantum chemistry codes have built-in TREXIO export:

#### Quantum Package

```bash
# After running calculation
qp run fci
qp set trexio trexio_file molecule_fci.trexio
qp run export_trexio
```

## Inspecting TREXIO Files

### Check File Contents

```bash
# List available groups
h5ls molecule.hdf5

# Detailed HDF5 structure
h5dump -n molecule.hdf5
```

## Converting TREXIO to CHAMP Text Files

If you need individual text files instead of using TREXIO directly, CHAMP provides a converter script.

### Using trex2champ.py

Located in `champ/tools/trex_tools/`, this Python script extracts data from TREXIO into separate CHAMP-format text files.

#### Basic Usage

```bash
python3 /path/to/champ/tools/trex_tools/trex2champ.py \
  --trex "molecule.hdf5" \
  --backend "HDF5" \
  --basis_prefix "BFD-aug-cc-pVDZ" \
  --lcao \
  --ecp \
  --sym \
  --geom \
  --basis \
  --det
```

#### Parameters

| Flag | Description | Output File |
|------|-------------|-------------|
| `--trex` | Input TREXIO filename | - |
| `--backend` | Backend format (HDF5 or TEXT) | - |
| `--basis_prefix` | Prefix for basis file names | - |
| `--geom` | Extract geometry | `molecule.xyz` |
| `--basis` | Extract basis sets | `basis.bas` |
| `--lcao` | Extract MO coefficients | `orbitals.lcao` |
| `--det` | Extract determinants | `determinants.det` |
| `--ecp` | Extract pseudopotentials | `ecp.dat` |
| `--sym` | Extract orbital symmetries | `symmetry.sym` |

#### Example: Complete Extraction

```bash
cd /path/to/calculation

python3 ~/champ/tools/trex_tools/trex2champ.py \
  --trex "h2o_ccpvtz.hdf5" \
  --backend "HDF5" \
  --basis_prefix "cc-pVTZ" \
  --geom \
  --basis \
  --lcao \
  --det \
  --ecp \
  --sym
```

This creates:
- `molecule.xyz` - Geometry file
- `cc-pVTZ.bas` - Basis set file
- `orbitals.lcao` - Molecular orbital coefficients
- `determinants.det` - Wavefunction determinants
- `ecp.dat` - Effective core potentials
- `symmetry.sym` - Orbital symmetries

#### Help and Options

```bash
python3 trex2champ.py --help
```

## TREXIO File Formats

### HDF5 Backend (Recommended)

```bash
# Create HDF5 TREXIO
trexio convert-from --back_end=HDF5 molecule.hdf5
```

### TEXT Backend

```bash
# Create TEXT TREXIO
trexio convert-from --back_end=TEXT molecule.trexio
```

## Additional Resources

- [TREXIO Official Documentation](https://trex-coe.github.io/trexio/)
- [trexio_tools Repository](https://github.com/TREX-CoE/trexio_tools)
- [TREX-CoE Project](https://trex-coe.eu/)
- [Quantum Package](https://quantumpackage.github.io/qp2/)

## Next Steps

- Review [Jastrow Factors](jastrow.md) for creating correlation factors
- See [Calculation Guide](../calculations/index.md) for running VMC/DMC
- Try [Tutorials](../tutorials/index.md) for complete examples
- Explore [Wavefunction Optimization](../calculations/optimization/index.md)

## Getting Help

- Check [Troubleshooting Guide](../troubleshooting/index.md)
- Visit [TREXIO GitHub Issues](https://github.com/TREX-CoE/trexio/issues)
- Open an issue on [CHAMP GitHub](https://github.com/filippi-claudia/champ)