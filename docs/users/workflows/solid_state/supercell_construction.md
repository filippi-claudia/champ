# Supercell Construction

QMC calculations for solids are typically performed on supercells to reduce finite-size effects.

## Overview

A supercell is a multiple of the primitive unit cell. For example, a $2 \times 2 \times 2$ supercell contains 8 times the number of atoms of the primitive cell.

## Input Configuration

The lattice vectors and atomic positions are defined in the input.

```python
load molecule
    ... (atomic positions in the supercell) ...
end

%module cell
    a1  10.0  0.0  0.0
    a2   0.0 10.0  0.0
    a3   0.0  0.0 10.0
%endmodule
```

## Tools

External tools (like `ase` or `pymatgen`) are often used to generate the supercell geometry and convert it to CHAMP format.
