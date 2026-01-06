# Basis & MO Generation

CHAMP requires a trial wavefunction consisting of a determinantal part (Orbitals) and a Jastrow factor.

## Recommended Workflow: TREXIO

The most robust way to generate input files is using the **TREXIO** format and the `trex2champ` converter.

1.  **Generate TREXIO File**:
    -   Run a calculation in a supported code (Quantum Package, GAMESS, PySCF).
    -   Export the result to a TREXIO file (e.g., `molecule.hdf5`).
    -   See [Using TREXIO Files](../preparation/trexio.md) for details.

2.  **Convert to CHAMP Format**:
    -   Use the `trex2champ` tool to extract the necessary files.
    ```bash
    trex2champ molecule.hdf5
    ```
    -   This will generate:
        -   `molecule.xyz` (Geometry)
        -   `<basis_name>.basis.<element>` (Basis sets for each element)
        -   `orbitals.lcao` (Molecular orbitals)
        -   `determinants.det` (Determinants)
        -   `ecp.dat` (Pseudopotentials, if applicable)

3.  **Organize Files**:
    -   Move the generated files to a `pool/` directory.
    -   CHAMP looks for basis files in the `pool/` directory by default.

## Manual File Creation

If you are not using TREXIO, you must ensure your files follow the specific CHAMP formats.

### Basis Sets
-   **Format**: Radial grid representation (not Gaussian exponents directly).
-   **Naming**: `<basis_name>.basis.<element>` (e.g., `BFD.basis.C`).
-   **Location**: Must be in the `pool/` directory.
-   See [Basis Sets](../preparation/basis.md) for the file format specification.

### Molecular Orbitals
-   **Format**: `.lcao` or `.orb` file containing LCAO coefficients.
-   **Ordering**: Atomic orbitals must follow the **TREXIO alphabetical convention** (e.g., X, Y, Z for p-orbitals).
-   **Spin**: For open-shell systems, you need separate files for alpha and beta orbitals.
-   See [Molecular Orbitals](../preparation/orbitals.md) for details.

## Input Configuration

In your CHAMP input file:

```python
%module general
    pool   './pool/'
    basis  'BFD'      # Prefix for basis files
%endmodule

load orbitals  $pool/orbitals.lcao
```
