# Geometry Conversion

CHAMP uses a specific format for defining molecular geometry in the input file.

!!! warning "Important: Units are in Bohr"
    All atomic coordinates in CHAMP must be specified in **Bohr** (atomic units).
    1 Ångström = 1.8897259886 Bohr.

## Input Formats

The geometry can be defined in the `load molecule` section, a separate file, or using an inline block.

### 1. Standard Format (Automatic Valence)

For all-electron calculations, CHAMP automatically assigns valence charges.

```python
load molecule
    3
    H2O molecule
    O  0.000000  0.000000  0.000000
    H  0.000000  1.430428  1.107230
    H  0.000000 -1.430428  1.107230
end
```

### 2. Explicit Valence (For ECPs)

When using Effective Core Potentials (ECPs), you **must** specify the number of valence electrons for each atom. This is the 5th column in the geometry file.

```python
load molecule
    3
    H2O with ECP on Oxygen (6 valence e-)
    O  0.000000  0.000000  0.000000  6.0
    H  0.000000  1.430428  1.107230  1.0
    H  0.000000 -1.430428  1.107230  1.0
end
```

### 3. Inline Block

For small molecules, you can use the `%block molecule` syntax directly in the input file.

```python
%block molecule
    3
    H2O inline
    O  0.000000  0.000000  0.000000
    H  0.000000  1.430428  1.107230
    H  0.000000 -1.430428  1.107230
%endblock
```

## Using TREXIO

If you are using TREXIO, the geometry is automatically read from the TREXIO file.

```python
load trexio  'my_molecule.trexio'
```

## Detailed Reference

For a complete description of geometry formats and options, see the [Molecular Geometry](../preparation/molecule.md) page.
