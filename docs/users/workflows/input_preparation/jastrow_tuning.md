# Jastrow Tuning

Before starting a full wavefunction optimization, it is often necessary to "tune" the initial Jastrow factor.

## Workflow

1.  **Define Jastrow Form**: Choose the expansion orders (`norda`, `nordb`, `nordc`).
2.  **Initial Parameters**: Set initial parameters.
3.  **Short VMC**: Run a short VMC calculation to check the energy and variance.

## Jastrow File Format

The Jastrow file (`jastrow.jas`) must follow a specific structure.

### Example: 3-Body Jastrow

To include a 3-body (e-e-n) Jastrow term (`nordc > 0`):

```python
jastrow_parameter   1
  5  5  5           norda,nordb,nordc
   0.60000000         scalek
   0.0 0.0 -0.4 -0.2 -0.04 0.08  (a parameters for atom type 1)
   0.0 0.0 -0.1 -0.01 0.01 0.01  (a parameters for atom type 2)
   0.5 0.37 0.07 0.01 -0.01 -0.01 (b parameters)
   (c parameters for atom type 1)
   (c parameters for atom type 2)
end
```

**Key Parameters**:
-   `norda`: Order of electron-nucleus terms (one set per atom type).
-   `nordb`: Order of electron-electron terms (one set for all).
-   `nordc`: Order of electron-electron-nucleus terms (one set per atom type).

See [Jastrow Factors](../preparation/jastrow.md) for the full specification of the file format and parameter counts.

## Tuning Strategy

1.  **Start Simple**: Use `nordc=0` (no 3-body) initially.
2.  **Optimize Jastrow**: Use `ioptjas=1` in VMC to optimize the parameters.
3.  **Add 3-Body**: Once 1- and 2-body terms are good, increase `nordc` and re-optimize.
