# H₂ Dimer Test

(Reference: `tests/CI_test/VMC-H2`)

This is the simplest test case, suitable for a quick sanity check.

## System

-   Molecule: H₂
-   Bond length: 1.4 a.u.
-   Basis: cc-pVDZ (or similar simple basis)

## Workflow

1.  **VMC**: Run a short VMC calculation.
2.  **Check Energy**: The total energy should be approximately -1.17 a.u. (depending on the quality of the trial wavefunction).

## Input Example

```python
load molecule
   2
   H  1.0  0.0  0.0  0.0
   H  1.0  0.0  0.0  1.4
end

%module blocking_vmc
    vmc_nstep     10
    vmc_nblk      100
%endmodule
```
