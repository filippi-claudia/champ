# Open-shell Systems

(Reference: `tests/CI_test/VMC-HNO-cipsi_322_dets_143_csfs_noopt`)

Open-shell systems (molecules with unpaired electrons) require special attention to spin states.

## Configuration

1.  **Spin Polarization**:
    -   Specify the number of up (`nup`) and total (`nelec`) electrons correctly in the input.
    -   $N_{dn} = N_{elec} - N_{up}$.

2.  **Multi-Determinant Wavefunctions**:
    -   Open-shell systems often require multi-determinant wavefunctions (e.g., from CIPSI or CASSCF) to capture static correlation.
    -   Ensure the expansion includes the necessary determinants for the target spin state.

## Example (HNO)

HNO radical has 7 up and 6 down electrons (with pseudo-potentials).

```python
load molecule
    ...
end

%module electrons
    nup         7
    nelec       13
%endmodule
```
