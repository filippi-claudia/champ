# VMC Optimization (Solid State)

Optimization in solids is similar to molecules but requires attention to k-point sampling and finite-size effects.

## Optimization Strategy

1.  **Jastrow Optimization**: Optimize the Jastrow factor first. The Jastrow factor is often short-ranged and transferable.
2.  **Orbital Optimization**: Optimizing orbitals in solids can be expensive and is sometimes skipped if the DFT orbitals are good enough.

## Example

```python
%module optwf
    ioptwf      1
    ioptjas     1
    ioptorb     0
    method      'linear'
%endmodule
```
