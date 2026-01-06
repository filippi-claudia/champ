# Linear Method

The Linear Method is a powerful optimization algorithm that expands the wavefunction to first order in the parameter changes and solves the resulting generalized eigenvalue problem.

## Overview

The Linear Method can be more stable and efficient than simple gradient descent, especially for optimizing many parameters simultaneously (e.g., Jastrow + Orbitals).

It is selected with `method = 'linear'` or `method = 'lin_d'`.

## Key Features

- **Robustness**: Can handle large parameter sets.
- **Efficiency**: Often requires fewer iterations than simple gradient methods.
- **Cost Function**: Can target energy, variance, or a mix.

## Usage Example

```python
%module optwf
    ioptwf      1
    method      'linear'
    nopt_iter   10
%endmodule
```
