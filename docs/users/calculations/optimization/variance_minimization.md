# Variance Minimization

Variance minimization aims to minimize the variance of the local energy:

$$ \sigma^2 = \frac{\int |\Psi_T|^2 (E_L - E_T)^2 d\mathbf{R}}{\int |\Psi_T|^2 d\mathbf{R}} $$

## Overview

Historically, variance minimization was widely used because it is more stable than energy minimization for simple algorithms. However, energy minimization (via SR or Linear Method) is now often preferred as it directly optimizes the quantity of interest.

Variance minimization is still useful for:
- Generating initial parameters.
- Optimizing Jastrow factors when energy optimization is difficult.

## Implementation

In CHAMP, variance minimization can be achieved by specific settings in the optimization module, often by using the Linear Method with a target function that emphasizes variance.

(Note: Specific keywords for pure variance minimization in CHAMP should be checked in the latest source code, as energy minimization is the default).
