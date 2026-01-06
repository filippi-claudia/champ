# Variance Minimization

Quantum Monte Carlo methods are some of the most accurate and efficient methods for treating many-body systems. The success of these methods differs largely on the flexibility of the trial wave functions and the capability to efficiently optimize their parameters for both Variational Monte Carlo (VMC) and Diffusion Monte Carlo (DMC).

The **variance minimization** method has become a frequently used method for optimizing many-body wave functions because it is often more efficient than straightforward energy minimization.

Mathematically, it aims to minimize the variance of the local energy:

$$ \sigma^2 = \frac{\int |\Psi_T|^2 (E_L - E_T)^2 d\mathbf{R}}{\int |\Psi_T|^2 d\mathbf{R}} $$

## Why use Variance Minimization?

For a sufficiently flexible variational wave function, straightforward energy minimization on a finite set of Monte Carlo (MC) configurations can be problematic. It is possible to lower the energy on the finite sample on which the optimization is performed, while in fact raising the true expectation value of the energy.

On the other hand, if the **variance** of the local energy is minimized, each term in the sum over MC configurations is bounded from below by zero. This makes the optimization problem far less severe and numerically more stable.

## Implementation in CHAMP

