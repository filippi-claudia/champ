---
title: Backflow Transformations
tags:
    - backflow
    - correlation
    - wavefunction
    - optimization
---

# Backflow Transformations

Backflow replaces the electron coordinates used in the Slater determinants by
collective **quasiparticle coordinates**. Unlike a Jastrow factor, which
multiplies the wavefunction, backflow changes the nodal surface of the
determinantal part of the trial wavefunction.

The current backflow implementation in CHAMP is inspired by Rios *et al.* (2006). It
can be used in VMC and DMC, and its parameters can be included in wavefunction
optimization. It is currently not compatible with force calculations, multiple states, and multiple Slater determinants.


## Transformation

For electron $i$, CHAMP evaluates the orbitals at the quasiparticle coordinate
$\mathbf q_i$:

$$
\begin{aligned}
\mathbf q_i = \mathbf r_i
&+ \sum_{j\ne i} \eta_{ij}(r_{ij})\mathbf r_{ij}
+ \sum_I \chi_{iI}(r_{iI})\mathbf r_{iI} \\
&+ \sum_{j\ne i}\sum_I
   \left[
      \phi_{ijI}\mathbf r_{ij}
      + \theta_{ijI}\mathbf r_{iI}
   \right],
\end{aligned}
$$

where

$$
\mathbf r_{ij}=\mathbf r_i-\mathbf r_j,
\qquad
\mathbf r_{iI}=\mathbf r_i-\mathbf R_I.
$$

The transformation contains three optional correlation terms:

1. **Electron-electron (e-e)**: $\eta_{ij}$, with either one shared spin
   channel or separate parallel- and antiparallel-spin channels.
2. **Electron-nucleus (e-n)**: $\chi_{iI}$, with a separate parameter set for
   each nuclear type.
3. **Electron-electron-nucleus (e-e-n)**: $\phi_{ijI}$ and $\theta_{ijI}$,
   with separate nuclear-type and, optionally, spin-channel parameters.

All terms use the compact cutoff function

$$
g(r;L,C)=
\begin{cases}
\left(1-r/L\right)^C, & r\le L,\\
0, & r>L,
\end{cases}
$$

where $L$ is the cutoff length and $C$ is the integer cutoff exponent stored in
the backflow file.

The two-body functions are polynomials multiplied by this cutoff:

$$
\eta_{ij}(r)=g(r;L_b,C)\sum_{p=0}^{N_b}b_p^{\sigma_{ij}}r^p,
$$

$$
\chi_{iI}(r)=g(r;L_{a,t(I)},C)\sum_{p=0}^{N_a}a_p^{t(I)}r^p.
$$

Here $t(I)$ is the nuclear type and $\sigma_{ij}$ is the e-e spin channel.
For the three-body terms, CHAMP uses

$$
\phi_{ijI}=
g(r_{iI};L_{c,t(I)},C)g(r_{jI};L_{c,t(I)},C)
\sum_{p+q+s\le N_c}c^\phi_{pqs,t(I),\sigma_{ij}}
r_{iI}^p r_{jI}^q r_{ij}^s,
$$

with the same expansion for $\theta_{ijI}$ using coefficients $c^\theta$.

## Enabling Backflow

Backflow parameters are read from a separate text file. Add a top-level load
command to the main CHAMP input:

```perl
load backflow  backflow.bf
```

Loading a backflow file enables backflow. To disable backflow, omit the
`load backflow` command.

Backflow is not stored in the TREXIO input read by CHAMP, so the external
backflow file is also required for TREXIO-based calculations.

## File Format

The recommended complete layout is:

```text
backflow  1
N_a  N_b  N_c  nparm  nspin
C
L_a(1) ... L_a(ntype)  cutoff-a
L_b                       cutoff-b
L_c(1) ... L_c(ntype)  cutoff-c

a coefficients for nuclear type 1
...
a coefficients for nuclear type ntype
b coefficients for parallel/shared spin
b coefficients for antiparallel spin       # only when nspin = 2

phi coefficients for parallel/shared spin and nuclear type 1
...
phi coefficients for all remaining spin/type combinations
theta coefficients in the same spin/type order
end
```

The cutoff lines are present only for enabled components. The literal labels
`cutoff-a`, `cutoff-b`, and `cutoff-c` are required on their corresponding
lines.

### Header and orders

The first line must be

```text
backflow  1
```

where `1` selects Rios *et al.* backflow.

The second line contains five integers:

```text
N_a  N_b  N_c  nparm  nspin
```

- `N_a` is the e-n polynomial order. Each nuclear type has `N_a + 1`
  coefficients, ordered as $a_0,a_1,\ldots,a_{N_a}$.
- `N_b` is the e-e polynomial order. Each spin channel has `N_b + 1`
  coefficients, ordered as $b_0,b_1,\ldots,b_{N_b}$.
- `N_c` is the e-e-n total polynomial order. Set it to `0` to omit the
  three-body term.
- `nparm` is the total number of stored cutoff and polynomial parameters. CHAMP
  derives its internal array sizes from the orders, nuclear-type count, and
  spin count; this values does not need to be updated manually. It is only for informational purposes in generated files.
- `nspin` is `1` for one shared e-e channel or `2` for separate parallel and
  antiparallel channels. Values other than `2` are treated as `1`.

For the current text reader, use positive `N_a` and `N_b`. A zero `N_c` is
supported and disables the three-body section.

The third line is the cutoff exponent `C`. Currently, only the value `3` is supported.

### Cutoff lengths

Cutoff lengths use CHAMP's internal length unit (bohr) and must be positive.
Their records are ordered as follows:

1. `ntype` e-n cutoffs followed by `cutoff-a`;
2. one e-e cutoff followed by `cutoff-b`;
3. `ntype` e-e-n cutoffs followed by `cutoff-c`, when `N_c > 0`.

`ntype` is the number of nuclear types in the molecule. Values and subsequent
coefficient lines follow the same nuclear-type order as the molecular input. When `nspin` is equal to `2`, both spin channels use the same cutoff factor.


## Complete Template

The following file is a valid one-nuclear-type template with first-order e-n,
e-e, and e-e-n expansions and separate spin channels. Its 25 parameters are
all zero except for the three cutoff lengths, so it initially gives the
identity transformation $\mathbf q_i=\mathbf r_i$.

```text
backflow    1
    1    1    1   25    2   norda,nordb,nordc,nparm,nspin_bf_ee
  3                             C
  3.00000000                    cutoff-a
  3.00000000                    cutoff-b
  3.00000000                    cutoff-c
  0.00000000  0.00000000       e-n, type 1
  0.00000000  0.00000000       e-e, parallel
  0.00000000  0.00000000       e-e, antiparallel
  0.00000000  0.00000000  0.00000000  0.00000000   phi, parallel, type 1
  0.00000000  0.00000000  0.00000000  0.00000000   phi, antiparallel, type 1
  0.00000000  0.00000000  0.00000000  0.00000000   theta, parallel, type 1
  0.00000000  0.00000000  0.00000000  0.00000000   theta, antiparallel, type 1
end
```

For a system with multiple nuclear types:

- put all `cutoff-a` values on the same line and all `cutoff-c` values on the
  same line;
- add one e-n line per nuclear type;
- add one $\phi$ and one $\theta$ line per nuclear type in every spin channel;


## Optimizing Backflow

Set `ioptbf = 1` in the `optwf` module to include backflow parameters in a
wavefunction optimization:

```perl
%module optwf
    ioptwf          1
    ioptbf          1
    ibf_opt_cutoff  true
    method          'sr_n'
    nopt_iter       20
%endmodule
```

`ibf_opt_cutoff` controls whether the cutoff lengths are optimized:

- `true` (default): optimize cutoffs and polynomial coefficients;
- `false`: keep all cutoffs fixed and optimize only the coefficients.

Backflow can be optimized together with Jastrow, orbital, or CI parameters by
also enabling `ioptjas`, `ioptorb`, or `ioptci`.

