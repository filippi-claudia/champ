---
title: Home
---

# Cornell-Holland Ab-initio Materials Package

<figure markdown="span">
  ![CHAMP logo](assets/logo.webp){ width="300" }
</figure>

## Overview

The Cornell-Holland Ab-initio Materials Package (CHAMP-EU) is a quantum Monte Carlo suite of programs for electronic structure calculations of atomic and molecular systems. The code is a sister code of the homonymous program originally developed by Cyrus Umrigar and Claudia Filippi of which it retains the accelerated Metropolis method and the efficient diffusion Monte Carlo algorithms.

The European branch of the code is currently developed by Claudia Filippi and Saverio Moroni,
with significant contributions by Ravindra Shinde, Emiel Slootman, Nicolas Renaud, Victor Azizi, Edgar Landinez, and Stuart Shepard.

## Key Features

CHAMP has three basic capabilities:

* Metropolis or variational Monte Carlo (VMC)
* Diffusion Monte Carlo (DMC)
* Optimization of many-body wave functions by energy minimization (VMC) for ground and excited states

Noteworthy features of CHAMP are:

* Efficient wave function optimization also in a state-average and a state-specific fashion for multiple states of the same symmetry (VMC)
* Efficient computation of analytical interatomic forces (VMC)
* Compact formulation for a fast evaluation of multi-determinant expansions and their derivatives (VMC and DMC)
* Multiscale VMC and DMC calculations in classical point charges (MM), polarizable continuum model (PCM), and polarizable force fields (MMpol)

## Documentation

For full documentation, including detailed installation guides and tutorials, visit the [Online User's Documentation](https://filippi-claudia.github.io/champ/). 

## Citation

If you use CHAMP-EU in your research, please cite this source:

Shinde, R., Landinez Borda, E. J., Shepard, S., Slootman, E., Cuzzocrea, A., Azizi, V., Lopez-Tarifa, P., Renaud, N., Umrigar, C., Moroni, S., & Filippi, C. (2024). Cornell-Holland Ab-initio Materials Package (CHAMP-EU) v2.3.0. Zenodo. [https://doi.org/10.5281/zenodo.11369538](https://doi.org/10.5281/zenodo.11369538)


## Contributing

We welcome contributions! Please see [Contributing Guidelines](https://github.com/filippi-claudia/champ/blob/main/Contributing.md) for details.

## License

Distributed under the GPL-3.0 License. See [LICENSE](https://github.com/filippi-claudia/champ/blob/main/LICENSE) for more information.

## Support

For support, please contact the PI [Claudia Filippi](mailto:c.filippi@utwente.nl) or the developer [Ravindra Shinde](mailto:r.l.shinde@utwente.nl).

## Discussion

For discussions or FAQs, please use the [GitHub Discussions](https://github.com/filippi-claudia/champ/discussions) page.

## Disclaimer

The authors make no claims about the correctness of the program suite and it is provided without warranty under GPL-3.0. 




