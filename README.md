# Cornell-Holland Ab-initio Materials Package (CHAMP)

<p align="center">
  <img src="docs/assets/logo.webp" alt="CHAMP logo" width=600/>
</p>


## CI/CD Status

| Category | Status |
|----------|--------|
| **Build - Release - Intel oneAPI** | [![CHAMP-EU CI - Release build - Intel oneAPI hpckit toolchain](https://github.com/filippi-claudia/champ/actions/workflows/build_champ_intel.yml/badge.svg)](https://github.com/filippi-claudia/champ/actions/workflows/build_champ_intel.yml) |
| **Build - Release - (VMC/DMC)** | [![CHAMP-EU CI - Release build - GNU without TREX (VMC/DMC Tests)](https://github.com/filippi-claudia/champ/actions/workflows/build_champ_vmc_dmc.yml/badge.svg)](https://github.com/filippi-claudia/champ/actions/workflows/build_champ_vmc_dmc.yml) |
| **Build - Release - TREX (TREXIO)** | [![CHAMP-EU CI - Release build - GNU with TREXIO (TREXIO Tests)](https://github.com/filippi-claudia/champ/actions/workflows/build_champ_trexio.yml/badge.svg)](https://github.com/filippi-claudia/champ/actions/workflows/build_champ_trexio.yml) |
| **Build - Release - TREX (QMCkl)** | [![CHAMP-EU CI - Release build - GNU with QMCkl (QMCKL Tests)](https://github.com/filippi-claudia/champ/actions/workflows/build_champ_qmckl.yml/badge.svg)](https://github.com/filippi-claudia/champ/actions/workflows/build_champ_qmckl.yml) |
| **Build - Debug - Intel + GNU** | [![CHAMP-EU CI - Debug build - Intel oneAPI and GNU toolchain](https://github.com/filippi-claudia/champ/actions/workflows/debug_champ_intel_and_gnu.yml/badge.svg)](https://github.com/filippi-claudia/champ/actions/workflows/debug_champ_intel_and_gnu.yml) |
| **Build - Docker Container** | [![CHAMP-EU CI - Docker image build and tagging for Docker Hub](https://github.com/filippi-claudia/champ/actions/workflows/docker-image.yml/badge.svg)](https://github.com/filippi-claudia/champ/actions/workflows/docker-image.yml) |
| **Publish - Docker Container** | [![CHAMP-EU CI - Docker image publish to Docker Hub Registry](https://github.com/filippi-claudia/champ/actions/workflows/docker-publish.yml/badge.svg)](https://github.com/filippi-claudia/champ/actions/workflows/docker-publish.yml) |
| **Tests - Unit Tests** | [![CHAMP-EU CI - Unit tests for CHAMP with GNU using FortUTF](https://github.com/filippi-claudia/champ/actions/workflows/build_champ_unit_test.yml/badge.svg)](https://github.com/filippi-claudia/champ/actions/workflows/build_champ_unit_test.yml) |
| **Tests - Python Interface** | [![CHAMP-EU CI - Python tests for the TREXIO tools integration](https://github.com/filippi-claudia/champ/actions/workflows/test_python.yml/badge.svg)](https://github.com/filippi-claudia/champ/actions/workflows/test_python.yml) |
| **Deployment - User Guide** | [![CHAMP-EU CI - Deploy Doxygen-generated content to Pages](https://github.com/filippi-claudia/champ/actions/workflows/build_documentation_devel.yml/badge.svg)](https://github.com/filippi-claudia/champ/actions/workflows/build_documentation_devel.yml) |

## Repository Statistics

| Category | Badges |
|----------|--------|
| **Activity** | ![Github Issues](https://img.shields.io/github/issues/filippi-claudia/champ) ![Github Pull Requests](https://img.shields.io/github/issues-pr/filippi-claudia/champ) ![Github Last Commit](https://img.shields.io/github/last-commit/filippi-claudia/champ) [![Commit Activity](https://img.shields.io/github/commit-activity/w/filippi-claudia/champ)](https://img.shields.io/github/commit-activity/t/filippi-claudia/champ) |
| **Release** | ![Last release tag](https://img.shields.io/github/v/tag/filippi-claudia/champ) |
| **Community** | ![Github forks](https://img.shields.io/github/forks/filippi-claudia/champ) ![Github stars](https://img.shields.io/github/stars/filippi-claudia/champ) |
| **Size** | ![Repo Size](https://img.shields.io/github/repo-size/filippi-claudia/champ) ![Code Size](https://img.shields.io/github/languages/code-size/filippi-claudia/champ) |
| **License** | ![Github license](https://img.shields.io/github/license/filippi-claudia/champ) |

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

We welcome contributions! Please see [CONTRIBUTING.md](CONTRIBUTING.md) for details.

## License

Distributed under the GPL-3.0 License. See [LICENSE](LICENSE) for more information.

## Support

For support, please contact the PI [Claudia Filippi](mailto:c.filippi@utwente.nl) or the developer [Ravindra Shinde](mailto:r.l.shinde@utwente.nl).

## Discussion

For discussions or FAQs, please use the [GitHub Discussions](https://github.com/filippi-claudia/champ/discussions) page.

## Disclaimer

The authors make no claims about the correctness of the program suite and it is provided without warranty under GPL-3.0. 


