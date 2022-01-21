```
 ____________________________________________________________________


  .d8888b.   888    888         d8888  888b     d888  8888888b.
 d88P  Y88b  888    888        d88888  8888b   d8888  888   Y88b
 888    888  888    888       d88P888  88888b.d88888  888    888
 888         8888888888      d88P 888  888Y88888P888  888   d88P
 888         888    888     d88P  888  888 Y888P 888  8888888P"
 888    888  888    888    d88P   888  888  Y8P  888  888
 Y88b  d88P  888    888   d8888888888  888   "   888  888
  "Y8888P"   888    888  d88P     888  888       888  888

 ____________________________________________________________________
```
------
![Logo](https://github.com/filippi-claudia/champ/blob/main//docs/logo_small.png?raw=true | width=200)


[![Self-hosted Intel OneAPI Build and Testing on ccpgate/ccp01](https://github.com/filippi-claudia/champ/actions/workflows/self_hosted_build_champ_intel_fdfparser.yml/badge.svg)](https://github.com/filippi-claudia/champ/actions/workflows/self_hosted_build_champ_intel_fdfparser.yml)

[![CI using Doxygen generated doc](https://github.com/filippi-claudia/champ/actions/workflows/CI.yml/badge.svg)](https://github.com/filippi-claudia/champ/actions/workflows/CI.yml)


The Cornell-Holland Ab-initio Materials Package (CHAMP) is a quantum Monte Carlo
suite of programs for electronic structure calculations of atomic and molecular systems.
The code is a sister code of the homonymous program originally developed by Cyrus Umrigar
and Claudia Filippi of which it retains the accelerated Metropolis method and the efficient
diffusion Monte Carlo algorithms.

The European branch of the code is currently developed by Claudia Filippi and Saverio Moroni,
with significant contributions by Claudio Amovilli and other collaborators.

CHAMP has three basic capabilities:

* Metropolis or variational Monte Carlo (VMC)
* Diffusion Monte Carlo (DMC)
* Optimization of many-body wave functions by energy minimization (VMC) for ground and excited states

Noteworthy features of CHAMP are:

* Efficient wave function optimization also in a state-average fashion for multiple states of the same symmetry (VMC)
* Efficient computation of analytical interatomic forces (VMC)
* Compact formulation for a fast evaluation of multi-determinant expansions and their derivatives (VMC and DMC)
* Multiscale VMC and DMC calculations in classical point charges (MM), polarizable continuum model (PCM), and polarizable force fields (MMpol)

**NOTE**

You should neither obtain this program from any other source nor should you distribute it
or any portion thereof to any person, including people in the same research group.

It is expected that users of the programs will do so in collaboration
with one of the principal authors.  This serves to ensure both that the
programs are used correctly and that the principal authors get adequate
scientific credit for the time invested in developing the programs.

**Usual disclaimer**

The authors make no claims about the correctness of
the program suite and people who use it do so at their own risk.

------------------------------------------------------------------------

CHAMP relies on various other program packages:

1. [Parser](https://github.com/neelravi/mpi-libfdf-parser):
   An easy-to-use and easy-to-extend keyword-value pair based input file parser written in Fortran 2008.  This parser uses a heavily modified libFDF library and is written by Ravindra Shinde. The parser is built in the code, however, the parser folder can be easily adapted by any other Fortran-based code.

2. GAMESS:
   For finite systems the starting wavefunction is obtained from the
   quantum chemistry program GAMESS, written by Mike Schmidt and
   collaborators at Iowa State University.

3. GAMESS_Interface:
   The wavefunction produced by GAMESS has to be cast in a form
   suitable for input to CHAMP.  This is a lot more work than first meets
   the eye. The Perl script was written by Friedemann Schautz.

4. MOLCAS_Interface: recently added thanks to Csaba Daday and Monika Dash

------------------------------------------------------------------------

### Requirements
1. cmake >= 3.21
2. gfortran/gcc >= 9.3.0 or Intel Fortran 2020 onwards
3. BLAS/LAPACK or Intel MKL
4. openMPI >= 3.0 or Intel MPI
5. [Optional] TREXIO library >= 2.0.0
6. [Optional] doxygen (for documentation)


### Installation Using CMake
To install **Champ** using [cmake](https://cmake.org/) you need to run the following commands:
```
cmake -H. -Bbuild
cmake --build build -- -j4
```
The first command is only required to set up the build directory and needs to be
executed only once. Compared to the previous Makefiles the dependencies for the
include files (e.g include/vmc.h) are correctly setup and no `--clean-first` is
required.

#### CMAKE Options

To select a given compiler, you can type:
```
cmake -H. -Bbuild -D CMAKE_Fortran_COMPILER=ifort
```
To use LAPACK and BLAS installed locally, include the path to the libraries:
```
cmake -H. -Bbuild -D CMAKE_Fortran_COMPILER=ifort -D BLAS_blas_LIBRARY=/home/user/lib/BLAS/blas_LINUX.a -D LAPACK_lapack_LIBRARY=/home/user/lib/LAPACK/liblapack.a
```
To compile only e.g. VMC serial:
```
cmake --build build --target vmc.mov1
```
Clean and build:
```
cmake --build build --clean-first
```
Compared to the previous Makefiles the dependencies for the include files
(e.g include/vmc.h) are correctly setup and no `--clean-first` is required.

#### CMAKE Recipes

Here are a couple of recipes for commonly used computing facilities, which can be easily adapted.
* **Cartesius** (cartesius.surfasara.nl):
	- To compile the code, first load the required modules:
		```bash
		module load 2021
		module load git
        module load CMake/3.20.1-GCCcore-10.3.0
        module load intel-compilers/2021.2.0
        module load imkl/2021.2.0-iimpi-2021a
        module load impi/2021.2.0-intel-compilers-2021.2.0
		```
		then set-up the build:
		```bash
		cmake -H. -Bbuild -DCMAKE_Fortran_COMPILER=mpiifort
		```
		and finally build:
		```bash
		cmake --build build -j8 --clean-first
		```
	- To run the code, you need to submit a job to the queue system:
		```bash
		sbatch job.cmd
		```
		where `job.cmd` is a SLURM script that looks like this:

		```bash
		#!/bin/bash
        #SBATCH -t 0-12:00:00           # time in (day-hours:min:sec)
        #SBATCH -N 60                   # number of nodes
        #SBATCH -n 1440                 # number of cores
        #SBATCH --ntasks-per-node 24    # tasks per node
        #SBATCH -J cn19-B1              # name of the job
        #SBATCH -o vmc.%j.out           # std output file name for slurm
        #SBATCH -e vmc.%j.err           # std error file name for slurm
        #SBATCH --constraint=haswell    # specific requirements about procs
		#SBATCH -p normal               # partition (queue)
        #
        module purge
		module load 2021
		module load imkl
        #
        export I_MPI_PMI_LIBRARY=/usr/lib64/libpmi2.so
        cd $PWD
		srun champ/bin/vmc.mov1 -i input.inp -o output.out -e error
		```
* **CCPGate**:
	- To build with ifort set the variables for the Intel Compiler and MPI:
		If you use CSH:
		```
		source /software/intel/intel_2019.0.117/compilers_and_libraries_2019.1.144/linux/bin/compilervars.csh -arch intel64 -platform linux
		source /software/intel/intel_2019.0.117/compilers_and_libraries_2019.0.117/linux/mpi/intel64/bin/mpivars.csh -arch intel64 -platform linux
		```
		If you use BASH:
		```
		. /software/intel/intel_2019.0.117/compilers_and_libraries_2019.1.144/linux/bin/compilervars.sh intel64
		. /software/intel/intel_2019.0.117/compilers_and_libraries_2019.0.117/linux/mpi/intel64/bin/mpivars.sh intel64
		```
		Setup the build:
		```
		cmake -H. -Bbuild -DCMAKE_Fortran_COMPILER=mpiifort
		```
 	- To build with gfortran:
		If you use CSH:
		```
		source /software/intel/intel_2019.0.117/impi/2019.0.117/intel64/bin/mpivars.sh -arch intel64 -platform linux
		```
		If you use BASH:
		```
		. /software/intel/intel_2019.0.117/impi/2019.0.117/intel64/bin/mpivars.sh intel64
		```
		Setup the build:
		```
		cmake -H. -Bbuild -DCMAKE_Fortran_COMPILER=mpif90
		```
		which will use LAPACK & BLAS from the Ubuntu repository. (Cmake should find them already if none of the Intel MKL variables are set.) Combining gfortran with the Intel MKL is possible but requires special care to work with the compiler flag `-mcmodel=large`.
	- To run the code:
		```
		mpirun -s all -np "n process" -machinefile "machinefile"
		```
* **Ubuntu desktop**:
	- Ubuntu 18.04:
		Install the required packages:
		```
		sudo apt install gfortran openmpi-bin gawk libblacs-mpi-dev liblapack-dev
		```
		Set-up the build:
		```
		cmake -H. -Bbuild -DCMAKE_Fortran_COMPILER=mpifort
		```
		Build:
		```
		cmake --build build -- -j2
		```
		To run in parallel:
		```
		mpirun --stdin all -n 2 path_to_CHAMP/bin/vmc.mov1 < vmc.inp > vmc.out
		```
	- Ubuntu 20.04:
	We are still working on the CHAMP built on the latest Unbuntu release using a OpenMPI v4.X version. The code compiles but fails to run the test in parallel. For the time being, we urge the user to use an older version of Ubuntu, as shown above.

------------------------------------------------------------------------

### Documentation
CHAMP developer documentation can be generated using [Doxygen](http://www.doxygen.nl/) tool. To install the package, we advise to follow the instructions at the Doxygen web page: <http://www.doxygen.nl/download.html>.

The Doxyfile file provided in CHAMP docs directory contains all the settings needed to generate the documentation. Once Doxygen is installed, at the docs folder of CHAMP simply run:
```
doxygen Doxyfile
```
Then two folders will be created, /docs/developers/html and ./docs/developers/latex, containing the documentation about modules, subroutines and functions in both html and latex formats.
