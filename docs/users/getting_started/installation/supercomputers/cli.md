---
title: CHAMP CLI
icon: material/console
tags:
    - usage
---

This document explains how to use the CHAMP program from the command line. It is based directly on the program’s source-code modules responsible for execution control, calculation modes (VMC/DMC), periodic settings, file handling, and command-line parsing.

## Overview

CHAMP is controlled at runtime using both:

- Command-line flags — specify input/output/error files, verbosity, debug mode, etc.

- Input files (.in, .inp, .dat) — define calculation parameters for:

    - VMC (vmc.inp)

    - DMC (dmc.inp)

The program supports both serial and MPI-aware execution.


## Basic Command Structure

The basic command line structure of CHAMP looks like:

```bash
/path/to/champ/bin/vmc.mov1 -i input.inp -o output.log -e error
```

Only the master MPI rank (rank 0) writes output; worker ranks redirect output to /dev/null for clean parallel runs.

A successfully built CHAMP provides two executables in the `bin` directory.

- `vmc.mov1` for VMC one-electron move
- `dmc.mov1` for DMC one-electron move

## Command-Line Options

CHAMP recognizes several command-line flags, handled inside the contrl_file module.

### Required / Most Common Flags

| Flag(s) | Description |
|---------|-------------|
| `-i`, `-in`, `-inp`, `-input`, `--input` | Specify input file (*.in, *.inp, *.dat) |
| `-o`, `-out`, `--output` | Output file (master rank only; workers redirect to /dev/null) |
| `-e`, `-err`, `--error` | Error and warnings log file |

### Optional Flags

| Flag(s) | Description |
|---------|-------------|
| `-v`, `--version` | Print version information and exit |
| `-h`, `--help` | Print help and exit |
| `-V`, `--verbose` | Enable verbose output |
| `-d`, `--debug` | Enable debug mode |
| `-p`, `--prefix` | Add prefix to generated files |

### Example: Get help

```bash
champ/bin/vmc.mov1 -h
```

## Execution Modes

Execution mode is controlled by the control module. You can set mode inside input files or via initialization routines:

- `vmc` — Variational Monte Carlo
- `dmc` — Diffusion Monte Carlo

Verbosity (`ipr`) controls output level:

| ipr value | Meaning |
|-----------|----------|
| -1 | minimal output |
| 0 | normal (default) |
| 1 | increasingly verbose, debug-level output |
| &gt;2 | All possible print output |


## Typical Usage Examples



### Run VMC

```bash
champ/bin/vmc.mov1 -i vmc.inp -o vmc.out -e error
```

### Run DMC

```bash
champ/bin/dmc.mov1 -i dmc.inp -o dmc.log -e error
```


### Using MPI

For production calculations, you should use MPI. 

```bash
mpirun -np 1024 champ/bin/vmc.mov1 -i dmc.inp -o dmc.out
```

