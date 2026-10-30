# VMC Output Analysis

## Standard Output

The standard output of a VMC run contains information about the progress of the simulation and the final results.

### Energy and Variance

Look for lines containing `total E` or `Energy` to find the estimated energy.

```text
Block    1  E = -15.8765 +/- 0.0012
...
Final Energy = -15.8760 +/- 0.0005
```

### Blocking Analysis

CHAMP performs blocking analysis to estimate the statistical error. The output will show the error estimate as a function of block size or number of blocks.

### Correlation Time

The output may also report the correlation time, which indicates how many steps are needed to generate an independent sample.

## Output Files

- **`mc_configs`**: Contains the current walker configurations.
- **`mc_configs_newX`**: Generated if `vmc_nconf_new > 0`. Contains walkers saved for subsequent runs (e.g., DMC).
- **`vmc.out`** (or similar): The main log file.

## Visualization

You can use the `analysis` tools (see [Analysis](../analysis/index.md)) to plot the energy trace and check for equilibration.
