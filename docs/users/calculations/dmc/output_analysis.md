# DMC Output Analysis

## Energy Trace

The DMC output reports the energy at each block.

```text
Block    1  E = -15.9200 +/- 0.0020
...
Final Energy = -15.9215 +/- 0.0008
```

## Population Control

DMC involves a fluctuating population of walkers. The code adjusts the trial energy $E_T$ slightly to keep the population stable around the target number (`dmc_nconf` * number of MPI tasks).

Check the output for population dynamics. Large fluctuations might indicate instabilities or a poor trial wavefunction.

## Equilibration

Ensure that the energy has equilibrated. The initial blocks might still be relaxing from the VMC distribution to the DMC distribution (mixed distribution). Use `dmc_nblkeq` to discard these initial blocks.

## Analysis Tools

Use the provided analysis scripts to plot the energy trace and perform blocking analysis to estimate errors accurately.
