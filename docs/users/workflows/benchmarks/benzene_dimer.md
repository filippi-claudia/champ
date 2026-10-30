# Benzene Dimer

This benchmark is used to test the code's performance on larger molecular systems and its parallel scaling efficiency.

## System

-   System: Two benzene molecules (C₆H₆) in a stacked configuration.
-   Electrons: 84 valence electrons (with ECPs).

## Workflow

1.  **MPI Setup**: Run with varying numbers of MPI tasks (e.g., 1, 2, 4, 8 nodes).
2.  **Performance**: Measure the time per block (`time/blk` in output).
3.  **Scaling**: The speedup should be nearly linear with the number of cores.

## Note

Ensure that `vmc_nconf` (walkers per core) is kept constant or adjusted appropriately for strong/weak scaling tests.
