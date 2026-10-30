# VMC Input Keywords

These keywords are typically found in the `blocking_vmc` module or the general `vmc` control section.

## `blocking_vmc` Module

| Keyword | Type | Default | Description |
| :--- | :--- | :--- | :--- |
| `vmc_nstep` | Integer | 1 | Number of Monte Carlo steps per block. |
| `vmc_nblk` | Integer | 1 | Number of blocks to accumulate for statistics. |
| `vmc_nblkeq` | Integer | 2 | Number of equilibration blocks (discarded). |
| `vmc_nconf_new` | Integer | 0 | Number of configurations (walkers) to save for future runs (e.g., DMC). If 0, no configurations are saved. |
| `vmc_nblk_max` | Integer | `vmc_nblk` | Maximum number of blocks (used in optimization to increase precision). |

## Other Keywords

| Keyword | Type | Default | Description |
| :--- | :--- | :--- | :--- |
| `etrial` | Real | 1.0 | Trial energy. Required for some calculations. |
| `vmc_nconf` | Integer | 1 | Number of configurations per MPI task (walkers). |

## Example

```python
%module blocking_vmc
    vmc_nstep     20      # 20 steps per block
    vmc_nblk      200     # 200 blocks
    vmc_nblkeq    10      # 10 blocks for equilibration
    vmc_nconf_new 100     # Save 100 walkers per core
%endmodule
```
