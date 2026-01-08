# DMC Input Keywords

## `blocking_dmc` Module

| Keyword | Type | Default | Description |
| :--- | :--- | :--- | :--- |
| `dmc_nstep` | Integer | 1 | Number of steps per block. |
| `dmc_nblk` | Integer | 1 | Number of blocks. |
| `dmc_nblkeq` | Integer | 0 | Number of equilibration blocks. |
| `dmc_nconf` | Integer | 1 | Number of walkers per MPI task. |

## `dmc` Module

| Keyword | Type | Default | Description |
| :--- | :--- | :--- | :--- |
| `tau` | Real | 0.01 | Time step (atomic units). |
| `etrial` | Real | Required | Trial energy (should be close to VMC energy). |
| `icasula` | Integer | 0 | Controls non-local pseudopotential moves. <br> 0: Locality approximation <br> -1: T-moves (Casula 2006) |

## Example

```python
%module blocking_dmc
    dmc_nstep     40
    dmc_nblk      100
    dmc_nconf     100
%endmodule

%module dmc
    tau           0.005
    etrial        -17.24
    icasula       -1
%endmodule
```
