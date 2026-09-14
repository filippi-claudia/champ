---
title: Outputs Keywords
tags:
    - input
    - outputs
    - keywords
    - forces
---

# Output File Keywords

The `%module outputs` block names the optional files a run writes. Every keyword
in it has a default, so the block itself is optional: leaving it out reproduces
the file names CHAMP has always used.

Its purpose is to keep concurrent runs from writing over each other. A VMC and a
DMC calculation started in the same directory both write their analytic forces
to `force_analytic` unless they are told otherwise, and whichever finishes last
wins.

```perl
%module outputs
    file_force_analytic  'force_vmc.dat'
%endmodule
```

## The `outputs` module

<div class="grid cards single-col" markdown>

-   __Analytic force file__

    ---

    Name of the file the analytic forces are written to at the end of a VMC or
    DMC run. Only written when `iforce_analy 1` is set in `%module general`.
    Include the extension you want; the name is used verbatim, relative to the
    working directory. Quotes are optional.

    Default: `force_analytic`

    ```perl
    file_force_analytic  'force_dmc.dat'
    ```

</div>

## Running VMC and DMC in the same directory

Give each run its own name and both sets of forces survive:

```perl
# vmc.inp
%module general
    title           'butadiene forces'
    mode            'vmc_one_mpi'
    iforce_analy     1
%endmodule

%module outputs
    file_force_analytic  'force_vmc.dat'
%endmodule
```

```perl
# dmc.inp
%module general
    title           'butadiene forces'
    mode            'dmc_one_mpi1'
    iforce_analy     1
%endmodule

%module outputs
    file_force_analytic  'force_dmc.dat'
%endmodule
```
