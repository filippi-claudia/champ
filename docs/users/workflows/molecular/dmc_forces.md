# DMC Forces

Calculating atomic forces in DMC allows for geometry optimization and vibrational analysis at the QMC level.

## Overview

Force calculations in DMC are more expensive than energy calculations alone because they require computing the gradient of the energy with respect to nuclear coordinates.

## Configuration

To enable force calculations, you typically need to:
1.  Enable force computation in the input.
2.  Ensure that the pseudopotentials (if used) support force calculations.

## Example

```python
%module opt_geo
    ...
    iforce_analy        1   # Analytic forces
    alfgeo              0.5 # Force step size
%endmodule
```

*Check `tests/CI_test/VMC-C4H6-ci44_pVQZ_20_dets_12_csf_optall/vmc_optall_ci44.inp` for VMC force examples, which often translate to DMC.*
