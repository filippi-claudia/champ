! constants, flags, common blocks for
! additional properties
! F. Schautz

! number of properties
      parameter(MAXPROP=6) 

! iprop: flag whether properties should be sampled
! nprop: number of properties sampled, print-out
      common /prp000/ iprop,nprop,ipropprt

! values
      common /prp001/ vprop(MAXPROP)
! sum, cum, square
      common /prp003/ vprop_sum(MAXPROP),vprop_cum(MAXPROP)
      common /prp003/ vprop_cm2(MAXPROP),cc_nuc(3)
