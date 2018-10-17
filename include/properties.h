C constants, flags, common blocks for
C additional properties
C F. Schautz

c number of properties
      parameter(MAXPROP=6) 

c iprop: flag whether properties should be sampled
c nprop: number of properties sampled, print-out
      common /prp000/ iprop,nprop,ipropprt

c values
      common /prp001/ vprop(MAXPROP)
c sum, cum, square
      common /prp003/ vprop_sum(MAXPROP),vprop_cum(MAXPROP)
     $     ,vprop_cm2(MAXPROP),cc_nuc(3)
