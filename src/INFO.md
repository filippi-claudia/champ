Info for 3d-grid

1) Use of real\*8, real\*4 

-> Choice in 3dgrid.h => store in \*4/\*8

-> Choice in spline_mo, setup_3dsplorb => calculate in \*4/\*8

Note: call to library subroutine need to be changed in the
code and must be compatible with the choice of storage \*4/\*8

In Lagrange interpolation, all operations are in real\*8 but 
you can still store in real\*4 <=> real\*4 in 3dgrid.h

2) settings_make

dyncom -> common blocks with 3dgrid

The memory is allocated only if you enter in the routines.
Therefore, if you do not use 3dgrid, you problem with size

3) read_input

`i3dsplorb  0, 1, 2`

2 -> always spline
1 -> only in determinante.f (call in hpsie not hpsi)
