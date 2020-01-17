The Makefile is (UNIX-) platform-independend, but you must have GNU-make
(gmake) to use it. (On Linux make is gmake)

The subdirectory make_config contains platform specific files with 
settings for the fortran compiler used by the Makefile.

You might have to create the settings for your platform. Choose the 
appropriate file and copy (or symlink) it to settings.make

`make clean` to get rid of old object files,
`make vmc` to build the vmc code
`make dmc` to build the dmc code

Shortcuts:

`make seq` to build the vmc,fit,dmc codes
`make mpi` to build the mpi version of the vmc,fit,dmc
`make all` to build seq,mpi

In the subdirectory `input` there are awk-scripts which are used
in context with the input reader (during compilation) and you might 
have to adjust the path to awk at the top of those scripts.
