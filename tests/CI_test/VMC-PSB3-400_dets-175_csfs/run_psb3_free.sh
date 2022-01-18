#!/bin/bash
set -e
if [ $# -eq 0 ]; then
  echo "No argument supplied, running with default of 1 processor"
  nproc=1
else
nproc=$1
fi 
filename=revised_free

mpirun -n $nproc ../../../bin/vmc.mov1 -i ${filename}.inp \
         > ${filename}.$nproc.out 2> >(tee ${filename}.$nproc.err >&2)
