#!/bin/bash
set -e
if [ $# -eq 0 ]; then
  echo "No argument supplied, running with default of 1 processor"
  nproc=1
else
nproc=$1
fi 
filename=vmc_optimization_15000_3body

mpirun -n $nproc ../../../../bin/vmc.mov1 -i ${filename}.inp \
         > ${filename}.$nproc.out 2> >(tee ${filename}.$nproc.err >&2)
#../../../../tools/compare_value.py ${filename}.$nproc.out \
#   "total E" -1.0137312  0.0070186
