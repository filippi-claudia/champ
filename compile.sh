#!/bin/bash

rm -rf build bin 
cmake -H. -Bbuild -DCMAKE_Fortran_COMPILER=mpifort
cmake --build build -- -j4 > output 2>&1
