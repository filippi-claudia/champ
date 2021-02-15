#!/bin/bash 
rm -rf bin build
OMPI=/home/plopez/Programs/openmpi-2.1.1_linux/build/bin
cmake -H. -Bbuild \
	-DCMAKE_CXX_COMPILER=/usr/bin/g++ \
        -DCMAKE_Fortran_COMPILER=$OMPI/mpifort \
        -DENABLE_TEST=ON  \
        -DLAPACK_lapack_LIBRARY=/home/plopez/Programs/lapack-3.8.0/build/lib/liblapack.a \
        -DBLAS_blas_LIBRARY=/home/plopez/Programs/lapack-3.8.0/build/lib/libblas.a

 
#cmake -H. -Bbuild -DCMAKE_CXX_COMPILER=/usr/bin/g++ -DCMAKE_Fortran_COMPILER=mpifort -DENABLE_TEST=ON
cmake --build build --target all -- -j2
#cmake --build build --target vmc.mov1 -- -j2
#cmake --build build --target dmc.mov1 -- -j2

#        -DENABLE_QMMM=ON \
#        -DENABLE_PERIODIC=ON \
