module load 2022
module load CMake/3.23.1-GCCcore-11.3.0
module load OpenMPI/4.1.4-NVHPC-22.7-CUDA-11.7.0

cmake -H. -B build \
    -DCMAKE_Fortran_COMPILER=mpif90 \
    -DCMAKE_C_COMPILER=mpicc \
    -DENABLE_GPU=ON \
    -DNVFORTRAN_PATH=/sw/arch/RHEL8/EB_production/2022/software/NVHPC/22.7-CUDA-11.7.0/Linux_x86_64/22.7/compilers/bin
cmake --build build -j 12
