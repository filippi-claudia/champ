echo "TREXIO water 3body Jastrow"

input="vmc.inp"
output="vmc_3body"

# Multicore test
N=2
ReferenceEnergy=-17.2252374
ReferenceError=0.0015519
mpirun -np $N ../../../bin/vmc.mov1 -i ${input} -o ${output}_core_${N}.out -e error
echo "Comparing energy with reference Core=${N}           (total E = $ReferenceEnergy +-  $ReferenceError ) "
../../../tools/compare_value.py ${output}_core_${N}.out     "total E"  $ReferenceEnergy     $ReferenceError
