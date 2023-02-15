echo "TREXIO water optimization"

input="vmc_h2o_dft_optall.inp"
output="vmc_h2o_dft_optall"

# Multicore test
N=2
ReferenceEnergy=-17.2163601
ReferenceError=0.0040069
mpirun -np $N ../../../bin/vmc.mov1 -i ${input} -o ${output}_core_${N}.out -e error
echo "Comparing energy with reference Core=${N}           (total E = $ReferenceEnergy +-  $ReferenceError ) "
../../../tools/compare_value.py ${output}_core_${N}.out     "total E"  $ReferenceEnergy     $ReferenceError
