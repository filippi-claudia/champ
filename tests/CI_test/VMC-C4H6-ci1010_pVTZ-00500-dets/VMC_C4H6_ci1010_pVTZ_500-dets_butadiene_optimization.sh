echo "VMC butadiene ci1010 pVTZ 500 determinants"

input="vmc_optimization_500.inp"
output="vmc_optimization_500"

# unicore test
N=1
ReferenceEnergy=-26.1993642
ReferenceError=0.0233401
mpirun -np $N ../../../bin/vmc.mov1 -i $input -o ${output}_core_${N}.out -e error
echo "Comparing energy with reference Core=$N           (total E = $ReferenceEnergy +-  $ReferenceError ) "
../../../tools/compare_value.py ${output}_core_${N}.out     "total E"  $ReferenceEnergy     $ReferenceError


# Multicore test
N=2
ReferenceEnergy=-26.181223
ReferenceError=0.0178963
mpirun -np $N ../../../bin/vmc.mov1 -i $input -o ${output}_core_${N}.out -e error
echo "Comparing energy with reference Core=$N           (total E = $ReferenceEnergy +-  $ReferenceError ) "
../../../tools/compare_value.py ${output}_core_${N}.out     "total E"  $ReferenceEnergy     $ReferenceError
