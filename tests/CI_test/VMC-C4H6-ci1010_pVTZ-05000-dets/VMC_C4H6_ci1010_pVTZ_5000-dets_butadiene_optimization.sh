echo "VMC butadiene ci1010 pVTZ 5000 determinants"

input="vmc_optimization_5000.inp"
output="vmc_optimization_5000"

# unicore test
N=1
ReferenceEnergy=-26.2363885
ReferenceError=0.0225517
mpirun -np $N ../../../bin/vmc.mov1 -i $input -o ${output}_core_${N}.out -e error
echo "Comparing energy with reference Core=$N           (total E = $ReferenceEnergy +-  $ReferenceError ) "
../../../tools/compare_value.py ${output}_core_${N}.out     "total E"  $ReferenceEnergy     $ReferenceError


# Multicore test
N=2
ReferenceEnergy=-26.2077357
ReferenceError=0.0154788
mpirun -np $N ../../../bin/vmc.mov1 -i $input -o ${output}_core_${N}.out -e error
echo "Comparing energy with reference Core=$N           (total E = $ReferenceEnergy +-  $ReferenceError ) "
../../../tools/compare_value.py ${output}_core_${N}.out     "total E"  $ReferenceEnergy     $ReferenceError
