echo "VMC butadiene ci1010 pVTZ 500 determinants"

input="vmc_optimization_500.inp"
output="vmc_optimization_500"

# unicore test
N=1
ReferenceEnergy=-26.1886146
ReferenceError=0.0219769
mpirun -np $N ../../../bin/vmc.mov1 -i $input -o ${output}_core_${N}.out -e error
echo "Comparing energy with reference Core=$N           (total E = $ReferenceEnergy +-  $ReferenceError ) "
../../../tools/compare_value.py ${output}_core_${N}.out     "total E"  $ReferenceEnergy     $ReferenceError


# Multicore test
N=2
ReferenceEnergy=-26.1363913
ReferenceError=0.0244361
mpirun -np $N ../../../bin/vmc.mov1 -i $input -o ${output}_core_${N}.out -e error
echo "Comparing energy with reference Core=$N           (total E = $ReferenceEnergy +-  $ReferenceError ) "
../../../tools/compare_value.py ${output}_core_${N}.out     "total E"  $ReferenceEnergy     $ReferenceError
