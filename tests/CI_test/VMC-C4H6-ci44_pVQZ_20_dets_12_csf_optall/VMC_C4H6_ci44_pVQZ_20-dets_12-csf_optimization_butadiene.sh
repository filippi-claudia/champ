echo "VMC butadiene ci44 pVQZ all optimizations "

input="vmc_optall_ci44.inp"
output="vmc_optall_ci44"

# unicore test
N=1
ReferenceEnergy=-26.2087347
ReferenceError=0.0167763
mpirun -np $N ../../../bin/vmc.mov1 -i $input -o ${output}_core_${N}.out -e error
echo "Comparing energy with reference Core=$N           (total E = $ReferenceEnergy +-  $ReferenceError ) "
../../../tools/compare_value.py ${output}_core_${N}.out     "total E"  $ReferenceEnergy     $ReferenceError


# Multicore test
N=2
ReferenceEnergy=-26.2110877
ReferenceError=0.0114863
mpirun -np $N ../../../bin/vmc.mov1 -i $input -o ${output}_core_${N}.out -e error
echo "Comparing energy with reference Core=$N           (total E = $ReferenceEnergy +-  $ReferenceError ) "
../../../tools/compare_value.py ${output}_core_${N}.out     "total E"  $ReferenceEnergy     $ReferenceError
