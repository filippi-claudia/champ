echo "VMC butadiene ci44 pVQZ energy only"

input="vmc_noopt_ci44.inp"
output="vmc_noopt_ci44"

# unicore test
N=1
ReferenceEnergy=-26.2428668
ReferenceError=0.0135012
mpirun -np $N ../../../bin/vmc.mov1 -i $input -o ${output}_core_${N}.out -e error
echo "Comparing energy with reference Core=$N           (total E = $ReferenceEnergy +-  $ReferenceError ) "
../../../tools/compare_value.py ${output}_core_${N}.out     "total E"  $ReferenceEnergy     $ReferenceError

# Multicore test
N=2
ReferenceEnergy=-26.2380698
ReferenceError=0.0098286
mpirun -np $N ../../../bin/vmc.mov1 -i $input -o ${output}_core_${N}.out -e error
echo "Comparing energy with reference Core=$N           (total E = $ReferenceEnergy +-  $ReferenceError ) "
../../../tools/compare_value.py ${output}_core_${N}.out     "total E"  $ReferenceEnergy     $ReferenceError
