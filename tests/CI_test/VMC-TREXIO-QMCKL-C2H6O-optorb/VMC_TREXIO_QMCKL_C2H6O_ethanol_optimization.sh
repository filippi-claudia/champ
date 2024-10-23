echo "QMCKL + TREXIO VMC ethanol orbital optimization "

input="vmc_ethanol_opt.inp"
output="vmc_ethanol_opt"

N=2
mpirun -np $N ../../../bin/vmc.mov1 -i $input -o ${output}_core_${N}.out -e error

ReferenceEnergy=-30.5921225
ReferenceError=0.0152373
echo "Comparing energy with reference Core=$N           (total E = $ReferenceEnergy +-  $ReferenceError ) "
../../../tools/compare_value.py ${output}_core_${N}.out  "total E" $ReferenceEnergy     $ReferenceError
