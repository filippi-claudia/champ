echo "TREXIO DMC water 2body jastrow tau=0.05  "

input="vmc_h2o_dft_jas2body.inp"
output="vmc_h2o_dft_jas2body"

N=8
\rm -f mc_configs
mpirun -np $N ../../../bin/vmc.mov1 -i $input -o ${output}_core_${N}.out -e error

cat mc_configs_new* >> mc_configs
\rm mc_configs_new*

input="dmc_h2o_dft_jas2body_tau0.05.inp"
output="dmc_h2o_dft_jas2body_tau0.05"
N=8
ReferenceEnergy=-17.2624300
ReferenceError=0.0009044
mpirun -np $N ../../../bin/dmc.mov1 -i $input -o ${output}_core_${N}.out -e error
echo "Comparing energy with reference Core=$N           (total E = $ReferenceEnergy +-  $ReferenceError ) "
../../../tools/compare_value.py ${output}_core_${N}.out  "total energy ( 100) "  $ReferenceEnergy     $ReferenceError
