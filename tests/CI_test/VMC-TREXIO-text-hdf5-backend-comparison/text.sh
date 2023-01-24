echo "TREXIO backend comparison : TEXT"
input="vmc_opt_ci1010_pVTZ_1522_text.inp"
output="vmc_opt_ci1010_pVTZ_1522_text"

# unicore test
N=1
mpirun -np $N ../../../bin/vmc.mov1 -i ${input} -o ${output}_core_${N}.out -e error

