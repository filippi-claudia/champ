mpirun -np 8 ../../../bin/vmc.mov1  -i vmc_h2o_dft_jas2body_optall.inp -o vmc_h2o_dft_jas2body_optall.out -e error

cat mc_configs_new* >> mc_configs
\rm mc_configs_new*

mpirun -np 8 ../../../bin/dmc.mov1  -i dmc_h2o_dft_jas2body_optall.inp -o dmc_h2o_dft_jas2body_optall.out -e error

\rm problem*
\rm mc_configs_new*
