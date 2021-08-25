cp -v mc_configs_start_save mc_configs_start
mpirun -np 1  ../../../../bin/vmc.mov1 -i vmc_only.inp  -o vmc_only_single.out  -e error_vmc_only 
cp -v mc_configs_start_save mc_configs_start
mpirun -np 2  ../../../../bin/vmc.mov1 -i vmc_only.inp  -o vmc_only_double.out  -e error_vmc_only
