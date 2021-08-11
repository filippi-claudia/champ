cp -v mc_configs_start_save mc_configs_start
mpirun -np 1  ../../../../bin/vmc.mov1 -i vmc_sr_opt_geom.inp  -o vmc_sr_opt_geom_single.out  -e error_vmc_optimization
mpirun -np 2  ../../../../bin/vmc.mov1 -i vmc_sr_opt_geom.inp  -o vmc_sr_opt_geom_double.out  -e error_vmc_optimization
