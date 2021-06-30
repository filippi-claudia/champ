cp -v mc_configs_start_save mc_configs_start
mpirun -np 1  /home/ravindra/softwares/neelravi-champ-modern/bin/vmc.mov1 -i vmc_optimization.inp  -o vmc_optimization.out  -e error_vmc_opt_only 
