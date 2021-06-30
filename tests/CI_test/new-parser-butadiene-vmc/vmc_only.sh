#cp -v mc_configs_start_save mc_configs_start
mpirun -np 8  /home/ravindra/softwares/neelravi-champ-modern/bin/vmc.mov1 -i vmc_only.inp  -o vmc_only.out  -e error_vmc_only > screen 
