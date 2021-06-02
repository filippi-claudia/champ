cp -v mc_configs_start_save mc_configs_start
mpirun -np 1  /home/ravindra/softwares/neelravi-champ-modern/bin/vmc.mov1 -i test-champ.inp  -o test.out  -e error > screen 
#mpirun -np 1  /home/ravindra/softwares/neelravi-champ-modern/bin/dmc.mov1 -i test-champ.inp  -o test.out  -e error 
