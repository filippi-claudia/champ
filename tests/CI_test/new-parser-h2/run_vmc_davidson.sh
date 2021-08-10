mpirun -np 1 /home/ravindra/softwares/neelravi-champ-modern/bin/vmc.mov1 -i revised_vmc_davidson_check.inp -o revised_vmc_davidson_check_single.out -e error 
mpirun -np 2 /home/ravindra/softwares/neelravi-champ-modern/bin/vmc.mov1 -i revised_vmc_davidson_check.inp -o revised_vmc_davidson_check_double.out -e error 
