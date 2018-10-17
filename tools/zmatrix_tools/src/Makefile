# Name of Fortran Compiler
FC=ifort

# Fortran Compiler Options - Normal options
FC_Options_64_normal=
			
Install_dir=~/bin


all: normal_version molden_version

normal_version:  cart2zmat.x zmat2cart.x zmat2cart_rc.x 


molden_version:  cart2zmat_molden.x zmat2cart_molden.x zmat2cart_rc_molden.x 


cart2zmat.x:  geo_cart2zmat.f90 convcartzmat.f90


	$(FC) $(FC_Options_64_normal) \
	-o cart2zmat.x \
	geo_cart2zmat.f90 convcartzmat.f90 
	 

cart2zmat_molden.x:  geo_cart2zmat_molden.f90 convcartzmat.f90


	$(FC) $(FC_Options_64_normal) \
	-o cart2zmat_molden.x \
	geo_cart2zmat_molden.f90 convcartzmat.f90 


zmat2cart.x:  geo_zmat2cart.f90 convcartzmat.f90


	$(FC) $(FC_Options_64_normal) \
	-o zmat2cart.x \
	geo_zmat2cart.f90 convcartzmat.f90 
	

zmat2cart_molden.x:  geo_zmat2cart_molden.f90 convcartzmat.f90


	$(FC) $(FC_Options_64_normal) \
	-o zmat2cart_molden.x \
	geo_zmat2cart_molden.f90 convcartzmat.f90 


zmat2cart_rc.x:  geo_zmat2cart_rc.f90 convcartzmat.f90


	$(FC) $(FC_Options_64_normal) \
	-o zmat2cart_rc.x \
	geo_zmat2cart_rc.f90 convcartzmat.f90


zmat2cart_rc_molden.x:  geo_zmat2cart_rc_molden.f90 convcartzmat.f90


	$(FC) $(FC_Options_64_normal) \
	-o zmat2cart_rc_molden.x \
	geo_zmat2cart_rc_molden.f90 convcartzmat.f90 


clean:
	
	rm *.mod


veryclean: 
	           
#	rm cart2zmat.x zmat2cart.x zmat2cart_rc.x
#	rm cart2zmat_molden.x zmat2cart_molden.x zmat2cart_rc_molden.x


install:  all

	cp cart2zmat.x $(Install_dir)
	cp zmat2cart.x $(Install_dir)
	cp zmat2cart_rc.x $(Install_dir)
	cp cart2zmat_molden.x $(Install_dir)
	cp zmat2cart_molden.x $(Install_dir)
	cp zmat2cart_rc_molden.x $(Install_dir)
	

