#python3 -m pip install .
#pip install .
trexio convert-from --type gamess --input ci44_pVTZ.out gamess_butadiene_ci44_cc-VTZ.hdf5 --back_end=HDF5
trexio convert2champ gamess_butadiene_ci44_cc-VTZ.hdf5 --input ci44_pVTZ.out --back_end=HDF5
