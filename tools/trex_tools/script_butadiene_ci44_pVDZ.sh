#!/bin/bash
python trex2champ.py 	--trex 		butadiene_ci44_pVDZ.hdf5 \
			--gamess 	butadiene_ci44_pVDZ.out \
			--motype 	"GUGA" \
			--backend	"HDF5" \
			--basis_prefix  "BFD-D" \
			--lcao \
			--geom \
			--basis \
			--ecp \
			--sym \
			--det

