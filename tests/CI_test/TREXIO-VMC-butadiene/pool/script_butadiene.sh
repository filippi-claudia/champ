#!/bin/bash
python3 ../../../../tools/trex_tools/trex2champ.py \
			--trex 	"butadiene_ci1010_pVTZ_1.hdf5" \
			--backend	"HDF5" \
			--basis_prefix  "BFD-T" \
			--lcao \
			--ecp \
			--sym \
			--geom \
			--basis

