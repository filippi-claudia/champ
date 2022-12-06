#!/bin/bash
python trex2champ.py 	--trex butadiene_ci1010_pVTZ.hdf5 \
			--motype 	"GUGA" \
			--backend	"HDF5" \
			--basis_prefix  "BFD-T" \
			--lcao \
			--geom \
			--basis \
			--ecp \
			--sym \
			--det \

