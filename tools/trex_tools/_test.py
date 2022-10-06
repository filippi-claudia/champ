#!/usr/bin/env python3
#   trex2champ is a tool which allows to read output files of quantum
#   chemistry codes (GAMESS and trexio files) and write input files for
#   CHAMP in V3.0 format.
#
# Copyright (c) 2021, TREX Center of Excellence
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, this
#    list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
#
# 3. Neither the name of the copyright holder nor the names of its
#    contributors may be used to endorse or promote products derived from
#    this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
#   Ravindra Shinde
#   University of Twente
#   Enschede, The Netherlands
#   r.l.shinde@utwente.nl



import sys
from trex2champ import Champ

# Before we do anything else, we need to check if trexio and resultsFile are installed
try:
    import trexio
except:
    print("Error: The TREXIO Python library is not installed")
    sys.exit(1)


def test_benzene_ground_state():
    champ = Champ()
    champ.filename="benzene.hdf5"
    champ.motype="RHF"
    champ.back_end='hdf5'
    champ.gamessfile=None
    champ.save_geometry=True
    champ.save_lcao = True
    champ.save_basis = True
    champ.save_eigenvalues = False
    champ.save_ecp = True
    champ.save_symmetry = False
    champ.save_determinants = False
    champ.save_csfs = False

    # Optional argument for controlling the names of the output files
    champ.basis_prefix = "TEST1"

    champ.run()
    assert champ is not None
    assert champ.nucleus_num == 12
    assert champ.ao_num == 114
    assert champ.mo_num == 108
    assert champ.shell_num == 48
    assert champ.prim_num == 186
    assert champ.ecp_num == 42



def test_formaldehyde_ground_state():
    champ = Champ()
    champ.filename="COH2_GS.trexio"
    champ.motype="RHF"
    champ.back_end='HDF5'
    champ.gamessfile=None
    champ.save_geometry=True
    champ.save_lcao = True
    champ.save_basis = True
    champ.save_eigenvalues = False
    champ.save_ecp = True
    champ.save_symmetry = False
    champ.save_determinants = True
    champ.save_csfs = False

    # Optional argument for controlling the names of the output files
    champ.basis_prefix = "TEST2"

    champ.run()
    assert champ is not None
    assert champ.nucleus_num == 4
    assert champ.ao_num == 66
    assert champ.mo_num == 66
    assert champ.shell_num == 26
    assert champ.prim_num == 62
    assert champ.ecp_num == 14
    assert champ.num_dets == 1862
    assert champ.num_states == 1


def test_butadiene():
    champ = Champ()
    champ.filename="butadiene_ci44_pVDZ.hdf5"
    champ.motype="GUGA"
    champ.back_end='HDF5'
    champ.gamessfile=None
    champ.save_geometry=True
    champ.save_lcao = True
    champ.save_basis = True
    champ.save_eigenvalues = False
    champ.save_ecp = True
    champ.save_symmetry = False
    champ.save_determinants = True
    champ.save_csfs = False

    # Optional argument for controlling the names of the output files
    champ.basis_prefix = "TEST3"

    champ.run()
    assert champ is not None
    assert champ.nucleus_num == 10
    assert champ.ao_num == 86
    assert champ.mo_num == 86
    assert champ.shell_num == 38
    assert champ.prim_num == 114
    assert champ.ecp_num == 34


if __name__ == "__main__":
    champ = Champ()
    champ.__main__()
    champ.run()
