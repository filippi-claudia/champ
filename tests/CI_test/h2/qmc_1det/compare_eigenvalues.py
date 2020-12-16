#!/usr/bin/env python
import numpy as np
import scipy.linalg as spla
import subprocess


def get_champ_values(outfile='vmc_davidson_check.out'):
    """get the values outputed by CHAMP."""
    out = subprocess.Popen(['grep', 'DAV: eigv', 'vmc_davidson_check.out'],
                           stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    stdout, stderr = out.communicate()
    vals = stdout.decode('utf-8').split('\n')[-2].split()[-2:]
    return [float(v) for v in vals]


def get_ref_values(hpath='H.dat', spath='S.dat'):
    """Compute eigenvalues using scipy."""
    h = np.loadtxt(hpath).reshape(10, 10)
    s = np.loadtxt(spath).reshape(10, 10)
    v, u = spla.eig(h, s)
    v = v.real
    v.sort()
    return v[:2]


# get all values
vchamp = get_champ_values()
vref = get_ref_values()

# compare the difference
delta = np.abs(vchamp-vref)
assert np.all(delta < 0.01)
