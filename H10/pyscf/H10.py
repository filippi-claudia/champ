#!/usr/bin/env python

import numpy as np
import sys
import os 
#sys.path.append('/p/home/jusers/landinez1/juwels/Codes/pyscf/') # if local installation


import pyscf
from pyscf import gto, scf, dft, cc

os.environ['OMP_NUM_THREADS'] = "72"







#passing coordinates in Borh
coords=np.zeros((10,3),dtype=np.float64)
coords[:,2]=(-4.5+np.arange(0,10,1))*1.5

basisset={'H': gto.parse('''
##NWchem
H S
      23.843185  0.00411490
      10.212443  0.01046440
       4.374164  0.02801110
       1.873529  0.07588620
       0.802465  0.18210620
       0.343709  0.34852140
       0.147217  0.37823130
       0.063055  0.11642410
H S
       0.139013  1.00000000
H P
       0.740212  1.00000000
'''
)}



mol = pyscf.M(atom=[('H', r) for r in coords],
              basis=basisset,
              max_memory=40000)

mol.ecp={'H': gto.basis.parse_ecp(
'''
H nelec 0
H ul
1 21.24359508259891 1.00000000000000
3 21.24359508259891 21.24359508259891
2 21.77696655044365 -10.85192405303825
H S
2 1.000000000000000 0.00000000000000
'''
)}

mol.unit     = 'B'
mol.symmetry = False
mol.verbose = 5
#cell.build(cart=True)
mol.build(cart=False)
print(mol.atom_coords())

#initial guess
mf = scf.RHF(mol)
mf.chkfile = 'H-diis.chk'
mf.max_cycle = 2
mf.kernel()



#fast newton
mf = scf.RHF(mol)
mf.__dict__.update(scf.chkfile.load('H-diis.chk', 'scf'))
scfdump= "RHF.dump"               
mf.chkfile = scfdump;
mf=scf.fast_newton(mf)
print('HF energy %.17g ' % mf.e_tot)
ehf=mf.e_tot

#Get kinetic Energy

dm = mf.make_rdm1()

#build h1e first
#pure kinetic part
t = mol.intor_symmetric('int1e_kin')
print("t shape", t.shape)
t = np.einsum('ij,ji', t, dm).real
print("Kinectic energy only", t)

#pseudopotential part/ pen
pp=0.0
pen1=0.0

if mol._pseudo:
    # Although mol._pseudo for GTH PP is only available in Cell, GTH PP
    # may exist if mol is converted from cell object.
    from pyscf.gto import pp_int
    pen1 = pp_int.get_gth_pp(mol)
else:
    pen1 = mol.intor_symmetric('int1e_nuc')
    
if len(mol._ecpbas) > 0:
    pp = mol.intor_symmetric('ECPscalar')

pen=pp+pen1

print("pp",pp)
print("pen1",pen1)
print("pen",pen)
    
pen = np.einsum('ij,ji', pen, dm).real
print("pen (energy)",pen)

h1e_t=t+pen
print('h1e_t',h1e_t)


# pure pyscf function for testing
h1e = mf.get_hcore()
h1e = np.einsum('ij,ji->', h1e, dm).real
print('h1e (energy)',h1e)

##Returns matrix Vhf = 2*J - K.  Vhf can be a list matrices, corresponding to the rdm
#pyscf function
vhf = mf.get_veff(mol, dm)
ee_coul = np.einsum('ij,ji->', vhf, dm).real * .5
print('ee_coul',ee_coul)

te_t=h1e_t+ee_coul
te=h1e+ee_coul

print("Total electron Energy test", te_t)
print("Total electron Energy pyscf", te)


#nuclear repulsion energy
p_nn = mf.energy_nuc()
print("p_nn",p_nn)


print("Total Energy test", te_t+p_nn)
print("Total Energy split pyscf", te+p_nn)
print("Total Energy pyscf", ehf)

