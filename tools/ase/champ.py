"""CHAMP ASE calculator

"""

from pathlib import Path
import os
from datetime import datetime
import numpy as np
from champio import write_jastrow, write_jastrow_der, write_determinant

from ase.calculators.genericfileio import (BaseProfile, CalculatorTemplate,
                                           GenericFileIOCalculator)
from ase.io import write

class CHAMPProfile(BaseProfile):
    def __init__(self, binarydir, **kwargs):
        super().__init__(**kwargs)
        self.binarydir = binarydir

    def get_calculator_command(self, mode):
        if mode == 'vmc_one_mpi1':
            append = '/bin/vmc.mov1'
        elif mode == 'dmc_one_mpi1':
            append = '/bin/dmc.mov1'
        return [self.binarydir + append, '-i', 'input.inp', '-o', 'output.out', '-e', 'error']
    
    def version(self):
        return 0

class CHAMPTemplate(CalculatorTemplate):
    _label = 'CHAMP'

    def __init__(self):
        super().__init__('CHAMP', ['energy', 'forces'])

    def write_input(self, profile, directory, atoms, parameters, properties):
        self.profile = profile
        self.directory = directory
        self.atoms = atoms
        self.parameters = parameters


    def execute(self, directory, profile):
        homedir = self.directory / f'champ_{datetime.now().strftime("%Y%m%d%H%M%S")}'
        homedir.mkdir(exist_ok=True)
        os.chdir(homedir)

        for i in range(len(self.parameters)):
            step = self.parameters[i]
            workdir = self.directory / str(i)
            workdir.mkdir(exist_ok=True, parents=True)
            os.chdir(workdir)

            if step[0] == 'CHAMP':
                if step[1].previous:
                    if self.parameters[i-1][0] == 'PySCF':
                        import trexio
                        trexiofile = trexio.File(f'../{i-1}/molecule.hdf5', mode="r", back_end=trexio.TREXIO_HDF5)
                        electron_num = trexio.read_electron_num(trexiofile)
                        trexiofile.close()
                        # step[1].trexio = f'../{i-1}/molecule.hdf5'
        
                        step[1].trexio = self.directory / '..' / str(i-1)/ 'molecule.hdf5'
                        write_jastrow(len(set(self.atoms.get_chemical_symbols())))
                        write_jastrow_der(len(set(self.atoms.get_chemical_symbols())))
                        write_determinant(electron_num)
                        step[1].jastrow = self.directory / '..' / str(i)/'jastrow.0'
                        step[1].jastrow_der = self.directory / '..' / str(i)/'jastrow_der'
                        step[1].determinants = self.directory / '..' / str(i)/'det.0'
                    if self.parameters[i-1][0] == 'qp2':
                        import trexio
                        trexiofile = trexio.File(f'../{i-1}/cipsi.hdf5', mode="r", back_end=trexio.TREXIO_HDF5)
                        electron_num = trexio.read_electron_num(trexiofile)
                        trexiofile.close()
                        # step[1].trexio = f'../{i-1}/molecule.hdf5'
        
                        step[1].trexio = self.directory / '..' / str(i-1)/ 'cipsi.hdf5'
                        write_jastrow(len(set(self.atoms.get_chemical_symbols())))
                        write_jastrow_der(len(set(self.atoms.get_chemical_symbols())))
                        write_determinant(electron_num)
                        step[1].jastrow = self.directory / '..' / str(i)/'jastrow.0'
                        step[1].jastrow_der = self.directory / '..' / str(i)/'jastrow_der'
                        step[1].determinants = self.directory / '..' / str(i-1)/'det.cipsi'
                    if self.parameters[i-1][0] == 'CHAMP':
                        step[1].jastrow_der = self.parameters[i-1][1].jastrow_der 
                        step[1].trexio = self.parameters[i-1][1].trexio
                        optsteps = self.parameters[i-1][1].optwf.nopt_iter

                        if self.parameters[i-1][1].optwf.ioptwf == 1 and self.parameters[i-1][1].optwf.ioptjas ==1:
                            step[1].jastrow = self.directory / '..' / str(i-1) / f'jastrow_optimal.1.iter{optsteps}'
                        else:
                            step[1].jastrow = self.parameters[i-1][1].jastrow

                        if self.parameters[i-1][1].optwf.ioptwf == 1 and self.parameters[i-1][1].optwf.ioptci ==1:
                            step[1].determinants = self.directory / '..' / str(i-1) / f'det_optimal.1.iter{optsteps}'
                        else:
                            step[1].determinants = self.parameters[i-1][1].determinants

                        if self.parameters[i-1][1].optwf.ioptwf == 1 and self.parameters[i-1][1].optwf.ioptorb ==1:
                            step[1].orbitals= self.directory / '..' / str(i-1) / f'orbitals_optimal.1.iter{optsteps}'
                        else:
                            step[1].orbitals = self.parameters[i-1][1].orbitals
                        if step[1].mode == 'dmc_one_mpi1':
                            shutil.copy(self.directory/'..'/str(i-1)/'mc_configs', self.directory/'mc_configs')

                step[1].write(self.directory / 'input.inp')
                profile.run(directory, inputfile=step[1].mode, outputfile='output.out')

                if step[1].mode == 'vmc_one_mpi1' and step[1].blocking_vmc.vmc_nconf_new > 0:
                    import glob
                    import shutil

                    # Get a list of all files that match the pattern
                    files = glob.glob('mc_configs_new*')

                    # Open the target file in append mode
                    with open('mc_configs', 'ab') as outfile:
                        for filename in files:
                            # Open each source file and append it to the target file
                            with open(filename, 'rb') as infile:
                                shutil.copyfileobj(infile, outfile)

                # execute
            elif  step[0] == 'PySCF':
                from pyscf import gto, scf
                from pyscf.lib import chkfile
                from trexio_tools.converters import pyscf_to_trexio
                from trexio_tools.converters import convert_to

                write('temp.xyz', self.atoms)
                mol = gto.M(atom='temp.xyz')
                mol.basis =  step[1]['basis']
                mol.ecp = step[1]['ecp']
                mol.symmetry=True
                mol.cart=False
                mol.build()
                mf = scf.RHF(mol)
                mf.chkfile = 'pyscf.chk'
                mf.init_guess = 'huckel'
                mf.run()

                pyscf_to_trexio.pyscf_to_trexio('pyscf.chk', 'pyscf.hdf5')
                convert_to.run('pyscf.hdf5', 'molecule.hdf5', 'cartesian', None)

            elif step[0] == 'qp2' and self.parameters[i-1][0] == 'PySCF':
                import subprocess
                from ezfio import ezfio_obj
                ezfio = ezfio_obj()
                this_env = os.environ.copy()
                subprocess.check_call(f'qp_import_trexio.py ../{i-1}/molecule.hdf5 -o cipsi', shell=True, env=this_env)
                # subprocess.check_call(f'ezfio set_file cipsi', shell=True, env=os.environ)
                ezfio.set_file('cipsi')
                for j in range(len(step[1])):
                    command = 'set_' + step[1][j][0] + '_' + step[1][j][1]
                    f = getattr(ezfio,command)
                    f(step[1][j][2])
                    # subprocess.check_call(f'echo {step[1][j][2]} | ezfio.py set {step[1][j][0]} {step[1][j][1]}', shell=True, env={**this_env, 'EZFIO_FILE': 'cipsi'})
                subprocess.check_call(f'qp_run fci cipsi > out', shell=True, env={**this_env, 'EZFIO_FILE': 'cipsi'})
                subprocess.check_call(f'qp_run save_for_champ cipsi', shell=True, env={**this_env, 'EZFIO_FILE': 'cipsi'})
                # subprocess.check_call(f'echo cipsi.hdf5 | ezfio.py set trexio trexio_file', shell=True, env={**this_env, 'EZFIO_FILE': 'cipsi'})
                ezfio.set_trexio_trexio_file('cipsi.hdf5')
                subprocess.check_call(f'qp_run export_trexio cipsi', shell=True, env={**this_env, 'EZFIO_FILE': 'cipsi'})
                with open('det.cipsi', 'r+') as file:
                        lines = file.readlines()
                        lines[1] = lines[1].strip() + ' 1\n'
                        file.seek(0)
                        file.writelines(lines)
            os.chdir('..')

    def read_results(self, directory):
        self.results = {}
        dir = len(self.parameters)-1
        os.chdir(self.directory / str(dir))
        atom_num = len(self.atoms.get_chemical_symbols())
        if self.parameters[-1][1].mode == 'vmc_one_mpi1':
            for line in reversed(open("output.out").readlines()):
                if line[:7] == 'total E':
                    self.results['energy'] = float(line.split()[3])
                    break
            self.results['forces'] = -1*np.loadtxt('force_analytic', dtype=float, usecols=[1,2,3])[-atom_num:].reshape((atom_num,3))
        elif self.parameters[-1][1].mode == 'dmc_one_mpi1':
            for line in reversed(open("output.out").readlines()):
                if line[:19] == 'total energy ( 100)':
                    self.results['energy'] = float(line.split()[5])
                    break
            self.results['forces'] = -1*np.loadtxt('force_analytic', dtype=float, usecols=[2,3,4]).reshape((atom_num,3))
        os.chdir('../..')
        return self.results


    def load_profile(self, cfg, **kwargs):
        return CHAMPProfile.from_config(cfg, self.name, **kwargs)


class CHAMP(GenericFileIOCalculator):
    def __init__(
        self,
        *,
        profile=None,
        directory='.',
        parallel_info=None,
        parallel=True,
        **kwargs):
        
        template = CHAMPTemplate()
        super().__init__(
            profile=profile,
            template=template,
            directory=directory,
            parallel_info=parallel_info,
            parallel=parallel,
            parameters=kwargs)
        


