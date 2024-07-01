"""Functions and classes for CHAMP I/O."""
from glob import glob
from os import remove
from os.path import exists
from pathlib import Path, PosixPath
from typing import Optional, Type, Union, Any

from pydantic import BaseModel, DirectoryPath, Field, FilePath, Extra, validator


class Settings(BaseModel, extra=Extra.allow):
    """Data class containing CHAMP configuration.

    This class can hold the neccessery configuration to run CHAMP.
    """
    # class General(BaseModel, extra=Extra.allow):
    #     """General module class.
    #     """

    title: str = 'champ'
    """:obj:`str`: title of the simulation."""

    mode: str = 'vmc_one_mpi1'
    """:obj:`str`, optional: QMC mode."""

    seed: int = 1837465927472523
    """int, optional: seed for CHAMP."""

    class Force(BaseModel, extra=Extra.allow):
        """Force module class.
        """

        alfgeo: float = 0.0
        iforce_analy: int = 1
        node_cutoff: int = 1
        enode_cutoff: float = 0.1

    class Optwf(BaseModel, extra=Extra.allow):
        """Wavefunction optimization module class.
        """

        ioptwf: int = 1
        ioptci: int = 1
        ioptjas: int = 1
        ioptorb: int = 1
        method: str = 'sr_n'
        nextorb: int = 100
        nblk_max: int = 100
        no_active: int = 0 
        nopt_iter: int = 10
        sr_tau: float = 0.05
        sr_eps: float = 0.005
        sr_adiag: float = 0.01

    class BlockingVmc(BaseModel, extra=Extra.allow):
        """VMC module class.
        """

        vmc_nstep: int = 20
        vmc_nblk: int = 400
        vmc_nblkeq: int = 1
        vmc_nconf_new: int = 0

    class BlockingDmc(BaseModel, extra=Extra.allow):
        """DMC module class.
        """

        dmc_nstep: int = 40
        dmc_nblk: int = 20
        dmc_nblkeq: int = 1
        dmc_nconf: int = 100
        dmc_ivd: int = 1
        ipathak: int = 1
        eps_max: float = 0.05
        deps: float = 0.005
        nwprod: int = 20
        tau: float = 0.05
        etrial: float = -31.0
        icasula: int = -1


    # general: General = General()

    trexio: Optional[Path] = Field(default=None, prefix="load ")
    """:obj:`path <pathlib.Path>`: path to trexio."""

    determinants: FilePath = Field(default=None, prefix="load ")
    """:obj:`path <pathlib.Path>`: path to the determinants file."""

    orbitals: FilePath = Field(default=None, prefix="load ")
    """:obj:`path <pathlib.Path>`: path to the orbitals file."""

    jastrow: FilePath = Field(default=None, prefix="load ")
    """:obj:`path <pathlib.Path>`: path to the jastrow file."""

    jastrow_der: FilePath = Field(default=None, prefix="load ")
    """:obj:`path <pathlib.Path>`: path to the jastrow derivatives file."""

    previous: bool = False


    optgeo: Optional[Force] = Force()
    optwf: Optional[Optwf] = Optwf()
    blocking_vmc: Optional[BlockingVmc] = BlockingVmc()
    blocking_dmc: Optional[BlockingDmc] = BlockingDmc()

    # @validator('blocking_dmc', always=True)
    # def mutually_exclusive(cls, v, values):
    #     if values['blocking_vmc'] is not None and v is not None:
    #         raise ValueError('blocking_vmc and blocking_dmc are mutually exclusive')
    #     if values['optwf'] is not None and v is not None:
    #         raise ValueError('optwf and blocking_dmc are mutually exclusive')

    def model_post_init(self, __context: Any) -> None:
        if self.mode == 'vmc_one_mpi1':
            self.blocking_dmc = None
        elif self.mode == 'dmc_one_mpi1':
            self.blocking_vmc = None
            self.optwf = None

    def write(self, filename='vmc.inp'):
        """Write a this dataclass containing the CHAMP configuration to an input file.

        Args:
            filename (:obj:`str`, optional): input file to write to.
        """

        input_file = open(filename, 'w', encoding="utf-8")
        exclude = ['previous']
        schema = self.schema()['properties']
        for _, item in enumerate(self):
            if item[0] in exclude:
                continue
            if isinstance(item[1], BaseModel):
                if callable(getattr(item[1], 'schema')):
                    schema2 = item[1].schema()['properties']
                else:
                    schema2 = []
                input_file.write("\n%module " + item[0] + "\n")
                for _, item2 in enumerate(item[1]):
                    if item2[0] in schema2 and 'postfix' in schema2[item2[0]]:
                        if item2[1] is not None:
                            input_file.write("\t" + item2[0] + " " + str(item2[1]) +\
                                 schema2[item2[0]]['postfix'] + "\n")
                    else:
                        if item2[1] is not None:
                            input_file.write("\t" + item2[0] + " " + str(item2[1]) + "\n")
                input_file.write("%endmodule\n\n")
            else:
                if item[0] in schema and 'prefix' in schema[item[0]]:
                    if item[1] is not None:
                        input_file.write(schema[item[0]]['prefix'] + item[0] + " " +\
                             str(item[1]) + '\n')
                else:
                    if item[1] is not None:
                        input_file.write(item[0] + " " + str(item[1]) + '\n')

        input_file.close()

    @classmethod
    def read(cls: Type['BaseModel'], filename: Union[str, Path]) -> 'BaseModel':
        """Read the CHAMP input file and convert it to a dictionary format

        Args:
            filename (:obj:`str`, optional): file of the input file to read and write to.

        Returns:
            :obj:`Settings`: Settings object containing CHAMP configuration.
        """
        # Check if the file exists
        if not exists(filename):
            raise FileNotFoundError(filename + " was not found!")

        path = PosixPath(filename).resolve().parent

        try:
            path = path.relative_to(Path.cwd()).as_posix() + "/"
        except:
            path = path.as_posix() + "/"

        output = {}

        input_file = open(filename, 'r', encoding='utf-8')

        # Set some variables to keep track of modules
        curmod = ""
        inmod = False

        for line in input_file:
            # We skip empty lines
            if not line.strip():
                continue
            # Beginning of a module
            if line.startswith("%module"):
                curmod = line[8:-1]
                inmod = True
                output[curmod] = {}
            # End of a module
            elif line.startswith("%endmodule"):
                inmod = False
            # Loading in a file
            elif line.startswith("load"):
                temp = line[5:-1]
                temp = temp.split()
                key = temp[0]
                temp.pop(0)
                value = ' '.join(temp)
                output[key] = path + value
            # All other tags
            else:
                temp = line[:-1].split()
                key = temp[0]
                temp.pop(0)
                value = ' '.join(temp)
                if inmod:
                    # If we are in a module, add the tag to that subdictionary.
                    if key != "pool":
                        output[curmod][key] = value
                    else:
                        output[curmod][key] = path + value
                else:
                    output[key] = value

        input_file.close()

        return cls(**output)


    def use_opt_wf(self):
        """Function to replace the orbitals, determinants and jastrow with the optimized files.
        """
        opt = ["det_optimal.1.iter*", "orbitals_optimal.1.iter*", "jastrow_optimal.1.iter*"]
        keys = ["determinants", "orbitals", "jastrow"]

        for ind, name in enumerate(opt):
            num = len(glob(name))
            # If there are no optimized WF files, we do not use them
            index = min(num, self.optwf.nopt_iter)
            if num >= index and index > 0:
                opt[ind] = opt[ind].strip('*') + str(index)
                setattr(self, keys[ind], opt[ind])

    def todict(self):
        """Return a dictionary of the class.

        Returns:
            :obj:`str`: json dictionary.
        """
        return self.json(exclude_none=True)
    
def write_jastrow(atoms, filename='jastrow.0'):
    with open(filename, 'w') as jastrow_file:
        jastrow_file.write('jastrow_parameter   1\n')
        jastrow_file.write('  5  5  0           norda,nordb,nordc\n')
        jastrow_file.write('   0.60000000         scalek\n')
        for i in range(atoms):
            jastrow_file.write('   0.00000000   0.00000000   0.00000000   0.00000000  0.00000000   0.00000000 (a(iparmj),iparmj=1,nparma)\n')
        jastrow_file.write('   0.50000000   0.00000000   0.00000000   0.00000000   0.00000000   0.00000000 (b(iparmj),iparmj=1,nparmb)\n')
        for i in range(atoms):
            jastrow_file.write(' (c(iparmj),iparmj=1,nparmc)\n')
        jastrow_file.write('end')

def write_jastrow_der(atoms, filename='jastrow_der'):
    with open(filename, 'w') as der_file:
        der_file.write('jasderiv\n')
        der_file.write('4 4 4  5 0 0 0 0 0 0 nparma,nparmb,nparmc,nparmf\n')
        for i in range(atoms):
            der_file.write('3 4 5 6   (iwjasa(iparm),iparm=1,nparma)\n')
        der_file.write('2 3 4 5 6 (iwjasb(iparm),iparm=1,nparmb)\n')
        for i in range(atoms):
            der_file.write('(c(iparmj),iparmj=1,nparmc)\n')
        der_file.write('end')


def write_determinant(electron_num, filename='det.0'):
    with open(filename, 'w') as file:
        file.write('determinants 1 1 \n')
        file.write('1.00000000\n')
        dets = [i+1 for i in range(electron_num//2)]
        file.write(' '.join(map(str, dets)))
        file.write(' ')
        file.write(' '.join(map(str, dets)))
        file.write('\n')
        file.write('end')        



def cleanup(*args):
    """Remove files created by CHAMP.

    This function can cleanup the directory of the simulation.
    Add files as arguments to include (if not in default list)
    or exclude (if in default list) files.

    Args:
        *args: filesnames to include (if not in default list) or exclude (if in default list).
    """
    # Default list of the files that will be removed
    standard = ["force_analytic", "restart_vmc", "output.log", "parser.log", "orbitals_optimal.1.*",
                "jastrow_optimal.1.*", "det_optimal.1.*", "geo_optimal.*", "geo_optimal_final",
                "mc_configs_start", "nohup.out", "vmc_temp.inp"]

    # Add or remove file depending on the input arguments
    for item in args:
        if item in standard:
            standard.remove(item)
        else:
            standard.append(item)

    for file in standard:
        # Take care of the wildcard files
        if file.find('*') != -1:
            for file2 in glob(file):
                remove(file2)
        else:
            if exists(file):
                remove(file)