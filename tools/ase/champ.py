"""CHAMP ASE calculator

"""

from pathlib import Path

from ase.calculators.genericfileio import (BaseProfile, CalculatorTemplate,
                                           GenericFileIOCalculator)

class CHAMPProfile(BaseProfile):
    def __init__(self, binarydir, **kwargs):
        super().__init__(**kwargs)
        self.binarydir = binarydir

    def get_calculator_command(self, inputfile):
        return [self.binarydir, str(inputfile)]

class CHAMPTemplate(CalculatorTemplate):
    _label = 'CHAMP'

    def __init__(self):
        super().__init('CHAMP', ['energy', 'forces'])

    def write_input(self, profile, directory, atoms, parameters, properties):
        self.profile = profile
        self.directory = directory
        self.atoms = atoms
        self.parameters = parameters


    def execute(self, directory, profile):
        for i in range(self.parameters):
            step = self.parameters[i]
            self.directory.mkdir(i)
            if step[0] == 'CHAMP':
                step[1].write(self.directory / 'input.inp')
                # execute
            elif  step[0] == 'PySCF':
                # coonect via ase
                # do the HF
                # convert to trexio
                # convert to cartesian
                # fix ECP
        # Here I can run the custom recipe

    def read_results(self, directory: PathLike) -> Mapping[str, Any]:
        return 0


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
        


