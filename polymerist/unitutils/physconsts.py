'''Unit-aware compendium of useful physical constants'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

from scipy.constants import physical_constants
from .unitstr import unit_from_unit_str


PHYS_CONSTS_WITH_UNITS = {}
for key, (value, unit_str, uncert) in physical_constants.items():
    if (unit_val := unit_from_unit_str(unit_str)) is not None:
        PHYS_CONSTS_WITH_UNITS[key] = value * unit_val

# variable-names of physical constants from scipy.constants, which were in turn taken from NIST CODATA (https://physics.nist.gov/cuu/Constants/)
named_phys_consts = { 
    'c'         : ('speed of light in vacuum', 'speed_of_light'),
    'mu_0'      : ('vacuum mag. permeability', None),
    'epsilon_0' : ('vacuum electric permittivity', None),
    'h'         : ('Planck constant', 'Planck'),
    'hbar'      : ('Planck constant over 2 pi', None),
    'G'         : ('Newtonian constant of gravitation', 'gravitational_constant'),
    'g'         : ('standard acceleration of gravity', None),
    'e'         : ('elementary charge', 'elementary_charge'),
    'R'         : ('molar gas constant', 'gas_constant'),
    'alpha'     : ('fine-structure constant', 'fine_structure'),
    'N_A'       : ('Avogadro constant', 'Avogadro'),
    'k'         : ('Boltzmann constant', 'Boltzmann'),
    'sigma'     : ('Stefan-Boltzmann constant', 'Stefan_Boltzmann'),
    'Wien'      : ('Wien wavelength displacement law constant', None),
    'Rydberg'   : ('Rydberg constant', None),
    'm_e'       : ('electron mass', 'electron_mass'),
    'm_p'       : ('proton mass', 'proton_mass'),
    'm_n'       : ('neutron mass', 'neutron_mass'),
    'm_u'       : ('atomic mass constant', 'atomic_mass'),
}

for var_name, (CODATA_key, alias) in named_phys_consts.items():
    const_with_units = PHYS_CONSTS_WITH_UNITS[CODATA_key]
    globals()[var_name] = const_with_units # inject into global namespace
    if alias is not None:
        globals()[alias] = const_with_units