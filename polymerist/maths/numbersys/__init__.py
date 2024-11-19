'''Implementations of various number systems and representations of numbers'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

from .commonbases import FixedRadixNumberSystem, FactorialNumberSystem

# initialization of some common bases
Factoradic = FactorialNumberSystem()

COMMON_BASES : dict[str, int] = {
    'binary'      : 2,
    'ternary'     : 3,
    'seximal'     : 6,
    'octal'       : 8,
    'decimal'     : 10,
    'duodecimal'  : 12,
    'hexadecimal' : 16,
}
for base_name, base in COMMON_BASES.items():
    base_sys = FixedRadixNumberSystem(base)
    globals()[base_name.capitalize()]     = base_sys # systematic name
    globals()[f'base{base}'.capitalize()] = base_sys # common name

    globals()[f'nega{base_name}'.capitalize()] = FixedRadixNumberSystem(-base) # also include negative base for kicks