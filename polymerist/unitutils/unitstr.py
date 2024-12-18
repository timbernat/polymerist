'''Utilities for looking up and producing Units from strings'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

from typing import Optional
from openmm import unit as openmm_unit
from openmm.unit import Unit, dimensionless

from . import extraunits
from ..genutils.iteration import product


# GENERATING REFERENCE TABLES OF DEFAULTS AND ADDITIONAL UNITS
UNITS_BY_SYMBOL : dict[str, Unit]= {}
UNITS_BY_NAME   : dict[str, Unit]= {}

for _module in (openmm_unit, extraunits):
    for attr in dir(_module):
        attr_val = getattr(_module, attr)
        if isinstance(attr_val, Unit):
            UNITS_BY_SYMBOL[attr_val.get_symbol()] = attr_val
            UNITS_BY_NAME[  attr_val.get_name()  ] = attr_val


# READING UNITS FROM STRINGS
def unit_from_unit_str(unit_str : str, unit_delimiter : str=' ', exp_delimiter : str='^') -> Optional[Unit]:
    '''Take a string describing units, with units and exponents delimited consistently (ex: "kg m^2 s^-2")
    and convert it into a proper OpenMM Unit if possible, returning None if not'''
    indiv_units = []
    for unit_part in unit_str.split(unit_delimiter):
        unit_and_exp = unit_part.split(exp_delimiter)
        unit_sym, exp = unit_and_exp if (len(unit_and_exp) == 2) else (unit_and_exp[0], 1)
        
        for unit_table in (UNITS_BY_SYMBOL, UNITS_BY_NAME):
            if (unit_val := unit_table.get(unit_sym)):
                indiv_units.append(unit_val**int(exp))
                break
        else:
            return
    
    overall_unit = product(indiv_units)
    if overall_unit == UNITS_BY_SYMBOL['']:
        overall_unit = dimensionless

    return overall_unit