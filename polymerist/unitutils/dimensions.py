'''For checking dimensionality and presence of units'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

from typing import Any, Union

from openmm.unit import Quantity, Unit, length_dimension
from pint import ( # this is also the base classes for all OpenFF-style units
    Unit as PintUnit,
    Quantity as PintQuantity, 
)


# CHECKING FOR AND REMOVING UNITS
class MissingUnitsError(Exception):
    pass

def hasunits(obj : Any) -> bool:
    '''Naive but effective way of checking for pint and openmm units'''
    return any(hasattr(obj, attr) for attr in ('unit', 'units')) 

def strip_units(coords : Union[tuple, PintQuantity, Quantity]) -> tuple[float]:
    '''
    Sanitize coordinate tuples for cases which require unitless quantities
    Specifically needed since OpenMM and pint each have their own Quantity and Units classes
    '''
    if isinstance(coords, PintQuantity):
        return coords.magnitude
    elif isinstance(coords, Quantity):
        return coords._value

    return coords

# CHECKING DIMENSIONALITY
def is_volume(unit_val : Union[Unit, Quantity]) -> bool:
    '''Return whether a unit corresponds to a volume'''
    if isinstance(unit_val, Quantity):
        unit_val = unit_val.unit # extract just the unit component if a Quantity is passed

    for i, (dim, exp) in enumerate(unit_val.iter_base_dimensions()):
        if i > 0:
            return False # immediate rule out if more than just one unit is present
        
    if (dim == length_dimension) and (exp == 3.0): # if monodimensional, check that the single dimension is L^3
        return True
    return False