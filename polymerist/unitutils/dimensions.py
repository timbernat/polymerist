'''For checking dimensionality and presence of units'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

from typing import Any, Union, TypeVar
T = TypeVar('T')

from numpy import ndarray
from openmm.unit import (
    Unit as OpenMMUnit, 
    Quantity as OpenMMQuantity,
    length_dimension,
)
OpenMMUnitLike = Union[OpenMMUnit, OpenMMQuantity] # TODO: add union type checkers

from pint import ( # this is also the base classes for all OpenFF-style units
    Unit as PintUnit,
    Quantity as PintQuantity, 
)
PintUnitLike = Union[PintUnit, PintQuantity] # TODO: add union type checkers

Unit     = Union[PintUnit    , OpenMMUnit]
Quantity = Union[PintQuantity, OpenMMQuantity]


# CHECKING FOR AND REMOVING UNITS
class MissingUnitsError(Exception):
    pass

def hasunits(obj : Any) -> bool:
    '''Naive but effective way of checking for pint and openmm units'''
    return any(hasattr(obj, attr) for attr in ('unit', 'units')) 

def strip_units(coords : Union[T, PintQuantity, OpenMMQuantity]) -> Union[T, ndarray[Any]]:
    '''
    Sanitize coordinate tuples for cases which require unitless quantities
    Specifically needed since OpenMM and pint each have their own Quantity and Units classes
    '''
    if isinstance(coords, PintQuantity):
        return coords.magnitude # for container-like values (e.g. tuples), will always return numpy array instead (not type-safe!)
    elif isinstance(coords, OpenMMQuantity):
        return coords._value

    return coords

# CHECKING DIMENSIONALITY
def _is_volume_openmm(unitlike : OpenMMUnitLike) -> bool:
    '''Check whether an OpenMM Unit/Quantity dimensionally corresponds to a volume'''
    if isinstance(unitlike, OpenMMQuantity):
        unitlike = unitlike.unit # extract just the unit component if a Quantity is passed

    for i, (dim, exp) in enumerate(unitlike.iter_base_dimensions()):
        if i > 0:
            return False # immediate rule out if more than just one dimension is present
        
    if (dim == length_dimension) and (exp == 3.0): # if monodimensional, check that the single dimension is L^3
        return True
    return False

def _is_volume_pint(unitlike : PintUnitLike) -> bool:
    '''Check whether an Pint Unit/Quantity dimensionally corresponds to a volume'''
    return unitlike.dimensionality == '[length]**3' # "dimensionality" attr is present on both the Unit and Quantity classes in Pint

def is_volume(unitlike : Union[Unit, Quantity]) -> bool:
    '''
    Check whether a Unit or Quantity dimensionally corresponds to a volume
    Accepts both OpenMM-style and Pint-style unit-like objects
    '''
    if isinstance(unitlike, OpenMMUnitLike):
        return _is_volume_openmm(unitlike)
    elif isinstance(unitlike, PintUnitLike):
        return _is_volume_pint(unitlike)
    else:
        # raise TypeError(f'Cannot interpret object of type "{type(unitlike).__name__}" as unit-like')
        return False # strictly speaking, anything which has no notion of units cannot be a volume