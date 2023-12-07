'''Tools for unit checking'''

from typing import Any, Callable, Union, TypeVar
R = TypeVar('R')   # for representing generic return values
Q = TypeVar('Q')   # for representing generic Quantity-like objects

from pint import Quantity as PintQuantity # this is also the base class for all OpenFF-style units
from openmm.unit import Quantity, Unit, length_dimension
from openff.units.openmm import (
    from_openmm as openmm_to_openff,
    to_openmm as openff_to_openmm,
)


# CHECKING FOR AN REMOVING UNITS
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


# DECORATORS FOR INTEROP BETWEEN OpenMM AND OpenFF UNIT SYSTEM
def allow_openmm_units(funct : Callable[[Q], R]) -> Callable[[Q], R]:
    '''Allow a Callable which expects ALL of its args to be OpenFF Quanitities to also accept equivalent OpenMM Quanitites'''
    def wrapper(*args, **kwargs) -> R:
        new_args = [
            openmm_to_openff(arg) if isinstance(arg, Quantity) else arg
                for arg in args
        ]

        new_kwargs = {
            key : openmm_to_openff(kwarg) if isinstance(kwarg, Quantity) else kwarg
                for key, kwarg in kwargs.items()
        }

        return funct(*new_args, **new_kwargs)
    return wrapper

def allow_openff_units(funct : Callable[[Q], R]) -> Callable[[Q], R]:
    '''Allow a Callable which expects ALL of its args to be OpenMM Quanitities to also accept equivalent OpenFF Quanitites'''
    def wrapper(*args, **kwargs) -> R:
        new_args = [
            openff_to_openmm(arg) if isinstance(arg, PintQuantity) else arg
                for arg in args
        ]

        new_kwargs = {
            key : openff_to_openmm(kwarg) if isinstance(kwarg, PintQuantity) else kwarg
                for key, kwarg in kwargs.items()
        }

        return funct(*new_args, **new_kwargs)
    return wrapper


# CHECKING DIMENSIONALITY
@allow_openff_units
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