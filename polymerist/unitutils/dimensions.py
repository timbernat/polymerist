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
    
# COMPARING QUANTITIES
def _quantities_approx_equal_openmm(quantity_expected : OpenMMQuantity, quantity_actual : OpenMMQuantity, rel_tol : float=1E-8) -> bool:
    '''Determine whether two OpenMM Quantities with compatible dimensions are equal to within some relative error'''
    # NOTE: OpenMM returns a dimensionless quantity as simply a float, bypassing the need for any further conversion
    # Also, the subtraction here raises TypeError when dimensions are incompatible
    rel_err = abs(quantity_expected - quantity_actual) / quantity_actual 
    assert isinstance(rel_err, float) # verify that ratio is in fact dimensionless just to be safe :)
    
    return rel_err < rel_tol

def _quantities_approx_equal_pint(quantity_expected : PintQuantity, quantity_actual : PintQuantity, rel_tol : float=1E-8) -> bool:
    '''Determine whether two Pint Quantities with compatible dimensions are equal to within some relative error'''
    # NOTE: Pint, on the other hand, returns dimensionless quantity with explicit "dimensionless" attached
    rel_err = abs(quantity_expected - quantity_actual) / quantity_actual
    assert rel_err.unitless 
    
    return rel_err.magnitude < rel_tol # compare the magnitude of the dimensionless quantity to the relative tolerance
    
def quantities_approx_equal(
    quantity_expected : Quantity,
    quantity_actual   : Quantity,
    rel_tol : float=1E-8
) -> bool:
    '''
    Check whether two Quantity objects with compatible dimensions
    are equal to within some set relative error (default 1E-8)
    
    Accepts both OpenMM-style and Pint-style quantity-like objects
    '''
    if not all([isinstance(quantity_expected, Quantity), isinstance(quantity_actual, Quantity)]):
        raise TypeError(f'Comparison only valid between two Quantity instances, not objects of type "{type(quantity_expected)}" and "{type(quantity_actual)}"')
    if type(quantity_expected) != type(quantity_actual):
        raise TypeError(f'Comparison between mixed Quantity types not currently supported: quantities must either both be {OpenMMQuantity.__name__} or both be {PintQuantity.__name__} instances')
    
    # NOTE: the below comparisons working relies on the above checks ensuring both quantity arguments have the same Quantity-like type
    if isinstance(quantity_expected, OpenMMQuantity):
        return _quantities_approx_equal_openmm(quantity_expected, quantity_actual, rel_tol)
    elif isinstance(quantity_expected, PintQuantity):
        return _quantities_approx_equal_pint(quantity_expected, quantity_actual, rel_tol)