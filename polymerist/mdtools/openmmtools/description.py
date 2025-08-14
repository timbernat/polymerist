'''Utilities from describing the various parameters and settings possessed by OpenMM objects'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

from typing import Any, Union
from enum import Enum

from openmm import Force, CustomNonbondedForce, NonbondedForce, System

from ...genutils.attrs import compile_argfree_getable_attrs
from ...genutils.textual.prettyprint import dict_to_indented_str


# DEVNOTE: this affords a mechanism for identifying which Forces in an OpenMM system were added by this library
# While not infallible, this value was chosen to afford as many other forces as possible, assuming sequential numbering
# (see https://docs.openmm.org/latest/api-python/generated/openmm.openmm.Force.html#openmm.openmm.Force.setForceGroup)
# and should UNDER NO CIRCUMSTANCES be altered; the global constant here is for synchrony between internal functions
_POLYMERIST_FORCE_GROUP : int = 31

# REFERENCE FOR NONBONDED METHODS
## referenced from "Attributes" section of docs page: https://docs.openmm.org/latest/api-python/generated/openmm.openmm.NonbondedForce.html#nonbondedforce
## necessary to assign meaningful names to the otherwise-inaccessible int enum for these methods in the OpenMM API
_NONBOND_METHOD_NAMES : tuple[str] = (
    'NoCutoff',
    'CutoffNonPeriodic',
    'CutoffPeriodic',
    'Ewald',
    'PME',
    'LJPME',
)
_NONBOND_METHOD_VALUES : dict[str, int] = {}
for method_name in _NONBOND_METHOD_NAMES:
    int_value = getattr(NonbondedForce, method_name)
    # NOTE: order matters here; want capitalized version AFTER the correct-case version, 
    # so the resulting Enum aliases duplicate values to the correct-case version of the entry
    _NONBOND_METHOD_VALUES[method_name] = int_value
    _NONBOND_METHOD_VALUES[method_name.upper()] = int_value
LongRangeNonbondedMethod = Enum('LongRangeNonbondedMethod', _NONBOND_METHOD_VALUES)


# DESCRIBING PARAMETERS CONTAINED IN FORCES
def describe_force(force : Force) -> dict[str, Any]:
    '''Provides a dictionary which summarizes the parameters of a Force'''
    force_attrs = compile_argfree_getable_attrs(force, getter_re='\Aget', repl_str='') # getter string here asserts that "get" is at the start of the attribute name
    force_attrs['Type'] = type(force).__name__
    
    # NOTE: the overlap in long range method enum values is coincidental for NonbondedForce, CustomNonbondedForce
    # and is !NOT! the case in general (e.g. see attributes of https://docs.openmm.org/latest/api-python/generated/openmm.openmm.AmoebaMultipoleForce.html)
    if isinstance(force, (NonbondedForce, CustomNonbondedForce)):
        log_range_method = LongRangeNonbondedMethod(force.getNonbondedMethod())
        force_attrs['NonbondedMethod'] = log_range_method.name

    return force_attrs

def describe_forces(ommsys : System, as_str : bool=False) -> Union[str, dict[str, dict[str, Any]]]:
    '''Provides a dictionary (keyed by force names) which summarizes the parameters of each Force in an OpenMM system'''
    force_desc_dict = {
        force.getName(): describe_force(force)
            for force in ommsys.getForces()
    }
    if as_str:
        return dict_to_indented_str(force_desc_dict)
    return force_desc_dict