'''For describing, labelling, and extracting information from OpenMM Forces'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

from typing import Any, Container, Union

from enum import Enum
from collections import defaultdict

from openmm import Force, CustomNonbondedForce, NonbondedForce, System

from ...genutils.attrs import compile_argfree_getable_attrs
from ...genutils.textual.prettyprint import dict_to_indented_str


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
        force_attrs['LongRangeNonbondedMethod'] = log_range_method.name

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


# ASSIGNING AND READING FORCE GROUPS
_POLYMERIST_FORCE_GROUP : int = 31
# DEVNOTE: this affords a mechanism for identifying which Forces in an OpenMM system were added by this library
# While not infallible, this value was chosen to afford as many other forces as possible, assuming sequential numbering
# (see https://docs.openmm.org/latest/api-python/generated/openmm.openmm.Force.html#openmm.openmm.Force.setForceGroup)
# and should UNDER NO CIRCUMSTANCES be altered; the global constant here is for synchrony between internal functions

def force_added_by_polymerist(force : Force) -> bool:
    '''Determines if a Force was added by this library'''
    return force.getForceGroup() == _POLYMERIST_FORCE_GROUP

def force_groups_are_unique(system : System) -> bool:
    '''Check whether the Forces in an OpenMM system share their force group ID with any other forces'''
    seen = set()
    for force in system.getForces():
        force_group = force.getForceGroup()
        if force_group in seen:
            return False
        seen.add(force_group)
    else:
        return True

def uniquify_force_groups(system : System, except_for : Container[int]=None) -> None:
    # DEVNOTE: docstring assigned below - dynamic docstrings cannot be assigned in-place
    if except_for is None:
        except_for = {_POLYMERIST_FORCE_GROUP} # don't squash polymerist-assigned forceGroups unless SPECIFICALLY requested

    force_grp_idx : int = 0
    for force in system.getForces():
        if force.getForceGroup() in except_for:
            continue

        force.setForceGroup(force_grp_idx)
        force_grp_idx += 1
    # TODO : add labelling (depends partially on Interchange's NonbondedForce separation)
uniquify_force_groups.__doc__ = f'''
Assigns each Force in an OpenMM System with a unique force group index

Parameters
----------
system : System
    The OpenMM System whose Forces are to be labelled
except_for : Container[int], optional, default={{{_POLYMERIST_FORCE_GROUP}}}
    An optional set of forceGroup which should NOT be relabelled,
    EVEN IF not relabelling then would prevent forces from being uniquified

    By default, only contains {_POLYMERIST_FORCE_GROUP} (_POLYMERIST_FORCE_GROUP), to allow tracking provenance 
    of creation of and modifications to forces done by this library
'''

def impose_unique_force_groups(ommsys : System, except_for : Container[int]=None) -> None:
    '''Impose unique labels on Forces in an OpenMM System'''
    if not force_groups_are_unique(ommsys):
        uniquify_force_groups(ommsys, except_for=except_for)

def forces_by_force_group(system : System, denest : bool=False) -> dict[int, Union[Force, list[Force]]]:
    '''Compile the Forces in an OpenMM system according to their assigned force group labels'''
    force_dict = defaultdict(list)
    for force in system.getForces():
        force_dict[force.getForceGroup()].append(force)
    force_dict = dict(force_dict) # convert to vanilla dict to avoid non-builtin covariant typing
    
    if denest:
        for group_id, force_list in force_dict.items():
            if len(force_list) == 1: # for single-force lists, extract item and make value directly
                force_dict[group_id] = force_list.pop() 
    
    return force_dict