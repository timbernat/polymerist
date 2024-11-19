'''Utilities from describing the various parameters and settings possessed by OpenMM objects'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

from typing import Any, Union
from openmm import NonbondedForce, System

from ...genutils.attrs import compile_argfree_getable_attrs
from ...genutils.textual.prettyprint import dict_to_indented_str


# REFERENCE TABLES FOR NONBONDED METHODS
NONBOND_METHOD_KEY : str = 'NonbondedMethod'
NONBOND_CUTOFF_METHOD_NAMES = (
    'NoCutoff',
    'CutoffNonPeriodic',
    'CutoffPeriodic',
    'Ewald',
    'PME',
    'LJPME',
)
NONBOND_CUTOFF_METHODS = {
    idx : method_name
        for idx, method_name in sorted( # sort in ascending order by integer code
            (getattr(NonbondedForce, method_name), method_name)
                for method_name in NONBOND_CUTOFF_METHOD_NAMES
        )
}

# DESCRIBING OpenMM PARAMETERS
def describe_forces(ommsys : System, as_str : bool=False) -> Union[str, dict[str, dict[str, Any]]]:
    '''Provides a dictionary (keyed by force names) which summarizes the parameters of each Force in an OpenMM system'''
    force_desc_dict = {}
    for force in ommsys.getForces():
        force_attrs = compile_argfree_getable_attrs(force, getter_re='\Aget', repl_str='') # getter string here asserts that "get" is at the start of the attribute name
        force_attrs['Type'] = type(force).__name__
        
        if (nonbond_id := force_attrs.get(NONBOND_METHOD_KEY)) is not None:
            force_attrs[nonbond_id] = NONBOND_CUTOFF_METHODS.get(nonbond_id, nonbond_id) # attempt a lookup of the nonbonded method name
        force_desc_dict[force.getName()] = force_attrs

    if as_str:
        return dict_to_indented_str(force_desc_dict)
    return force_desc_dict