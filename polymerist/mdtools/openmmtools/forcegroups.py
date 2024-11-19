'''Tools for labelling and extracting force groups for Forces within an OpenMM System'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

from typing import Union
from collections import defaultdict
from openmm import Force, NonbondedForce, System


# HANDLING FORCE GROUPS
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

def uniquify_force_groups(system : System) -> None:
    '''Assigns each Force in an OpenMM System with a unique force group index'''
    for i, force in enumerate(system.getForces()):
        force.setForceGroup(i)
    # TODO : add labelling (depends partially on Interchange's NonbondedForce separation)

def impose_unique_force_groups(ommsys : System) -> None:
    '''Impose unique labels on Forces in an OpenMM System'''
    if not force_groups_are_unique(ommsys):
        uniquify_force_groups(ommsys)

def forces_by_force_group(system : System, denest : bool=False) -> dict[int, Union[Force, list[Force]]]:
    '''Compile the Forces in an'''
    force_dict = defaultdict(list)
    for force in system.getForces():
        force_dict[force.getForceGroup()].append(force)
    force_dict = dict(force_dict) # convert to vanilla dict to avoid non-builtin covariant typing
    
    if denest:
        for group_id, force_list in force_dict.items():
            if len(force_list) == 1: # for single-force lists, extract item and make value directly
                force_dict[group_id] = force_list.pop() 
    
    return force_dict