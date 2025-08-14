'''Tools for labelling and extracting force groups for Forces within an OpenMM System'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

from typing import Container, Union
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

def uniquify_force_groups(system : System, except_for : Container[int]=None) -> None:
    '''Assigns each Force in an OpenMM System with a unique force group index'''
    if except_for is None:
        except_for = set()

    force_grp_idx : int = 0
    for force in system.getForces():
        if force.getForceGroup() in except_for:
            continue

        force.setForceGroup(force_grp_idx)
        force_grp_idx += 1
    # TODO : add labelling (depends partially on Interchange's NonbondedForce separation)

def impose_unique_force_groups(ommsys : System, except_for : Container[int]=None) -> None:
    '''Impose unique labels on Forces in an OpenMM System'''
    if not force_groups_are_unique(ommsys):
        uniquify_force_groups(ommsys, except_for=except_for)

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