'''For extracting properties from OpenMM Contexts (e.g. positions, energies, etc)'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

from typing import Optional

from openmm import Context
from openmm.unit import Unit, Quantity

from .forcegroups import impose_unique_force_groups, forces_by_force_group


# POSITIONS
def get_context_positions(context : Context) -> Quantity:
    '''Extract coordinates from the current state of a simulation''' 
    # NOTE : forcing numpy output for now, as OpenMM Vec3's don't seem particularly useful; may change this in the future
    return context.getState(getPositions=True).getPositions(asNumpy=True) 

# ENERGIES
def get_openmm_energies(context : Context, preferred_unit : Optional[Unit]=None, force_group_names : Optional[dict[int, str]]=None) -> tuple[dict[str, Quantity], ...]:
    '''Evaluate energies of an OpenMM Context, both total and by individual force group
    Returns two dicts, the first of potential energies and the second by kinetic'''
    if force_group_names is None:
        force_group_names = {}

    energy_types = {
        'potential' : 'getPotentialEnergy',
        'kinetic'   : 'getKineticEnergy',
    }

    energies = []
    global_state = context.getState(getEnergy=True) # initialize shared global state
    for energy_label, energy_funct_attr in energy_types.items():
        states = {
            'Total' : global_state
        }
        for group_id in forces_by_force_group(context.getSystem()):
            group_label = force_group_names.get(group_id, f'Group {group_id}')
            states[group_label] = context.getState(getEnergy=True, groups={group_id})

        energy_dict : dict[str, Quantity] = {}
        for state_label, state in states.items():
            energy = getattr(state, energy_funct_attr)() # look up energy evaluation function dynamically, then evaluate (assumes no args need be passed) 
            if preferred_unit is not None:
                energy = energy.in_units_of(preferred_unit)
            energy_dict[f'{state_label} {energy_label} energy'] = energy
        energies.append(energy_dict)

    return tuple(energies)
eval_openmm_energies = get_openmm_energies

def get_openmm_energies_separated(context : Context, preferred_unit : Optional[Unit]=None) -> tuple[dict[str, Quantity], ...]:
    '''Evaluate energies of an OpenMM Context, both total and by individual force group
    Enforces separation of each Force's contribution to the energy (i.e. by imposing unique force groups)
    Returns two dicts, the first of potential energies and the second by kinetic'''
    ommsys = context.getSystem()
    impose_unique_force_groups(ommsys)
    force_group_names = {
        group_id : force.getName()
            for group_id, force in forces_by_force_group(ommsys, denest=True).items()
    }

    return get_openmm_energies(context, preferred_unit=preferred_unit, force_group_names=force_group_names)
eval_openmm_energies_separated = get_openmm_energies_separated