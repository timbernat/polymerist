'''For extracting properties from OpenMM Contexts (e.g. positions, energies, etc)'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

from typing import Optional

from openmm import Context, State
from openmm.unit import Unit, Quantity

from .forces import impose_unique_force_groups, forces_by_force_group


# POSITIONS
def get_context_positions(context : Context) -> Quantity:
    '''Extract coordinates from the current state of a simulation''' 
    # NOTE : forcing numpy output for now, as OpenMM Vec3's don't seem particularly useful; may change this in the future
    return context.getState(getPositions=True).getPositions(asNumpy=True) 

# ENERGIES
def get_openmm_energies(
        context : Context,
        preferred_unit : Optional[Unit]=None,
        force_group_names : Optional[dict[int, str]]=None,
    ) -> dict[str, Quantity]:
    '''
    Evaluate energies of an OpenMM Context
    
    Returns dict, keyed by contribution name, containing the total kinetic and potential energies,
    as well as the individual potential energies contributed by each distinct force group
    '''
    if force_group_names is None:
        force_group_names = {}

    energies : dict[str, Quantity] = {}
    ## global
    global_state : State = context.getState(getEnergy=True) # initialize shared global state
    energies['Total potential energy'] = global_state.getPotentialEnergy()
    energies['Total kinetic energy'  ] = global_state.getKineticEnergy()
    
    ## by contribution
    for group_id in forces_by_force_group(context.getSystem()):
        local_state : State = context.getState(getEnergy=True, groups={group_id})
        group_label = force_group_names.get(group_id, f'Group {group_id}')
        energies[f'{group_label} potential energy'] = local_state.getPotentialEnergy()

    if preferred_unit is not None:
        for key, energy in energies.items():
            energies[key] = energy.in_units_of(preferred_unit)

    return energies
eval_openmm_energies = get_openmm_energies

def get_openmm_energies_separated(context : Context, preferred_unit : Optional[Unit]=None) -> dict[str, Quantity]:
    '''
    Evaluate energies of an OpenMM Context
    Enforces separation of each Force's contribution to the energy (i.e. by imposing unique force groups)
    
    Returns dict, keyed by contribution name, containing the total kinetic and potential energies,
    as well as the individual potential energies contributed by each distinct force group
    '''
    ommsys = context.getSystem()
    impose_unique_force_groups(ommsys)
    
    return get_openmm_energies(
        context,
        preferred_unit=preferred_unit,
        force_group_names={
            group_id : force.getName()
                for group_id, force in forces_by_force_group(ommsys, denest=True).items()
        },
    )
eval_openmm_energies_separated = get_openmm_energies_separated