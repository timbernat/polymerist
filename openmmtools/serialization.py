'''For reading and writing OpenMM components to files'''

from pathlib import Path

from openmm import System, Context, State
from openmm.app import Simulation, PDBFile
from openmm import XmlSerializer


# GENERIC PARAMETER SETS
DEFAULT_STATE_PARAMS : dict[str, bool] = {
    'getPositions'  : True,
    'getVelocities' : True,
    'getForces'     : False,
    'getEnergy'     : True,
    'getParameters' : True,
    'getParameterDerivatives' : False,
    'getIntegratorParameters' : False
}

# SERIALIZATION FUNCTIONS
def serialize_state_and_sys(sim : Simulation, out_dir : Path, out_name : str, state_params : dict[str, bool]=DEFAULT_STATE_PARAMS) -> None:
    '''For saving State and System info of a Simulation to disc'''
    sim_dict = {
        'system' : sim.system,
        'state' : sim.context.getState(**state_params)
    }
    
    for affix, save_data in sim_dict.items():
        save_path = out_dir / f'{out_name}_{affix}.xml'
        save_path.touch()

        with save_path.open('w') as file:
            file.write( XmlSerializer.serialize(save_data) )

def apply_state_to_sim(context : Context, state : State) -> None:
    '''For applying saved State data to an existing OpenMM Simulation'''
    context.setPeriodicBoxVectors(*state.getPeriodicBoxVectors())
    context.setPositions(state.getPositions())
    context.setVelocities(state.getVelocities())
    context.setTime(state.getTime())

    context.reinitialize(preserveState=True)

def save_sim_snapshot(sim : Simulation, pdb_path : Path, keep_ids : bool=True) -> None:
    '''Saves a PDB of the current state of a simulation's Topology'''
    curr_state = sim.context.getState(getPositions=True, getEnergy=True)
    with pdb_path.open('w') as output:
        PDBFile.writeFile(sim.topology, curr_state.getPositions(), output, keepIds=keep_ids)