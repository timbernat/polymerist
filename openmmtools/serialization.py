'''For reading and writing OpenMM components to files'''

import logging
LOGGER = logging.getLogger(__name__)

from typing import Any, Optional
from dataclasses import dataclass
from pathlib import Path

from openmm import System, Context, State
from openmm.app import Simulation, PDBFile
from openmm import XmlSerializer

from ..genutils.decorators.functional import allow_string_paths
from ..genutils.fileutils.jsonio import JSONifiable, JSONSerializable


# DEFINING AND STORING SIMULATION PATHS
@allow_string_paths
def assemble_sim_file_path(out_dir : Path, out_name : str, extension : str, affix : str='') -> Path:
    '''Combine output, naming, descriptive, and filetype info to generate a complete Simulation-related file Path'''
    if extension[0] == '.':
        extension = extension[1:] # remove leading dots if included
    path_name = f'{out_name}{"_" if affix else ""}{affix}.{extension}'

    return out_dir / path_name

@dataclass
class SimulationPaths(JSONifiable):
    '''Encapsulates Paths to various files associated with an OpenMM Simulation'''
    sim_params : Optional[Path] = None

    topology   : Optional[Path] = None
    system     : Optional[Path] = None
    state      : Optional[Path] = None

    checkpoint : Optional[Path] = None
    trajectory : Optional[Path] = None

    state_data   : Optional[Path] = None
    time_data    : Optional[Path] = None
    spatial_data : Optional[Path] = None

    @staticmethod
    def serialize_json_dict(unser_jdict : dict[Any, Any]) -> dict[str, JSONSerializable]:
        '''Convert all Paths to strings'''
        ser_jdict = {}
        for key, value in unser_jdict.items():
            if isinstance(value, Path):
                ser_jdict[key] = str(value)
            else:
                ser_jdict[key] = value

        return ser_jdict
    
    @staticmethod
    def unserialize_json_dict(ser_jdict : dict[str, JSONSerializable]) -> dict[Any, Any]:
        '''For de-serializing JSON-compatible data into a form that the __init__method can accept'''
        unser_jdict = {}
        for key, value in ser_jdict.items():
            if value is not None:
                unser_jdict[key] = Path(value)
            else:
                unser_jdict[key] = value

        return unser_jdict


# SERIALIZATION AND DESERIALIZATION FUNCTIONS
DEFAULT_STATE_PROPS : dict[str, bool] = {
    'getPositions'  : True,
    'getVelocities' : True,
    'getForces'     : True,
    'getEnergy'     : True,
    'getParameters' : True,
    'getParameterDerivatives' : False,
    'getIntegratorParameters' : False
}

def serialize_system(system : System, out_dir : Path, out_name : str) -> Path:
    '''For saving an existing OpenMM System to file'''
    sys_path = assemble_sim_file_path(out_dir, out_name, extension='xml', affix='state')
    sys_path.touch()

    with sys_path.open('w') as file:
        file.write(XmlSerializer.serialize(system))

    return sys_path

def serialize_context_state(context : Context, out_dir : Path, out_name : str, state_params : dict[str, bool]=DEFAULT_STATE_PROPS) -> Path:
    '''For saving State data within an existing OpenMM Context to file'''
    state = context.getState(**state_params)
    state_path = assemble_sim_file_path(out_dir, out_name, extension='xml', affix='state')
    state_path.touch()

    with state_path.open('w') as file:
        file.write(XmlSerializer.serialize(state))

    return state_path

def apply_state_to_context(context : Context, state : State) -> None: # TOSELF : this might be replaced with Context.getState() followed by context.reinitialize(preserveState=True)
    '''For applying saved State data to an existing OpenMM Context'''
    context.setPeriodicBoxVectors(*state.getPeriodicBoxVectors())
    context.setPositions(state.getPositions())
    context.setVelocities(state.getVelocities())
    context.setTime(state.getTime())

    context.reinitialize(preserveState=True)

def save_sim_snapshot(sim : Simulation, pdb_path : Path, keep_ids : bool=True) -> None:
    '''Saves a PDB of the current state of a simulation's Topology'''
    curr_state = sim.context.getState(getPositions=True, getEnergy=True)
    with pdb_path.open('w') as output:
        PDBFile.writeFile(sim.topology, curr_state.getPositions(), output, keepIds=keep_ids) # TODO : generalize file output beyond just PDB