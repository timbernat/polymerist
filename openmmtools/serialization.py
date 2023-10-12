'''For reading and writing OpenMM components to files'''

import logging
LOGGER = logging.getLogger(__name__)

from typing import Any, Optional
from dataclasses import dataclass
from pathlib import Path

from openmm import System, Context, State
from openmm.app import Simulation, PDBFile
from openmm import XmlSerializer

from .parameters import SimulationParameters
from ..genutils.decorators.functional import allow_string_paths
from ..genutils.fileutils.pathutils import assemble_path
from ..genutils.fileutils.jsonio.jsonify import make_jsonifiable
from ..genutils.fileutils.jsonio.serialize import PathSerializer


# DEFINING AND STORING SIMULATION PATHS
@make_jsonifiable(type_serializer=PathSerializer)
@dataclass
class SimulationPaths:
    '''Encapsulates Paths to various files associated with an OpenMM Simulation'''
    parameters_path   : Optional[Path] = None
    paths_path        : Optional[Path] = None

    system_path       : Optional[Path] = None
    topology_path     : Optional[Path] = None
    state_path        : Optional[Path] = None
    checkpoint_path   : Optional[Path] = None
    trajectory_path   : Optional[Path] = None
    
    state_data_path   : Optional[Path] = None
    time_data_path    : Optional[Path] = None
    spatial_data_path : Optional[Path] = None

    @allow_string_paths
    def init_top_and_sys_paths(self, out_dir : Path, prefix : str, record : bool=True) -> tuple[Path, Path]:
        '''Initialize Topology and System output paths for a given directory'''
        topology_path = assemble_path(out_dir, prefix, extension='pdb', postfix='topology')
        system_path   = assemble_path(out_dir, prefix, extension='xml', postfix='system')
        state_path    = assemble_path(out_dir, prefix, extension='xml', postfix='state')

        if record:
            self.topology_path = topology_path
            self.system_path = system_path
            self.state_path = state_path

        return topology_path, system_path, state_path
    
    @classmethod
    def from_dir_and_parameters(cls, out_dir : Path, prefix : str, sim_params : SimulationParameters, touch : bool=True) -> 'SimulationPaths':
        '''Create file directory and initialize simulationPaths object from a set of SimulationParameters'''
        path_obj = cls() # create empty path instance

        path_obj.parameters_path = assemble_path(out_dir, prefix, extension='json', postfix='parameters')
        path_obj.init_top_and_sys_paths(out_dir, prefix, record=True)  # initialize simulation saving files

        sim_params.reporter_params.init_reporter_paths(out_dir, prefix) # initialize simulation reporter files 
        for rep_label, path in sim_params.reporter_params.reporter_paths.items():
            setattr(path_obj, f'{rep_label}_path', path) # register reporter output paths; NOTE :relies on ReporterConfigs having matching labels

        if touch: # ensure files exist 
            out_dir.mkdir(exist_ok=True) # ensure the directory exists before attempting output
            sim_params.to_file(path_obj.parameters_path) # record simulation parameters
        
            sim_paths_path = assemble_path(out_dir, prefix, extension='json', postfix='paths')
            path_obj.paths_path = sim_paths_path # record path the SimulationsPaths object is being stored at within the object itself
            path_obj.to_file(sim_paths_path)
            
            for val in path_obj.__dict__.values():
                if isinstance(val, Path):
                    val.touch()

        return path_obj


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

@allow_string_paths
def serialize_system(sys_path : Path, system : System) -> None:
    '''For saving an existing OpenMM System to file'''
    with sys_path.open('w') as file:
        file.write(XmlSerializer.serialize(system))

@allow_string_paths
def serialize_topology_from_simulation(pdb_path : Path, sim : Simulation, keep_ids : bool=False) -> None:
    '''Saves a PDB of the current state of a simulation's Topology'''
    curr_state = sim.context.getState(getPositions=True)
    with pdb_path.open('w') as output:
        PDBFile.writeFile(sim.topology, curr_state.getPositions(), output, keepIds=keep_ids) # TODO : generalize file output beyond just PDB

# TODO : add assertions for file extensions
@allow_string_paths
def serialize_state_from_context(state_path : Path, context : Context, state_params : dict[str, bool]=DEFAULT_STATE_PROPS) -> None:
    '''For saving State data within an existing OpenMM Context to file'''
    state = context.getState(**state_params)
    with state_path.open('w') as file:
        file.write(XmlSerializer.serialize(state))

def apply_state_to_context(context : Context, state : State) -> None: # TOSELF : this might be replaced with Context.getState() followed by context.reinitialize(preserveState=True)
    '''For applying saved State data to an existing OpenMM Context'''
    context.setPeriodicBoxVectors(*state.getPeriodicBoxVectors())
    context.setPositions(state.getPositions())
    context.setVelocities(state.getVelocities())
    context.setTime(state.getTime())

    context.reinitialize(preserveState=True)