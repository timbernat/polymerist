'''For reading and writing OpenMM components to files'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

import logging
LOGGER = logging.getLogger(__name__)

from typing import Any, Optional, Union
from numpy.typing import NDArray
from dataclasses import dataclass

from pathlib import Path
from collections import Counter

from openmm import System, Context, State
from openmm import XmlSerializer, Vec3
from openmm.app import Simulation, PDBFile
from openmm.app import Topology as OpenMMTopology

from .parameters import SimulationParameters
from ...genutils.decorators.functional import allow_string_paths
from ...genutils.fileutils.pathutils import assemble_path
from ...genutils.fileutils.jsonio.jsonify import make_jsonifiable
from ...genutils.fileutils.jsonio.serialize import PathSerializer
from ...molfiles.pdb import SerialAtomLabeller


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


# SERIALIZATION AND DESERIALIZATION FUNCTIONS # TODO : add assertions for file extensions
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


@allow_string_paths
def serialize_system(sys_path : Path, system : System) -> None:
    '''For saving an existing OpenMM System to file'''
    with sys_path.open('w') as file:
        file.write(XmlSerializer.serialize(system))

@allow_string_paths
def serialize_openmm_pdb(
        pdb_path : Path,
        topology : OpenMMTopology,
        positions : Union[NDArray, list[Vec3]],
        keep_chain_and_res_ids : bool=True,
        atom_labeller : Optional[SerialAtomLabeller]=SerialAtomLabeller(),
        resname_map : Optional[dict[str, str]]=None,
    ) -> None:
    '''Configure and write an Protein DataBank File from an OpenMM Topology and array of positions
    Provides options to configure atom ID numbering, residue numbering, and residue naming'''
    if resname_map is None:
        resname_map = {} # avoids mutable default

    # chain config
    for chain in topology.chains():
        chain.id = str(chain.id)

    # residue config
    for residue in topology.residues():
        residue.id = str(residue.id) # avoids TypeError when specifying keepIds during PDB write
        repl_res_name = resname_map.get(residue.name, None) # lookup current residue name to see if a replacement is called for
        if repl_res_name is not None:
            residue.name = repl_res_name

    # individual atom config
    if atom_labeller: # implicitly, preserves extant atom names if a labeller is not given
        for atom in topology.atoms():
            atom.name = atom_labeller.get_atom_label(atom.element.symbol)

    # file write
    with pdb_path.open('w') as file:
        PDBFile.writeFile(topology, positions, file, keepIds=keep_chain_and_res_ids)

@allow_string_paths
def serialize_topology_from_simulation(pdb_path : Path, sim : Simulation, keep_ids : bool=False) -> None:
    '''Saves a PDB of the current state of a simulation's Topology'''
    serialize_openmm_pdb(
        pdb_path,
        topology=sim.topology,
        positions=sim.context.getState(getPositions=True).getPositions(),
        keep_chain_and_res_ids=keep_ids
    )