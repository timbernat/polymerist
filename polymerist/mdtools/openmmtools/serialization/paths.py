'''For managing directories and files associated with OpenMM Simulations'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

from typing import Optional
from dataclasses import dataclass

from pathlib import Path

from ..parameters import SimulationParameters
from ....genutils.fileutils.pathutils import assemble_path, allow_string_paths
from ....genutils.fileutils.jsonio.jsonify import make_jsonifiable
from ....genutils.fileutils.jsonio.serialize import PathSerializer


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
    def init_top_and_sys_paths(
            self,
            out_dir : Path,
            prefix : str,
            record : bool=True,
        ) -> tuple[Path, Path]:
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
    def from_dir_and_parameters(
            cls,
            out_dir : Path,
            prefix : str,
            sim_params : SimulationParameters,
            touch : bool=True,
        ) -> 'SimulationPaths':
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