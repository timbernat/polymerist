'''Utilities for handlling setup of Simulation Reporters'''

from functools import partial
from dataclasses import dataclass, field
from typing import Any, Callable, Optional, TypeAlias, Union

from pathlib import Path

from openmm.app import PDBReporter, PDBxReporter, DCDReporter
from openmm.app import StateDataReporter, CheckpointReporter

from ..genutils.typetools import Args, KWArgs
from ..genutils.decorators.functional import allow_string_paths
from ..genutils.fileutils.jsonio import JSONifiable, JSONSerializable


# REPORTER-SPECIFIC TYPEHINTS
TrajectoryReporter  : TypeAlias = Union[PDBReporter, PDBxReporter, DCDReporter]
Reporter            : TypeAlias = Union[StateDataReporter, CheckpointReporter, TrajectoryReporter]
ReporterInitializer : TypeAlias = Callable[[Args, KWArgs], Reporter]


# REPORTER FILE EXTENSIONS AND INITIALIZERS
@dataclass(frozen=True)
class ReporterConfig:
    '''Interface for generating a Reporter and associated Path in a uniform and standardized way'''
    label     : str
    extension : str
    reporter_dict : dict[str, ReporterInitializer] # dictionary of intializers keyed by extension
    extra_kwargs  : dict[str, Any] = field(default_factory=dict)

    @property
    def initializer(self) -> ReporterInitializer:
        '''Returns the particular class or initializer'''
        return self.reporter_dict[self.extension] # explicitly called with getitem to raise KeyError if extension is invalid

TRAJ_REPORTERS = { # index output formats of reporters by file extension
    'pdb'  : PDBReporter,
    'pdbx' : PDBxReporter,
    'dcd'  : DCDReporter,
}

CHK_REPORTERS = {
    'chk' : partial(CheckpointReporter, writeState=False),
    'xml' : partial(CheckpointReporter, writeState=True ),
}

STATE_DAT_REPORTERS = {
    'csv' : StateDataReporter
}


# SERIALIZABLE REPORTER PARAMETER STORAGE AND GENERATION
DEFAULT_STATE_DATA_PROPS = {
    'step'            : True,
    'time'            : True,
    'potentialEnergy' : True,
    'kineticEnergy'   : True,
    'totalEnergy'     : True,
    'temperature'     : True,
    'volume'          : True,
    'density'         : True,
    'speed'           : True,
    'progress'        : False,
    'remainingTime'   : False,
    'elapsedTime'     : False
}

@allow_string_paths
def assemble_sim_file_path(out_dir : Path, out_name : str, extension : str, affix : str='') -> Path:
    '''Combine output, naming, descriptive, and filetype info to generate a complete Simulation-related file Path'''
    if extension[0] == '.':
        extension = extension[1:] # remove leading dots if included
    path_name = f'{out_name}{"_" if affix else ""}{affix}.{extension}'

    return out_dir / path_name

@dataclass
class ReporterParameters(JSONifiable):
    '''Parameters for specifying Simulation reporters'''
    report_checkpoint : bool = True
    report_state      : bool = True
    report_trajectory : bool = True
    report_state_data : bool = True

    traj_ext   : Optional[str] = 'dcd'
    num_steps  : Optional[int] = None 
    state_data : Optional[dict[str, bool]] = field(default_factory=lambda : DEFAULT_STATE_DATA_PROPS)

    @property
    def rep_configs(self) -> list[ReporterConfig]:
        '''Generate Reporter configurations for currently specified reporters'''
        rep_configs = []
        
        if self.report_trajectory:
            rep_configs.append(ReporterConfig(label='trajectory', extension=self.traj_ext, reporter_dict=TRAJ_REPORTERS))

        if self.report_checkpoint:
            rep_configs.append(ReporterConfig(label='checkpoint', extension='chk', reporter_dict=CHK_REPORTERS))

        if self.report_state:
            rep_configs.append(ReporterConfig(label='state', extension='xml', reporter_dict=CHK_REPORTERS))

        if self.report_state_data:
            rep_configs.append(ReporterConfig(label='state_data', extension='csv', reporter_dict=STATE_DAT_REPORTERS, extra_kwargs={**self.state_data, 'totalSteps' : self.num_steps}))
                    
        return rep_configs

    def prepare_reporters(self, out_dir : Path, out_name : str, record_freq : int) -> tuple[dict[str, Path], list[Reporter]]:
        '''Prepare all reporters specified by internal parameters and return Paths linked to said reporters'''
        rep_paths, reps = {}, []

        for rep_config in self.rep_configs:
            rep_path = assemble_sim_file_path(out_dir, out_name, extension=rep_config.extension, affix=rep_config.label) 
            rep = rep_config.initializer(str(rep_path), reportInterval=record_freq, **rep_config.extra_kwargs)

            rep_paths[rep_config.label] = rep_path
            reps.append(rep)

        return rep_paths, reps