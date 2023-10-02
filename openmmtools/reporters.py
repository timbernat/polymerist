'''Utilities for handlling setup of Simulation Reporters'''

import logging
LOGGER = logging.getLogger(__name__)

from typing import Any, Callable, Optional, TypeAlias, Union
from dataclasses import dataclass, field

from functools import partial
from pathlib import Path

from openmm.app import PDBReporter, PDBxReporter, DCDReporter
from openmm.app import StateDataReporter, CheckpointReporter

from .serialization import assemble_sim_file_path
from ..genutils.typetools import Args, KWArgs
from ..genutils.fileutils.jsonio.jsonify import make_jsonifiable


# REPORTER-SPECIFIC TYPES AND TYPEHINTS
class StateReporter(CheckpointReporter):
    '''CheckpointReporter which is forced to report States; greatly simplifies interface, documentation, and logging'''
    def __init__(self, file, reportInterval):
        super().__init__(file, reportInterval, writeState=True)

TrajectoryReporter  : TypeAlias = Union[PDBReporter, PDBxReporter, DCDReporter]
Reporter            : TypeAlias = Union[CheckpointReporter, StateReporter, StateDataReporter, TrajectoryReporter]
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
    'chk' : CheckpointReporter, # TODO : consider adding checks to guarantee that this can't ever be called with writeState=True
    'xml' : StateReporter
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

@make_jsonifiable
@dataclass
class ReporterParameters:
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

    def prepare_reporters(self, out_dir : Path, out_name : str, report_interval : int) -> tuple[dict[str, Path], list[Reporter]]:
        '''Prepare all reporters specified by internal parameters and return Paths linked to said reporters'''
        rep_paths, reps = {}, []

        for rep_config in self.rep_configs:
            rep_path = assemble_sim_file_path(out_dir, out_name, extension=rep_config.extension, affix=rep_config.label)
            rep_path.touch() # ensure the reporter output file actually exists

            rep = rep_config.initializer(str(rep_path), reportInterval=report_interval, **rep_config.extra_kwargs)
            LOGGER.info(f'Prepared {rep.__class__.__name__} which reports to {rep_path}')

            rep_paths[rep_config.label] = rep_path
            reps.append(rep)

        return rep_paths, reps