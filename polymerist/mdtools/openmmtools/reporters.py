'''Utilities for handlling setup of Simulation Reporters'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

import logging
LOGGER = logging.getLogger(__name__)

from typing import Any, Callable, Optional, TypeAlias, Union
from dataclasses import dataclass, field
from pathlib import Path

from openmm.app import PDBReporter, PDBxReporter, DCDReporter
from openmm.app import StateDataReporter, CheckpointReporter

from ...genutils.typetools.parametric import Args, KWArgs
from ...genutils.fileutils.pathutils import assemble_path
from ...genutils.fileutils.jsonio.jsonify import make_jsonifiable
from ...genutils.fileutils.jsonio.serialize import PathSerializer


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

@make_jsonifiable(type_serializer=PathSerializer)
@dataclass
class ReporterParameters:
    '''Parameters for specifying Simulation reporters'''
    report_checkpoint : bool = True
    report_state      : bool = True
    report_trajectory : bool = True
    report_state_data : bool = True

    traj_ext   : Optional[str] = 'dcd'
    num_steps  : Optional[int] = None 
    state_data     : Optional[dict[str, bool]] = field(default_factory=lambda : DEFAULT_STATE_DATA_PROPS)
    reporter_paths : Optional[dict[str, Path]] = field(default=None)

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

    def init_reporter_paths(self, out_dir : Path, prefix : str) -> None: # TODO : add mechanism to automatically re-initialize when bool params change (might onvolve caching dir info)
        '''Generate paths for all configured and enabled reporters'''
        self.reporter_paths = {
            rep_config.label : assemble_path(out_dir, prefix, extension=rep_config.extension, postfix=rep_config.label)
                for rep_config in self.rep_configs
        }

    def prepare_reporters(self, report_interval : int) -> list[Reporter]:
        '''Prepare all reporters specified by internal parameters and return Paths linked to said reporters'''
        if self.reporter_paths is None:
            raise ValueError('Cannot create Simulation reporters without initializing reporter output paths')
        
        reporters = []
        for rep_config in self.rep_configs:
            rep_path = self.reporter_paths[rep_config.label] # will raise KeyErro if config changes between Path init
            rep_path.touch() # ensure the reporter output file actually exists

            rep = rep_config.initializer(str(rep_path), reportInterval=report_interval, **rep_config.extra_kwargs)
            LOGGER.info(f'Prepared {rep.__class__.__name__} which reports to {rep_path}')
            reporters.append(rep)

        return reporters