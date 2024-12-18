'''For recording, storing, and organizing parameters associated wtih a Simulation'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

from dataclasses import dataclass, field
from typing import Any

import numpy as np
from openmm.unit import Quantity

from .thermo import ThermoParameters
from .reporters import ReporterParameters
from ...genutils.fileutils.jsonio.jsonify import make_jsonifiable, dataclass_serializer_factory
from ...genutils.fileutils.jsonio.serialize import PathSerializer, QuantitySerializer, MultiTypeSerializer


# PARAMETER SET HANDLERS
@make_jsonifiable(type_serializer=QuantitySerializer)
@dataclass
class IntegratorParameters:
    '''For recording total, time step, recoridng frequency, and other integration time parameters'''
    time_step   : Quantity
    total_time  : Quantity
    num_samples : int

    @property
    def num_steps(self) -> int:
        '''Total number of steps in the simulation'''
        return round(self.total_time / self.time_step)
    
    @property
    def report_interval(self) -> int:
        '''Number of steps between successive samples'''
        return round(self.num_steps / self.num_samples)
    
    @property
    def report_duration(self) -> float:
        '''Length of time between successive samples'''
        return self.report_interval * self.time_step
    
    @property
    def time_points(self) -> np.ndarray[int]:
        '''An array of the time data points represented by the given sampling rate and sim time'''
        return (np.arange(0, self.num_steps, step=self.report_interval) + self.report_interval)* self.time_step # extra offset by recording frequency need to align indices (not 0-indexed)

# UNIFIED SIMULATION PARAMETER SETS
@make_jsonifiable
@dataclass
class SimulationParameters:
    '''Unified class for storing simulation parameters'''
    integ_params : IntegratorParameters
    thermo_params : ThermoParameters
    reporter_params : ReporterParameters
    
    def __post_init__(self) -> None:
        self.reporter_params.num_steps = self.integ_params.num_steps