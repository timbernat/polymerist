'''For recording, storing, and organizing parameters associated wtih a Simulation'''

from dataclasses import dataclass, field
from typing import Any

import numpy as np
from openmm.unit import Quantity

from .thermo import ThermoParameters
from .reporters import ReporterParameters
from ..genutils.fileutils import jsonio


# PARAMETER SET HANDLERS
@dataclass
class IntegratorParameters(jsonio.JSONifiable):
    '''For recording the parameters used to run an OpenMM Simulation'''
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

    # JSON serialization
    @staticmethod
    def serialize_json_dict(unser_jdict : dict[Any, Any]) -> dict[str, jsonio.JSONSerializable]:
        '''Serialize unit-ful Quantity attrs in a way the JSON can digest'''
        return jsonio.json_serialize_quantities(unser_jdict)
    
    @staticmethod
    def unserialize_json_dict(ser_jdict : dict[str, jsonio.JSONSerializable]) -> dict[Any, Any]:
        '''For unserializing unit-ful Quantities upon load from json file'''
        return jsonio.json_unserialize_quantities(ser_jdict)
    

@dataclass
class SimulationParameters(jsonio.JSONifiable):
    '''Unified class for storing simulation parameters'''
    integ_params : IntegratorParameters
    thermo_params : ThermoParameters
    reporter_params : ReporterParameters

    @staticmethod
    def serialize_json_dict(unser_jdict : dict[Any, Any]) -> dict[str, jsonio.JSONSerializable]:
        '''Convert all Paths to strings'''
        return {
            field_name : params.serialize_json_dict(params.__dict__)
                for field_name, params in unser_jdict.items()
        }
    
    @staticmethod
    def unserialize_json_dict(ser_jdict : dict[str, jsonio.JSONSerializable]) -> dict[Any, Any]:
        '''For de-serializing JSON-compatible data into a form that the __init__method can accept'''
        unser_jdict = {}
        for key, value in ser_jdict.items():
            subparam_field = SimulationParameters.__dataclass_fields__.get(key) 
            if subparam_field is not None:
                subparam_dict = subparam_field.type.unserialize_json_dict(value)
                unser_jdict[key] = subparam_field.type(**subparam_dict)# Deserialize sub-parameter set at top-level
            else:
                unser_jdict[key] = value # leaves sub-parameter values untouched

        return unser_jdict
    