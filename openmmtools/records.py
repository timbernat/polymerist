'''For recording, storing, and organizing parameters and files attached to a Simulation'''

import logging
LOGGER = logging.getLogger(__name__)

from dataclasses import dataclass, field
from typing import Any

import re
import numpy as np
from pathlib import Path

# from .. import OPENFF_DIR
from openforcefields.openforcefields import get_forcefield_dirs_paths
OPENFF_DIR = Path(get_forcefield_dirs_paths()[0])

import openmm.unit
from openmm.unit import Quantity
from openmm.unit import picosecond, femtosecond
from openmm.unit import atmosphere, kelvin

from ..genutils.fileutils.jsonio import JSONSerializable, JSONifiable


@dataclass
class SimulationParameters(JSONifiable):
    '''For recording the parameters used to run an OpenMM Simulation'''
    total_time  : Quantity
    num_samples : int

    ensemble : str
    periodic : bool = True
    forcefield_name : str = 'openff-2.0.0.offxml' # will be looked up in openforcefields resource module - uses Sage by default

    affix : str = '' # optional descriptive string for simulation
    binary_traj : bool = True  # whether to save trajectory as compact binary (.dcd) or human-readable (.pdb) format
    save_state  : bool = False # whether to save State or Checkpoint when updating simulation checkpoints
    reported_state_data : dict[str, bool] = field(default_factory=dict)

    timestep       : Quantity = field(default_factory=lambda : (2 * femtosecond)) # just specifying Quantities as default doesn't cut it, since these are (evidently) mutable defaults which dataclasses can't digest
    temperature    : Quantity = field(default_factory=lambda : (300 * kelvin))
    pressure       : Quantity = field(default_factory=lambda : (1 * atmosphere))
    friction_coeff : Quantity = field(default_factory=lambda : (1 / picosecond))
    barostat_freq : int = 25

    @property
    def num_steps(self) -> int:
        '''Total number of steps in the simulation'''
        return round(self.total_time / self.timestep)
    
    @property
    def record_freq(self) -> int:
        '''Number of steps between each taken sample'''
        return round(self.num_steps / self.num_samples)
    
    @property
    def record_interval(self) -> float:
        '''Length of time between successive samples'''
        return self.record_freq * self.timestep
    
    @property
    def time_points(self) -> np.ndarray[int]:
        '''An array of the time data points represented by the given sampling rate and sim time'''
        return (np.arange(0, self.num_steps, step=self.record_freq) + self.record_freq)* self.timestep # extra offset by recording frequency need to align indices (not 0-indexed)

    @property
    def forcefield_path(self) -> Path:
        '''Returns the path to the official OpenFF Forcefield named in the parameter set'''
        ff_path = OPENFF_DIR / self.forcefield_name
        assert(ff_path.exists()) # make sure the forcefield requested genuinely exists
        
        return ff_path

    # JSON serialization
    @staticmethod
    def serialize_json_dict(unser_jdict : dict[Any, Any]) -> dict[str, JSONSerializable]:
        '''Serialize unit-ful Quantity attrs in a way the JSON can digest'''
        ser_jdict = {}
        for attr_name, attr_val in unser_jdict.items():
            if not isinstance(attr_val, Quantity):
                ser_jdict[attr_name] = attr_val # nominally, copy all non-Quantities directly
            else:
                ser_jdict[f'{attr_name}_value'] = attr_val._value
                ser_jdict[f'{attr_name}_unit' ] = str(attr_val.unit)

        return ser_jdict
    
    @staticmethod
    def unserialize_json_dict(ser_jdict : dict[str, JSONSerializable]) -> dict[Any, Any]:
        '''For unserializing unit-ful Quantities upon load from json file'''
        unser_jdict = {}
        for attr_name, attr_val in ser_jdict.items():
            if attr_name.endswith('_unit'): # skip pure units
                continue 
            elif attr_name.endswith('_value'): # reconstitute Quantities from associated units
                quant_name = re.match('(.*)_value', attr_name).groups()[0] # can't use builtin str.strip, as it removes extra characters
                unit_name = ser_jdict[f'{quant_name}_unit']
                
                if unit_name.startswith('/'): # special case needed to handle inverse units
                    unit = 1 / getattr(openmm.unit, unit_name.strip('/'))
                else:
                    unit = getattr(openmm.unit, unit_name)

                unser_jdict[quant_name] = Quantity(attr_val, unit)
            else: # for non-Quantity entries, load as-is
                unser_jdict[attr_name] = attr_val
        
        return unser_jdict
    
@dataclass
class SimulationPaths(JSONifiable):
    '''Stores paths to various files associated with a completed MD simulation'''
    sim_params : Path
    trajectory : Path
    state_data : Path = None
    checkpoint : Path = None
    
    time_data    : Path = None
    spatial_data : Path = None

    # JSON serialization
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