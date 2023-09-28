'''For recording, storing, and organizing parameters and files attached to a Simulation'''

from dataclasses import dataclass, field
from typing import Any

import re
import numpy as np

import openmm.unit
from openmm.unit import Quantity

from polysaccharide2.genutils.fileutils.jsonio import JSONSerializable, JSONifiable


def json_serialize_quantities(unser_jdict : dict[Any, Any]) -> dict[str, JSONSerializable]:
    '''Serialize unit-ful Quantity attrs in a way the JSON can digest'''
    ser_jdict = {}
    for attr_name, attr_val in unser_jdict.items():
        if not isinstance(attr_val, Quantity):
            ser_jdict[attr_name] = attr_val # nominally, copy all non-Quantities directly
        else:
            ser_jdict[f'{attr_name}_value'] = attr_val._value
            ser_jdict[f'{attr_name}_unit' ] = str(attr_val.unit)

    return ser_jdict

def json_unserialize_quantities(ser_jdict : dict[str, JSONSerializable]) -> dict[Any, Any]:
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
class IntegratorParameters(JSONifiable):
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
    def serialize_json_dict(unser_jdict : dict[Any, Any]) -> dict[str, JSONSerializable]:
        '''Serialize unit-ful Quantity attrs in a way the JSON can digest'''
        return json_serialize_quantities(unser_jdict)
    
    @staticmethod
    def unserialize_json_dict(ser_jdict : dict[str, JSONSerializable]) -> dict[Any, Any]:
        '''For unserializing unit-ful Quantities upon load from json file'''
        return json_unserialize_quantities(ser_jdict)
    