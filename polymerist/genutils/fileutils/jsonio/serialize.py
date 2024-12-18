'''Interfaces for extending what types of objects can be serialized to JSON'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

from typing import Any, ClassVar, Optional, Type, TypeVar, Union
from abc import ABC, abstractstaticmethod
from inspect import isclass
T = TypeVar('T') # generic type

from pathlib import Path
import numpy as np
import openmm.unit

from ...decorators.classmod import register_subclasses, register_abstract_class_attrs


# CHECKING IF AN OBJECT IS SERIALIZABLE TO JSON BY DEFAULT
JSONSerializable = Union[str, bool, int, float, tuple, list, dict] 

# ABSTRACT INTERFACE FOR DEFINING CUSTOM SERIALIZERS (ENCODER + DECODER)
@register_subclasses(key_attr='python_type')
@register_abstract_class_attrs('python_type')
class TypeSerializer(ABC):
    '''Interface for defining how types which are not JSON serializable by default should be encoded and decoded'''
    python_type : ClassVar[Type[T]] # NOTE: this is kept here purely for static typehinting purposes

    @abstractstaticmethod
    def encode(python_obj : T) -> JSONSerializable:
        pass

    @abstractstaticmethod
    def decode(json_obj : JSONSerializable) -> T:
        pass

    @classmethod
    def encoder_default(cls, python_obj : Any) -> JSONSerializable: # NOTE : this is only called on objects which cannot be JSON serialized by default (i.e. don't need base case)
        '''Augmented Encoder for encoding registered objects along with type info for decoding'''
        if isinstance(python_obj, cls.python_type):
            return {
                '__class__' : cls.python_type.__name__,
                # '__class__' : python_obj.__class__.__name__, # supports subclasses of the base Python type
                '__values__' : cls.encode(python_obj),
            }
        else:
            raise TypeError(f'Object of type {python_obj.__class__.__name__} is not JSON serializable')

    @classmethod
    def decoder_hook(cls, json_dict : dict[JSONSerializable, JSONSerializable]) -> Union[dict, T]:
        type_name : Optional[str] = json_dict.get('__class__', None) # remove __class__ attr if present, returning None if not
        if type_name is None:
            return json_dict # return unmodified dict for untyped entries
        
        # raise Exception when attempting to decode typed entry of the wrong type
        if type_name != cls.python_type.__name__:
            raise TypeError(f'{cls.python_type.__name__} decoder cannot decode JSON-serialized object of type {type_name}')

        return cls.decode(json_dict['__values__']) # extract and decode values if entry has the right type

class MultiTypeSerializer:
    '''For dynamically merging multiple TypeSerializer encoders and decoders'''
    def __init__(self, *type_sers : tuple[Type[TypeSerializer]]) -> None:
        self._type_sers = []
        for ts in type_sers:
            self.add_type_serializer(ts)

    @property
    def type_sers(self) -> list[Type[TypeSerializer]]:
        '''Read-only wrapper for the internal registry of TypeSerializers'''
        return self._type_sers
    
    def add_type_serializer(self, obj : Union[TypeSerializer, 'MultiTypeSerializer']) -> None: # TODO: add uniqueness check
        '''For type, instance, and uniqueness checking of Type '''
        if isinstance(obj, TypeSerializer):
            self._type_sers.append(obj)
        elif isclass(obj) and issubclass(obj, TypeSerializer):
            self._type_sers.append(obj()) # instantiate, then add to registered serializers
        elif isinstance(obj, MultiTypeSerializer):
            for ts in obj.type_sers:
                self.add_type_serializer(ts)
        else:
            raise TypeError(f'Object of type "{obj.__name__ if isclass(obj) else type(obj).__name__}" is not a valid TypeSerializer')

    def encoder_default(self, python_obj : Any) -> JSONSerializable:
        for type_ser in self.type_sers:
            try:
                return type_ser.encoder_default(python_obj)
            except:
                pass # keep trying rest of encoders (don't immediately raise error) - TODO : make this less redundant-looking?
        else:
            raise TypeError(f'Object of type {python_obj.__class__.__name__} is not JSON serializable')

    def decoder_hook(self, json_dict : dict[JSONSerializable, JSONSerializable]) -> Any:
        for type_ser in self.type_sers:
            try:
                return type_ser.decoder_hook(json_dict)
            except:
                pass # keep trying rest of encoders (don't immediately raise error) - TODO : make this less redundant-looking?
        else: # only raised if no return occurs in any iteration - each decoder works for default-decodable values
            raise TypeError(f'No registered decoders for dict : {json_dict}')


# CONCRETE IMPLEMENTATIONS
class PathSerializer(TypeSerializer, python_type=Path):
    '''For JSON-serializing OpenMM Quantities'''
    @staticmethod
    def encode(python_obj : Path) -> str:
        '''Separate openmm.unit.Quantity's value and units to serialize as a single dict'''
        return str(python_obj)

    @staticmethod
    def decode(json_obj : str) -> Path:
        '''Unpack a value-unit string dict back into a usable openmm.unit.Quantity'''
        return Path(json_obj)
        
class QuantitySerializer(TypeSerializer, python_type=openmm.unit.Quantity):
    '''For JSON-serializing OpenMM Quantities'''
    @staticmethod
    def encode(python_obj : openmm.unit.Quantity) -> dict[str, Union[str, float]]:
        '''Separate openmm.unit.Quantity's value and units to serialize as a single dict'''
        value = python_obj._value
        if isinstance(value, np.ndarray): # supports numpy array serialization
            value = value.tolist()

        return {
            'value' : value,
            'unit'  : str(python_obj.unit),
        }

    @staticmethod
    def decode(json_obj : dict[str, Union[str, float]]) -> openmm.unit.Quantity:
        '''Unpack a value-unit string dict back into a usable openmm.unit.Quantity'''
        unit_name = json_obj['unit']
        if unit_name.startswith('/'): # special case needed to handle inverse units
            unit = 1 / getattr(openmm.unit, unit_name.strip('/'))
        else:
            unit = getattr(openmm.unit, unit_name)

        value = json_obj['value']
        if isinstance(value, list):
            value = np.array(value) # de-serialize numpy arrays; TOSELF : is there ever a case where this should remain a list (technically a valid Quanitity value)

        return openmm.unit.Quantity(value, unit)