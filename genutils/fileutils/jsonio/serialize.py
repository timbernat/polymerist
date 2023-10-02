'''Interfaces for extending what types of objects can be serialized to JSON'''

from typing import Any, ClassVar, Optional, Type, TypeVar, Union
from abc import ABC, abstractstaticmethod
T = TypeVar('T') # generic type

from pathlib import Path
import openmm.unit

from ...decorators.classmod import register_subclasses


# CHECKING IF AN OBJECT IS SERIALIZABLE TO JSON BY DEFAULT
JSONSerializable = Union[str, bool, int, float, tuple, list, dict] 

# ABSTRACT INTERFACE FOR DEFINING CUSTOM SERIALIZERS (ENCODER + DECODER)

@register_subclasses(key_attr='python_type')
class TypeSerializer(ABC):
    '''Interface for defining how types which are not JSON serializable by default should be encoded and decoded'''
    python_type : ClassVar[Type[T]]

    @abstractstaticmethod
    def encode(python_obj : T) -> JSONSerializable:
        pass

    @abstractstaticmethod
    def decode(json_obj : JSONSerializable) -> T:
        pass

    @classmethod
    def encoder_default(cls, python_obj : Any) -> JSONSerializable: # NOTE : this is only called on objects which cannot be JSON serialized by default (i.e. don;t need base case)
        '''Augmented Encoder for encoding registered objects along with type info for decoding'''
        if isinstance(python_obj, cls.python_type):
            return {
                '__class__' : cls.python_type.__name__,
                '__values__' : cls.encode(python_obj),
            }
        else:
            raise TypeError(f'Object of type {python_obj.__class__.__name__} is not JSON serializable')

    @classmethod
    def decoder_hook(cls, json_dict : dict[JSONSerializable, JSONSerializable]) -> Union[dict, T]:
        type_name : Optional[str] = json_dict.pop('__class__', None) # remove __class__ attr if present, returning None if not
        if type_name is None:
            return json_dict
        
        if type_name == cls.python_type.__name__:
            return cls.decode(json_dict['__values__']) # extract and decode values

class MultiTypeSerializer:
    '''For dynamically merging multiple TypeSerializer encoders and decoders'''
    def __init__(self, *type_sers) -> None:
        self.type_sers = type_sers

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
            return type_ser.decoder_hook(json_dict)
        else: # only raised if no return occurs in any iteration - each decoder works for default-decodable values
            raise TypeError(f'No registered decoders for dict : {json_dict}')


# CONCRETE IMPLEMENTATIONS
class PathSerializer(TypeSerializer):
    '''For JSON-serializing OpenMM Quantities'''
    python_type = Path

    @staticmethod
    def encode(python_obj : Path) -> str:
        '''Separate openmm.unit.Quantity's value and units to serialize as a single dict'''
        return str(python_obj)

    @staticmethod
    def decode(json_obj : str) -> Path:
        '''Unpack a value-unit string dict back into a usable openmm.unit.Quantity'''
        return Path(json_obj)
        
class QuantitySerializer(TypeSerializer):
    '''For JSON-serializing OpenMM Quantities'''
    python_type = openmm.unit.Quantity

    @staticmethod
    def encode(python_obj : openmm.unit.Quantity) -> dict[str, Union[str, float]]:
        '''Separate openmm.unit.Quantity's value and units to serialize as a single dict'''
        return {
            'value' : python_obj._value,
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

        return openmm.unit.Quantity(json_obj['value'], unit)