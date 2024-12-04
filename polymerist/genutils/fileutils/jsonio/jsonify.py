'''Tools for making existing classes easily readable/writable to JSON'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

from typing import Any, Callable, ClassVar, Optional, Type, TypeVar, Union
C = TypeVar('C') # generic type for classes

from dataclasses import dataclass, is_dataclass
from functools import update_wrapper, wraps
from inspect import signature, isclass

import json
from pathlib import Path

from ...decorators.functional import allow_string_paths
from .serialize import TypeSerializer, MultiTypeSerializer


# TYPEHINTING AND SERIALIZATION FOR JSONIFIABLE CLASSES 
class JSONifiable: # documentation class, makes type-checking for jsonifiable-modified classes easier
    '''For type-hinting classes which are jsonifiable'''
    pass

def dataclass_serializer_factory(cls : Type[C]) -> TypeSerializer:
    '''For generating a custom TypeSerializer for a JSONifiable dataclass'''
    assert(is_dataclass(cls)) # can enforce modification only to dataclasses (makes behavior a little more natural)

    class DataclassSerializer(TypeSerializer, python_type=cls):
        f'''JSON encoder and decoder for the {cls.__name__} dataclass'''
        @staticmethod
        def encode(python_obj : Path) -> dict[str, Any]:
            '''Extract dictionary of attributes (may need other external converters to be fully serialized)'''
            return python_obj.__dict__

        @staticmethod
        def decode(json_obj : dict[str, Any]) -> C:
            '''Load registered JSONifiable class from a (presumed to be a dictionary during __values__ parse)'''
            json_obj = {
                attr : value
                    for attr, value in json_obj.items()
                        if cls.__dataclass_fields__[attr].init # only pass fields which are allowed to be initialized with values
            }
            return cls(**json_obj)
        
    # dynamically update signatures for readability
    non_generic_name = f'{cls.__name__}Serializer' # dynamically set the name of the serializer to match with the wrapped class
    setattr(DataclassSerializer, '__name__', non_generic_name)
    setattr(DataclassSerializer, '__qualname__', non_generic_name)
    setattr(DataclassSerializer, '__module__', cls.__module__)
        
    return DataclassSerializer


# DECORATOR METHOD FOR MAKING CLASSES JSON-SERIALIZABLE
def make_jsonifiable(cls : Optional[C]=None, type_serializer : Optional[Union[TypeSerializer, MultiTypeSerializer]]=None) -> C:
    '''
    Modify a dataclass to make its attributes writeable-to and readable-from JSON files
    Can optionally specify additional TypeSerializers to support objects with attributes whose types are, by default, not JSON-serializable
    '''
    def jsonifiable_factory(cls : C) -> C:
        '''Factory method, defines new class which inherits from original class and adds serialization methods'''
        assert(is_dataclass(cls)) # can enforce modification only to dataclasses (makes behavior a little more natural)

        multi_serializer = MultiTypeSerializer()
        if type_serializer is not None:
            multi_serializer.add_type_serializer(type_serializer)

        # check if any of the init fields of the dataclass are also JSONifiable, register respective serializers if they are
        for init_param in signature(cls).parameters.values(): # TODO: find a way to have this recognize serializable classes in default containers (e.g. list[JSONifiable])
            if isclass(init_param.annotation) and issubclass(init_param.annotation, JSONifiable): # check if the init field is itself a JSONifiable class
                multi_serializer.add_type_serializer(init_param.annotation.serializer)

        # generate serializable wrapper class
        @wraps(cls, updated=()) # copy over docstring, module, etc; set updated to empty so as to not attempt __dict__updates (classes don;t have these)
        @dataclass
        class WrappedClass(cls, JSONifiable):
            '''Class which inherits from modified class and adds JSON serialization capability'''
            serializer : ClassVar[MultiTypeSerializer] = multi_serializer

            @allow_string_paths
            def to_file(self, save_path : Path) -> None:
                '''Store parameters in a JSON file on disc'''
                assert(save_path.suffix == '.json')
                with save_path.open('w') as dumpfile:
                    json.dump(self, dumpfile, default=self.serializer.encoder_default, indent=4)

            @classmethod
            @allow_string_paths
            def from_file(cls, load_path : Path) -> cls.__class__:
                assert(load_path.suffix == '.json')
                with load_path.open('r') as loadfile:
                    return json.load(loadfile, object_hook=cls.serializer.decoder_hook)

        # !CRITICAL! that the custom serializer be registered for WrappedClass and NOT cls; otherwise, decoded instances will have different type to the parent class
        CustomSerializer = dataclass_serializer_factory(WrappedClass)
        multi_serializer.add_type_serializer(CustomSerializer)
        
        return WrappedClass

    if cls is None:
        return jsonifiable_factory
    return jsonifiable_factory(cls)