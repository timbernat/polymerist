'''Tools for making existing classes easily readable/writable to JSON'''

from typing import Any, Optional, Type, TypeVar, Union
C = TypeVar('C') # generic type for classes
from dataclasses import dataclass, is_dataclass
from functools import update_wrapper, wraps

import json
from pathlib import Path

from ...decorators.functional import allow_string_paths
from .serialize import TypeSerializer, MultiTypeSerializer


# TYPEHINTING AND SERIALIZATION FOR JSONIFIABLE CLASSES 
class JSONifiable: # documentation class, makes type-checking for jsonifiable-modified classes easier
    '''For type-hinting classes which are jsonifiable'''
    pass

def jsonifiable_serializer_factory(cls : Type[C]) -> TypeSerializer:
    '''For generating a custom TypeSerializer for a JSONifiable dataclass'''
    assert(is_dataclass(cls)) # can enforce modification only to dataclasses (makes behavior a little more natural)

    class JSONifiableSerializer(TypeSerializer):
        f'''JSON encoder and decoder for the {cls.__name__} dataclass'''
        python_type = cls

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
    nongeneric_name = f'{cls.__name__}Serializer' # dynamically set the name of the serializer to match with the wrapped class
    setattr(JSONifiableSerializer, '__name__', nongeneric_name)
    setattr(JSONifiableSerializer, '__qualname__', nongeneric_name)
    setattr(JSONifiableSerializer, '__module__', cls.__module__)
        
    return JSONifiableSerializer


# DECORATOR METHOD FOR MAKING CLASSES JSON-SERIALIZABLE
def make_jsonifiable(cls : Optional[C]=None, type_serializer : Optional[Union[TypeSerializer, MultiTypeSerializer]]=None) -> C:
    '''
    Modify a dataclass to make its attributes writeable-to and readable-from JSON files
    Can optionally specify additional TypeSerializers to support objects with attributes whose types are, by default, not JSON-serializable
    '''
    def jsonifiable_factory(cls : C) -> C:
        '''Factory method, defines new class which inherits from original class and adds serialization methods'''
        assert(is_dataclass(cls)) # can enforce modification only to dataclasses (makes behavior a little more natural)

        # assign encoder(s) and decoder(s)
        if type_serializer is None:
            Encoder     = None
            object_hook = None
        else:
            class Encoder(json.JSONEncoder): # define new class with custom encoder
                '''Encoder with overriden default'''
                def default(self, python_obj : Any) -> Any:
                    return type_serializer.encoder_default(python_obj)
            object_hook = type_serializer.decoder_hook

        # generate serializable child classmake_jsonifiable @t attempt __dict__updates (classes don't have these)
        @wraps(cls, updated=()) # copy over docstring, module, etc; set updated to empty so as to not attempt __dict__updates (classes don;t have these)
        @dataclass
        class WrappedClass(cls, JSONifiable):
            '''Class which inherits from modified class and adds JSON serialization capability'''
            _Encoder = Encoder
            _object_hook = object_hook

            @allow_string_paths
            def to_file(self, save_path : Path) -> None:
                '''Store parameters in a JSON file on disc'''
                assert(save_path.suffix == '.json')
                with save_path.open('w') as dumpfile:
                    json.dump(self.__dict__, dumpfile, cls=Encoder, indent=4)
            setattr(cls, 'to_file', to_file) # bind new attribute to class

            @classmethod
            @allow_string_paths
            def from_file(cls, load_path : Path) -> cls.__class__:
                assert(load_path.suffix == '.json')
                with load_path.open('r') as loadfile:
                    params = json.load(loadfile, object_hook=object_hook)
                    params = {
                        attr : value
                            for attr, value in params.items()
                                if cls.__dataclass_fields__[attr].init # only pass fields which are allowed to be initialized with values
                    }
                return cls(**params)
        return WrappedClass

    if cls is None:
        return jsonifiable_factory
    return jsonifiable_factory(cls)