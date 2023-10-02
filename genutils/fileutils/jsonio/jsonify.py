'''Tools for making existing classes easily readable/writable to JSON'''

from typing import Any, Optional, TypeVar, Union
C = TypeVar('C') # generic type for classes
from dataclasses import is_dataclass

import json
from pathlib import Path

from ...decorators.functional import allow_string_paths
from .serialize import TypeSerializer, MultiTypeSerializer


def make_jsonifiable(cls : Optional[C]=None, type_serializer : Optional[Union[TypeSerializer, MultiTypeSerializer]]=None, dataclasses_only : bool=True) -> C:
    '''
    Modify a class to make its attributes writeable-to and readable-from JSON files
    Can optionally specify additional TypeSerializers to support objects with attributes whose types are, by default, not JSON-serializable
    '''
    def jsonifiable_factory(cls : C) -> C:
        '''Factory method, defines new class which inherits from original class and adds serialization methods'''
        if dataclasses_only:
            assert(is_dataclass(cls)) # can enforce modification only to dataclasses (makes behavior a little more natural)

        # assign encoder(s) and decoder(s)
        if type_serializer is None:
            Encoder = None
            object_hook = None
        else:
            class Encoder(json.JSONEncoder): # define new class with custom encoder
                '''Encoder with overriden default'''
                def default(self, python_obj : Any) -> Any:
                    return type_serializer.encoder_default(python_obj)
            object_hook = type_serializer.decoder_hook

        setattr(cls, '_Encoder', Encoder) # 
        setattr(cls, '_decoder_hook', object_hook)

        # define additional methods
        @allow_string_paths
        def to_file(self, save_path : Path) -> None:
            '''Store parameters in a JSON file on disc'''
            assert(save_path.suffix == '.json')
            with save_path.open('w') as dumpfile:
                json.dump(self.__dict__, dumpfile, cls=Encoder, indent=4)
        setattr(cls, 'to_file', to_file) # bind new attribute to class

        @classmethod
        @allow_string_paths
        def from_file(cls, load_path : Path) -> cls:
            assert(load_path.suffix == '.json')
            with load_path.open('r') as loadfile:
                params = json.load(loadfile, object_hook=object_hook)

            return cls(**params)
        setattr(cls, 'from_file', from_file) # bind new attribute to class

        # @staticmethod
        # def update_checkpoint(funct : Callable[[Any], T]) -> Callable[[Any, Args, KWArgs], T]: # NOTE : this deliberately doesn't have a "self" arg!
        #     '''Decorator for updating the on-disc checkpoint file after a function updates a Polymer attribute'''
        #     def update_fn(self, *args, **kwargs) -> Optional[Any]:
        #         ret_val = funct(self, *args, **kwargs) # need temporary value so update call can be made before returning
        #         self.to_file()
        #         return ret_val
        #     return update_fn

        return cls

    if cls is None:
        return jsonifiable_factory
    return jsonifiable_factory(cls)