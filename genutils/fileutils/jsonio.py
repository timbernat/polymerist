'''Utilities for simplifying the loading and writing of various classes to JSON'''

from typing import Any, Callable, Optional, Union
from abc import ABC, abstractstaticmethod

import json
from pathlib import Path

from ..typetools import T, Args, KWArgs
from ..decorators.functional import allow_string_paths


# JSON-SPECIFIC FUNCTIONS
JSONSerializable = Union[str, bool, int, float, tuple, list, dict] 

@allow_string_paths
def append_to_json(json_path : Path, **kwargs) -> None:
    '''Add an entry to an existing JSON file'''
    with json_path.open('r') as json_file:
        jdat = json.load(json_file)

    jdat.update(**kwargs)

    with json_path.open('w') as json_file:
        jdat = json.checkpoint(jdat, json_file, indent=4)

# JSON CLASSES
class JSONifiable: # TODO : implement encode/decode methods as identity return by default, ditch abstract class behavior
    '''Base class which allows a child class to have its attributes written to and from a JSON file on-disc between interpreter sessions
    Children must implement how dict data (i.e. self.__dict__) is encoded to and decoded from JSON formatted dict'''

    # JSON encoding and decoding
    @staticmethod
    def serialize_json_dict(unser_jdict : dict[Any, Any]) -> dict[str, JSONSerializable]:
        '''For converting selfs __dict__ data into a form that can be serialized to JSON'''
        return unser_jdict # by default, perform no encoding/decoding - child classes can override this 
    
    @staticmethod
    def unserialize_json_dict(ser_jdict : dict[str, JSONSerializable]) -> dict[Any, Any]:
        '''For de-serializing JSON-compatible data into a form that the __init__method can accept'''
        return ser_jdict # by default, perform no encoding/decoding - child classes can override this

    # File I/O
    @allow_string_paths
    def to_file(self, savepath : Path) -> None:
        '''Store parameters in a JSON file on disc'''
        assert(savepath.suffix == '.json')
        with savepath.open('w') as dumpfile:
            json.dump(self.serialize_json_dict(self.__dict__), dumpfile, indent=4)

    @classmethod
    @allow_string_paths
    def from_file(cls, loadpath : Path) -> 'JSONifiable':
        assert(loadpath.suffix == '.json')
        with loadpath.open('r') as loadfile:
            params = json.load(loadfile, object_hook=cls.unserialize_json_dict)

        return cls(**params)
    
    @staticmethod
    def update_checkpoint(funct : Callable[[Any], T]) -> Callable[[Any, Args, KWArgs], T]: # NOTE : this deliberately doesn't have a "self" arg!
        '''Decorator for updating the on-disc checkpoint file after a function updates a Polymer attribute'''
        def update_fn(self, *args, **kwargs) -> Optional[Any]:
            ret_val = funct(self, *args, **kwargs) # need temporary value so update call can be made before returning
            self.to_file()
            return ret_val
        return update_fn

class JSONDict(dict):
    '''Dict subclass which also updates an underlying JSON file - effectively and on-disc dict
    !NOTE! - JSON doesn't support non-string keys, so all keys given will be stringified - plan accordingly!'''
    def __init__(self, json_path : Path, *args, **kwargs):
        if isinstance(json_path, str):
            json_path = Path(json_path) # make input arg a bit more flexible to str input from user end

        if json_path.suffix != '.json':
            raise ValueError(f'The path "{json_path}" does not point to a .json file')
        
        self.json_path : Path = json_path
        if self.json_path.exists():
            try:
                kwargs.update(self._read_file(json_path))
            except json.JSONDecodeError: # catches Paths which point to incorrectly formatted JSONs - TODO: revise terrible except-pass structure
                pass
        else:
            self.json_path.touch()

        super().__init__(*args, **kwargs)
        self._update_file() # ensure file contains current contents post-init

    @staticmethod
    def _read_file(json_path : Path) -> dict:
        with json_path.open('r') as file:
            return json.load(file)        

    def _update_file(self, indent : int=4):
        '''Save current dict contents to JSON file'''
        with self.json_path.open('w') as file:
            json.dump(self, file, indent=indent)

    def __setitem__(self, __key: str, __value: JSONSerializable) -> None:
        super().__setitem__(__key, __value)
        self._update_file()

    def __delitem__(self, __key: str) -> None:
        super().__delitem__(__key)
        self._update_file()