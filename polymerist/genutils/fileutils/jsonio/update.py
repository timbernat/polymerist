'''Tools for statically or dynamically updating JSON files'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

import json
from pathlib import Path

from .serialize import JSONSerializable
from ...decorators.functional import allow_string_paths


@allow_string_paths
def append_to_json(json_path : Path, **kwargs) -> None:
    '''Add an entry to an existing JSON file'''
    with json_path.open('r') as json_file:
        jdat = json.load(json_file)

    jdat.update(**kwargs)

    with json_path.open('w') as json_file:
        jdat = json.dump(jdat, json_file, indent=4)

class JSONDict(dict): # TODO : add allow_str_path decorators
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