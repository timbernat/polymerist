'''Utilities for categorizing and representing file extensions/suffixes'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

import mimetypes
from collections import defaultdict


# CUSTOM EXCEPTIONS
class FileTypeError(Exception):
    '''Raise when file extension is not valid for a particular application'''
    pass

# GLOBAL REGISTRY OF COMMON EXTENSIONS
EXT_REG = defaultdict(dict) # registry of builtin extensions for common file types
for ext, desc in mimetypes.types_map.items():
    ext_type, app = desc.split('/')
    EXT_REG[ext_type][ext] = app