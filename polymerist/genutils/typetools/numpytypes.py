'''Type aliases specific to numpy arrays'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

from typing import Type, TypeAlias, TypeVar, Union
import numpy.typing as npt

import re
import numpy as np


 # dynamically registering Unions of numeric numpy types
for num_type in (int, float):
    num_type_name = num_type.__name__
    num_type_re = re.compile(rf'\A{num_type_name}\d+$') # attribute name contains type name and numbers and nothing more (e.g. "int64")
    type_variants_present : list[Type] = [
        getattr(np, match.group(0))
            for attr in dir(np)
                if (match := re.match(num_type_re, attr)) is not None
    ]
    globals()[f'NP{num_type_name.capitalize()}'] = Union[tuple(type_variants_present)] # NOTE: can't star-unpack here

# Array shape hinting - TODO: implement support for arithmetic in numeric generics to allow syntax like ndarray[Shape[N//2, M], int], etc.
D = TypeVar('D', bound=int) # for typing arbitrary array dimensionality
N = TypeVar('N', bound=int) # for typing arbitrary array dimension size
M = TypeVar('M', bound=int) # for typing arbitrary array dimension size

DType = TypeVar('DType', bound=np.generic)
Shape : TypeAlias = tuple
